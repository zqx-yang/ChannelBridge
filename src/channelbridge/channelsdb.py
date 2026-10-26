"""Tunnels from ChannelsDB 2.0, instead of running CAVER yourself.

If the structure you are modelling is a deposited PDB entry, its tunnels may
already have been computed, annotated by **what they lead to**, and published.
Fetching them is faster than running CAVER, needs no CAVER installation, and
comes with a second opinion: every entry carries tunnels from both MOLE and
CAVER, computed independently.

What you get, per tunnel, is a sphere chain -- ``X, Y, Z, Radius,
FreeRadius, Distance`` at each point -- which is the same object CAVER
produces, so it drops straight into :class:`channelbridge.tunnel.Tunnel` and
through the rest of the pipeline.

**The fields are not all what they look like.** ``Distance`` is cumulative
arc length in the MOLE entries and something else entirely in the
CAVER-derived ones, where it is non-monotonic and saturates near 1. So
:attr:`DbTunnel.length` measures the coordinates instead, and
:attr:`DbTunnel.distance_is_arc_length` says whether the reported series
could be used. Read that before trusting any per-point quantity.

What you do NOT get is energetics. ChannelsDB stores geometry and the
physicochemistry of the lining residues; it does not serve the CaverDock
energy profiles from the large-scale annotation. Checked, not assumed: there
is no energy field anywhere in a report.

**The categories are the useful part**, because they say what a tunnel is
*for*:

``ProcognateTunnels``
    to the site of the enzyme's cognate ligand -- the substrate channel, and
    usually the one this pipeline wants
``CofactorTunnels``
    to a bound cofactor. For a pair like GDH/FNR, where the shuttled species
    IS the cofactor, this is the relevant one
``CSATunnels``
    to the catalytic site as annotated in the Catalytic Site Atlas
``ReviewedChannels``
    manually curated

Two cautions worth knowing before relying on it:

* **Coverage is per ENTRY, not per protein.** 4-coumarate:CoA ligase has
  thirteen entries in the PDB and two of them carry tunnels; 5BST has none
  while 5BSR has 25. Search the enzyme, then pick an annotated entry.
* **The coordinates are that entry's.** Using these tunnels means modelling
  on that structure, since the sphere chain is in its frame and nothing
  re-fits it for you.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np

from .tunnel import Tunnel

__all__ = [
    "API",
    "CATEGORIES",
    "DbTunnel",
    "Match",
    "best_match",
    "compare",
    "fetch",
    "read_report",
    "superpose",
    "to_tunnel",
    "tunnels",
]

API = "https://channelsdb2.biodata.ceitec.cz/api/channels/pdb/%s"

#: In the order this pipeline prefers them: the substrate channel first.
CATEGORIES = ("ProcognateTunnels", "CofactorTunnels", "ReviewedChannels",
              "CSATunnels")


@dataclass
class DbTunnel:
    """One tunnel as ChannelsDB reports it."""

    category: str            #: ProcognateTunnels, CofactorTunnels, ...
    method: str              #: "Caver" or "MOLE"
    tunnel_id: str
    xyz: np.ndarray
    radii: np.ndarray        #: the ``Radius`` series
    free_radii: np.ndarray   #: the ``FreeRadius`` series
    distance: np.ndarray     #: arc length from the start point
    properties: Dict[str, float]

    @property
    def label(self) -> str:
        return "%s_%s_%s" % (self.category.replace("Tunnels", "")
                             .replace("Channels", ""), self.method, self.tunnel_id)

    @property
    def length(self) -> float:
        """Arc length along the sphere centres, measured from COORDINATES.

        Not from the report's own ``Distance`` series, which means different
        things in the two methods. MOLE writes cumulative arc length, as you
        would expect: 0.0, 0.086, 0.172, ... 17.354. The CAVER-derived
        entries write something else -- 0.0, 0.827, 0.729, ... 1.000, 1.0002
        -- which is not monotonic and saturates near 1, so it is a per-step
        or normalised quantity and certainly not a length.

        Taking ``Distance[-1]`` therefore reported every ChannelsDB CAVER
        tunnel as about 1 A long. Nothing complained: a 1 A tunnel is a
        plausible-looking number, and the chains themselves were fine, so
        the error only showed up as a selection rule picking two-sphere
        "tunnels" in preference to real ones. Coordinates are the ground
        truth here and are always present, so this measures them.
        """
        if len(self.xyz) < 2:
            return 0.0
        return float(np.linalg.norm(np.diff(self.xyz, axis=0), axis=1).sum())

    @property
    def reported_length(self) -> float:
        """The last value of the report's own ``Distance`` series, as given.

        Kept so the two can be compared rather than one quietly replacing
        the other; see :attr:`distance_is_arc_length`.
        """
        return float(self.distance[-1]) if len(self.distance) else 0.0

    @property
    def distance_is_arc_length(self) -> bool:
        """Whether this tunnel's ``Distance`` series is a cumulative length.

        True for MOLE, false for the CAVER-derived entries. Worth checking
        before using ``Distance`` for anything: the convention is per
        method, not per database, and the report does not say which it is.
        """
        d = np.asarray(self.distance, dtype=float)
        if len(d) < 2:
            return False
        if np.any(np.diff(d) < -1e-9):
            return False
        # bool(), not the numpy comparison: a numpy bool_ fails `is True`,
        # which is how a caller would naturally test this.
        return bool(abs(d[-1] - self.length) <= max(0.5, 0.05 * self.length))

    @property
    def bottleneck(self) -> float:
        """Narrowest point, using ``Radius``.

        Not ``FreeRadius``: the two differ by about 0.1-0.2 A and only one of
        them is what CAVER's own bottleneck means. Both are kept so a caller
        can check rather than take this on trust.
        """
        return float(self.radii.min()) if len(self.radii) else 0.0


def fetch(pdb_id: str, cache_dir: Optional[str] = None, timeout: float = 30.0) -> Dict:
    """Download one entry's report, caching it if asked.

    Cached by PDB id. The database is updated periodically, so a cache is a
    convenience for repeated runs rather than a permanent record; delete the
    file to refresh.
    """
    pdb_id = str(pdb_id).strip().lower()
    if not pdb_id:
        raise ValueError("no PDB id given")
    path = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        path = os.path.join(cache_dir, "channelsdb_%s.json" % pdb_id)
        if os.path.exists(path):
            return read_report(path)
    try:
        with urllib.request.urlopen(API % pdb_id, timeout=timeout) as fh:
            raw = fh.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise LookupError(
                "ChannelsDB has no entry for %s. Coverage is per PDB ENTRY, "
                "not per protein: another entry of the same enzyme may well "
                "be annotated (4CL's 5BST has none, 5BSR has 25)."
                % pdb_id.upper()) from exc
        raise
    if path:
        with open(path, "w", encoding="utf-8") as out:
            out.write(raw)
    return json.loads(raw)


def read_report(path: str) -> Dict:
    """A report already on disk, from the site's download button or a cache."""
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def tunnels(report: Dict, categories: Sequence[str] = CATEGORIES,
            method: Optional[str] = None) -> List[DbTunnel]:
    """Every tunnel in a report, in category order.

    ``method`` filters to ``"Caver"`` or ``"MOLE"``. Leaving it unset gives
    both, which is the point: two independent computations of the same
    channel, and where they disagree that is a measure of how determined the
    channel really is.
    """
    out: List[DbTunnel] = []
    channels = report.get("Channels", report)
    for cat in categories:
        for key, entries in channels.items():
            if not key.startswith(cat) or not entries:
                continue
            got_method = key.split("_")[-1]
            if method and got_method.lower() != method.lower():
                continue
            for entry in entries:
                profile = entry.get("Profile") or []
                if len(profile) < 2:
                    continue
                out.append(DbTunnel(
                    category=cat, method=got_method,
                    tunnel_id=str(entry.get("Id", "?")),
                    xyz=np.array([[p["X"], p["Y"], p["Z"]] for p in profile], float),
                    radii=np.array([p.get("Radius", 0.0) for p in profile], float),
                    free_radii=np.array([p.get("FreeRadius", 0.0) for p in profile], float),
                    distance=np.array([p.get("Distance", 0.0) for p in profile], float),
                    properties=dict(entry.get("Properties") or {})))
    return out


def superpose(mobile, target):
    """Kabsch fit of *mobile* onto *target*, on coordinates ALREADY PAIRED.

    Returns ``(rotation, translation, rms)``. The pairing is the caller's
    job and is the part that matters: a whole-chain fit answers "is this the
    same fold", while what a tunnel comparison needs is "is this the same
    SITE". Those differ by a few angstroms on a homolog, and a tunnel
    comparison inherits the difference -- so pair on the site's residues
    when you have them.
    """
    P = np.atleast_2d(np.asarray(mobile, float))
    Q = np.atleast_2d(np.asarray(target, float))
    if P.shape != Q.shape or P.shape[-1] != 3:
        raise ValueError("need matching (N, 3) arrays, got %r and %r"
                         % (P.shape, Q.shape))
    if len(P) < 3:
        raise ValueError("need at least three paired atoms, got %d" % len(P))
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = qc - R @ pc
    rms = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
    return R, t, rms


@dataclass
class Match:
    """How well one of our tunnels corresponds to one of theirs.

    **Not "closest approach".** Every CAVER cluster of a run starts from the
    same starting point, so all of them share their first spheres; the
    minimum distance between any two is therefore set by that shared origin
    and is the same for cluster 1 and cluster 9. Using it to pick a cluster
    was wrong and said nothing.

    What distinguishes clusters is where they come OUT. So the comparison is
    exit-led, with a shape check behind it:

    ``exit_gap``
        distance between the two exits. The discriminating number.
    ``hausdorff``
        the worst case over both directions: the furthest either chain
        strays from the other. Catches two routes that share an exit but
        take different paths to it.
    ``coverage``
        fraction of OUR spheres within ``tol`` of their chain.
    ``db_coverage``
        fraction of THEIR spheres within ``tol`` of ours. The other
        direction, and not redundant -- see :meth:`recovery`.
    """

    db_label: str
    exit_gap: float
    hausdorff: float
    coverage: float
    db_coverage: float = 0.0

    def verdict(self, exit_tol: float = 3.0, coverage_min: float = 0.5) -> str:
        """Are these the same tunnel? Exit-led, for comparing like with like.

        Use this when the two chains have comparable extent -- two clusters
        of one CAVER run, or our tunnel against a database one when both
        reach the surface. It is the WRONG test when one chain is
        deliberately longer than the other; :meth:`recovery` is for that.
        """
        if self.exit_gap <= exit_tol and self.coverage >= coverage_min:
            return "same tunnel"
        if self.exit_gap <= exit_tol:
            return "same exit, different route"
        if self.coverage >= coverage_min:
            return "overlapping, different exit"
        return "different tunnel"

    def recovery(self, good: float = 0.8, partial: float = 0.5) -> str:
        """Did OUR chain find THEIR route, allowing ours to run further?

        The question when a database chain is used to seed a CAVER run
        because it stopped before the surface. There, ours is longer by
        design, so both of :meth:`verdict`'s inputs are biased against a
        correct answer: the exits cannot coincide, and ``coverage`` counts
        our extra spheres as misses. Judging those seeded runs by
        ``verdict`` called six recovered routes "different tunnel" --
        including one the chains agree on over 100 % of the database's
        length.

        So this reads ``db_coverage`` only: how much of the chain they
        published lies on the chain we computed.
        """
        if self.db_coverage >= good:
            return "route recovered"
        if self.db_coverage >= partial:
            return "route partly recovered"
        return "different route"


def compare(ours, db: "DbTunnel", rotation=None, translation=None,
            tol: float = 2.0) -> Match:
    """Compare one of our tunnels with one of theirs, in OUR frame.

    ``rotation``/``translation`` move THEIR chain into our frame -- the
    superposition of their structure onto ours. Without it the two are in
    different coordinate systems and every number is meaningless, which is
    not detectable from the numbers themselves: a 19-40 A "disagreement" is
    what an unsuperposed comparison looks like.
    """
    a = np.asarray(ours.xyz, float)
    b = np.asarray(db.xyz, float)
    if rotation is not None:
        b = b @ np.asarray(rotation, float).T
    if translation is not None:
        b = b + np.asarray(translation, float)
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return Match(db_label=db.label,
                 exit_gap=float(np.linalg.norm(a[-1] - b[-1])),
                 hausdorff=float(max(d.min(axis=1).max(), d.min(axis=0).max())),
                 coverage=float((d.min(axis=1) <= tol).mean()),
                 db_coverage=float((d.min(axis=0) <= tol).mean()))


def best_match(ours, db_tunnels: Sequence["DbTunnel"], rotation=None,
               translation=None, tol: float = 2.0) -> Optional[Match]:
    """The database tunnel whose EXIT is nearest ours.

    Exit-led because that is what separates one cluster from another; ties
    on the exit are broken by the shape term.
    """
    scored = [compare(ours, t, rotation, translation, tol) for t in db_tunnels]
    if not scored:
        return None
    return min(scored, key=lambda m: (round(m.exit_gap, 1), m.hausdorff))


def to_tunnel(db: DbTunnel, *, use_free_radius: bool = False) -> Tunnel:
    """Convert to the pipeline's own :class:`Tunnel`.

    The chain is assumed to run **inner to outer**, as both MOLE and CAVER
    write it: point 0 is at the start point and distance increases outward.
    Everything downstream depends on that order -- the exit is the last
    sphere -- so :meth:`Tunnel.oriented` should still be applied if there is
    any doubt, exactly as it is for a CAVER chain.
    """
    radii = db.free_radii if use_free_radius else db.radii
    return Tunnel(xyz=np.asarray(db.xyz, float), radii=np.asarray(radii, float),
                  source="ChannelsDB %s" % db.label)
