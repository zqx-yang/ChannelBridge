"""CAVER tunnels as geometry, made comparable between CAVER builds.

Why this module is not a thin reader
------------------------------------

A CAVER tunnel cluster is a chain of probe spheres. Three properties the
pipeline depends on -- the exit position, the tunnel span, and the direction
the tunnel points at its mouth -- are **not** intrinsic to the protein. They
depend on how finely the build discretised the path and on how far past the
protein surface it kept going. Measured on the same GDH structure, the same
config and the same seed, with two CAVER builds:

================================  =================  ===================
property                          build A (earlier)  build B (cluster)
================================  =================  ===================
spheres in cluster 1              16                 145
mean sphere spacing               0.98 A             0.17 A
span (first to last sphere)       12.1 A             19.1 A
burial at the outermost sphere    21                 0
bottleneck radius                 2.915 A            2.910 A
O_IC arm = last sphere step       0.91 A             0.26 A
================================  =================  ===================

It is the same tunnel: build B's chain passes 0.078 A from build A's exit,
with the same burial there, and then runs a further ~9 A into bulk solvent.
The bottleneck agrees, because it is found in the interior where the shell
parameters have no influence.

So two normalisations make the geometry portable, and both are on by default:

:meth:`Tunnel.truncate_at_surface`
    Scan outward from the buried inner end and stop where the chain leaves
    the solvent-excluded surface, judged by burial. On the data above this
    lands 1.4 A from the independently derived exit of the other build,
    instead of 9 A out in solvent.

:meth:`Tunnel.resample`
    Re-space the chain at a fixed arc length (1.0 A by default, the density
    the published iMARS construction was calibrated on).

and one definition removes the remaining build dependence:

:meth:`Tunnel.axis_point`
    P1 for the O_IC angle is taken a fixed arc length inward from the exit,
    interpolated along the path, instead of "the second-to-last sphere". On a
    densely sampled chain the sphere step is a fraction of an angstrom and
    sphere-placement jitter dominates the direction.

Burial
------

``burial`` is the number of protein atoms within 8 A of a point, counted
against the atoms CAVER actually searched -- the input **minus**
``exclude_residue_names``. Counting a cofactor CAVER never saw inflates it.
Calibrated on GDH and FNR with a 1.5 A grid:

=============  ==============================================
burial         meaning
=============  ==============================================
< 15           past the surface: convex, fully exposed
15 - 40        a real mouth on the surface (the wanted value)
40 - 60        recessed: a groove or shallow pocket
> 60           still inside the protein, not an exit at all
=============  ==============================================
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, replace
from typing import Sequence, Tuple

import numpy as np

from .io.pdb import read_pdb

__all__ = [
    "BURIAL_BANDS",
    "BURIAL_RADIUS",
    "BURIAL_SURFACE_MIN",
    "Tunnel",
    "burial_at",
    "burial_verdict",
    "fitted_axis",
]

#: Radius of the sphere burial is counted in, in angstroms.
BURIAL_RADIUS = 8.0

#: Lower edge of the "a real mouth on the surface" band. Used as the default
#: truncation threshold: below this the chain is past the surface.
BURIAL_SURFACE_MIN = 15

#: ``(upper_exclusive, verdict)`` in order; the last entry is open-ended.
BURIAL_BANDS = (
    (15, "past the surface"),
    (40, "surface mouth"),
    (60, "recessed"),
    (None, "buried"),
)


def burial_verdict(burial: float) -> str:
    """The calibrated band a burial count falls in."""
    for upper, verdict in BURIAL_BANDS:
        if upper is None or burial < upper:
            return verdict
    raise AssertionError("unreachable: BURIAL_BANDS must end with an open band")


def burial_at(points, protein_xyz, radius: float = BURIAL_RADIUS) -> np.ndarray:
    """Protein atoms within *radius* of each point.

    *protein_xyz* must already have the CAVER ``exclude_residue_names``
    residues removed; see :meth:`channelbridge.io.pdb.AtomArray.exclude_resnames`.
    """
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    prot = np.asarray(protein_xyz, dtype=float).reshape(-1, 3)
    if pts.shape[-1] != 3:
        raise ValueError("points must be (N, 3), got %r" % (pts.shape,))
    if len(prot) == 0:
        return np.zeros(len(pts), dtype=int)
    out = np.empty(len(pts), dtype=int)
    # Chunked to keep the pairwise block bounded for a tetramer (~8k atoms)
    # against a dense chain (~150 spheres).
    for i in range(0, len(pts), 64):
        block = pts[i : i + 64]
        d = np.linalg.norm(prot[None, :, :] - block[:, None, :], axis=2)
        out[i : i + 64] = (d <= radius).sum(axis=1)
    return out


def fitted_axis(xyz, n_spheres: int = 5, *, source: str = "<points>") -> np.ndarray:
    """Outward mouth direction fitted over the last *n_spheres* of a chain.

    A free function because two callers need it and neither should own a
    second copy: :meth:`Tunnel.fitted_axis`, and the stage-3 shim, which has
    the sphere coordinates from PyMOL rather than a Tunnel. It WAS duplicated,
    and the inline copy was the one nothing tested.

    NOT the same vector as the O_IC arm, and the difference is not noise: on
    these tunnels the two disagree by 9-12 degrees. The arm is a short
    baseline because O_IC is defined on one; this is a fit over several
    spheres, because consecutive CAVER spheres sit ~1 A apart and a single
    pair is dominated by placement jitter. Steady enough to classify a whole
    surface against, which is what the residue classes need -- using the pair
    instead made the ported stage 3 disagree with its reference.
    """
    pts = np.asarray(xyz, dtype=float)
    n = max(2, int(n_spheres))
    tail = pts[-min(len(pts), n):]
    if len(tail) < 2:
        raise ValueError("tunnel %s: need at least two spheres to fit an axis" % source)
    u = np.linalg.svd(tail - tail.mean(axis=0))[2][0]
    if float(u @ (pts[-1] - pts[0])) < 0:
        u = -u
    return u / float(np.linalg.norm(u))


@dataclass(frozen=True)
class Tunnel:
    """A chain of CAVER probe spheres, ordered **inner to outer**.

    Attributes
    ----------
    xyz:
        ``(N, 3)`` sphere centres.
    radii:
        ``(N,)`` probe radii, in angstroms.
    source:
        Where it came from, for provenance in output tables.
    history:
        The normalisations already applied, in order. Carried into reports so
        a D_IC can be traced to the processing that produced it.
    """

    xyz: np.ndarray
    radii: np.ndarray
    source: str = "<memory>"
    history: Tuple[str, ...] = ()

    # ----------------------------------------------------------- construction

    def __post_init__(self):
        xyz = np.asarray(self.xyz, dtype=float).reshape(-1, 3)
        radii = np.asarray(self.radii, dtype=float).reshape(-1)
        if len(xyz) != len(radii):
            raise ValueError("got %d centres but %d radii" % (len(xyz), len(radii)))
        if len(xyz) < 2:
            raise ValueError(
                "a tunnel needs at least 2 spheres to have a direction, got %d" % len(xyz)
            )
        if not np.all(np.isfinite(xyz)) or not np.all(np.isfinite(radii)):
            raise ValueError("tunnel contains non-finite values")
        object.__setattr__(self, "xyz", xyz)
        object.__setattr__(self, "radii", radii)

    @classmethod
    def from_caver_cluster(cls, path, *, model="first") -> "Tunnel":
        """Read ``data/clusters/tun_cl_*.pdb``.

        Each MODEL in such a file is one tunnel of the cluster; with
        ``one_tunnel_in_snapshot cheapest`` and a single snapshot there is
        exactly one. The radius comes from the B-factor column.

        The sphere order as written is **not** trusted -- call
        :meth:`oriented` with the protein before using the ends.
        """
        atoms = read_pdb(path, models="first" if model == "first" else "all")
        if len(atoms) < 2:
            raise ValueError(
                "%s holds %d sphere(s); a tunnel needs at least 2. An empty or "
                "near-empty cluster file usually means CAVER found no pathway "
                "(check the starting point is not inside an atom)." % (path, len(atoms))
            )
        radii = atoms.radii
        if np.allclose(radii, 0.0):
            raise ValueError(
                "%s has no radii in the B-factor column. CAVER writes probe "
                "radii there; a file of zeros is not a CAVER tunnel." % (path,)
            )
        return cls(atoms.xyz.copy(), radii.copy(), source=str(path))

    # --------------------------------------------------------------- geometry

    def __len__(self) -> int:
        return len(self.xyz)

    @property
    def inner(self) -> np.ndarray:
        """The innermost sphere centre (deepest in the protein)."""
        return self.xyz[0]

    @property
    def exit(self) -> np.ndarray:
        """The outermost sphere centre: P2 of the iMARS construction."""
        return self.xyz[-1]

    @property
    def bottleneck(self) -> float:
        """The narrowest probe radius along the chain, in angstroms.

        Quotable across builds: the bottleneck sits in the interior, where
        the shell parameters have no influence.
        """
        return float(self.radii.min())

    @property
    def span(self) -> float:
        """Straight-line distance from the inner end to the exit."""
        return float(np.linalg.norm(self.exit - self.inner))

    @property
    def path_length(self) -> float:
        """Arc length along the sphere chain."""
        return float(self._arc()[-1])

    @property
    def mean_spacing(self) -> float:
        """Mean distance between consecutive spheres: the sampling density."""
        return float(self.path_length / (len(self) - 1))

    def _arc(self) -> np.ndarray:
        """Cumulative arc length, ``(N,)``, starting at 0."""
        steps = np.linalg.norm(np.diff(self.xyz, axis=0), axis=1)
        return np.concatenate([[0.0], np.cumsum(steps)])

    # ---------------------------------------------------------- normalisation

    def oriented(self, protein_xyz=None, *, radius: float = BURIAL_RADIUS) -> "Tunnel":
        """Return the chain ordered inner to outer.

        With *protein_xyz* the ends are told apart by **burial**: the inner
        end is the more buried one. That is the reliable test. Without it,
        the fallback is distance from the protein centroid, which fails on a
        tunnel that runs tangentially or on a multimer whose centroid sits in
        the middle of the assembly rather than inside the protomer.
        """
        if protein_xyz is not None:
            prot = np.asarray(protein_xyz, dtype=float).reshape(-1, 3)
            b = burial_at(np.vstack([self.xyz[0], self.xyz[-1]]), prot, radius=radius)
            flip = b[0] < b[-1]
            how = "burial"
        else:
            centre = self.xyz.mean(axis=0)
            flip = np.linalg.norm(self.xyz[0] - centre) > np.linalg.norm(self.xyz[-1] - centre)
            how = "centroid-distance"
        if not flip:
            return replace(self, history=(*self.history, "oriented(%s, kept)" % how))
        return replace(
            self,
            xyz=self.xyz[::-1].copy(),
            radii=self.radii[::-1].copy(),
            history=(*self.history, "oriented(%s, reversed)" % how),
        )

    def burial_profile(self, protein_xyz, *, radius: float = BURIAL_RADIUS) -> np.ndarray:
        """Burial at every sphere, ``(N,)``."""
        return burial_at(self.xyz, protein_xyz, radius=radius)

    def burial_exit(self, protein_xyz, *, radius: float = BURIAL_RADIUS) -> int:
        """Burial at the exit sphere."""
        return int(burial_at(self.exit, protein_xyz, radius=radius)[0])

    def truncate_at_surface(
        self,
        protein_xyz,
        *,
        burial_min: int = BURIAL_SURFACE_MIN,
        radius: float = BURIAL_RADIUS,
        persist: float = 2.0,
        keep_min_spheres: int = 2,
    ) -> "Tunnel":
        """Drop the outer spheres that have left the protein surface.

        The exit is taken where the tunnel **first** leaves the surface,
        scanning outward from the buried inner end: the last sphere before
        burial drops below *burial_min* and stays there.

        Scanning outward rather than inward from the far end matters. A
        neighbouring subunit, a symmetry copy or a crystallisation additive
        sitting in the solvent beyond the mouth also raises burial, and
        scanning inward from the end would accept it and place the exit out
        past the real mouth. The tunnel starts buried by construction, so the
        first crossing from the inside is the mouth of the protein the tunnel
        belongs to.

        Parameters
        ----------
        persist:
            A drop must hold for this many angstroms of arc length to count
            as the exit, so a single noisy sphere cannot truncate the chain
            early. A dip shorter than this is stepped over.

        A chain whose inner end is already below *burial_min*, or which never
        reaches it, is returned unchanged with a warning: that means the
        tunnel does not start inside the protein, which is a result about the
        structure and not something to clip silently.
        """
        prot = np.asarray(protein_xyz, dtype=float).reshape(-1, 3)
        b = self.burial_profile(prot, radius=radius)
        n = len(self)
        s = self._arc()

        if int(b.max()) < burial_min:
            warnings.warn(
                "tunnel %s never reaches burial %d (max %d); not truncating. The "
                "whole chain lies past the solvent-excluded surface."
                % (self.source, burial_min, int(b.max())),
                stacklevel=2,
            )
            return replace(self, history=(*self.history, "truncate_at_surface(no-crossing)"))

        if b[0] < burial_min:
            warnings.warn(
                "tunnel %s starts at burial %d, below %d: the inner end is not "
                "inside the protein, so there is no first crossing to find. Not "
                "truncating -- check the sphere order and the exclude list."
                % (self.source, int(b[0]), burial_min),
                stacklevel=2,
            )
            return replace(self, history=(*self.history, "truncate_at_surface(inner-not-buried)"))

        last_good = 0
        i = 1
        while i < n:
            if b[i] >= burial_min:
                last_good = i
                i += 1
                continue
            # A drop at i. Does it hold for `persist` angstroms?
            j = i
            resumed = False
            while j < n and (s[j] - s[i]) <= persist:
                if b[j] >= burial_min:
                    resumed = True
                    break
                j += 1
            if not resumed:
                break
            last_good = j
            i = j + 1

        if last_good + 1 < keep_min_spheres:
            warnings.warn(
                "tunnel %s: the surface crossing at sphere %d would leave fewer "
                "than %d spheres; not truncating."
                % (self.source, last_good, keep_min_spheres),
                stacklevel=2,
            )
            return replace(self, history=(*self.history, "truncate_at_surface(too-short)"))
        if last_good == n - 1:
            return replace(
                self,
                history=(*self.history, "truncate_at_surface(b>=%d, nothing to drop)" % burial_min),
            )
        return replace(
            self,
            xyz=self.xyz[: last_good + 1].copy(),
            radii=self.radii[: last_good + 1].copy(),
            history=(*self.history, "truncate_at_surface(b>=%d, dropped %d outer)" % (burial_min, n - 1 - last_good)),
        )

    def resample(self, spacing: float = 1.0) -> "Tunnel":
        """Re-space the chain at approximately *spacing* angstroms.

        Both endpoints are preserved exactly, so the exit does not move; the
        achieved spacing is ``path_length / (n - 1)``, the nearest value that
        divides the path into whole steps. Radii are interpolated linearly
        along arc length.

        The default 1.0 A is the sampling density the published iMARS
        construction was calibrated on.
        """
        if spacing <= 0:
            raise ValueError("spacing must be positive, got %r" % (spacing,))
        s = self._arc()
        total = float(s[-1])
        if total == 0.0:
            raise ValueError("tunnel %s has zero path length" % (self.source,))
        n = max(2, round(total / spacing) + 1)
        targets = np.linspace(0.0, total, n)
        xyz = np.column_stack([np.interp(targets, s, self.xyz[:, k]) for k in range(3)])
        radii = np.interp(targets, s, self.radii)
        # Endpoints exactly, not within floating-point interpolation error.
        xyz[0] = self.xyz[0]
        xyz[-1] = self.xyz[-1]
        radii[0] = self.radii[0]
        radii[-1] = self.radii[-1]
        return replace(
            self,
            xyz=xyz,
            radii=radii,
            history=(*self.history, "resample(%.3g A -> %d spheres at %.3g A)" % (spacing, n, total / (n - 1))),
        )

    def normalised(
        self,
        protein_xyz,
        *,
        burial_min: int = BURIAL_SURFACE_MIN,
        spacing: float = 1.0,
        radius: float = BURIAL_RADIUS,
    ) -> "Tunnel":
        """Orient, truncate at the surface, then resample: the default chain.

        This is what makes D_IC, O_IC, span and burial_exit comparable
        between CAVER builds, and between a static structure and MD frames.
        """
        return (
            self.oriented(protein_xyz, radius=radius)
            .truncate_at_surface(protein_xyz, burial_min=burial_min, radius=radius)
            .resample(spacing)
        )

    # ------------------------------------------------------- the O_IC arm

    def axis_point(self, arm: float = 1.0, *, strict: bool = True) -> np.ndarray:
        """P1: the point *arm* angstroms inward from the exit, along the path.

        *arm* is **arc length** along the sphere chain, not straight-line
        distance, so on a curved tunnel the chord ``|exit - axis_point|`` is
        slightly shorter than *arm* (see :meth:`arm_chord`). Arc length is
        used because it is monotone and well defined whatever the sampling
        density; the two agree to about 1% over a 1 A arm on real tunnels.

        Interpolated along the polyline, so it does not depend on the sphere
        spacing. This replaces "the second-to-last sphere", which on a
        densely sampled chain gives an arm of a fraction of an angstrom where
        sphere-placement jitter dominates the direction.

        Parameters
        ----------
        strict:
            If the chain is shorter than *arm*, raise (default) rather than
            quietly fall back to the inner end and report an angle measured
            over an unknown baseline.
        """
        if arm <= 0:
            raise ValueError("arm must be positive, got %r" % (arm,))
        s = self._arc()
        total = float(s[-1])
        if arm > total:
            if strict:
                raise ValueError(
                    "tunnel %s is %.2f A long, shorter than the requested %.2f A "
                    "O_IC arm. Lower the arm, or use a longer tunnel: an angle "
                    "measured over the whole of a short chain is not comparable "
                    "with one measured over an arm." % (self.source, total, arm)
                )
            return self.xyz[0].copy()
        target = total - arm
        return np.array([np.interp(target, s, self.xyz[:, k]) for k in range(3)])

    def arm_chord(self, arm: float = 1.0, *, strict: bool = True) -> float:
        """Straight-line length of the O_IC arm, ``<= arm`` on a curved tunnel."""
        return float(np.linalg.norm(self.exit - self.axis_point(arm, strict=strict)))

    def axis(self, arm: float = 1.0, *, strict: bool = True) -> np.ndarray:
        """Unit vector pointing outward along the mouth, over *arm* angstroms."""
        v = self.exit - self.axis_point(arm, strict=strict)
        n = float(np.linalg.norm(v))
        if n == 0.0:
            raise ValueError("tunnel %s: zero-length axis over a %.2f A arm" % (self.source, arm))
        return v / n

    def fitted_axis(self, n_spheres: int = 5) -> np.ndarray:
        """Outward mouth direction, fitted over the last *n_spheres*.

        NOT the same vector as :meth:`axis`, and the difference is not noise:
        on these tunnels the two disagree by 9-12 degrees.

        :meth:`axis` is the O_IC arm -- exit minus a point one arm-length
        back -- and it is what the iMARS angle is defined on, so it must stay
        a short baseline however jittery that makes it. This one is the first
        principal component of the outer spheres, sign-corrected to point out
        of the protein. Consecutive CAVER spheres sit ~1 A apart, so a single
        pair is noisy; a fit over several is steady enough to classify a
        whole surface against, which is what the residue classes need.

        Use this for geometry ABOUT the mouth (is a residue in the shadow,
        how far off the axis does it sit). Use :meth:`axis` for O_IC itself.
        """
        return fitted_axis(self.xyz, n_spheres, source=self.source)

    def axis_spread(self, arms: Sequence[float] = (0.5, 1.0, 2.0), *, strict: bool = False) -> float:
        """Largest angle, in degrees, between the mouth axes over *arms*.

        A quality flag, not a correction. The mouth direction of a tunnel
        that curves as it opens genuinely depends on the length scale it is
        measured over, so O_IC inherits that ambiguity. Measured on the GDH
        primary route, two CAVER builds normalised identically still disagree
        by 15 deg at a 1 A arm and 36 deg at 2 A, while agreeing to 5 deg at
        8 A. A large spread means the reported O_IC should be quoted with
        that uncertainty: at the published coefficients, 6.6 deg of O_IC is
        worth 1 A of D_IC.
        """
        axes = []
        for a in arms:
            try:
                axes.append(self.axis(a, strict=strict))
            except ValueError:
                continue
        if len(axes) < 2:
            return 0.0
        worst = 0.0
        for i in range(len(axes)):
            for j in range(i + 1, len(axes)):
                cos = float(np.clip(np.dot(axes[i], axes[j]), -1.0, 1.0))
                worst = max(worst, float(np.degrees(np.arccos(cos))))
        return worst

    # ----------------------------------------------------------------- output

    def write_caver_pdb(self, path, *, resname: str = "FIL", chain: str = "T") -> int:
        """Write the chain in CAVER's tunnel-cluster format, inner end first.

        One sphere per ATOM record, radius in the **B-factor** column and the
        occupancy column left blank, which is what CAVER itself writes and
        what downstream readers expect. LF line endings unconditionally.

        The point of this is to hand a *normalised* chain to a tool that reads
        CAVER output directly, without that tool having to change. A chain
        resampled to 1.0 A also makes the published "second-to-last sphere"
        O_IC construction equal a 1.0 A arm, so a reader implementing the
        paper literally gets the normalised answer.

        Returns the number of spheres written.
        """
        with open(path, "w", encoding="ascii", newline="\n") as fh:
            fh.write("MODEL        1\n")
            for i, (xyz, r) in enumerate(zip(self.xyz, self.radii), start=1):
                fh.write(
                    "ATOM  %5d  H   %-3s %1s%4d    %8.3f%8.3f%8.3f      %6.2f\n"
                    % (i, resname, chain, 1, xyz[0], xyz[1], xyz[2], r)
                )
            fh.write("ENDMDL\nEND\n")
        return len(self)

    # ------------------------------------------------------------- reporting

    def describe(self, protein_xyz=None) -> dict:
        """A flat dict of the numbers worth putting in a table."""
        out = {
            "source": self.source,
            "spheres": len(self),
            "span_A": round(self.span, 3),
            "path_length_A": round(self.path_length, 3),
            "mean_spacing_A": round(self.mean_spacing, 3),
            "bottleneck_A": round(self.bottleneck, 3),
            "exit_x": round(float(self.exit[0]), 3),
            "exit_y": round(float(self.exit[1]), 3),
            "exit_z": round(float(self.exit[2]), 3),
            "oic_arm_spread_deg": round(self.axis_spread(), 2),
            "history": "; ".join(self.history) or "raw",
        }
        if protein_xyz is not None:
            b = self.burial_exit(protein_xyz)
            out["burial_exit"] = b
            out["verdict"] = burial_verdict(b)
        return out
