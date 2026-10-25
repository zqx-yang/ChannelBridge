"""Which poses get a .pse, an SVG and a PDB -- and what goes in each.

Figures are deliberately a **separate step from the DO sampling**. Sampling
produces a ranked list of poses; only after seeing the scores can anyone say
which are worth drawing, or whether any are. So this module takes the pose
table and a few short strings and works out exactly what to build.

The selection language
----------------------

Each of ``--pse``, ``--svg`` and ``--pdb`` takes a comma-separated list. Terms
may be mixed freely, duplicates collapse, and the result keeps rank order.

===================  =====================================================
``none`` / empty     nothing
``all``              every pose
``5`` or ``top5``    the best 5 by DO
``rank3`` / ``#3``   the single pose at rank 3
``rank2-4``          ranks 2, 3 and 4
``geom007``          that pose by name
===================  =====================================================

So ``top3,geom010,rank7`` is three best, plus one by name, plus one by rank.

The dependency
--------------

An SVG is traced from renders of a session, so **an SVG needs its .pse**. Ask
for an SVG of a pose whose session was neither requested nor already on disk
and the session is added for you, and said so in the plan -- rather than
failing, or silently producing nothing.

Scene variants
--------------

``render`` is optional and takes
``none | pair | assembly | both``:

* ``pair``  -- the two posed protomers alone (the ``_complex`` structure)
* ``assembly`` -- the full assembly with symmetry copies (``_occupancy``)
* ``both``     -- both, as separate scenes sharing one camera
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence

__all__ = [
    "VARIANTS",
    "FigurePlan",
    "Pose",
    "parse_spec",
    "parse_variant",
    "rank_poses",
    "resolve",
]

#: Accepted values for the optional ``render`` scene.
VARIANTS = ("none", "pair", "assembly", "both")

_RANK_RANGE = re.compile(r"^(?:rank|#)(\d+)\s*-\s*(\d+)$", re.I)
_RANK_ONE = re.compile(r"^(?:rank|#)(\d+)$", re.I)
_TOP = re.compile(r"^top[:\s]*(\d+)$", re.I)
_INT = re.compile(r"^\d+$")


@dataclass(frozen=True)
class Pose:
    """One sampled pose: its id, its DO, and its 1-based rank by DO."""

    pose_id: str
    do: float
    rank: int


def rank_poses(poses: Iterable) -> List[Pose]:
    """Rank ``(pose_id, do)`` pairs by DO descending.

    Ties break on pose id so the ranking is stable between runs -- a figure
    called "rank 3" must mean the same pose tomorrow.
    """
    items = []
    for p in poses:
        if isinstance(p, Pose):
            items.append((p.pose_id, float(p.do)))
        else:
            pid, do = p
            items.append((str(pid), float(do)))
    seen = {}
    for pid, _ in items:
        if pid in seen:
            raise ValueError("duplicate pose id %r in the pose table" % pid)
        seen[pid] = True
    ordered = sorted(items, key=lambda t: (-t[1], t[0]))
    return [Pose(pid, do, i) for i, (pid, do) in enumerate(ordered, start=1)]


def _by_id(ranked: Sequence[Pose]) -> Dict[str, Pose]:
    return {p.pose_id: p for p in ranked}


def _suggest(term: str, known: Sequence[str]) -> str:
    """A 'did you mean' for an unrecognised id, since these are hand-typed."""
    t = term.lower()
    close = [k for k in known if k.lower().startswith(t[:4])] or \
            [k for k in known if t in k.lower()]
    return ("; did you mean %s?" % ", ".join(close[:3])) if close else ""


def parse_spec(spec, ranked: Sequence[Pose], what: str = "selection") -> List[str]:
    """Resolve one selection string to pose ids, in rank order.

    Returns ids, de-duplicated, ordered by rank (not by the order typed) so
    that outputs are produced best-first whatever the user wrote.
    """
    if spec is None:
        return []
    if isinstance(spec, int) and not isinstance(spec, bool):
        spec = str(spec)
    s = str(spec).strip()
    if not s or s.lower() == "none":
        return []

    ids = _by_id(ranked)
    n = len(ranked)
    chosen = set()

    for raw in s.split(","):
        term = raw.strip()
        if not term:
            continue
        low = term.lower()
        if low == "none":
            continue
        if low == "all":
            chosen.update(p.pose_id for p in ranked)
            continue

        m = _RANK_RANGE.match(term)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            if lo < 1 or hi < lo:
                raise ValueError(
                    "%s: %r is not a usable rank range (ranks start at 1 and "
                    "the first number must not exceed the second)" % (what, term)
                )
            if lo > n:
                raise ValueError(
                    "%s: rank %d was asked for but there are only %d pose(s)"
                    % (what, lo, n)
                )
            chosen.update(p.pose_id for p in ranked[lo - 1 : min(hi, n)])
            continue

        m = _RANK_ONE.match(term)
        if m:
            r = int(m.group(1))
            if r < 1:
                raise ValueError("%s: ranks start at 1, got %r" % (what, term))
            if r > n:
                raise ValueError(
                    "%s: rank %d was asked for but there are only %d pose(s)"
                    % (what, r, n)
                )
            chosen.add(ranked[r - 1].pose_id)
            continue

        m = _TOP.match(term) or (_INT.match(term) and re.match(r"^(\d+)$", term))
        if m:
            k = int(m.group(1))
            chosen.update(p.pose_id for p in ranked[:k])
            continue

        if term in ids:
            chosen.add(term)
            continue

        raise ValueError(
            "%s: cannot read %r. Use a count (5, top5), a rank (rank3, #3, "
            "rank2-4), a pose id, 'all' or 'none'%s"
            % (what, term, _suggest(term, [p.pose_id for p in ranked]))
        )

    return [p.pose_id for p in ranked if p.pose_id in chosen]


def parse_variant(spec, what: str = "scene") -> str:
    """Validate a ``none|pair|assembly|both`` option.

    ``pair`` is one protomer of each enzyme -- the two chains the DO geometry
    is defined on. ``assembly`` is every protomer of both. The earlier names
    ``monomer``/``multimer`` are accepted as aliases but were misleading: the
    "monomer" scene was never a monomer.
    """
    if spec is None:
        return "none"
    s = str(spec).strip().lower()
    if s == "":
        return "none"
    # The old spellings stay as aliases: "monomer" was a bad name -- the
    # scene holds TWO chains, one protomer of each enzyme -- but anything
    # already written against it should keep working.
    aliases = {"mono": "pair", "monomer": "pair", "protomers": "pair",
               "multi": "assembly", "multimer": "assembly", "complex": "pair",
               "occupancy": "assembly", "full": "assembly", "all": "both"}
    s = aliases.get(s, s)
    if s not in VARIANTS:
        raise ValueError(
            "%s: %r is not one of %s" % (what, spec, ", ".join(VARIANTS))
        )
    return s


@dataclass
class FigurePlan:
    """What to build, after dependencies and existing files are accounted for."""

    ranked: List[Pose] = field(default_factory=list)
    pse: List[str] = field(default_factory=list)
    svg: List[str] = field(default_factory=list)
    pdb: List[str] = field(default_factory=list)
    pse_added_for_svg: List[str] = field(default_factory=list)
    pse_reused: List[str] = field(default_factory=list)
    render: str = "none"
    notes: List[str] = field(default_factory=list)

    @property
    def nothing_to_do(self) -> bool:
        return not (self.pse or self.svg or self.pdb)

    def rank_of(self, pose_id: str) -> int:
        for p in self.ranked:
            if p.pose_id == pose_id:
                return p.rank
        raise KeyError(pose_id)

    def summary(self) -> str:
        lines = []
        if self.nothing_to_do:
            return ("nothing selected. Pass --pse / --svg / --pdb with a count "
                    "(5), a rank (rank3, rank2-4) or a pose id (geom007).")
        for kind, items in (("pse", self.pse), ("svg", self.svg), ("pdb", self.pdb)):
            if not items:
                continue
            lines.append("%-4s %d: %s" % (
                kind, len(items),
                ", ".join("%s (rank %d)" % (i, self.rank_of(i)) for i in items)))
        if self.pse_added_for_svg:
            lines.append("     + session(s) added because an SVG needs one: %s"
                         % ", ".join(self.pse_added_for_svg))
        if self.pse_reused:
            lines.append("     + session(s) already on disk, not rebuilt: %s"
                         % ", ".join(self.pse_reused))
        lines.append("scenes   classes, occupancy, render=%s" % self.render)
        lines.extend("note     %s" % n for n in self.notes)
        return "\n".join(lines)


def resolve(
    poses: Iterable,
    pse=None,
    svg=None,
    pdb=None,
    render="none",
    existing_pse: Sequence[str] = (),
    rebuild: bool = False,
) -> FigurePlan:
    """Work out what to build.

    Parameters
    ----------
    poses:
        ``(pose_id, do)`` pairs, or :class:`Pose` objects.
    existing_pse:
        Pose ids whose session is already on disk. An SVG for one of these
        needs no rebuild, so it is not added to the work list.
    rebuild:
        Build sessions even when they already exist.
    """
    ranked = rank_poses(poses)
    plan = FigurePlan(ranked=ranked)
    plan.render = parse_variant(render, "render")

    plan.pse = parse_spec(pse, ranked, "--pse")
    plan.svg = parse_spec(svg, ranked, "--svg")
    plan.pdb = parse_spec(pdb, ranked, "--pdb")

    have = set(existing_pse)
    want_pse = set(plan.pse)
    for pid in plan.svg:
        if pid in want_pse:
            continue
        if pid in have and not rebuild:
            plan.pse_reused.append(pid)
            continue
        want_pse.add(pid)
        plan.pse_added_for_svg.append(pid)
    if rebuild:
        for pid in list(plan.pse):
            if pid in have:
                plan.notes.append("%s: session exists and is being rebuilt (--rebuild)" % pid)
    else:
        already = [p for p in plan.pse if p in have]
        if already:
            plan.notes.append(
                "session(s) already on disk and will be overwritten: %s"
                % ", ".join(already))

    plan.pse = [p.pose_id for p in ranked if p.pose_id in want_pse]
    plan.pse_added_for_svg = [p.pose_id for p in ranked if p.pose_id in set(plan.pse_added_for_svg)]
    plan.pse_reused = [p.pose_id for p in ranked if p.pose_id in set(plan.pse_reused)]

    if plan.svg and plan.render == "none":
        plan.notes.append(
            "an SVG was asked for but render is 'none'; the SVG is traced from "
            "that scene, so set at least one variant")
    return plan
