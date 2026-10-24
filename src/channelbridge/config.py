"""The pair file: everything about one enzyme pair, in one document.

The P5 prototype carried its choices as ~40 command-line flags and a few
hard-coded chain ids. A pair file replaces them, so a run is reproducible
from one artefact that can sit next to the results in version control.

YAML is used when PyYAML is installed; JSON always works, so the package has
no hard YAML dependency. Unknown keys are an error, not a warning: a typo in
``channel_ligand`` that silently disables the channel is the kind of bug this
pipeline cannot afford.

A minimal pair file::

    name: gdh_fnr
    channelled: NADPH
    upstream: 1
    targets:
      - id: 1
        name: GDH
        file: targets/gdh.pdb
        chain: B
        site: "resn NAP and name C4N"
        vent: "resn NAP and name N1A+C8A"
      - id: 2
        name: FNR
        file: targets/fnr.pdb
        chain: N
        site: "resn FDA and name N5"
        vent: "resn A2P"
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .geometry import DO_EXACT, DO_PUBLISHED, DO_ST2A, DOCoefficients
from .render.style import RenderStyle

__all__ = [
    "DO_SETS",
    "CaverSpec",
    "ConfigError",
    "PairFile",
    "PoseSpec",
    "Target",
    "TunnelSpec",
    "load_pair_file",
]

#: DO coefficient sets selectable by name in a pair file. See
#: :mod:`channelbridge.geometry` and docs/imars-verification.md for which
#: iMARS table each reproduces.
DO_SETS = {
    "imars-fig2c": DO_PUBLISHED,
    "imars-exact": DO_EXACT,
    "imars-st2a": DO_ST2A,
}


class ConfigError(ValueError):
    """A pair file is malformed. The message names the offending path."""


def _reject_unknown(where: str, got: Dict[str, Any], allowed: Sequence[str]) -> None:
    extra = sorted(set(got) - set(allowed))
    if extra:
        raise ConfigError(
            "%s: unknown key(s) %s. Allowed: %s"
            % (where, ", ".join(repr(e) for e in extra), ", ".join(sorted(allowed)))
        )


def _need(where: str, got: Dict[str, Any], key: str) -> Any:
    if key not in got or got[key] is None:
        raise ConfigError("%s: missing required key %r" % (where, key))
    return got[key]


def _as_str_list(where: str, value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)):
        out = []
        for v in value:
            if not isinstance(v, str):
                raise ConfigError("%s: expected a string, got %r" % (where, v))
            out.append(v)
        return out
    raise ConfigError("%s: expected a string or list of strings, got %r" % (where, value))


@dataclass
class CaverSpec:
    """How to run CAVER for one target's opening."""

    exclude: List[str] = field(default_factory=lambda: ["HOH", "WAT"])
    start: Optional[Tuple[float, float, float]] = None
    probe_radius: float = 1.4
    shell_depth: float = 6.0
    clustering_threshold: float = 3.5
    approximating_balls: int = 12
    seed: int = 1
    shell_radius: float = 3.5
    shell_radius_sweep: List[float] = field(default_factory=list)
    cluster: Optional[int] = None

    _ALLOWED = (
        "exclude",
        "start",
        "probe_radius",
        "shell_depth",
        "clustering_threshold",
        "approximating_balls",
        "seed",
        "shell_radius",
        "shell_radius_sweep",
        "cluster",
    )

    @classmethod
    def from_dict(cls, where: str, raw: Optional[Dict[str, Any]]) -> "CaverSpec":
        raw = dict(raw or {})
        _reject_unknown(where, raw, cls._ALLOWED)
        start = raw.get("start")
        if start is not None:
            if not isinstance(start, (list, tuple)) or len(start) != 3:
                raise ConfigError("%s.start: expected [x, y, z], got %r" % (where, start))
            start = tuple(float(v) for v in start)
        sweep = raw.get("shell_radius_sweep") or []
        if not isinstance(sweep, (list, tuple)):
            raise ConfigError("%s.shell_radius_sweep: expected a list, got %r" % (where, sweep))
        spec = cls(
            exclude=[s.strip().upper() for s in _as_str_list(where + ".exclude", raw.get("exclude", ["HOH", "WAT"]))],
            start=start,
            probe_radius=float(raw.get("probe_radius", 1.4)),
            shell_depth=float(raw.get("shell_depth", 6.0)),
            clustering_threshold=float(raw.get("clustering_threshold", 3.5)),
            approximating_balls=int(raw.get("approximating_balls", 12)),
            seed=int(raw.get("seed", 1)),
            shell_radius=float(raw.get("shell_radius", 3.5)),
            shell_radius_sweep=[float(v) for v in sweep],
            cluster=None if raw.get("cluster") is None else int(raw["cluster"]),
        )
        if spec.probe_radius <= 0:
            raise ConfigError("%s.probe_radius must be positive" % where)
        if spec.cluster is not None and spec.cluster < 1:
            raise ConfigError("%s.cluster is 1-based; got %r" % (where, spec.cluster))
        return spec


@dataclass
class Target:
    """One enzyme of the pair."""

    id: int
    name: str
    file: Path
    chain: str
    site: str
    vent: str
    chains: Optional[str] = None
    symmetry: Optional[str] = None
    ligands: List[str] = field(default_factory=list)
    channel_ligand: Optional[str] = None
    prosthetic: List[str] = field(default_factory=list)
    extra_openings: List[str] = field(default_factory=list)
    avoid: List[str] = field(default_factory=list)
    avoid_from: List[str] = field(default_factory=list)
    caver: CaverSpec = field(default_factory=CaverSpec)

    _ALLOWED = (
        "id",
        "name",
        "file",
        "chain",
        "site",
        "vent",
        "chains",
        "symmetry",
        "ligands",
        "channel_ligand",
        "prosthetic",
        "extra_openings",
        "avoid",
        "avoid_from",
        "caver",
    )

    @classmethod
    def from_dict(cls, where: str, raw: Dict[str, Any], base: Path) -> "Target":
        if not isinstance(raw, dict):
            raise ConfigError("%s: expected a mapping, got %r" % (where, raw))
        _reject_unknown(where, raw, cls._ALLOWED)
        tid = int(_need(where, raw, "id"))
        if tid not in (1, 2):
            raise ConfigError("%s.id must be 1 or 2 (a pair has two targets); got %r" % (where, tid))
        name = str(_need(where, raw, "name"))
        if not name.isidentifier():
            raise ConfigError(
                "%s.name %r is used in file names and chain labels; use letters, "
                "digits and underscores only" % (where, name)
            )
        path = Path(str(_need(where, raw, "file")))
        if not path.is_absolute():
            path = (base / path).resolve()
        chain = str(_need(where, raw, "chain"))
        if len(chain) != 1:
            raise ConfigError("%s.chain is a single PDB chain id; got %r" % (where, chain))
        return cls(
            id=tid,
            name=name,
            file=path,
            chain=chain,
            site=str(_need(where, raw, "site")),
            vent=str(_need(where, raw, "vent")),
            chains=None if raw.get("chains") is None else str(raw["chains"]),
            symmetry=None if raw.get("symmetry") is None else str(raw["symmetry"]),
            ligands=[s.strip().upper() for s in _as_str_list(where + ".ligands", raw.get("ligands"))],
            channel_ligand=(
                None if raw.get("channel_ligand") is None else str(raw["channel_ligand"]).strip().upper()
            ),
            prosthetic=[s.strip().upper() for s in _as_str_list(where + ".prosthetic", raw.get("prosthetic"))],
            extra_openings=_as_str_list(where + ".extra_openings", raw.get("extra_openings")),
            avoid=_as_str_list(where + ".avoid", raw.get("avoid")),
            avoid_from=_as_str_list(where + ".avoid_from", raw.get("avoid_from")),
            caver=CaverSpec.from_dict(where + ".caver", raw.get("caver")),
        )

    def validate_against_disk(self) -> List[str]:
        """Problems that need the files present; returns human-readable lines."""
        problems = []
        if not self.file.is_file():
            problems.append("%s: file not found: %s" % (self.name, self.file))
            return problems
        if self.file.suffix.lower() in (".cif", ".mmcif"):
            problems.append(
                "%s: %s is mmCIF. CAVER reads PDB only; write .pdb at the prep stage."
                % (self.name, self.file.name)
            )
        return problems


@dataclass
class TunnelSpec:
    """How a CAVER sphere chain is turned into geometry.

    The defaults are the normalisation that makes the geometry independent of
    which CAVER build produced it; see :mod:`channelbridge.tunnel`.
    """

    spacing: float = 1.0
    burial_min: int = 15
    burial_radius: float = 8.0
    persist: float = 2.0
    oic_arm: float = 1.0
    normalise: bool = True

    _ALLOWED = ("spacing", "burial_min", "burial_radius", "persist", "oic_arm", "normalise")

    @classmethod
    def from_dict(cls, where: str, raw: Optional[Dict[str, Any]]) -> "TunnelSpec":
        raw = dict(raw or {})
        _reject_unknown(where, raw, cls._ALLOWED)
        spec = cls(
            spacing=float(raw.get("spacing", 1.0)),
            burial_min=int(raw.get("burial_min", 15)),
            burial_radius=float(raw.get("burial_radius", 8.0)),
            persist=float(raw.get("persist", 2.0)),
            oic_arm=float(raw.get("oic_arm", 1.0)),
            normalise=bool(raw.get("normalise", True)),
        )
        if spec.spacing <= 0:
            raise ConfigError("%s.spacing must be positive" % where)
        if spec.oic_arm <= 0:
            raise ConfigError("%s.oic_arm must be positive" % where)
        if spec.burial_min < 0:
            raise ConfigError("%s.burial_min must be >= 0" % where)
        return spec


@dataclass
class PoseSpec:
    """The rigid-body pose grid."""

    dic_min: float = 6.0
    dic_max: float = 36.0
    dic_step: float = 1.0
    oic: List[float] = field(default_factory=lambda: [180, 165, 150, 135, 120, 105, 90])
    azimuth: List[float] = field(default_factory=lambda: list(range(0, 360, 30)))
    spin: List[float] = field(default_factory=lambda: list(range(0, 360, 30)))
    tilt: List[float] = field(default_factory=lambda: [0, 30, 60, 90])
    binder_gap: float = 6.0
    max_ok: int = 50

    _ALLOWED = (
        "dic_min",
        "dic_max",
        "dic_step",
        "oic",
        "azimuth",
        "spin",
        "tilt",
        "binder_gap",
        "max_ok",
    )

    @classmethod
    def from_dict(cls, where: str, raw: Optional[Dict[str, Any]]) -> "PoseSpec":
        raw = dict(raw or {})
        _reject_unknown(where, raw, cls._ALLOWED)

        def nums(key, default):
            v = raw.get(key)
            if v is None:
                return list(default)
            if not isinstance(v, (list, tuple)):
                raise ConfigError("%s.%s: expected a list, got %r" % (where, key, v))
            return [float(x) for x in v]

        spec = cls(
            dic_min=float(raw.get("dic_min", 6.0)),
            dic_max=float(raw.get("dic_max", 36.0)),
            dic_step=float(raw.get("dic_step", 1.0)),
            oic=nums("oic", [180, 165, 150, 135, 120, 105, 90]),
            azimuth=nums("azimuth", range(0, 360, 30)),
            spin=nums("spin", range(0, 360, 30)),
            tilt=nums("tilt", [0, 30, 60, 90]),
            binder_gap=float(raw.get("binder_gap", 6.0)),
            max_ok=int(raw.get("max_ok", 50)),
        )
        if spec.dic_step <= 0:
            raise ConfigError("%s.dic_step must be positive" % where)
        if spec.dic_min > spec.dic_max:
            raise ConfigError(
                "%s: dic_min %g is above dic_max %g" % (where, spec.dic_min, spec.dic_max)
            )
        for o in spec.oic:
            if not (0.0 <= o <= 180.0):
                raise ConfigError("%s.oic: %g is not in [0, 180] degrees" % (where, o))
        return spec

    @property
    def n_poses(self) -> int:
        """Grid size, for a budget estimate before anything is run."""
        import math

        n_d = math.floor((self.dic_max - self.dic_min) / self.dic_step) + 1
        return n_d * len(self.oic) * len(self.azimuth) * len(self.spin) * len(self.tilt)


def _render_style(raw):
    """Figure appearance from the pair file's ``render:`` section.

    The style layer raises plain ValueError because it is usable on its own;
    a problem that arrived through a pair file has to surface as a
    ConfigError like every other pair-file problem, naming the section.
    """
    try:
        return RenderStyle.from_dict(raw, where="render")
    except ConfigError:
        raise
    except (ValueError, TypeError, KeyError) as exc:
        raise ConfigError("render: %s" % exc) from None


@dataclass
class PairFile:
    """A validated pair file."""

    name: str
    channelled: str
    upstream: int
    targets: List[Target]
    tunnel: TunnelSpec = field(default_factory=TunnelSpec)
    pose: PoseSpec = field(default_factory=PoseSpec)
    do_coefficients: DOCoefficients = DO_PUBLISHED
    render: RenderStyle = field(default_factory=RenderStyle)
    source: Optional[Path] = None

    _ALLOWED = (
        "name",
        "channelled",
        "upstream",
        "targets",
        "tunnel",
        "pose",
        "do_coefficients",
        "render",
    )

    @property
    def upstream_target(self) -> Target:
        return self.target(self.upstream)

    @property
    def downstream_target(self) -> Target:
        return self.target(2 if self.upstream == 1 else 1)

    def target(self, tid: int) -> Target:
        for t in self.targets:
            if t.id == tid:
                return t
        raise KeyError("no target with id %r" % (tid,))

    @classmethod
    def from_dict(cls, raw: Dict[str, Any], *, base: Path, source=None) -> "PairFile":
        if not isinstance(raw, dict):
            raise ConfigError("the pair file must be a mapping at the top level")
        _reject_unknown("<root>", raw, cls._ALLOWED)
        name = str(_need("<root>", raw, "name"))
        raw_targets = _need("<root>", raw, "targets")
        if not isinstance(raw_targets, (list, tuple)) or len(raw_targets) != 2:
            raise ConfigError(
                "<root>.targets: a pair needs exactly 2 targets, got %s"
                % (len(raw_targets) if isinstance(raw_targets, (list, tuple)) else type(raw_targets).__name__)
            )
        targets = [
            Target.from_dict("targets[%d]" % i, t, base) for i, t in enumerate(raw_targets)
        ]
        if {t.id for t in targets} != {1, 2}:
            raise ConfigError(
                "targets must have ids 1 and 2, got %s" % sorted(t.id for t in targets)
            )
        if targets[0].name == targets[1].name:
            raise ConfigError("the two targets have the same name %r" % targets[0].name)
        if targets[0].chain == targets[1].chain:
            raise ConfigError(
                "both targets use chain %r; the pose stage needs them distinguishable"
                % targets[0].chain
            )
        upstream = int(_need("<root>", raw, "upstream"))
        if upstream not in (1, 2):
            raise ConfigError(
                "<root>.upstream must be 1 or 2 (O_IC is one-sided, so the order "
                "matters); got %r" % (upstream,)
            )
        do_name = str(raw.get("do_coefficients", "imars-fig2c"))
        if do_name not in DO_SETS:
            raise ConfigError(
                "<root>.do_coefficients: unknown set %r. Known: %s"
                % (do_name, ", ".join(sorted(DO_SETS)))
            )
        return cls(
            name=name,
            channelled=str(_need("<root>", raw, "channelled")),
            upstream=upstream,
            targets=sorted(targets, key=lambda t: t.id),
            tunnel=TunnelSpec.from_dict("tunnel", raw.get("tunnel")),
            pose=PoseSpec.from_dict("pose", raw.get("pose")),
            render=_render_style(raw.get("render")),
            do_coefficients=DO_SETS[do_name],
            source=None if source is None else Path(source),
        )

    def validate_against_disk(self) -> List[str]:
        problems = []
        for t in self.targets:
            problems.extend(t.validate_against_disk())
        return problems

    def summary(self) -> str:
        up, dn = self.upstream_target, self.downstream_target
        lines = [
            "pair            %s" % self.name,
            "channelled      %s" % self.channelled,
            "upstream        %s (target %d)" % (up.name, up.id),
            "downstream      %s (target %d)" % (dn.name, dn.id),
            "DO coefficients %s  (1 A of D_IC = %.1f deg of O_IC)"
            % (self.do_coefficients.name, self.do_coefficients.degrees_per_angstrom),
            "tunnel          %s"
            % (
                "normalised: resample %.2f A, exit at burial >= %d, O_IC arm %.2f A"
                % (self.tunnel.spacing, self.tunnel.burial_min, self.tunnel.oic_arm)
                if self.tunnel.normalise
                else "RAW sphere chain (build-dependent exit and O_IC)"
            ),
            "pose grid       %s poses" % format(self.pose.n_poses, ","),
            "figures         background %s, %s"
            % (self.render.background,
               "an opaque role is set (overlay needs occlusion clipping)"
               if self.render.needs_occlusion_clipping
               else "nothing opaque, tunnels are the subject"),
        ]
        for t in self.targets:
            lines.append(
                "target %d        %s chain %s  %s"
                % (t.id, t.name, t.chain, t.file.name)
            )
        return "\n".join(lines)


def load_pair_file(path) -> PairFile:
    """Read and validate a pair file (``.yaml``/``.yml`` or ``.json``).

    Relative ``file:`` paths resolve against the pair file's own directory,
    so a pair file plus its targets can be moved as a unit.
    """
    path = Path(os.fspath(path))
    if not path.is_file():
        raise ConfigError("pair file not found: %s" % path)
    text = path.read_text(encoding="utf-8")
    suffix = path.suffix.lower()
    if suffix in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise ConfigError(
                "%s is YAML but PyYAML is not installed. Install it "
                "(pip install 'channelbridge[yaml]') or write the pair file as JSON."
                % path
            ) from exc
        try:
            raw = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            raise ConfigError("%s is not valid YAML: %s" % (path, exc)) from exc
    elif suffix == ".json":
        try:
            raw = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ConfigError("%s is not valid JSON: %s" % (path, exc)) from exc
    else:
        raise ConfigError(
            "%s: expected a .yaml, .yml or .json pair file, got %r" % (path, suffix)
        )
    if raw is None:
        raise ConfigError("%s is empty" % path)
    return PairFile.from_dict(raw, base=path.parent, source=path)
