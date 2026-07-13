# pycsamt/zonge/heads.py
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Header facade that aggregates survey-level property containers
(Hardware, Annotation, Configuration, Tx, Rx, Skip) and exposes
a uniform read/write API based on $keyword=value lines.

This module is intentionally *key/value* oriented, not tabular.
Use AVGComponentBase elsewhere for row/column data blocks.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import (
    Any,
)

from .property import (
    Hardware,
    Receiver,
    SkipFlag,
    SurveyAnnotation,
    SurveyConfiguration,
    Transmitter,
)

__all__ = ["Header", "Head"]

# ------------------------------------------------------------------ #
# hardware banner parsing (commented lines from legacy files)        #
# ------------------------------------------------------------------ #
_AMTAVG_RX = re.compile(
    r"""^\\\s*AMTAVG\s*(?P<ver>[\d.]+)\s*:\s*      # version
        "(?P<src>[^"]+)"\s*,\s*                   # source file
        Dated\s*(?P<dated>[^,]+)\s*,\s*           # dated
        Processed\s*(?P<proc>.+)$                 # processed
    """,
    re.IGNORECASE | re.VERBOSE,
)

_AST_RX = re.compile(
    r"""^\\\s*ASTATIC\s*
        (?P<ver>v?[\d.]+[a-z]?)\s*updated\s*data\s*on\s*
        (?P<updated>.+)$
    """,
    re.IGNORECASE | re.VERBOSE,
)

_TMA_RX = re.compile(
    r"""^\\\s*(?P<pts>\d+)[- ]point\s*TMA\s*Filter\s*at\s*
        (?P<freq>[\d.]+)\s*hertz
    """,
    re.IGNORECASE | re.VERBOSE,
)

_KV_RX = re.compile(r"^\s*\$(?P<key>[^=\s]+)\s*=\s*(?P<val>.*)\s*$")


class HeaderComponentBase:
    """
    Minimal base for *header* components that serialize to
    $keyword=value lines (no CSV tables).

    Subclasses implement:
      - update_from_keywords(meta: dict[str, Any]) -> None
      - to_keywords() -> dict[str, Any]

    This base provides a consistent read/from_avg/write surface
    and utility helpers but keeps the format very lightweight.
    """

    def read(
        self,
        df_unused: Any = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        """Populate from a keyword mapping (ignore df)."""
        if meta:
            self.update_from_keywords(dict(meta))

    def write(self) -> Sequence[str]:
        """Serialize to $keyword=value lines."""
        return _dict_to_kv_lines(self.to_keywords())

    # Subclasses should provide these two:
    def update_from_keywords(self, meta: dict[str, Any]) -> None:  # noqa: D401
        raise NotImplementedError

    def to_keywords(self) -> dict[str, Any]:  # noqa: D401
        raise NotImplementedError


def _kv_lines_to_dict(lines: Sequence[str]) -> dict[str, str]:
    """Parse $key=value lines into a dict (last key wins)."""
    out: dict[str, str] = {}
    for ln in lines:
        m = _KV_RX.match(ln)
        if not m:
            continue
        k = m.group("key").strip()
        v = m.group("val").strip().strip('"')
        out[k] = v
    return out


def _dict_to_kv_lines(d: Mapping[str, Any]) -> list[str]:
    """Format mapping as $key=value lines (stable key order)."""
    lines: list[str] = []
    for k in sorted(d.keys()):
        v = d[k]
        v_str = str(v).replace("\n", " ").strip()
        lines.append(f"${k}={v_str}")
    return lines


def _now_utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _parse_hardware_banner(lines: Sequence[str]) -> dict[str, Any]:
    r"""
    Extract a few useful fields from legacy banner lines starting
    with backslash, e.g.:

      \ AMTAVG 7.76: "LCS01.fld", Dated 99-01-01, Processed 22 Jul 16
      \ ASTATIC v3.60d updated data on 22/07/16
      \ 5-point TMA Filter at 1024 hertz
    """
    hw: dict[str, Any] = {}
    for ln in lines:
        if not ln.strip().startswith("\\"):
            continue
        m1 = _AMTAVG_RX.match(ln)
        if m1:
            hw["version"] = m1.group("ver")
            hw["source_file"] = Path(m1.group("src"))
            hw["dated"] = m1.group("dated").strip()
            hw["processed"] = m1.group("proc").strip()
            continue
        m2 = _AST_RX.match(ln)
        if m2:
            hw["astatic_ver"] = m2.group("ver")
            hw["updated"] = m2.group("updated").strip()
            continue
        m3 = _TMA_RX.match(ln)
        if m3:
            hw["tma_points"] = int(m3.group("pts"))
            hw["tma_freq"] = float(m3.group("freq"))
            continue
    return hw


def _emit_hardware_banner(hw: Hardware) -> list[str]:
    """
    Emit a human-friendly banner (comment lines) at the top of
    the header. Values are best-effort; missing fields omitted.
    """
    out: list[str] = []
    # AMTAVG line
    parts: list[str] = []  # noqa
    ver = hw.get("version")
    src = hw.get("source_file")
    dat = hw.get("dated")
    pro = hw.get("processed")
    if ver or src or dat or pro:
        p = []
        if ver:
            p.append(f"AMTAVG {ver}:")
        if src:
            name = Path(src).name if isinstance(src, str) else src.name
            p.append(f'"{name}",')
        if dat:
            p.append(f"Dated {dat},")
        if pro:
            p.append(f"Processed {pro}")
        if p:
            out.append("\\ " + " ".join(p).rstrip(","))

    # ASTATIC line
    aver = hw.get("astatic_ver")
    upd = hw.get("updated")
    if aver or upd:
        p = ["ASTATIC"]
        if aver:
            p.append(str(aver))
        if upd:
            p.append("updated data on " + str(upd))
        out.append("\\ " + " ".join(p))

    # TMA filter line
    pts = hw.get("tma_points")
    frq = hw.get("tma_freq")
    if pts or frq:
        segs = []
        if pts:
            segs.append(f"{int(pts)}-point")
        segs.append("TMA Filter")
        if frq:
            segs.append(f"at {float(frq):g} hertz")
        out.append("\\ " + " ".join(segs))

    return out


class Header(HeaderComponentBase):
    """
    Aggregate facade for header-level metadata components.
    It reads/writes $keyword=value lines and keeps a compact,
    human-friendly banner for legacy context.

    Composition:
      - hardware   : Hardware
      - annotation : SurveyAnnotation
      - config     : SurveyConfiguration
      - tx         : Transmitter
      - rx         : Receiver
      - skip       : SkipFlag

    Typical usage:
      >>> hdr = Header()
      >>> hdr.read(meta=kv_mapping)
      >>> lines = hdr.write()
      >>> print("\\n".join(lines))
    """

    def __init__(self) -> None:
        super().__init__()
        self.hardware = Hardware()
        self.annotation = SurveyAnnotation()
        self.config = SurveyConfiguration()
        self.tx = Transmitter()
        self.rx = Receiver()
        self.skip = SkipFlag()

        # Keep a copy of raw header lines if provided via
        # from_lines; useful for diagnostics or provenance.
        self._raw_banner: list[str] = []

    def read(
        self,
        df_unused: Any = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        """
        Populate all sub-components from a meta mapping.
        - hardware is populated best-effort if banner fields
          are present in meta (rare) or provided via
          `from_lines(...)`.
        - annotation/config/tx/rx rely on their KEYMAP/ALIASES
          to consume $Job.*, $Survey.*, $Tx.*, $Rx.* keys.
        """
        if not meta:
            return

        # Split mapping into source domains by prefix.
        # We accept keys with or without leading '$'.
        def norm_key(k: str) -> str:
            return k.lstrip("$").strip()

        ann: dict[str, Any] = {}
        cfg: dict[str, Any] = {}
        txm: dict[str, Any] = {}
        rxc: dict[str, Any] = {}

        for k, v in meta.items():
            nk = norm_key(k)
            lk = nk.lower()
            if lk.startswith("job."):
                ann[nk] = v
            elif (
                lk.startswith("survey.")
                or lk.startswith("stn.")
                or lk.startswith("unit.")
                or lk.startswith("gps.")
            ):
                cfg[nk] = v
            elif lk.startswith("tx.") or lk == "xmtr":
                txm[nk] = v
            elif lk.startswith("rx.") or lk == "aspace":
                rxc[nk] = v

        # Delegate to property classes (they resolve aliases).
        if ann:
            self.annotation.update_from_keywords(ann)
        if cfg:
            self.config.update_from_keywords(cfg)
        if txm:
            self.tx.update_from_keywords(txm)
        if rxc:
            self.rx.update_from_keywords(rxc)

        # Skip flag not part of modern headers; keep default.
        # If you pass a code in meta, accept it:
        sf = str(meta.get("skp", self.skip.code))
        try:
            self.skip.set(sf)
        except Exception:
            pass

    @classmethod
    def from_lines(
        cls,
        lines: Sequence[str],
    ) -> Header:
        """
        Parse raw header lines.  This is primarily for legacy
        files where the banner comes as comment lines plus a
        $keyword section.  We parse both.
        """
        # 1) Collect banner lines and keyword lines.
        banner: list[str] = []
        kv: list[str] = []
        for ln in lines:
            s = ln.strip()
            if not s:
                continue
            if s.startswith("\\"):
                banner.append(ln)
                continue
            if s.startswith("$"):
                kv.append(ln)

        # 2) Build meta map from $keyword lines.
        meta = _kv_lines_to_dict(kv)

        # 3) Create header and feed both pieces.
        obj = cls()
        obj._raw_banner = list(banner)

        # Hardware: parse banner explicitly.
        hw_fields = _parse_hardware_banner(banner)
        if hw_fields:
            obj.hardware.set(**hw_fields)

        # The rest: rely on keyword mapping.
        obj.read(meta=meta)
        return obj

    def write(self) -> Sequence[str]:
        """
        Emit a tidy header section:
          - legacy banner (comment lines) synthesized from
            Hardware (if any fields are known),
          - $Job.* lines (annotation),
          - $Survey.* + $Unit.* + $Stn.* lines (config),
          - $Tx.* and $Rx.* lines,
          - a $Written UTC stamp.
        """
        lines: list[str] = []

        # Banner
        banner = _emit_hardware_banner(self.hardware)
        if banner:
            lines.extend(banner)
            lines.append("")  # visual gap before $keywords

        # 2. Get all keywords and format them
        # keywords = self.to_keywords()
        # lines.extend(_dict_to_kv_lines(keywords))

        # $keywords by section
        lines.extend(_dict_to_kv_lines(self.annotation.to_keywords()))
        lines.extend(_dict_to_kv_lines(self.config.to_keywords()))
        lines.extend(_dict_to_kv_lines(self.tx.to_keywords()))
        lines.extend(_dict_to_kv_lines(self.rx.to_keywords()))

        # Stamp at the end for reproducibility.
        lines.append(f"$Written={_now_utc_stamp()}")
        return lines

    def update_from_keywords(self, meta: dict[str, Any]) -> None:
        """
        Populate all sub-components from a keyword mapping.
        This method dispatches keys to the appropriate component
        (annotation, config, tx, rx) based on their prefix.
        """
        # This method's logic is identical to read(), so we
        # can simply call it.
        self.read(meta=meta)

    def to_keywords(self) -> dict[str, Any]:
        """
        Aggregate all $keyword=value pairs from sub-components.

        This collects keywords from annotation, config, tx, and rx
        into a single dictionary suitable for writing to a file.
        """
        # Consolidate keywords from all child components
        # into a single dictionary.
        all_keywords: dict[str, Any] = {}
        all_keywords.update(self.annotation.to_keywords())
        all_keywords.update(self.config.to_keywords())
        all_keywords.update(self.tx.to_keywords())
        all_keywords.update(self.rx.to_keywords())

        return all_keywords

    # ---------------- diagnostics ---------------- #
    def asdict(self) -> dict[str, Any]:
        """Plain dict view of all sub-components."""
        return {
            "hardware": asdict(self.hardware),
            "annotation": asdict(self.annotation),
            "config": asdict(self.config),
            "transmitter": asdict(self.tx),
            "receiver": asdict(self.rx),
            "skip_flag": self.skip.code,
        }

    def __str__(self) -> str:
        line = self.config.line_name or "-"
        stn = (self.config.stn_left, self.config.stn_right)
        return (
            f"Header(project={self.annotation.project_name}, "
            f"line={line}, stn={stn})"
        )

    __repr__ = __str__


# Backward-compatible alias
Head = Header
