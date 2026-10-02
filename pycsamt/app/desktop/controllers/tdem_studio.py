# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TDEM studio (Qt-free): data sources, views, and conversion to EDI.

Sources
    A TEMAVG survey folder (``.AVG`` / ``.Z`` / ``.LOG``, optional
    coordinate table) or one file from a Geosoft ``.dat``, AMIRA, Zonge,
    WalkTEM or XYZ export -- the readers of :mod:`pycsamt.tdem.io`, with
    the acquisition they need (current, loop, receiver area, units).
Selection
    Soundings are grouped by profile (a TEMAVG file is one line:
    ``TEM100_120`` is point 120 of line TEM100).  Curves are drawn for the
    selected soundings only: drawing all 2 790 of a survey took minutes.
Views
    :data:`VIEWS` -- decay curves and apparent resistivity (selected
    soundings), TEMAVG / Z pseudo-sections and gate profiles (one profile
    or the survey), survey map, elevation profile, overview, dashboard --
    each with the options of its :mod:`pycsamt.tdem.plot` class.
Conversion
    :func:`convert` runs :class:`~pycsamt.tdem.transform.TEMtoEDI` (late-
    time or Fourier; skin-depth or diffusion frequency; homogeneous or
    Weidelt phase; loop-geometry correction) with an optional transmitter
    waveform, and returns :class:`~pycsamt.site.base.Sites` ready for the
    rest of pyCSAMT.
"""

from __future__ import annotations

from pycsamt.api.colormaps import COLORMAPS, colormap_label

from dataclasses import dataclass, field
from typing import Any, Callable

from pycsamt.app.desktop.controllers.inversion_engines import Field

__all__ = [
    "CONVERT_FIELDS",
    "SOURCES",
    "TDEMData",
    "VIEWS",
    "convert",
    "load",
    "render",
    "view",
]


# ── sources ───────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class Source:
    key: str
    label: str
    folder: bool
    filter: str
    fields: tuple


_ACQ = (
    Field("current", "Current", "float", 10.0, 0.001, 1e4, 1.0, 3,
          unit="A", help="Transmitter current"),
    Field("loop_side", "Loop side", "float", 100.0, 0.1, 1e5, 10.0, 1,
          unit="m", help="Square transmitter loop side"),
    Field("rx_area", "Receiver area", "float", 1.0, 0.0001, 1e6, 1.0, 2,
          unit="m²"),
    Field("rx_turns", "Receiver turns", "int", 1, 1, 10000, 1),
    Field("data_unit", "Data unit", "choice", "nV/Am2",
          choices=(("nV/Am2", "nV/(A·m²)"), ("uV/A", "µV/A"),
                   ("V/A", "V/A"), ("SI", "SI"))),
    Field("data_type", "Data", "choice", "dBdt",
          choices=(("dBdt", "dB/dt"), ("voltage", "Voltage"),
                   ("Bz", "B-field"))),
)
_TEMAVG = (
    Field("component", "Component", "choice", "Hz",
          choices=(("Hz", "Hz"), ("Hx", "Hx"), ("Hy", "Hy"))),
    Field("magnitude_unit", "Magnitude unit", "choice", "uV/A",
          choices=(("uV/A", "µV/A"), ("nV/Am2", "nV/(A·m²)"),
                   ("V/A", "V/A"))),
    Field("min_gates", "Min. gates", "int", 1, 1, 200, 1,
          help="Drop soundings with fewer valid gates"),
)

SOURCES: tuple[Source, ...] = (
    Source("temavg", "TEMAVG survey folder (.AVG / .Z / .LOG)", True, "",
           _TEMAVG),
    Source("geosoft", "Geosoft XYZ / .dat export", False,
           "Geosoft (*.dat *.xyz *.gdb);;All files (*)", _ACQ),
    Source("amira", "AMIRA TEM file", False,
           "AMIRA (*.amr *.txt *.dat);;All files (*)", _ACQ),
    Source("zonge", "Zonge TEM file", False,
           "Zonge (*.avg *.dat *.txt);;All files (*)", _ACQ),
    Source("walktem", "WalkTEM export", False,
           "WalkTEM (*.usf *.txt *.csv *.dat);;All files (*)",
           _ACQ + (Field("low_moment", "Low moment", "bool", True),)),
    Source("xyz", "Single sounding (time, data columns)", False,
           "Text (*.txt *.csv *.dat *.xyz);;All files (*)",
           _ACQ + (Field("time_unit", "Time unit", "choice", "s",
                         choices=(("s", "s"), ("ms", "ms"),
                                  ("us", "µs"))),)),
    Source("xyz_multi", "Several soundings in one XYZ (site column)", False,
           "Text (*.txt *.csv *.dat *.xyz);;All files (*)", ()),
)


def source(key: str) -> Source:
    return next(s for s in SOURCES if s.key == key)


@dataclass
class TDEMData:
    source: str = ""
    path: str = ""
    survey: Any = None  # TEMSurvey (TEMAVG folders only)
    soundings: list = field(default_factory=list)
    profiles: dict = field(default_factory=dict)  # name -> [indices]

    @property
    def n(self) -> int:
        return len(self.soundings)


def _profile_of(name: str) -> str:
    name = str(name)
    return name.rsplit("_", 1)[0] if "_" in name else "survey"


def load(source_key: str, path: str, values: dict | None = None, *,
         coordinates: str | None = None) -> TDEMData:
    """Read *path* with the reader of *source_key*."""
    import pycsamt.tdem as tdem
    from pycsamt.tdem import io

    v = dict(values or {})
    data = TDEMData(source=source_key, path=str(path))
    if source_key == "temavg":
        data.survey = tdem.read_temavg_survey(path,
                                              coordinate_file=coordinates)
        data.soundings = data.survey.to_soundings(
            component=v.get("component", "Hz"),
            magnitude_unit=v.get("magnitude_unit", "uV/A"),
            min_gates=int(v.get("min_gates", 1)))
    else:
        acq = {k: v[k] for k in ("current", "rx_area", "rx_turns",
                                 "data_unit", "data_type") if k in v}
        if "loop_side" in v:
            acq["loop_side"] = float(v["loop_side"])
        if source_key == "geosoft":
            out = io.read_geosoft_dat(path, **acq)
        elif source_key == "amira":
            out = io.read_amira(path, **acq)
        elif source_key == "zonge":
            out = io.read_zonge(path, **acq)
        elif source_key == "walktem":
            out = io.read_walkttem(path, low_moment=bool(
                v.get("low_moment", True)), **acq)
        elif source_key == "xyz":
            acq.pop("data_unit", None)
            out = io.read_xyz(path, time_unit=v.get("time_unit", "s"),
                              **acq)
        elif source_key == "xyz_multi":
            out = io.read_xyz_multisite(path)
        else:
            raise ValueError(f"unknown source {source_key!r}")
        data.soundings = list(out) if isinstance(out, (list, tuple)) else (
            list(out.values()) if isinstance(out, dict) else [out])
    if not data.soundings and data.survey is None:
        raise ValueError(f"no TEM soundings found in {path}")
    for i, s in enumerate(data.soundings):
        name = getattr(s, "station_name", "") or f"S{i + 1}"
        data.profiles.setdefault(_profile_of(name), []).append(i)
    return data


# ── views ─────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class TDEMView:
    key: str
    label: str
    group: str  # Soundings | Profiles | Survey
    needs: str  # "soundings" | "survey"
    help: str
    fields: tuple = ()
    per_profile: bool = False  # draws one profile (or the survey)


_CMAPS = (("", "Default"),) + tuple(
    (c, colormap_label(c)) for c in COLORMAPS)
_PHASE = Field("phase_mode", "Phase", "choice", "weidelt",
               choices=(("weidelt", "Weidelt (from ρa slope)"),
                        ("homogeneous", "Homogeneous (45°)")))
_FCONV = Field("freq_convention", "Frequency", "choice", "skin_depth",
               choices=(("skin_depth", "Skin-depth equivalence"),
                        ("diffusion", "Diffusion depth")))
_SECTION_Y = Field("y", "Vertical axis", "choice", "depth",
                   choices=(("depth", "Depth"), ("time_ms", "Time (ms)"),
                            ("window", "Gate number")))

VIEWS: tuple[TDEMView, ...] = (
    TDEMView("decay", "Decay curves", "Soundings", "soundings",
             "Transient decay of the selected soundings.",
             (Field("y_mode", "Show", "choice", "dBdt",
                    choices=(("dBdt", "dB/dt"), ("data", "As recorded"))),
              Field("show_error", "Error bars", "bool", True))),
    TDEMView("rho", "Apparent resistivity & phase", "Soundings",
             "soundings",
             "Late-time apparent resistivity (and phase) of the selected "
             "soundings.",
             (Field("show_phase", "Phase panel", "bool", True), _FCONV,
              _PHASE)),
    TDEMView("avg_section", "Resistivity pseudo-section", "Profiles",
             "survey", "TEMAVG values along the profile, against depth or "
             "time.",
             (Field("value", "Quantity", "choice", "ramp_app_res",
                    choices=(("ramp_app_res", "Apparent resistivity"),
                             ("magnitude", "Response magnitude"),
                             ("percent_magnitude", "Magnitude (%)"))),
              _SECTION_Y, Field("log_value", "Log scale", "bool", True),
              Field("cmap", "Colour map", "choice", "", choices=_CMAPS),
              Field("max_depth", "Max depth", "float", 0.0, 0.0, 1e5, 100.0,
                    0, unit="m", auto_zero=True,
                    help="Cut the section at this depth (0 = full)")),
             per_profile=True),
    TDEMView("z_section", "Z-file pseudo-section", "Profiles", "survey",
             "The .Z response along the profile against time.",
             (Field("log_value", "Log scale", "bool", True),
              Field("absolute", "Absolute value", "bool", True),
              Field("cmap", "Colour map", "choice", "", choices=_CMAPS)),
             per_profile=True),
    TDEMView("gates", "Gate profiles", "Profiles", "survey",
             "Selected gates plotted along the profile.",
             (Field("value", "Quantity", "choice", "magnitude",
                    choices=(("magnitude", "Magnitude"),
                             ("ramp_app_res", "Apparent resistivity"))),
              Field("log_y", "Log scale", "bool", True),
              Field("absolute", "Absolute value", "bool", True)),
             per_profile=True),
    TDEMView("dashboard", "Profile dashboard", "Profiles", "survey",
             "Section, Z section and the selected soundings together.",
             (), per_profile=True),
    TDEMView("map", "Station map", "Survey", "survey",
             "Plan view of every sounding.",
             (Field("color_by", "Colour by", "choice", "elevation",
                    choices=(("elevation", "Elevation"),
                             ("profile", "Profile"))),
              Field("annotate", "Label points", "bool", False),
              Field("contour", "Contours", "bool", False),
              Field("cmap", "Colour map", "choice", "", choices=_CMAPS))),
    TDEMView("elevation", "Elevation profile", "Survey", "survey",
             "Ground elevation along each profile.", ()),
    TDEMView("overview", "Survey overview", "Survey", "survey",
             "Map and elevation profile together.", ()),
)


def view(key: str) -> TDEMView:
    return next(v for v in VIEWS if v.key == key)


def _clean(values: dict) -> dict:
    return {k: v for k, v in values.items()
            if not (k == "cmap" and not v)}


def render(key: str, data: TDEMData, *, selected=(), profile: str = "",
           values: dict | None = None):
    """Draw view *key*; returns a matplotlib Figure."""
    import matplotlib.pyplot as plt

    import pycsamt.tdem as tdem

    v = view(key)
    kw = _clean(dict(values or {}))
    max_depth = float(kw.pop("max_depth", 0) or 0)
    snds = [data.soundings[i] for i in selected] if selected else []
    if v.needs == "soundings":
        if not snds:
            raise ValueError("select one or more soundings (Soundings list)")
    elif data.survey is None:
        raise ValueError("this view needs a TEMAVG survey folder "
                         "(.AVG / .Z files)")
    sv = data.survey
    target = sv.get(profile) if (v.per_profile and profile) else sv
    if key == "decay":
        out = tdem.PlotDecayCurve(snds, **kw).plot()
    elif key == "rho":
        out = tdem.PlotTransformedRho(snds, **kw).plot()
    elif key == "avg_section":
        out = tdem.PlotTEMAVGSection(target, **kw).plot()
    elif key == "z_section":
        z = sv.get_z(profile) if profile else sv
        out = tdem.PlotTEMZSection(z, **kw).plot()
    elif key == "gates":
        out = tdem.PlotGateProfile(target, **kw).plot()
    elif key == "dashboard":
        if not profile:
            raise ValueError("the dashboard shows one profile: pick it")
        chosen = snds or [data.soundings[i]
                          for i in data.profiles.get(profile, [])[:5]]
        out = tdem.PlotTEMDashboard(sv.get(profile), sv.get_z(profile),
                                    chosen).plot()
    else:
        cls = {"map": tdem.PlotSurveyMap,
               "elevation": tdem.PlotElevationProfile,
               "overview": tdem.PlotSurveyOverview}[key]
        try:
            out = cls(sv, **kw).plot()
        except ValueError as exc:
            if "coordinate" not in str(exc).lower():
                raise
            raise ValueError(
                "no station coordinates: put the coordinate table (e.g. "
                "'Coordinate of measuring point.xls') in the survey "
                "folder") from exc
    if hasattr(out, "savefig"):
        fig = out
    else:
        fig = getattr(out, "figure", None)
        if fig is None and isinstance(out, (list, tuple)) and out:
            fig = getattr(out[0], "figure", None)
        fig = fig or plt.gcf()
    if max_depth > 0 and kw.get("y", "depth") == "depth":
        for ax in fig.axes:
            if "depth" in ax.get_ylabel().lower():
                top = min(ax.get_ylim())
                ax.set_ylim(max_depth, top)  # depth axis points down
    return fig


# ── conversion ────────────────────────────────────────────────────────────
CONVERT_FIELDS: tuple[Field, ...] = (
    Field("method", "Method", "choice", "late_time",
          choices=(("late_time", "Late-time asymptote"),
                   ("fourier", "Fourier transform"))),
    _FCONV,
    Field("phase_mode", "Phase", "choice", "homogeneous",
          choices=(("homogeneous", "Homogeneous (45°)"),
                   ("weidelt", "Weidelt (from ρa slope)"))),
    Field("loop_geometry_correction", "Loop geometry correction", "bool",
          True, help="Correct for the finite loop size"),
    Field("waveform", "Transmitter waveform", "choice", "none",
          choices=(("none", "Ideal step (none)"),
                   ("square", "Square (= ideal step)"),
                   ("ramp", "Ramp (turn-off ramp)"),
                   ("halfsine", "Half-sine")),
          help="First-order waveform deconvolution (Fitterman & Stewart "
               "1986); a square wave has an instant turn-off, so it "
               "needs no correction"),
    Field("base_frequency", "Base frequency", "float", 25.0, 0.01, 1e4, 1.0,
          2, unit="Hz", advanced=True),
    Field("duty_cycle", "Duty cycle", "float", 0.5, 0.01, 1.0, 0.05, 2,
          advanced=True),
    Field("ramp_off_us", "Turn-off ramp", "float", 100.0, 0.0, 1e5, 10.0,
          1, unit="µs", advanced=True),
)


def _waveform(values: dict):
    import pycsamt.tdem as tdem

    kind = values.get("waveform", "none")
    bf = float(values.get("base_frequency", 25.0))
    if kind == "square":
        return tdem.SquareWaveform(base_frequency=bf, duty_cycle=float(
            values.get("duty_cycle", 0.5)))
    if kind == "ramp":
        return tdem.RampWaveform(base_frequency=bf, ramp_off=float(
            values.get("ramp_off_us", 100.0)) * 1e-6,
            duty_cycle=float(values.get("duty_cycle", 0.5)))
    if kind == "halfsine":
        return tdem.HalfSineWaveform(base_frequency=bf)
    return None


def convert(soundings, values: dict | None = None):
    """TEM soundings -> :class:`~pycsamt.site.base.Sites` (impedance)."""
    import pycsamt.tdem as tdem
    from pycsamt.site.base import to_sites

    v = dict(values or {})
    snds = list(soundings)
    if not snds:
        raise ValueError("no soundings to convert")
    wf = _waveform(v)
    if wf is not None:
        snds = [s.clone() for s in snds]
        for s in snds:
            s.waveform = wf
    conv = tdem.TEMtoEDI(
        method=v.get("method", "late_time"),
        freq_convention=v.get("freq_convention", "skin_depth"),
        phase_mode=v.get("phase_mode", "homogeneous"),
        loop_geometry_correction=bool(v.get("loop_geometry_correction",
                                            True)))
    return to_sites(conv.transform_many(snds))


def selection_every(indices, step: int) -> list[int]:
    """Every *step*-th of *indices* (first and last kept)."""
    idx = list(indices)
    if step <= 1 or len(idx) <= 2:
        return idx
    out = idx[::step]
    if idx[-1] not in out:
        out.append(idx[-1])
    return out


_: Callable = selection_every  # public helper (window + tests)
