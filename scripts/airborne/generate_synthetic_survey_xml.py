"""Generate synthetic ZTEM, MobileMT, and AFMAG EMTF-XML sample data.

SYNTHETIC DATA NOTICE
----------------------
No vendor delivery sample has been made available for any of these
technologies (see :mod:`pycsamt.airborne.mobilemt`'s module docstring
and ``pycsamt.airborne.io``'s for why that is a permanent, documented
gap rather than an oversight). Every file this script writes is
forward-modeled from a simplified, physically-motivated analytic
approximation -- a skin-depth-weighted, target-shaped anomaly on a
1-D background -- **not** a rigorous 2-D/3-D solution, and **not** a
reproduction of any paper's actual field data. Each generated EMTF
document's ``Description`` and ``ProcessingInfo`` say so explicitly.

Two synthetic survey locations are generated per technology, loosely
inspired by the geographic setting and published target parameters of
one real field example each (coordinates are representative, not the
paper's exact station coordinates, which are not published):

* ZTEM

  - ``gold_springs_nv`` -- low-sulphidation epithermal Au-Ag, Gold
    Springs, NV/UT, USA (Legault et al. 2012). Resistive mushroom-
    shaped silica cap under a topographic high.
  - ``forrestania_wa`` -- shallow massive-sulphide conductor IR2,
    Forrestania test range, WA, Australia (Sattel and Witherly 2012).

* MobileMT

  - ``flammefjeld_greenland`` -- Climax-style porphyry Mo-Cu breccia
    pipe, East Greenland (Zhdanov et al. 2024). Conductive alteration
    halo over a resistive core.
  - ``timiskaming_kimberlite_on`` -- KL-22 kimberlite pipe, Lake
    Timiskaming Kimberlite Field, ON, Canada (Prikhodko et al. 2022).

* AFMAG -- the two response families
  :mod:`pycsamt.airborne.afmag` keeps genuinely separate (a real
  scalar tilt versus a complex 3x2 interstation tensor -- *not* two
  configurations of one response, and *not* the same TF shape ZTEM's
  tipper uses, despite :mod:`pycsamt.emtools.afmag`'s docstring
  drawing a tipper equivalence for its own, narrower purpose):

  - ``abitibi_on`` -- original comparator AFMAG (Ward 1959), a real
    scalar tilt-angle response in degrees, no input/output channels
    at all. Historical Canadian-Shield-style massive-sulphide VMS
    exploration setting, using the two typical historical
    frequencies already published in
    :data:`pycsamt.airborne.afmag.AFMAG_ORIGINAL_TYPICAL_FREQUENCIES_HZ`
    (150, 510 Hz).
  - ``yulong_belt_cn`` -- tensor AFMAG/AirMt (motivated by Liu et al.
    2018's motion-induced-noise treatment of a modern digital
    airborne system -- :mod:`pycsamt.emtools.afmag`'s own
    "Motion-coupling physics" section already implements that
    paper's Eq. 1-14 directly), a complex ``(nf, 3, 2)`` interstation
    magnetic transfer function. Porphyry-Cu-style disseminated
    sulphide setting, Yulong belt, Sichuan/Yunnan border, China,
    within
    :data:`pycsamt.airborne.afmag.AFMAG_TENSOR_PRACTICAL_FREQUENCY_RANGE_HZ`
    (20-800 Hz).

Output layout (one EMTF-XML file per flight-line sample, matching the
one-station-per-file convention already used throughout
``data/xml-samples-generated``)::

    data/ZTEM/<slug>/<SAMPLE_ID>.xml
    data/mobileMT/<slug>/<SAMPLE_ID>.xml
    data/AFMAG/<slug>/<SAMPLE_ID>.xml

Run
---
python scripts/airborne/generate_synthetic_survey_xml.py

These files feed the future ``user_guide/airborne/`` documentation
for :mod:`pycsamt.emtools.ztem`, :mod:`pycsamt.emtools.mobilemt`, and
:mod:`pycsamt.emtools.afmag`.

STATUS FOR THE FUTURE DOCS TASK -- ``Sites`` loading
-----------------------------------------------------
:mod:`pycsamt.emtools.ztem` is ``Sites``-based (it reads
:attr:`~pycsamt.site.base.Site.tipper`). Every numeric accessor on
:class:`~pycsamt.site.base.Site` (``.tipper``, ``.z``, ``.rho``,
``.phase``) routes through :attr:`~pycsamt.site.base.Site.edi`, i.e.
:func:`~pycsamt.emtf.converters.edi.emtf_to_edi`, which deliberately
refuses to build an EDI from a tipper-only document (no impedance
transfer function) -- enforced by real, intentional tests for both
ZTEM and MobileMT. So pure ZTEM/AFMAG EMTF-XML never reached
``Sites`` through that bridge. **Resolved for ZTEM** (2026-08-19, same
session as this script) not by patching that bridge but by a new,
parallel, non-EDI container family,
:class:`~pycsamt.airborne.site.AirborneSites`, plus
:func:`~pycsamt.emtools._core.ensure_any_sites`, which
:mod:`pycsamt.emtools.ztem`'s ``sites=`` parameter now accepts
transparently (including a bare path/directory, auto-detected). The
``data/ZTEM/<slug>/*.xml`` files above are directly usable through
:mod:`pycsamt.emtools.ztem` today.

**Not yet done for AFMAG or MobileMT.** MobileMT was never meant to
use this path at all -- its admittance is a different tensor shape
entirely (see :mod:`pycsamt.emtools.mobilemt`, already built directly
on :class:`~pycsamt.airborne.AirborneEMDataset`). AFMAG is the real
remaining gap: :class:`~pycsamt.airborne.site.AirborneSite` currently
exposes ``.tipper``/``.z``/``.admittance`` but nothing for AFMAG's
own two shapes (a bare scalar ``afmag_tilt`` transfer function, or
the ``interstation_transfer_functions`` 3x2 tensor) -- wiring
:mod:`pycsamt.emtools.afmag` the same way ``ztem.py`` was wired needs
that extended first. Until then, every file below is directly usable
via :meth:`~pycsamt.emtf.EMTF.from_xml` (as
:func:`generate_mobilemt_site`'s reconstruction already does for
MobileMT), just not yet through ``AirborneSites``/``emtools.afmag``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

import pycsamt
from pycsamt.airborne import AirborneEMDataset, NavigationTrack
from pycsamt.airborne.afmag import (
    AFMAGReferenceStation,
    AirMtSystemSpec,
    OriginalAFMAGSystemSpec,
    build_airmt_dataset,
    build_airmt_line,
    build_original_afmag_dataset,
    build_original_afmag_line,
)
from pycsamt.airborne.mobilemt import (
    MobileMTReferenceStation,
    MobileMTSystemSpec,
    build_mobilemt_dataset,
    build_mobilemt_line,
)
from pycsamt.airborne.ztem import (
    ZTEMReferenceStation,
    ZTEMSystemSpec,
    build_ztem_dataset,
    build_ztem_line,
)
from pycsamt.metadata import (
    LocationMeta,
    OrientationMeta,
    ProcessingMeta,
    RemoteReferenceMeta,
    SiteMeta,
    Software,
    SurveyMeta,
)

DATA_ROOT = Path(__file__).resolve().parents[2] / "data"
RNG = np.random.default_rng(20260819)  # fixed seed: reproducible output

_DEG_PER_M_LAT = 1.0 / 111_320.0

_GENERATOR_SOFTWARE = Software(
    name="pycsamt.scripts.airborne.generate_synthetic_survey_xml",
    version=pycsamt.__version__,
)
_SYNTHETIC_NOTICE = (
    "SYNTHETIC DATA. Forward-modeled with a simplified skin-depth "
    "weighted analytic approximation for pyCSAMT documentation and "
    "testing only. Not a vendor delivery sample and not a "
    "reproduction of the cited paper's actual field data."
)


def _deg_per_m_lon(lat_deg: float) -> float:
    return 1.0 / (111_320.0 * max(float(np.cos(np.radians(lat_deg))), 1e-6))


def _offset_latlon(
    lat0: float, lon0: float, east_m: float, north_m: float
) -> tuple[float, float]:
    return (
        lat0 + north_m * _DEG_PER_M_LAT,
        lon0 + east_m * _deg_per_m_lon(lat0),
    )


def _skin_depth_m(rho_ohm_m: float, freq_hz: float) -> float:
    """Vozoff (1972) MT skin depth, already cited in emtools/ztem.py."""
    return 503.0 * np.sqrt(rho_ohm_m / freq_hz)


def _freq_weight(
    freq_hz: float,
    rho_host_ohm_m: float,
    target_depth_m: float,
    sigma_log: float = 1.3,
) -> float:
    """Peak near the frequency whose skin depth matches the target depth."""
    delta = _skin_depth_m(rho_host_ohm_m, freq_hz)
    log_ratio = (np.log10(delta) - np.log10(target_depth_m)) / sigma_log
    return float(np.exp(-0.5 * log_ratio**2))


def _crossover_shape(
    x_m: np.ndarray, x0_m: float, half_width_m: float
) -> np.ndarray:
    """Odd, bell-derivative crossover: 0 at the target, extrema nearby."""
    u = (x_m - x0_m) / half_width_m
    return u * np.exp(-(u**2))


def _gaussian_bump(
    x_m: np.ndarray, x0_m: float, half_width_m: float
) -> np.ndarray:
    u = (x_m - x0_m) / half_width_m
    return np.exp(-(u**2))


def _sample_ids(prefix: str, n: int) -> tuple[str, ...]:
    return tuple(f"{prefix}{i + 1:03d}" for i in range(n))


def _line_positions(
    n_stations: int, spacing_m: float, azimuth_deg: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (x_m along line, east offset, north offset), in metres."""
    x_m = np.arange(n_stations, dtype=float) * spacing_m
    az = np.radians(azimuth_deg)
    east = x_m * np.sin(az)
    north = x_m * np.cos(az)
    return x_m, east, north


# ─────────────────────────────────────────────────────────────────────────
# ZTEM
# ─────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ZtemSiteConfig:
    slug: str
    survey_name: str
    country: str
    year: int
    paper_citation: str
    origin_lat: float
    origin_lon: float
    origin_elev_m: float
    line_azimuth_deg: float
    n_stations: int
    spacing_m: float
    rho_host_ohm_m: float
    target_depth_m: float
    target_half_width_m: float
    polarity: float
    amp_real: float
    amp_imag_ratio: float
    tzy_ratio: float
    noise_std: float
    reference_offset_m: tuple[float, float]
    geology_notes: str
    # Multi-line survey fields (n_lines=1 preserves the original
    # single-line behaviour). strike_half_width_m is the cross-line
    # Gaussian half-width tapering the target's amplitude away from
    # the central line -- simulates a body of finite strike length
    # rather than an infinite 2-D structure; defaults to
    # 2.5 * line_spacing_m when None and n_lines > 1.
    n_lines: int = 1
    line_spacing_m: float = 150.0
    strike_half_width_m: float | None = None


ZTEM_FREQS_HZ = np.array([30.0, 45.0, 90.0, 180.0, 360.0, 720.0])

ZTEM_SITES = (
    ZtemSiteConfig(
        slug="gold_springs_nv",
        survey_name="Gold Springs ZTEM survey (synthetic)",
        country="USA",
        year=2011,
        paper_citation=(
            "Legault, J.M., Zhao, S., and Fitch, R. (2012), ASEG "
            "22nd Conference, Brisbane"
        ),
        origin_lat=37.98,
        origin_lon=-114.15,
        origin_elev_m=1750.0,
        line_azimuth_deg=90.0,  # west-east, matching the paper's flights
        n_stations=15,
        spacing_m=50.0,
        rho_host_ohm_m=150.0,
        target_depth_m=50.0,
        target_half_width_m=125.0,
        polarity=1.0,
        amp_real=0.18,
        amp_imag_ratio=0.4,
        tzy_ratio=0.2,
        noise_std=0.01,
        reference_offset_m=(-500.0, 300.0),
        geology_notes=(
            "Low-sulphidation epithermal Au-Ag; mushroom-shaped "
            "resistive silica cap below a topographic high, ~50 m "
            "thick x 250 m wide at ~50 m depth, in a 60-200 ohm-m "
            "weathered/fresh bedrock host (Legault et al. 2012, "
            "Fig. 6). A 7-line, 150 m-spaced block survey (loosely "
            "inspired by the real 200 m line spacing of Fig. 1) "
            "rather than a single flight line, so the target's "
            "finite along-strike extent and the resulting map-view "
            "'cluster of resistive areas' (Legault et al. 2012, "
            "Fig. 7) can be demonstrated honestly."
        ),
        n_lines=7,
        line_spacing_m=150.0,
    ),
    ZtemSiteConfig(
        slug="forrestania_wa",
        survey_name="Forrestania ZTEM test-range survey (synthetic)",
        country="Australia",
        year=2010,
        paper_citation=(
            "Sattel, D., and Witherly, K. (2012), 2012 NFEM Forum"
        ),
        origin_lat=-32.50,
        origin_lon=119.70,
        origin_elev_m=420.0,
        line_azimuth_deg=0.0,  # north-south, matching line 1075
        n_stations=15,
        spacing_m=40.0,
        rho_host_ohm_m=1000.0,
        target_depth_m=80.0,
        target_half_width_m=40.0,
        polarity=-1.0,
        amp_real=0.18,
        amp_imag_ratio=0.5,
        tzy_ratio=0.25,
        noise_std=0.01,
        reference_offset_m=(400.0, -600.0),
        geology_notes=(
            "Shallow massive-sulphide conductor IR2 (<100 m depth, "
            ">7000 S conductance, ~75x75 m, dipping 30-40 deg north) "
            "in resistive bedrock under conductive overburden "
            "(10-20 S) (Sattel and Witherly 2012)."
        ),
    ),
)


def _ztem_tipper(
    cfg: ZtemSiteConfig,
    x_m: np.ndarray,
    cross_line_gain: float = 1.0,
    *,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    n = x_m.size
    x0 = float(np.median(x_m))
    freqs = ZTEM_FREQS_HZ
    gen = rng if rng is not None else RNG
    tipper = np.zeros((n, freqs.size, 2), dtype=complex)
    for k, freq in enumerate(freqs):
        weight = _freq_weight(freq, cfg.rho_host_ohm_m, cfg.target_depth_m)
        shape = _crossover_shape(x_m, x0, cfg.target_half_width_m)
        real = cfg.polarity * cfg.amp_real * cross_line_gain * weight * shape
        imag = cfg.amp_imag_ratio * real
        tzx = real + 1j * imag
        tzy = cfg.tzy_ratio * tzx
        noise_x = cfg.noise_std * (
            gen.normal(size=n) + 1j * gen.normal(size=n)
        )
        noise_y = cfg.noise_std * (
            gen.normal(size=n) + 1j * gen.normal(size=n)
        )
        tipper[:, k, 0] = tzx + noise_x
        tipper[:, k, 1] = tzy + noise_y
    return tipper


def generate_ztem_site(cfg: ZtemSiteConfig) -> Path:
    n_lines = max(int(cfg.n_lines), 1)
    strike_half_width = cfg.strike_half_width_m or (
        2.5 * cfg.line_spacing_m if n_lines > 1 else 1.0
    )
    # Cross-line unit direction: perpendicular to the flight-line
    # azimuth, so parallel lines are offset across-strike from the
    # central line rather than along it.
    cross_az = np.radians(cfg.line_azimuth_deg + 90.0)
    cross_east_unit = float(np.sin(cross_az))
    cross_north_unit = float(np.cos(cross_az))

    ref_lat, ref_lon = _offset_latlon(
        cfg.origin_lat, cfg.origin_lon, *cfg.reference_offset_m
    )
    reference = ZTEMReferenceStation(
        station_id=f"{cfg.slug.upper()}_BASE",
        site=SiteMeta(
            site_id=f"{cfg.slug.upper()}_BASE",
            name="ZTEM ground reference station (synthetic)",
            location=LocationMeta(
                latitude=ref_lat,
                longitude=ref_lon,
                elevation=cfg.origin_elev_m - 30.0,
            ),
        ),
    )

    lines = []
    line_meta = []  # (line, sample_ids, lat, lon, x_m) per line, for XML pass
    total_stations = 0
    for li in range(n_lines):
        y_off = (li - (n_lines - 1) / 2.0) * cfg.line_spacing_m
        cross_line_gain = (
            _gaussian_bump(y_off, 0.0, strike_half_width)
            if n_lines > 1
            else 1.0
        )

        x_m, east, north = _line_positions(
            cfg.n_stations, cfg.spacing_m, cfg.line_azimuth_deg
        )
        east = east + y_off * cross_east_unit
        north = north + y_off * cross_north_unit
        lat = cfg.origin_lat + north * _DEG_PER_M_LAT
        lon = cfg.origin_lon + east * _deg_per_m_lon(cfg.origin_lat)
        elev = np.full(cfg.n_stations, cfg.origin_elev_m)
        prefix = cfg.slug.upper()[:2] + (
            f"_L{li + 1}_" if n_lines > 1 else "_"
        )
        sample_ids = _sample_ids(prefix, cfg.n_stations)
        line_id = f"{cfg.slug}_L{li + 1}" if n_lines > 1 else cfg.slug

        nav = NavigationTrack(
            sample_ids=sample_ids,
            latitude=lat,
            longitude=lon,
            terrain_elevation=elev,
            platform_elevation=elev + 70.0,
        )
        # Line 1 keeps drawing from the shared module-level RNG, in
        # exactly the same call pattern as the original single-line
        # generator -- so its output, and every downstream site's
        # (forrestania_wa, MobileMT, AFMAG) RNG-draw sequence in
        # main(), stay byte-identical to before this multi-line
        # extension. Extra lines (2..n_lines) draw from their own
        # independently seeded generator instead of consuming further
        # from the shared RNG, so they cannot perturb anything after
        # them in main()'s fixed call order.
        line_rng = (
            RNG if li == 0 else np.random.default_rng(20260819_000 + li)
        )
        tipper = _ztem_tipper(
            cfg, x_m, cross_line_gain=cross_line_gain, rng=line_rng,
        )
        variance = np.full(
            (cfg.n_stations, ZTEM_FREQS_HZ.size, 1, 2), cfg.noise_std**2
        )
        line = build_ztem_line(
            line_id,
            nav,
            tipper,
            frequency=ZTEM_FREQS_HZ,
            variance=variance,
            reference_station=reference,
            system_spec=ZTEMSystemSpec(),
            orientation=OrientationMeta(
                mode="orthogonal",
                angle_to_geographic_north=cfg.line_azimuth_deg,
            ),
            attrs={"paper_reference": cfg.paper_citation},
        )
        lines.append(line)
        line_meta.append((line, sample_ids, lat, lon, x_m))
        total_stations += cfg.n_stations

    dataset = build_ztem_dataset(
        cfg.slug,
        lines,
        survey=SurveyMeta(
            name=cfg.survey_name,
            project="pycsamt synthetic airborne demo",
            operator="pycsamt.examples",
            method="AEM",
            date_start=None,
            crs="WGS84",
            n_stations=total_stations,
            notes=_SYNTHETIC_NOTICE,
        ),
    )
    assert dataset.n_records == total_stations, (
        f"expected {total_stations} records, got {dataset.n_records}"
    )

    out_dir = DATA_ROOT / "ZTEM" / cfg.slug
    out_dir.mkdir(parents=True, exist_ok=True)
    for line, sample_ids, lat, lon, x_m in line_meta:
        for i, sample_id in enumerate(sample_ids):
            record = line.get_record(sample_id)
            doc = record.emtf
            doc.description = (
                f"{_SYNTHETIC_NOTICE} Survey: {cfg.survey_name}. Loosely "
                f"inspired by {cfg.paper_citation}. Flight line "
                f"{line.line_id}."
            )
            doc.site = SiteMeta(
                project="pycsamt synthetic airborne demo",
                survey=cfg.survey_name,
                year_collected=cfg.year,
                country=cfg.country,
                site_id=sample_id,
                name=sample_id,
                location=LocationMeta(
                    latitude=float(lat[i]),
                    longitude=float(lon[i]),
                    elevation=cfg.origin_elev_m,
                ),
                acquired_by="pycsamt.examples (synthetic)",
                extra={
                    "paper_reference": cfg.paper_citation,
                    "target_model": cfg.geology_notes,
                    "line_id": line.line_id,
                    "line_position_m": float(x_m[i]),
                },
            )
            remote = None
            if doc.processing is not None:
                remote = doc.processing.remote_reference
            doc.processing = ProcessingMeta(
                sign_convention="exp(-i ω t)",
                processed_by="pycsamt synthetic data generator",
                software=_GENERATOR_SOFTWARE,
                remote_reference=(
                    remote
                    if remote is not None
                    else RemoteReferenceMeta(
                        reference_type="fixed_ground_horizontal_magnetic",
                        site=reference.preferred_id,
                    )
                ),
                processing_tag=f"{cfg.slug}_synthetic_v1",
                extra={
                    "notice": _SYNTHETIC_NOTICE,
                    "line_id": line.line_id,
                },
            )
            # SiteMeta.extra/ProcessingMeta.extra do not round-trip
            # through EMTF-XML (confirmed empirically); stash the
            # flight-line id where it does survive -- the same
            # metadata["notes"] presentation block build_ztem_line
            # already populates -- so a bare directory re-read via
            # ensure_asites can recover per-station line membership.
            doc.metadata.setdefault("notes", {}).setdefault(
                "ZTEM", {}
            )["LineId"] = line.line_id
            # Re-run the site->station/lat/lon/elev legacy-alias bridge
            # (see EMTF._synchronize_site_bridge) now that .site was
            # set post-construction rather than through the
            # constructor.
            doc.validate()
            doc.write_xml(str(out_dir / f"{sample_id}.xml"))

    return out_dir


# ─────────────────────────────────────────────────────────────────────────
# MobileMT
# ─────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class MobilemtSiteConfig:
    slug: str
    survey_name: str
    country: str
    year: int
    paper_citation: str
    origin_lat: float
    origin_lon: float
    origin_elev_m: float
    line_azimuth_deg: float
    n_stations: int
    spacing_m: float
    rho_host_ohm_m: float
    anomaly_half_width_m: float
    conductive_contrast: float
    noise_frac: float
    reference_offset_m: tuple[float, float]
    geology_notes: str


MOBILEMT_FREQS_HZ = np.logspace(np.log10(25.0), np.log10(12_000.0), 10)

MOBILEMT_SITES = (
    MobilemtSiteConfig(
        slug="flammefjeld_greenland",
        survey_name="Flammefjeld MobileMT survey (synthetic)",
        country="Greenland",
        year=2023,
        paper_citation=(
            "Zhdanov, M.S., Gribenko, A., Prikhodko, A., Sabra, H.E., "
            "Jorgensen, M., and Cox, L.H. (2024), 1st ASEG DISCOVER "
            "Symposium"
        ),
        origin_lat=67.50,
        origin_lon=-32.00,
        origin_elev_m=650.0,
        line_azimuth_deg=90.0,
        n_stations=12,
        spacing_m=100.0,
        rho_host_ohm_m=800.0,
        anomaly_half_width_m=350.0,
        conductive_contrast=0.85,
        noise_frac=0.03,
        reference_offset_m=(-1200.0, 800.0),
        geology_notes=(
            "Climax-style porphyry Mo-Cu breccia pipe (500x800 m) "
            "at 400-600 m depth; conductive phyllic/potassic/"
            "argillic alteration halo over a resistive core "
            "(Zhdanov et al. 2024)."
        ),
    ),
    MobilemtSiteConfig(
        slug="timiskaming_kimberlite_on",
        survey_name=(
            "Lake Timiskaming Kimberlite Field MobileMT survey "
            "(synthetic)"
        ),
        country="Canada",
        year=2020,
        paper_citation="Prikhodko et al. (2022), Minerals, 12(5), 583",
        origin_lat=47.35,
        origin_lon=-79.60,
        origin_elev_m=300.0,
        line_azimuth_deg=0.0,
        n_stations=12,
        spacing_m=60.0,
        rho_host_ohm_m=300.0,
        anomaly_half_width_m=90.0,
        conductive_contrast=0.7,
        noise_frac=0.03,
        reference_offset_m=(700.0, -400.0),
        geology_notes=(
            "KL-22 kimberlite pipe, Lake Timiskaming Kimberlite "
            "Field, NE Ontario; near-surface conductivity increase "
            "delineates the pipe within a moderately resistive host "
            "(Prikhodko et al. 2022, Sec. 4.3)."
        ),
    ),
)


def _mobilemt_resistivity(
    cfg: MobilemtSiteConfig, x_m: np.ndarray
) -> np.ndarray:
    """Position-dependent resistivity: a conductive bump over the target."""
    x0 = float(np.median(x_m))
    bump = _gaussian_bump(x_m, x0, cfg.anomaly_half_width_m)
    return cfg.rho_host_ohm_m / (1.0 + cfg.conductive_contrast * bump)


def _mobilemt_admittance(
    cfg: MobilemtSiteConfig, x_m: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return (admittance (n, nf, 3, 2), apparent_conductivity (n, nf))."""
    n = x_m.size
    freqs = MOBILEMT_FREQS_HZ
    rho = _mobilemt_resistivity(cfg, x_m)  # (n,)
    admittance = np.zeros((n, freqs.size, 3, 2), dtype=complex)
    sigma = np.zeros((n, freqs.size), dtype=float)
    for i in range(n):
        for k, f in enumerate(freqs):
            amp = np.sqrt(5.0 * f * rho[i])
            off_diag = amp * (1 + 1j) / np.sqrt(2)
            # A small diagonal term represents the minor 3-D/near-
            # surface distortion real admittance tensors carry (a
            # perfectly zero diagonal would make the skew diagnostic
            # trivially zero everywhere, which is not representative).
            diag = 0.05 * off_diag * (
                RNG.normal() + 1j * RNG.normal()
            ) / np.sqrt(2)
            z2 = np.array(
                [[diag, off_diag], [-off_diag, diag]],
            )
            y2 = np.linalg.inv(z2)
            noise = cfg.noise_frac * (
                RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2))
            )
            y2 = y2 * (1.0 + noise)
            admittance[i, k, :2, :] = y2
            admittance[i, k, 2, 0] = 0.01 * y2[0, 1]
            admittance[i, k, 2, 1] = 0.01 * y2[1, 0]
            sigma[i, k] = (1.0 / rho[i]) * (
                1.0 + cfg.noise_frac * RNG.normal()
            )
    return admittance, sigma


def generate_mobilemt_site(cfg: MobilemtSiteConfig) -> Path:
    x_m, east, north = _line_positions(
        cfg.n_stations, cfg.spacing_m, cfg.line_azimuth_deg
    )
    lat = cfg.origin_lat + north * _DEG_PER_M_LAT
    lon = cfg.origin_lon + east * _deg_per_m_lon(cfg.origin_lat)
    elev = np.full(cfg.n_stations, cfg.origin_elev_m)
    sample_ids = _sample_ids(cfg.slug.upper()[:2] + "_", cfg.n_stations)

    nav = NavigationTrack(
        sample_ids=sample_ids,
        latitude=lat,
        longitude=lon,
        terrain_elevation=elev,
        platform_elevation=elev + 80.0,
    )
    admittance, sigma = _mobilemt_admittance(cfg, x_m)

    ref_lat, ref_lon = _offset_latlon(
        cfg.origin_lat, cfg.origin_lon, *cfg.reference_offset_m
    )
    reference = MobileMTReferenceStation(
        station_id=f"{cfg.slug.upper()}_BASE",
        site=SiteMeta(
            site_id=f"{cfg.slug.upper()}_BASE",
            name="MobileMT ground electric reference station (synthetic)",
            location=LocationMeta(
                latitude=ref_lat,
                longitude=ref_lon,
                elevation=cfg.origin_elev_m - 20.0,
            ),
        ),
    )

    line = build_mobilemt_line(
        cfg.slug,
        nav,
        admittance,
        frequency=MOBILEMT_FREQS_HZ,
        apparent_conductivity=sigma,
        reference_station=reference,
        system_spec=MobileMTSystemSpec(),
        attrs={"paper_reference": cfg.paper_citation},
    )
    dataset = build_mobilemt_dataset(
        cfg.slug,
        [line],
        survey=SurveyMeta(
            name=cfg.survey_name,
            project="pycsamt synthetic airborne demo",
            operator="pycsamt.examples",
            method="AEM",
            crs="WGS84",
            n_stations=cfg.n_stations,
            notes=_SYNTHETIC_NOTICE,
        ),
    )
    assert dataset.n_records == cfg.n_stations, (
        f"expected {cfg.n_stations} records, got {dataset.n_records}"
    )

    out_dir = DATA_ROOT / "mobileMT" / cfg.slug
    out_dir.mkdir(parents=True, exist_ok=True)
    for i, sample_id in enumerate(sample_ids):
        record = line.get_record(sample_id)
        doc = record.emtf
        # AirborneEMRecord.fields (where build_mobilemt_record put the
        # native apparent-conductivity vector) has no EMTF-XML
        # representation, so it would otherwise be lost on write_xml.
        # Stash it in the one part of the document that *does*
        # round-trip through arbitrary text -- metadata["notes"] --
        # so a future loader can recover it losslessly; see
        # generate_ztem_site's docstring for the ZTEM analogue that
        # does not need this (ZTEM has no comparable extra field).
        doc.metadata.setdefault("notes", {}).setdefault(
            "MobileMT", {}
        )["ApparentConductivitySm"] = ",".join(
            f"{v:.10g}" for v in sigma[i]
        )
        doc.description = (
            f"{_SYNTHETIC_NOTICE} Survey: {cfg.survey_name}. Loosely "
            f"inspired by {cfg.paper_citation}."
        )
        doc.site = SiteMeta(
            project="pycsamt synthetic airborne demo",
            survey=cfg.survey_name,
            year_collected=cfg.year,
            country=cfg.country,
            site_id=sample_id,
            name=sample_id,
            location=LocationMeta(
                latitude=float(lat[i]),
                longitude=float(lon[i]),
                elevation=cfg.origin_elev_m,
            ),
            acquired_by="pycsamt.examples (synthetic)",
            extra={
                "paper_reference": cfg.paper_citation,
                "target_model": cfg.geology_notes,
                "line_position_m": float(x_m[i]),
            },
        )
        remote = None
        if doc.processing is not None:
            remote = doc.processing.remote_reference
        doc.processing = ProcessingMeta(
            sign_convention="exp(-i ω t)",
            processed_by="pycsamt synthetic data generator",
            software=_GENERATOR_SOFTWARE,
            remote_reference=(
                remote
                if remote is not None
                else RemoteReferenceMeta(
                    reference_type="fixed_ground_horizontal_electric",
                    site=reference.preferred_id,
                )
            ),
            processing_tag=f"{cfg.slug}_synthetic_v1",
            extra={
                "notice": _SYNTHETIC_NOTICE,
                "paper_reference": cfg.paper_citation,
                "target_model": cfg.geology_notes,
                "line_position_m": float(x_m[i]),
                "country": cfg.country,
                "year_collected": cfg.year,
            },
        )
        doc.validate()
        doc.write_xml(str(out_dir / f"{sample_id}.xml"))

    return out_dir


# ─────────────────────────────────────────────────────────────────────────
# AFMAG (two response families, kept genuinely separate -- see the
# module docstring for why this does *not* reuse the ZTEM tipper
# generator, even though both are passive-source magnetic methods)
# ─────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class OriginalAfmagSiteConfig:
    slug: str
    survey_name: str
    country: str
    year: int
    paper_citation: str
    origin_lat: float
    origin_lon: float
    origin_elev_m: float
    line_azimuth_deg: float
    n_stations: int
    spacing_m: float
    rho_host_ohm_m: float
    target_depth_m: float
    target_half_width_m: float
    polarity: float
    amp_deg: float
    noise_std_deg: float
    geology_notes: str


@dataclass(frozen=True)
class AirmtSiteConfig:
    slug: str
    survey_name: str
    country: str
    year: int
    paper_citation: str
    origin_lat: float
    origin_lon: float
    origin_elev_m: float
    line_azimuth_deg: float
    n_stations: int
    spacing_m: float
    rho_host_ohm_m: float
    target_depth_m: float
    target_half_width_m: float
    polarity: float
    amp_real: float
    amp_imag_ratio: float
    tzy_ratio: float
    noise_std: float
    horizontal_noise_frac: float
    reference_offset_m: tuple[float, float]
    geology_notes: str


# Ward (1959)'s own typical comparator operating frequencies, already
# published as pycsamt.airborne.afmag.AFMAG_ORIGINAL_TYPICAL_FREQUENCIES_HZ.
AFMAG_ORIGINAL_FREQS_HZ = np.array([150.0, 510.0])

# Six log-spaced frequencies within the published AirMt practical band
# (pycsamt.airborne.afmag.AFMAG_TENSOR_PRACTICAL_FREQUENCY_RANGE_HZ,
# 20-800 Hz), mirroring ZTEM's own 6-frequency convention.
AFMAG_TENSOR_FREQS_HZ = np.logspace(np.log10(25.0), np.log10(600.0), 6)

ORIGINAL_AFMAG_SITES = (
    OriginalAfmagSiteConfig(
        slug="abitibi_on",
        survey_name="Abitibi comparator-AFMAG survey (synthetic)",
        country="Canada",
        year=1960,
        paper_citation="Ward, S.H. (1959), Geophysics, 24(4), 761-787",
        origin_lat=48.25,
        origin_lon=-79.30,
        origin_elev_m=305.0,
        line_azimuth_deg=90.0,
        n_stations=13,
        spacing_m=60.0,
        rho_host_ohm_m=800.0,
        target_depth_m=100.0,
        target_half_width_m=60.0,
        polarity=-1.0,
        amp_deg=12.0,
        noise_std_deg=0.6,
        geology_notes=(
            "Historical Canadian-Shield-style massive-sulphide VMS "
            "target in resistive volcanic host rock, Abitibi "
            "greenstone belt setting; the classic original-AFMAG "
            "comparator exploration ground (Ward 1959)."
        ),
    ),
)

AIRMT_SITES = (
    AirmtSiteConfig(
        slug="yulong_belt_cn",
        survey_name="Yulong belt AirMt survey (synthetic)",
        country="China",
        year=2019,
        paper_citation=(
            "Liu, F., Huang, L., Pang, Y., Shi, Z., Xiao, P., and "
            "Fang, G. (2018), Journal of Applied Geophysics, 158, "
            "129-138"
        ),
        origin_lat=31.30,
        origin_lon=97.20,
        origin_elev_m=4200.0,
        line_azimuth_deg=0.0,
        n_stations=11,
        spacing_m=80.0,
        rho_host_ohm_m=500.0,
        target_depth_m=250.0,
        target_half_width_m=220.0,
        polarity=1.0,
        amp_real=0.14,
        amp_imag_ratio=0.45,
        tzy_ratio=0.2,
        noise_std=0.012,
        horizontal_noise_frac=0.03,
        reference_offset_m=(-900.0, 500.0),
        geology_notes=(
            "Porphyry-Cu-style disseminated sulphide alteration "
            "halo, Yulong porphyry copper belt setting (Sichuan/"
            "Yunnan/Tibet border region, China); a modern digital "
            "tensor-AFMAG/AirMt system, motivated by Liu et al. "
            "(2018)'s motion-induced-noise treatment of exactly this "
            "class of instrument (see "
            "pycsamt.emtools.afmag's own motion-coupling-physics "
            "section, which implements that paper's equations "
            "directly)."
        ),
    ),
)


def _original_afmag_tilt(
    cfg: OriginalAfmagSiteConfig, x_m: np.ndarray
) -> np.ndarray:
    n = x_m.size
    x0 = float(np.median(x_m))
    freqs = AFMAG_ORIGINAL_FREQS_HZ
    tilt = np.zeros((n, freqs.size), dtype=float)
    shape = _crossover_shape(x_m, x0, cfg.target_half_width_m)
    for k, freq in enumerate(freqs):
        weight = _freq_weight(freq, cfg.rho_host_ohm_m, cfg.target_depth_m)
        tilt[:, k] = cfg.polarity * cfg.amp_deg * weight * shape
    tilt = tilt + cfg.noise_std_deg * RNG.normal(size=(n, freqs.size))
    return tilt


def generate_original_afmag_site(cfg: OriginalAfmagSiteConfig) -> Path:
    x_m, east, north = _line_positions(
        cfg.n_stations, cfg.spacing_m, cfg.line_azimuth_deg
    )
    lat = cfg.origin_lat + north * _DEG_PER_M_LAT
    lon = cfg.origin_lon + east * _deg_per_m_lon(cfg.origin_lat)
    elev = np.full(cfg.n_stations, cfg.origin_elev_m)
    sample_ids = _sample_ids(cfg.slug.upper()[:2] + "_", cfg.n_stations)

    nav = NavigationTrack(
        sample_ids=sample_ids,
        latitude=lat,
        longitude=lon,
        terrain_elevation=elev,
        platform_elevation=elev + 60.0,
    )
    tilt = _original_afmag_tilt(cfg, x_m)
    variance = np.full(
        (cfg.n_stations, AFMAG_ORIGINAL_FREQS_HZ.size), cfg.noise_std_deg**2
    )

    line = build_original_afmag_line(
        cfg.slug,
        nav,
        tilt,
        frequency=AFMAG_ORIGINAL_FREQS_HZ,
        response_kind="tilt_angle",
        variance=variance,
        system_spec=OriginalAFMAGSystemSpec(),
        orientation=OrientationMeta(
            mode="orthogonal",
            angle_to_geographic_north=cfg.line_azimuth_deg,
        ),
        attrs={"paper_reference": cfg.paper_citation},
    )
    dataset = build_original_afmag_dataset(
        cfg.slug,
        [line],
        survey=SurveyMeta(
            name=cfg.survey_name,
            project="pycsamt synthetic airborne demo",
            operator="pycsamt.examples",
            method="AEM",
            crs="WGS84",
            n_stations=cfg.n_stations,
            notes=_SYNTHETIC_NOTICE,
        ),
    )
    assert dataset.n_records == cfg.n_stations, (
        f"expected {cfg.n_stations} records, got {dataset.n_records}"
    )

    out_dir = DATA_ROOT / "AFMAG" / cfg.slug
    out_dir.mkdir(parents=True, exist_ok=True)
    for i, sample_id in enumerate(sample_ids):
        record = line.get_record(sample_id)
        doc = record.emtf
        doc.description = (
            f"{_SYNTHETIC_NOTICE} Survey: {cfg.survey_name}. Loosely "
            f"inspired by {cfg.paper_citation}."
        )
        doc.site = SiteMeta(
            project="pycsamt synthetic airborne demo",
            survey=cfg.survey_name,
            year_collected=cfg.year,
            country=cfg.country,
            site_id=sample_id,
            name=sample_id,
            location=LocationMeta(
                latitude=float(lat[i]),
                longitude=float(lon[i]),
                elevation=cfg.origin_elev_m,
            ),
            acquired_by="pycsamt.examples (synthetic)",
            extra={
                "paper_reference": cfg.paper_citation,
                "target_model": cfg.geology_notes,
                "line_position_m": float(x_m[i]),
            },
        )
        # No reference-station geometry for the original comparator
        # generation (see build_original_afmag_emtf's own docstring:
        # it accepts no reference_station at all), so this is a
        # plain ProcessingMeta with no remote_reference block.
        doc.processing = ProcessingMeta(
            processed_by="pycsamt synthetic data generator",
            software=_GENERATOR_SOFTWARE,
            processing_tag=f"{cfg.slug}_synthetic_v1",
            extra={"notice": _SYNTHETIC_NOTICE},
        )
        doc.validate()
        doc.write_xml(str(out_dir / f"{sample_id}.xml"))

    return out_dir


def _airmt_tensor(cfg: AirmtSiteConfig, x_m: np.ndarray) -> np.ndarray:
    n = x_m.size
    x0 = float(np.median(x_m))
    freqs = AFMAG_TENSOR_FREQS_HZ
    tensor = np.zeros((n, freqs.size, 3, 2), dtype=complex)
    shape = _crossover_shape(x_m, x0, cfg.target_half_width_m)
    for i in range(n):
        for k, f in enumerate(freqs):
            weight = _freq_weight(f, cfg.rho_host_ohm_m, cfg.target_depth_m)
            # Horizontal sub-block (Hxx,Hxy,Hyx,Hyy): near-identity
            # gain relating airborne to ground-reference horizontal
            # fields over a magnetically quiet background, plus
            # small instrument/processing noise.
            noise_h = cfg.horizontal_noise_frac * (
                RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2))
            )
            horiz = np.eye(2, dtype=complex) + noise_h
            # Vertical sub-block (Hzx,Hzy): the same skin-depth
            # weighted crossover model used for ZTEM's tipper, since
            # this is physically the same Hz-to-horizontal relation.
            real = cfg.polarity * cfg.amp_real * weight * shape[i]
            imag = cfg.amp_imag_ratio * real
            tzx = real + 1j * imag
            tzy = cfg.tzy_ratio * tzx
            noise_v = cfg.noise_std * (
                RNG.normal(size=2) + 1j * RNG.normal(size=2)
            )
            tensor[i, k, :2, :] = horiz
            tensor[i, k, 2, 0] = tzx + noise_v[0]
            tensor[i, k, 2, 1] = tzy + noise_v[1]
    return tensor


def generate_airmt_site(cfg: AirmtSiteConfig) -> Path:
    x_m, east, north = _line_positions(
        cfg.n_stations, cfg.spacing_m, cfg.line_azimuth_deg
    )
    lat = cfg.origin_lat + north * _DEG_PER_M_LAT
    lon = cfg.origin_lon + east * _deg_per_m_lon(cfg.origin_lat)
    elev = np.full(cfg.n_stations, cfg.origin_elev_m)
    sample_ids = _sample_ids(cfg.slug.upper()[:2] + "_", cfg.n_stations)

    nav = NavigationTrack(
        sample_ids=sample_ids,
        latitude=lat,
        longitude=lon,
        terrain_elevation=elev,
        platform_elevation=elev + 70.0,
    )
    tensor = _airmt_tensor(cfg, x_m)

    ref_lat, ref_lon = _offset_latlon(
        cfg.origin_lat, cfg.origin_lon, *cfg.reference_offset_m
    )
    reference = AFMAGReferenceStation(
        station_id=f"{cfg.slug.upper()}_BASE",
        site=SiteMeta(
            site_id=f"{cfg.slug.upper()}_BASE",
            name="AirMt ground magnetic reference station (synthetic)",
            location=LocationMeta(
                latitude=ref_lat,
                longitude=ref_lon,
                elevation=cfg.origin_elev_m - 30.0,
            ),
        ),
    )

    line = build_airmt_line(
        cfg.slug,
        nav,
        tensor,
        frequency=AFMAG_TENSOR_FREQS_HZ,
        reference_station=reference,
        system_spec=AirMtSystemSpec(),
        orientation=OrientationMeta(
            mode="orthogonal",
            angle_to_geographic_north=cfg.line_azimuth_deg,
        ),
        attrs={"paper_reference": cfg.paper_citation},
    )
    dataset = build_airmt_dataset(
        cfg.slug,
        [line],
        survey=SurveyMeta(
            name=cfg.survey_name,
            project="pycsamt synthetic airborne demo",
            operator="pycsamt.examples",
            method="AEM",
            crs="WGS84",
            n_stations=cfg.n_stations,
            notes=_SYNTHETIC_NOTICE,
        ),
    )
    assert dataset.n_records == cfg.n_stations, (
        f"expected {cfg.n_stations} records, got {dataset.n_records}"
    )

    out_dir = DATA_ROOT / "AFMAG" / cfg.slug
    out_dir.mkdir(parents=True, exist_ok=True)
    for i, sample_id in enumerate(sample_ids):
        record = line.get_record(sample_id)
        doc = record.emtf
        doc.description = (
            f"{_SYNTHETIC_NOTICE} Survey: {cfg.survey_name}. Loosely "
            f"inspired by {cfg.paper_citation}."
        )
        doc.site = SiteMeta(
            project="pycsamt synthetic airborne demo",
            survey=cfg.survey_name,
            year_collected=cfg.year,
            country=cfg.country,
            site_id=sample_id,
            name=sample_id,
            location=LocationMeta(
                latitude=float(lat[i]),
                longitude=float(lon[i]),
                elevation=cfg.origin_elev_m,
            ),
            acquired_by="pycsamt.examples (synthetic)",
            extra={
                "paper_reference": cfg.paper_citation,
                "target_model": cfg.geology_notes,
                "line_position_m": float(x_m[i]),
            },
        )
        remote = None
        if doc.processing is not None:
            remote = doc.processing.remote_reference
        doc.processing = ProcessingMeta(
            sign_convention="exp(-i ω t)",
            processed_by="pycsamt synthetic data generator",
            software=_GENERATOR_SOFTWARE,
            remote_reference=(
                remote
                if remote is not None
                else RemoteReferenceMeta(
                    reference_type="fixed_ground_magnetic",
                    site=reference.preferred_id,
                )
            ),
            processing_tag=f"{cfg.slug}_synthetic_v1",
            extra={"notice": _SYNTHETIC_NOTICE},
        )
        doc.validate()
        doc.write_xml(str(out_dir / f"{sample_id}.xml"))

    return out_dir


# ─────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────


def main() -> None:
    for cfg in ZTEM_SITES:
        out_dir = generate_ztem_site(cfg)
        n = len(list(out_dir.glob("*.xml")))
        print(f"[ZTEM]     {cfg.slug:28s} -> {n:3d} files in {out_dir}")

    for cfg in MOBILEMT_SITES:
        out_dir = generate_mobilemt_site(cfg)
        n = len(list(out_dir.glob("*.xml")))
        print(f"[MobileMT] {cfg.slug:28s} -> {n:3d} files in {out_dir}")

    for cfg in ORIGINAL_AFMAG_SITES:
        out_dir = generate_original_afmag_site(cfg)
        n = len(list(out_dir.glob("*.xml")))
        print(f"[AFMAG]    {cfg.slug:28s} -> {n:3d} files in {out_dir}")

    for cfg in AIRMT_SITES:
        out_dir = generate_airmt_site(cfg)
        n = len(list(out_dir.glob("*.xml")))
        print(f"[AirMt]    {cfg.slug:28s} -> {n:3d} files in {out_dir}")


if __name__ == "__main__":
    main()
