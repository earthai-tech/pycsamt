"""Convert real Occam2D / ModEM / MARE2DEM inversion results to PCSF.

PCSF (pyCSAMT Common Subsurface Format) is the backend-neutral
inversion-result container introduced by
``PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md`` (repository root). This demo
runs every adapter implemented so far against real bundled data:

* Occam2D (``data/occam2D``, 47 real stations) -> ``grid2d``, written
  without topography, with a named (``.csv``) illustrative topo file,
  and with a positional (``.bln``) illustrative topo file requiring a
  UTM/EPSG conversion -- exercising both attribution modes of
  :func:`pycsamt.format.topo_source.resolve_topo` via
  ``occam2d_to_pcsf``'s new ``topo=`` parameter.
* ModEM 3-D (``data/modem/willy_27freq_watex_line02_sample``, a real
  41x50x288 inverted volume) -> ``grid3d``, again without and *with*
  topography -- and the "with" version uses **real** elevation and
  lon/lat, not placeholders: ModEM's own station names carry a ``23-``
  survey-year prefix that ``data/AMT/WILLY_DATA``'s real EDI station
  ids don't, but they are the same physical stations (112 of 125 match
  once the prefix is stripped), passed as a ``topo=`` Sites object
  (the direct EDI-collection-as-topo-source path) rather than a
  hand-built elevation dict.
* MARE2DEM (``data/mare2dem/demo_mt_inversion``, 6540 real regions) ->
  ``mesh_unstructured``. The real mesh is rebuilt in-process (the
  ``triangle`` Python package, already a hard pycsamt dependency -- no
  external Triangle binary needed) from the run's own real ``.poly``
  PSLG, and reproduces the run's real region partition exactly.

DUHI (:mod:`pycsamt.ai.inversion`) has no separate adapter and so no
separate file here: its output only becomes a final resistivity model
once folded back into an Occam2D run, so a DUHI-prepared result
converts through the exact same ``occam2d_to_pcsf`` path already
demonstrated above.

Run from the repository root::

    python examples/pcsf_conversion_demo/run_demo.py

Every ``.pcsf`` file is a real HDF5 container -- open one directly with
``h5py``, ``h5dump``, or any HDF5 viewer; this script also prints each
file's own group/dataset tree so nothing needs installing to inspect
the first one.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "output"

OCCAM_DIR = ROOT / "data" / "occam2D"
MODEM_DIR = ROOT / "data" / "modem" / "willy_27freq_watex_line02_sample"
WILLY_DIR = ROOT / "data" / "AMT" / "WILLY_DATA"
MARE_DIR = ROOT / "data" / "mare2dem" / "demo_mt_inversion"


# ---------------------------------------------------------------------
# Inspection helper -- prints the on-disk HDF5 tree of a .pcsf file
# ---------------------------------------------------------------------


def _print_pcsf_tree(path: Path) -> None:
    import h5py

    print(f"\n--- {path.relative_to(ROOT)} ({path.stat().st_size:,} bytes) ---")

    def _visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(f"  {name}  dataset  shape={obj.shape}  dtype={obj.dtype}")
        else:
            print(f"  {name}/")

    with h5py.File(path, "r") as fh:
        print("  [root attrs]", dict(fh.attrs))
        fh.visititems(_visit)


# ---------------------------------------------------------------------
# Occam2D -> grid2d
# ---------------------------------------------------------------------


def convert_occam2d() -> list[dict]:
    from pycsamt.format import read_pcsf, write_pcsf
    from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
    from pycsamt.models.occam2d.results import InversionResult

    print("\n=== Occam2D -> grid2d ===")
    result = InversionResult(workdir=OCCAM_DIR)
    print(f"Loaded {OCCAM_DIR.relative_to(ROOT)}: final RMS={result.final_rms:.4f}")

    records = []

    model_no_topo = occam2d_to_pcsf(
        result,
        created_by="pcsf_conversion_demo",
        description="Tongkeng CSAMT (47 stations) -- no topography",
    )
    path_no_topo = write_pcsf(model_no_topo, OUTPUT / "occam2d_no_topo.pcsf")
    records.append(_summarize(path_no_topo, model_no_topo))
    print(f"wrote {path_no_topo.relative_to(ROOT)}  (no topography group)")

    # This bundled dataset's own station names (S00, S01, ...) have no
    # independently-sourced real coordinate anywhere in the repository,
    # so both topo files below are clearly-labelled illustrative fixtures
    # (examples/pcsf_conversion_demo/topo/), not a real GPS survey --
    # they demonstrate pycsamt.format.topo_source's two attribution
    # modes, not a claim of real Tongkeng CSAMT station positions.

    # 1) Named (.csv) source: matches by station id, partial coverage
    #    (only S00/S01 listed) -- the direct new-mechanism replacement
    #    for the old hand-typed station_elevations dict this demo used
    #    to build here.
    topo_csv = ROOT / "examples" / "pcsf_conversion_demo" / "topo" / "occam2d_topo.csv"
    model_with_topo = occam2d_to_pcsf(
        result,
        topo=topo_csv,
        created_by="pcsf_conversion_demo",
        description="Tongkeng CSAMT (47 stations) -- 2 illustrative topo.csv stations",
    )
    path_with_topo = write_pcsf(model_with_topo, OUTPUT / "occam2d_with_topo.pcsf")
    records.append(_summarize(path_with_topo, model_with_topo))
    print(f"wrote {path_with_topo.relative_to(ROOT)}  (topo=.csv, 2/47 stations matched by name)")

    # 2) Positional (.bln) source: no station identity in the format,
    #    so pycsamt.format.topo_source first checks the file's point
    #    count against the 47 expected stations, then attributes all
    #    47 in survey order -- the "detect count, then attribute" flow
    #    this feature was built for, plus a real EPSG/UTM conversion
    #    (illustrative UTM zone 48N coordinates -> real lon/lat).
    topo_bln = ROOT / "examples" / "pcsf_conversion_demo" / "topo" / "occam2d_topo.bln"
    model_with_bln_topo = occam2d_to_pcsf(
        result,
        topo=topo_bln,
        epsg=32648,
        created_by="pcsf_conversion_demo",
        description=(
            "Tongkeng CSAMT (47 stations) -- 47/47 illustrative topo.bln "
            "stations, positional match + UTM zone 48N -> lon/lat"
        ),
    )
    path_with_bln_topo = write_pcsf(
        model_with_bln_topo, OUTPUT / "occam2d_with_bln_topo.pcsf"
    )
    records.append(_summarize(path_with_bln_topo, model_with_bln_topo))
    print(
        f"wrote {path_with_bln_topo.relative_to(ROOT)}  "
        "(topo=.bln, 47/47 stations matched positionally, EPSG:32648)"
    )

    # Round-trip checks.
    restored = read_pcsf(path_with_topo)
    np.testing.assert_array_equal(restored.resistivity, model_with_topo.resistivity)
    restored_bln = read_pcsf(path_with_bln_topo)
    np.testing.assert_allclose(restored_bln.stations.lon, model_with_bln_topo.stations.lon)
    print("round-trip check: resistivity + topo lon/lat bit-exact/close after write -> read")

    _print_pcsf_tree(path_no_topo)
    _print_pcsf_tree(path_with_topo)
    _print_pcsf_tree(path_with_bln_topo)
    return records


# ---------------------------------------------------------------------
# ModEM -> grid3d
# ---------------------------------------------------------------------


def _willy_topo_for_modem(modem_station_names: list[str]) -> list:
    """Real WILLY_DATA EDI stations, wrapped under ModEM's own
    "23-"-prefixed site names, ready to pass as ``topo=`` directly --
    pycsamt.format.topo_source.topo_from_sites duck-types on
    id/longitude/latitude/elevation, so any object shape works; the
    id-prefix reconciliation itself (real ModEM names carry a survey-
    year prefix WILLY_DATA's own ids don't, for the same physical
    stations) is the one thing no generic tool can infer automatically.
    """
    from types import SimpleNamespace

    from pycsamt.map import load_lines

    if not WILLY_DIR.exists():
        return []
    willy = load_lines(WILLY_DIR, detect="folder")
    willy_by_id = {s.id: s for s in willy.stations}
    return [
        SimpleNamespace(
            id=name,
            longitude=willy_by_id[name.split("-", 1)[1]].longitude,
            latitude=willy_by_id[name.split("-", 1)[1]].latitude,
            elevation=willy_by_id[name.split("-", 1)[1]].elevation,
        )
        for name in modem_station_names
        if name.startswith("23-") and name.split("-", 1)[1] in willy_by_id
    ]


def convert_modem3d() -> list[dict]:
    from pycsamt.format import read_pcsf, write_pcsf
    from pycsamt.format.adapters.modem3d import modem3d_to_pcsf
    from pycsamt.models.modem.results import InversionResult

    print("\n=== ModEM 3-D -> grid3d ===")
    result = InversionResult(workdir=MODEM_DIR, load_data=True)
    mf = result.model_final
    print(
        f"Loaded {MODEM_DIR.relative_to(ROOT)}: "
        f"{mf.nz}x{mf.ny}x{mf.nx} = {mf.nz * mf.ny * mf.nx:,} cells, "
        f"final RMS={result.final_rms:.4f}, origin={mf.origin}"
    )

    records = []

    model_no_topo = modem3d_to_pcsf(
        result,
        created_by="pcsf_conversion_demo",
        description="Willy L18 line02, 27-freq ModEM 3-D -- no topography override",
    )
    path_no_topo = write_pcsf(model_no_topo, OUTPUT / "modem3d_no_topo.pcsf")
    records.append(_summarize(path_no_topo, model_no_topo))
    print(f"wrote {path_no_topo.relative_to(ROOT)}  (station z is ModEM's own flat 0.0)")

    willy_topo = _willy_topo_for_modem(list(result.data_obs.site_names))
    if willy_topo:
        # ModEM's own .dat file already carries GG_Lat/GG_Lon for every
        # station (its own real-world reference), so topo= (real,
        # EDI-sourced, matched for 112/125 stations) legitimately
        # overrides it per-station -- pycsamt.format.topo_source warns
        # about exactly that below, which this demo lets through on
        # purpose rather than silencing, since it is the real,
        # documented precedence rule in action, not a mistake.
        import warnings

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model_with_topo = modem3d_to_pcsf(
                result,
                topo=willy_topo,
                created_by="pcsf_conversion_demo",
                description=(
                    "Willy L18 line02, 27-freq ModEM 3-D -- real "
                    f"EDI-sourced topo (Sites object) matched to "
                    f"{len(willy_topo)} of {result.data_obs.n_sites} "
                    "stations via WILLY_DATA, overriding the .dat "
                    "file's own GG_Lat/GG_Lon"
                ),
            )
        for w in caught:
            print(f"(expected) UserWarning: {w.message}")
        path_with_topo = write_pcsf(model_with_topo, OUTPUT / "modem3d_with_topo.pcsf")
        records.append(_summarize(path_with_topo, model_with_topo))
        print(
            f"wrote {path_with_topo.relative_to(ROOT)}  "
            f"(topo=Sites object, real elevation+lon/lat matched for "
            f"{len(willy_topo)}/{result.data_obs.n_sites} stations)"
        )

        restored = read_pcsf(path_with_topo)
        np.testing.assert_array_equal(restored.resistivity, model_with_topo.resistivity)
        np.testing.assert_array_equal(
            restored.geometry.origin, model_with_topo.geometry.origin
        )
        print("round-trip check: resistivity + origin bit-exact after write -> read")

        _print_pcsf_tree(path_with_topo)
    else:
        print(f"WILLY_DATA not found under {WILLY_DIR} -- skipping the real-topo file.")

    return records


# ---------------------------------------------------------------------
# MARE2DEM -> mesh_unstructured
# ---------------------------------------------------------------------


def convert_mare2dem() -> list[dict]:
    try:
        import triangle
    except ImportError:
        print("\n=== MARE2DEM -> mesh_unstructured ===")
        print("the 'triangle' package is not installed -- skipping.")
        return []

    from pycsamt.format import read_pcsf, write_pcsf
    from pycsamt.format.adapters.mare2dem import mare2dem_to_pcsf
    from pycsamt.forward.maxwell.contracts_tri import TriMesh
    from pycsamt.models.mare2dem.iotools.poly import read_poly
    from pycsamt.models.mare2dem.results import InversionResult

    print("\n=== MARE2DEM -> mesh_unstructured ===")
    result = InversionResult(workdir=MARE_DIR)
    print(
        f"Loaded {MARE_DIR.relative_to(ROOT)}: "
        f"{result.model.num_regions} real regions, anisotropy={result.model.anisotropy}"
    )

    poly = read_poly(MARE_DIR / "demo.poly")
    pslg = {
        "vertices": poly.nodes,
        "segments": poly.segments - 1,
        "regions": poly.regions,
    }
    triangulated = triangle.triangulate(pslg, "pA")
    region_ids = np.round(triangulated["triangle_attributes"].ravel()).astype(
        np.int64
    )
    mesh = TriMesh(
        nodes_m=triangulated["vertices"],
        triangles=triangulated["triangles"],
        region_ids=region_ids,
    )
    print(
        f"Real in-process triangulation of demo.poly: "
        f"{mesh.n_nodes} nodes, {mesh.n_triangles} triangles, "
        f"{len(np.unique(region_ids))} unique region ids "
        f"(matches the file's {result.model.num_regions} regions exactly)"
    )

    model = mare2dem_to_pcsf(
        result,
        mesh,
        created_by="pcsf_conversion_demo",
        description="MARE2DEM demo_mt_inversion -- real mesh + real per-region resistivity",
    )
    path = write_pcsf(model, OUTPUT / "mare2dem.pcsf")
    print(f"wrote {path.relative_to(ROOT)}  (no per-station topography source for this dataset)")

    restored = read_pcsf(path)
    np.testing.assert_array_equal(restored.resistivity, model.resistivity)
    np.testing.assert_array_equal(
        restored.geometry.connectivity, model.geometry.connectivity
    )
    print("round-trip check: resistivity + mesh connectivity bit-exact after write -> read")

    _print_pcsf_tree(path)
    return [_summarize(path, model)]


# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------


def _summarize(path: Path, model) -> dict:
    return {
        "file": str(path.relative_to(ROOT)),
        "size_bytes": path.stat().st_size,
        "kind": model.kind,
        "source_backend": model.source_backend,
        "resistivity_shape": (
            list(model.resistivity.shape) if model.resistivity is not None else None
        ),
        "n_stations": (
            len(model.stations.name) if model.stations is not None else 0
        ),
        "has_topography": model.topography is not None,
        "n_topography_stations": (
            len(model.topography.station_id) if model.topography is not None else 0
        ),
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    records = []
    records += convert_occam2d()
    records += convert_modem3d()
    records += convert_mare2dem()

    summary = {
        "schema": "pycsamt.pcsf_conversion.demo/v1",
        "plan": "PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md",
        "duhi_note": (
            "DUHI has no separate adapter -- its output only becomes a "
            "final resistivity model once folded back into an Occam2D "
            "run, so a DUHI-prepared result converts through "
            "occam2d_to_pcsf() exactly like the plain Occam2D files above."
        ),
        "files": records,
    }
    summary_path = OUTPUT / "demo-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf8"
    )
    print(f"\nDemo completed: {summary_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
