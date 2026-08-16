"""Read a Zonge AVG survey line + its .stn topography, convert to EDI.

L14 is a real legacy (kind-1) CSAMT line: ``L14.avg`` (AMTAVG 7.40) plus
``L14.stn`` station coordinates in UTM zone 49N (Hunan, China). This script
is the reference pipeline used to build ``L14/edi/`` -- reading the AVG,
attaching topography, converting UTM to lat/lon, and exporting one EDI per
station via :class:`~pycsamt.transformers.jedi.AVGtoEDI`.

Usage (any cwd)::

    python data/avg/L14/avg_to_edi.py
"""

from __future__ import annotations

from pathlib import Path

from pycsamt.transformers.jedi import AVGtoEDI
from pycsamt.zonge import AVG

HERE = Path(__file__).resolve().parent
AVG_PATH = HERE / "L14.avg"
STN_PATH = HERE / "L14.stn"
EDI_DIR = HERE / "edi"

# L14 is a real CSAMT survey in Hunan, China -- longitude ~111 deg E falls
# in UTM zone 49N.
UTM_ZONE = "49N"


def main() -> None:
    avg = AVG.from_file(AVG_PATH, verbose=True)
    avg.add_topography(STN_PATH, utm_zone=UTM_ZONE)
    avg.topo.convert_coords(to="ll", inplace=True)

    print(avg.summary)

    collection = AVGtoEDI().transform(avg)
    result = collection.export(EDI_DIR)

    print(f"EDIs written: {len(result['successful'])}")
    if result["failed"]:
        print(f"EDIs failed: {len(result['failed'])}")
        for station, exc in result["failed"]:
            print(f"  {station}: {exc}")


if __name__ == "__main__":
    main()
