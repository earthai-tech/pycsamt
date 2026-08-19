"""Read an EMTF XML file and write a scientifically equivalent copy."""

from __future__ import annotations

import argparse

from pycsamt.emtf import EMTF


parser = argparse.ArgumentParser()
parser.add_argument("source")
parser.add_argument("target")
args = parser.parse_args()

doc = EMTF.from_xml(args.source)
print("station:", doc.station)
print("frequencies:", doc.n_periods)
print("transfer functions:", sorted(doc.transfer_functions))
doc.write(args.target)
