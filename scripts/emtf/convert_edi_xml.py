"""Convert EDI to EMTF XML through the neutral scientific model."""

from __future__ import annotations

import argparse

from pycsamt.emtf import EMTF


parser = argparse.ArgumentParser()
parser.add_argument("source")
parser.add_argument("target")
args = parser.parse_args()

doc = EMTF.from_edi(args.source)
doc.write(args.target)
print(f"wrote {args.target}")
