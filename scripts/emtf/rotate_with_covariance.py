"""Rotate an EMTF document while preserving full covariance when available."""

from __future__ import annotations

import argparse

from pycsamt.emtf import EMTF


parser = argparse.ArgumentParser()
parser.add_argument("source")
parser.add_argument("target")
parser.add_argument("angle", type=float)
args = parser.parse_args()

doc = EMTF.from_xml(args.source)
rotated = doc.rotate(args.angle, variance_policy="drop")
rotated.write(args.target)
print(f"rotated to {args.angle:g} degrees and wrote {args.target}")
