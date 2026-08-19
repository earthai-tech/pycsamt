"""Recover EMTFs and covariance from EDI SPECTRA, then archive as XML."""

from __future__ import annotations

import argparse

from pycsamt.emtf import EMTF


parser = argparse.ArgumentParser()
parser.add_argument("source")
parser.add_argument("target")
parser.add_argument("--angle", type=float, default=None)
args = parser.parse_args()

kwargs = {}
if args.angle is not None:
    kwargs["target_angle"] = args.angle

doc = EMTF.from_edi_spectra(args.source, **kwargs)
doc.write(args.target)
print("full covariance recovered:", doc.metadata["edi_spectra"]["full_covariance_recovered"])
