# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pyCSAMT Format Studio — standalone format-conversion desktop app.

A small, freezable (PyInstaller-able) PySide6 application dedicated to
converting between pyCSAMT's file formats:

- EDI ⇄ EMTF-XML transfer functions (:mod:`pycsamt.emtf`);
- ModEM 3-D / Occam2D / MARE2DEM inversion results, or any AI/DL array
  bundle, → PCSF/PCSM (:mod:`pycsamt.format.convert_engine`);
- PCSF ⇄ PCSM transcoding, validation, and inspection;
- building PCBH boreholes, PCGL geology legends, PCGS structural
  evidence, and PCPT points-of-interest documents from CSV/XLSX/LAS.

Unlike :mod:`pycsamt.app.desktop` (the full QC/inversion/interpretation
suite), this app has no survey-processing features -- it only converts
files, so it stays small enough to hand to a field/office user as a
single frozen binary. Run it with ``pycsamt-converter``,
``python -m pycsamt.app.converter``, or the packaged executable built
from ``packaging/pyinstaller/pycsamt_converter.spec``.

All conversion logic lives in :mod:`pycsamt.app.converter.jobs` as
plain, Qt-free functions -- widgets only collect input and run those
functions inside :class:`pycsamt.app.converter.workers.ConversionWorker`.
"""
