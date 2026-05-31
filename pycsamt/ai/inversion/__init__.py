# -*- coding: utf-8 -*-
"""
pycsamt.ai.inversion
====================

High-level AI-based EM inversion workflows.

Planned modules (Phase 2 / Phase 4)
-------------------------------------
inv1d.py  (Phase 2)
    :class:`EMInverter1D` — end-to-end 1-D inversion for MT, CSAMT,
    and TEM data.  Trains a network from a
    :class:`~pycsamt.forward.batch.ForwardDataset`, then maps
    :class:`~pycsamt.z.z.Z` objects to :class:`LayeredModel`
    predictions.

inv2d.py  (Phase 4)
    :class:`EMInverter2D` — 2-D section inversion using a U-Net.

joint.py  (Phase 4)
    :class:`JointInverter` — joint EM + seismic inversion using
    the DRCNN approach (Guo et al. 2021).

ensemble.py  (Phase 5)
    :class:`EnsembleInverter` — uncertainty quantification via
    ensemble predictions from multiple network snapshots.
"""
__all__: list = []
