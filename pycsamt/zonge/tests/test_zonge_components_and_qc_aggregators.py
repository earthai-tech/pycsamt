# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Coverage for the pure re-export aggregator modules
``pycsamt.zonge.components`` and ``pycsamt.zonge.qc``.
"""

from __future__ import annotations


def test_components_reexports_all_expected_symbols():
    import pycsamt.zonge.components as components

    for name in components.__all__:
        assert hasattr(components, name), f"missing re-export: {name}"


def test_components_symbols_match_their_source_modules():
    from pycsamt.zonge import components
    from pycsamt.zonge.meas import Amps, CompMeas, Frequency
    from pycsamt.zonge.qc import PcEmag, PcHmag, PcRho, SEphz, SHphz, SPhz
    from pycsamt.zonge.resphase import Phase, Resistivity
    from pycsamt.zonge.survey import Station
    from pycsamt.zonge.z import Z

    assert components.Frequency is Frequency
    assert components.Amps is Amps
    assert components.CompMeas is CompMeas
    assert components.Station is Station
    assert components.Resistivity is Resistivity
    assert components.Phase is Phase
    assert components.Z is Z
    assert components.PcEmag is PcEmag
    assert components.PcHmag is PcHmag
    assert components.PcRho is PcRho
    assert components.SEphz is SEphz
    assert components.SHphz is SHphz
    assert components.SPhz is SPhz


def test_qc_reexports_all_expected_symbols():
    import pycsamt.zonge.qc as qc

    for name in qc.__all__:
        assert hasattr(qc, name), f"missing re-export: {name}"


def test_qc_symbols_match_their_source_modules():
    from pycsamt.zonge import qc
    from pycsamt.zonge.var_pc import (
        EmagPctErr,
        HmagPctErr,
        PcEmag,
        PcHmag,
        PcRho,
        RhoPctErr,
    )
    from pycsamt.zonge.var_std import PhaseSigma, SEphz, SHphz, SPhz

    assert qc.PcEmag is PcEmag
    assert qc.PcHmag is PcHmag
    assert qc.PcRho is PcRho
    assert qc.EmagPctErr is EmagPctErr
    assert qc.HmagPctErr is HmagPctErr
    assert qc.RhoPctErr is RhoPctErr
    assert qc.SEphz is SEphz
    assert qc.SHphz is SHphz
    assert qc.SPhz is SPhz
    assert qc.PhaseSigma is PhaseSigma
