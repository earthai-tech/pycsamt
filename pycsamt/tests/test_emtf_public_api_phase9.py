# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 9 public API stabilization contracts for transfer functions."""

from __future__ import annotations

import io

import pycsamt
from pycsamt import io as tf_io


_MINIMAL_XML = b"""\
<EM_TF>
  <ProductId>P9.API.2026</ProductId>
  <Tags>impedance</Tags>
  <Data count="1">
    <Period value="1" units="secs">
      <Z type="complex" size="2 2" units="[mV/km]/[nT]">
        <value output="Ex" input="Hx">1 0</value>
        <value output="Ex" input="Hy">2 0</value>
        <value output="Ey" input="Hx">3 0</value>
        <value output="Ey" input="Hy">4 0</value>
      </Z>
    </Period>
  </Data>
</EM_TF>
"""


def test_io_short_aliases_match_long_names():
    assert tf_io.read_tf is tf_io.read_transfer_function
    assert tf_io.write_tf is tf_io.write_transfer_function


def test_top_level_short_alias_is_lazy_and_usable():
    doc = pycsamt.read_tf(_MINIMAL_XML)
    assert doc.product_id == "P9.API.2026"
    assert doc.z.shape == (1, 2, 2)


def test_top_level_write_tf_to_stream():
    doc = pycsamt.read_tf(_MINIMAL_XML)
    stream = io.BytesIO()
    pycsamt.write_tf(doc, stream, format="emtf_xml")
    assert b"<EM_TF>" in stream.getvalue()
