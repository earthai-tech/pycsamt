"""Check the record boundaries consumed by the native ModEM 2-D reader."""
from pathlib import Path
from uuid import uuid4

import numpy as np

from pycsamt.models.modem.data import ModEmData


def test_every_block_has_native_headers_and_roundtrips():
    root=Path(__file__).resolve().parents[4]
    work=root/'.tmp/modem_multiblock'/uuid4().hex
    work.mkdir(parents=True)
    data=ModEmData();data.site_names=['S1']
    for mode in ['TE','TM']:
        data.blocks.append({'component_type':mode+'_Impedance','units':'Ohm',
                            'sign_convention':'exp(+i\\omega t)',
                            'rows':[[.01,0,0.,500.,0.,mode,1.,2.,.1]]})
    path=data.write(work/'data.dat')
    lines=path.read_text().splitlines()
    starts=[i for i,line in enumerate(lines) if line in ['> TE_Impedance','> TM_Impedance']]
    assert len(starts)==2
    for i in starts:
        assert lines[i-2].startswith('# ') and lines[i-1].startswith('# ')
    back=ModEmData.read(path)
    assert back.component_types==['TE_Impedance','TM_Impedance']
    for original,restored in zip(data.blocks,back.blocks):
        assert original['rows'][0][5]==restored['rows'][0][5]
        np.testing.assert_allclose(original['rows'][0][6:],restored['rows'][0][6:])
