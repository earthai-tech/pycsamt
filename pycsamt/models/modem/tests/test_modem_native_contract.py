"""Native WS layout and actual forward-output contracts."""
from pathlib import Path
import numpy as np
import pytest
from unittest.mock import patch
from pycsamt.models.modem.model3d import ModEmModel3D
from pycsamt.models.modem.runner import ModEmRunner
from pycsamt.models.modem.config import ModEmConfig

@pytest.mark.parametrize('encoding', ['LOGE','LOG10','LINEAR'])
def test_ws_native_order_encoding_and_georeference(tmp_path, encoding):
    model=ModEmModel3D()
    model.x_widths=np.array([10.,20.,30.]);model.y_widths=np.array([40.,50.]);model.z_widths=np.array([60.])
    model.rho_loge=np.log(np.array([[[10.,20.,30.],[40.,50.,60.]]]))
    model.origin=np.array([-17.,-23.,0.]);model.rotation=12.;model.log_type=encoding
    path=model.write(tmp_path/'model.ws')
    lines=path.read_text().splitlines()
    row=np.array([float(x) for x in lines[6].split()])
    # Native read loop is i=nx,1,-1; physical first row is [10,20,30].
    expected=model.rho_loge[0,0,::-1]
    if encoding=='LOG10':expected=expected/np.log(10.)
    if encoding=='LINEAR':expected=np.exp(expected)
    np.testing.assert_allclose(row,expected,rtol=1e-5)
    read=ModEmModel3D.read(path)
    np.testing.assert_allclose(read.rho_loge,model.rho_loge,rtol=1e-5)
    np.testing.assert_array_equal(read.origin,model.origin)
    assert read.rotation==12.


def test_forward_zero_exit_without_output_fails(tmp_path):
    from types import SimpleNamespace
    runner=ModEmRunner(tmp_path,config=ModEmConfig(use_mpi=False))
    with patch('pycsamt.models.modem.runner._resolve_binary',return_value=Path('solver.exe')), \
         patch('pycsamt.models.modem.runner.subprocess.run',return_value=SimpleNamespace(check_returncode=lambda:None)):
        with pytest.raises(RuntimeError,match='without a nonempty'):
            runner.run_forward('model.ws','data.dat',load_result=False)


def test_existing_forward_response_is_rejected_before_launch(tmp_path):
    (tmp_path/'ModEM_Forward.dat').write_text('stale')
    runner=ModEmRunner(tmp_path,config=ModEmConfig(use_mpi=False))
    with patch('pycsamt.models.modem.runner._resolve_binary',return_value=Path('solver.exe')), \
         patch('pycsamt.models.modem.runner.subprocess.run') as run:
        with pytest.raises(FileExistsError):runner.run_forward('model.ws','data.dat',load_result=False)
    run.assert_not_called()



def test_impedance_builder_cannot_invent_tipper_observations():
    from pycsamt.models.modem.data import ModEmData
    with pytest.raises(NotImplementedError,match='impedance-only'):
        ModEmData.from_edi([],config=ModEmConfig(component_type='Full_Vertical_Components'))



def test_vertical_aliases_write_native_component_names(tmp_path):
    from pycsamt.models.modem.data import ModEmData
    data=ModEmData();data.site_names=['S0']
    data.blocks=[dict(component_type='Full_Vertical_Components',units='[]',
        rows=[(1.,0,0.,0.,0.,'HZX',.1,.2,.01),(1.,0,0.,0.,0.,'HZY',.3,.4,.01)])]
    path=data.write(tmp_path/'tipper.dat')
    read=ModEmData.read(path)
    assert [r[5] for r in read.blocks[0]['rows']]==['TX','TY']
    assert [r[6] for r in read.blocks[0]['rows']]==[.1,.3]



def test_missing_ws_origin_matches_native_centred_grid(tmp_path):
    path=tmp_path/'no_origin.ws'
    path.write_text('# model\n3 1 1 0 LOGE\n10 20 30\n40\n50\n3 2 1\n')
    model=ModEmModel3D.read(path)
    np.testing.assert_array_equal(model.origin,[-30.,-20.,0.])
    np.testing.assert_array_equal(model.rho_loge,[[[1.,2.,3.]]])



@pytest.mark.parametrize('lat,lon', [(1000.,0.), (0.,1000.), (float('nan'),0.)])
def test_invalid_geographic_inputs_fail_before_grid_construction(lat,lon):
    from pycsamt.models.modem.data import _latlon_to_xy
    with pytest.raises(ValueError,match='Invalid geographic'):
        _latlon_to_xy(lat,lon,0.,0.)
