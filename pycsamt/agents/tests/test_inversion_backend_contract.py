"""A requested scientific backend must never be silently replaced."""
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from pycsamt.agents.inversion_backend import InversionBackendAgent


@pytest.mark.parametrize('backend', ['unknown-solver', ''])
def test_unregistered_backend_does_not_run(backend):
    with patch('pycsamt.emtools._core.ensure_sites', return_value=object()), \
         patch('pycsamt.inversion.run_inversion') as run:
        result = InversionBackendAgent().execute({'path': 'unused', 'backend': backend})
    assert result.status == 'failed'
    assert 'not registered' in result.error
    run.assert_not_called()


@pytest.mark.parametrize('backend', ['simpeg', 'builtin'])
def test_unavailable_registered_backend_does_not_run(backend):
    with patch('pycsamt.emtools._core.ensure_sites', return_value=object()), \
         patch('pycsamt.inversion.backends.get_backend', side_effect=ImportError('missing dependency')), \
         patch('pycsamt.inversion.run_inversion') as run:
        result = InversionBackendAgent().execute({'path': 'unused', 'backend': backend})
    assert result.status == 'failed'
    assert backend in result.error and 'unavailable' in result.error
    run.assert_not_called()


def test_available_backend_is_preserved_and_runtime_failure_does_not_retry():
    with patch('pycsamt.emtools._core.ensure_sites', return_value=object()), \
         patch('pycsamt.inversion.backends.get_backend', return_value=object()), \
         patch('pycsamt.inversion.InversionConfig', side_effect=lambda **kw: SimpleNamespace(**kw)), \
         patch('pycsamt.inversion.run_inversion', side_effect=RuntimeError('solver failed')) as run:
        result = InversionBackendAgent().execute({'path': 'unused', 'backend': ' Occam2D '})
    assert result.status == 'failed'
    assert 'solver failed' in result.error
    run.assert_called_once()
    assert run.call_args.args[0].backend == 'occam2d'


@pytest.mark.parametrize('status,rms,metadata,expected', [
    ('needs_review', 0.1, {}, 'needs_review'),
    ('failed', 0.1, {'converged': True}, 'needs_review'),
    ('loaded', 0.1, {}, 'needs_review'),
    ('executed', 0.1, {}, 'needs_review'),
    ('success', 0.1, {}, 'needs_review'),
    ('converged', float('inf'), {}, 'needs_review'),
    ('converged', float('nan'), {}, 'needs_review'),
    ('converged', 0.1, {}, 'success'),
    ('success', 0.1, {'converged': True}, 'success'),
])
def test_completion_requires_convergence_and_finite_rms(status, rms, metadata, expected):
    native = SimpleNamespace(status=status, rms=rms, n_iter=2,
                             metadata=metadata, warnings=['native warning'])
    with patch('pycsamt.emtools._core.ensure_sites', return_value=object()), \
         patch('pycsamt.inversion.backends.get_backend', return_value=object()), \
         patch('pycsamt.inversion.run_inversion', return_value=native):
        result = InversionBackendAgent().execute({'path': 'unused'})
    assert result.status == expected
    assert 'native warning' in result.warnings
    assert result.data['scientific_acceptance'] is None


def test_explicit_execution_and_scientific_settings_are_preserved():
    request = dict(path='unused', backend='occam2d', dimension='2d',
                   run_external=True, phase_error=1.5, include_phase=False,
                   bounds={'log10_rho': [0, 5]}, tol=1e-7,
                   starting_model={'resistivities': [100, 100]},
                   reference_model={'resistivities': [30, 30]},
                   backend_options={'runner': {'timeout': 180, 'binary_path': 'candidate.exe'}})
    native = SimpleNamespace(status='loaded', rms=1.0, n_iter=7)
    with patch('pycsamt.emtools._core.ensure_sites', return_value=object()), \
         patch('pycsamt.inversion.backends.get_backend', return_value=object()), \
         patch('pycsamt.inversion.run_inversion', return_value=native) as run:
        result = InversionBackendAgent().execute(request)
    cfg = run.call_args.args[0]
    for key in ('run_external', 'phase_error', 'include_phase', 'bounds', 'tol',
                'starting_model', 'reference_model', 'backend_options'):
        assert getattr(cfg, key) == request[key]
    assert result.status == 'needs_review'


def test_text_execution_flag_is_rejected():
    with patch('pycsamt.emtools._core.ensure_sites', return_value=object()), \
         patch('pycsamt.inversion.backends.get_backend', return_value=object()), \
         patch('pycsamt.inversion.run_inversion') as run:
        result = InversionBackendAgent().execute({'path': 'unused', 'run_external': 'false'})
    assert result.status == 'failed'
    run.assert_not_called()


@pytest.mark.parametrize('failure', ['nonzero', 'exception'])
def test_failed_occam_execution_cannot_load_stale_model(failure):
    from pathlib import Path
    from uuid import uuid4
    from pycsamt.inversion import InversionConfig
    from pycsamt.inversion.backends.occam2d import Occam2DBackend
    workdir = Path('.tmp') / ('occam_failure_' + uuid4().hex)
    workdir.mkdir(parents=True)
    names = dict(data='data', mesh='mesh', model='model', startup='Startup')
    for name in names.values():
        (workdir / name).write_text('fixture')
    cfg = InversionConfig(backend='occam2d', dimension='2d', workdir=str(workdir),
                          run_external=True, backend_options={
                              'files': names, 'runner': {'timeout': 12}})
    kwargs = {'return_value': 2} if failure == 'nonzero' else {
        'side_effect': RuntimeError('process failed')}
    with patch('pycsamt.models.occam2d.OccamRunner.run', **kwargs) as run, \
         patch('pycsamt.models.occam2d.InversionResult') as load:
        result = Occam2DBackend(cfg).run()
    assert result.status == 'needs_review'
    assert result.rms != result.rms
    load.assert_not_called()
    assert run.call_args.kwargs['timeout'] == 12


@pytest.mark.parametrize('mutation,match', [
    ('binary', 'Executable SHA-256'), ('input', 'Input SHA-256'),
    ('budget', 'Iteration budget'), ('reuse', 'fresh output'),
    ('timeout', 'finite timeout'),
])
def test_prepared_project_rejects_changed_identity_or_unsafe_reuse(mutation, match):
    from pathlib import Path
    from uuid import uuid4
    import hashlib
    from pycsamt.inversion import InversionConfig
    from pycsamt.inversion.backends.occam2d import Occam2DBackend
    root = Path(__file__).resolve().parents[3]
    source = root/'PAPERS/CAGEO/R1/audit/heterogeneous_recovery/c7f5229227a64d93a33fae5642ac9192/background_100'
    binary = root/'PAPERS/CAGEO/R1/audit/fixed_boundary_endpoint/4dd438bfdb8a403cbdddeae1e1895162/build/solver.exe'
    if not binary.is_file() or not source.is_dir():
        pytest.skip('Local native-development archive unavailable')
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    names = ['Startup','Occam2DMesh','Occam2DModel','OccamDataFile.dat','FixedParameters.dat']
    contract = dict(project_dir=str(source), binary_path=str(binary), binary_sha256=digest(binary),
                    input_sha256={n:digest(source/n) for n in names}, timeout_seconds=180)
    out = root/'.tmp'/('prepared_negative_' + uuid4().hex)
    cfg = InversionConfig(backend='occam2d', dimension='2d', max_iter=8,
        workdir=str(out), run_external=True, backend_options={'prepared_project':contract})
    if mutation == 'binary': contract['binary_sha256'] = '0'*64
    if mutation == 'input': contract['input_sha256']['Startup'] = '0'*64
    if mutation == 'budget': cfg.max_iter = 9
    if mutation == 'reuse': out.mkdir(parents=True)
    if mutation == 'timeout': contract['timeout_seconds'] = float('nan')
    with patch('pycsamt.inversion.backends._occam_project.subprocess.run') as run:
        with pytest.raises(ValueError, match=match):
            Occam2DBackend(cfg).run()
    run.assert_not_called()


@pytest.mark.parametrize('damage', ['missing_response', 'wrong_identity', 'wrong_residual', 'changed_fixed', 'nonzero'])
def test_prepared_project_rejects_invalid_current_outputs(damage):
    from pathlib import Path
    from uuid import uuid4
    import hashlib
    import shutil
    import numpy as np
    from pycsamt.inversion import InversionConfig
    from pycsamt.inversion.backends.occam2d import Occam2DBackend
    from pycsamt.models.occam2d.startup import OccamIter
    root = Path(__file__).resolve().parents[3]
    source = root/'PAPERS/CAGEO/R1/audit/heterogeneous_recovery/c7f5229227a64d93a33fae5642ac9192/background_100'
    binary = root/'PAPERS/CAGEO/R1/audit/fixed_boundary_endpoint/4dd438bfdb8a403cbdddeae1e1895162/build/solver.exe'
    if not binary.is_file() or not source.is_dir():
        pytest.skip('Local native-development archive unavailable')
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    names = ['Startup','Occam2DMesh','Occam2DModel','OccamDataFile.dat','FixedParameters.dat']
    contract = dict(project_dir=str(source), binary_path=str(binary), binary_sha256=digest(binary),
                    input_sha256={n:digest(source/n) for n in names}, timeout_seconds=180)
    out = root/'.tmp'/('prepared_outputs_' + uuid4().hex)
    cfg = InversionConfig(backend='occam2d', dimension='2d', max_iter=8,
        workdir=str(out), run_external=True, backend_options={'prepared_project':contract})
    def execute(*args, **kwargs):
        for name in ['ITER00.iter', 'ITER07.iter']:
            shutil.copyfile(source/name, out/name)
        if damage != 'missing_response':
            response = np.loadtxt(source/'RESP07.resp')
            if damage == 'wrong_identity': response[0, 0] += 1
            if damage == 'wrong_residual': response[0, 6] += 0.2
            np.savetxt(out/'RESP07.resp', response)
        if damage == 'changed_fixed':
            model = OccamIter.read(out/'ITER07.iter')
            model.param_values[0] += 0.1
            model.write(out/'ITER07.iter')
        return SimpleNamespace(returncode=2 if damage == 'nonzero' else 0)
    with patch('pycsamt.inversion.backends._occam_project.subprocess.run', side_effect=execute):
        result = Occam2DBackend(cfg).run()
    assert result.status == 'needs_review'
    assert result.metadata['numerical_acceptance'] is False
    assert result.warnings
