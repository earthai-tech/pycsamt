# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Explicit, hash-bound Occam prepared-project execution (opt-in)."""
from pathlib import Path
import hashlib
import os
import shutil
import subprocess
import numpy as np
from ..results import InversionResult


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_project(cfg):
    """Copy frozen inputs into a new directory; validate current-run outputs.

    Numerical acceptance is restricted to this explicit fixed-boundary contract.
    It provides no geological or field acceptance.
    """
    from ...models.occam2d import OccamStartup
    from ...models.occam2d.startup import OccamIter
    from ...models.occam2d.data import OccamData
    contract = cfg.backend_options['prepared_project']
    if cfg.data is not None or cfg.method != 'mt' or cfg.dimension != '2d':
        raise ValueError('Prepared project requires data=None and explicit MT/2d.')
    if set(cfg.backend_options) != {'prepared_project'}:
        raise ValueError('Prepared project cannot be combined with implicit builder/runner overrides.')
    source = Path(contract['project_dir']).resolve()
    binary = Path(contract['binary_path']).resolve()
    names = {'Startup', 'Occam2DMesh', 'Occam2DModel', 'OccamDataFile.dat', 'FixedParameters.dat'}
    hashes = contract['input_sha256']
    if set(hashes) != names:
        raise ValueError('Exactly five declared native inputs are required.')
    if digest(binary) != contract['binary_sha256']:
        raise ValueError('Executable SHA-256 mismatch.')
    for name in names:
        if digest(source / name) != hashes[name]:
            raise ValueError('Input SHA-256 mismatch: ' + name)
    startup = OccamStartup.read(source / 'Startup')
    if startup.model_file != 'Occam2DModel' or startup.data_file != 'OccamDataFile.dat':
        raise ValueError('Unexpected startup references.')
    model_text = (source / 'Occam2DModel').read_text()
    fields = {line.split(':', 1)[0].strip().lower(): line.split(':', 1)[1].strip()
              for line in model_text.splitlines() if ':' in line}
    if (fields.get('mesh file') != 'Occam2DMesh'
            or fields.get('statics file', 'none').lower() != 'none'
            or fields.get('prejudice file', 'none').lower() != 'none'):
        raise ValueError('Unexpected auxiliary model references.')
    startup_fields = {line.split(':', 1)[0].strip().lower(): line.split(':', 1)[1].strip()
                      for line in (source/'Startup').read_text().splitlines() if ':' in line}
    limits = [float(x) for x in startup_fields.get('model limits', '').split(',')]
    if limits != [0., 5.]:
        raise ValueError('Tested prepared contract requires log10 model limits [0,5].')
    if startup.max_iterations != cfg.max_iter:
        raise ValueError('Iteration budget conflicts with frozen startup.')
    timeout = float(contract['timeout_seconds'])
    if not np.isfinite(timeout) or not 0 < timeout <= 3600:
        raise ValueError('Explicit finite timeout in (0,3600] required.')
    work = Path(cfg.workdir).resolve()
    # Atomic directory creation prevents reuse, including concurrent launches.
    if work.exists():
        raise ValueError('Prepared execution requires a fresh output directory.')
    if not cfg.run_external:
        return InversionResult(method=cfg.method, dimension=cfg.dimension,
            backend='occam2d', status='ready', metadata={'executed': False,
            'scientific_acceptance': None, 'binary_sha256': digest(binary)})
    work.mkdir(parents=True, exist_ok=False)
    for name in names:
        shutil.copyfile(source / name, work / name)
        if digest(work / name) != hashes[name]:
            raise ValueError('Copied input SHA-256 mismatch: ' + name)
    env = os.environ.copy()
    env['PYCSAMT_OCCAM_FIXED_PARAMETERS'] = 'FixedParameters.dat'
    warnings = []
    code = None
    timed_out = False
    try:
        with (work/'occam_stdout.log').open('w') as out, (work/'occam_stderr.log').open('w') as err:
            code = subprocess.run([str(binary), 'Startup'], cwd=work, env=env,
                                  stdout=out, stderr=err, timeout=timeout).returncode
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        warnings.append(str(exc))
    except OSError as exc:
        warnings.append(str(exc))
    accepted = False
    rms = float('nan')
    iteration = 0
    try:
        if code != 0:
            raise ValueError('Native execution failed: ' + str(code))
        files = sorted(work.glob('ITER*.iter'), key=lambda p: int(p.stem[4:]))
        if not files:
            raise ValueError('No native iterations.')
        last = OccamIter.read(files[-1])
        iteration = int(last.iteration)
        if not 0 < iteration <= cfg.max_iter or not last.misfit_reached:
            raise ValueError('Native target convergence not confirmed.')
        data = np.asarray(OccamData.read(work/'OccamDataFile.dat').data_blocks, dtype=float)
        response = np.loadtxt(work / f'RESP{iteration:02d}.resp', ndmin=2)
        if response.shape != (len(data), 7) or not np.isfinite(response).all():
            raise ValueError('Invalid paired response shape or nonfinite outputs.')
        if not np.array_equal(response[:, :3], data[:, :3]):
            raise ValueError('Response observation identities do not match.')
        residual = (data[:, 3] - response[:, 5]) / data[:, 4]
        # Native response text rounds predictions; allow its documented precision.
        if not np.allclose(residual, response[:, 6], atol=0.005, rtol=0.005):
            raise ValueError('Response residuals do not match declared observations/errors.')
        rms = float(np.sqrt(np.mean(response[:, 6] ** 2)))
        if not np.isfinite(rms) or rms > startup.target_misfit * 1.05:
            raise ValueError('Final paired-response RMS exceeds target tolerance.')
        locks = np.loadtxt(work/'FixedParameters.dat', skiprows=1, ndmin=2)
        header = (work/'FixedParameters.dat').read_text().splitlines()[0].split()
        if [int(v) for v in header] != [len(last.param_values), len(locks)]:
            raise ValueError('Fixed-parameter dimensions do not match.')
        idx = locks[:, 0].astype(int) - 1
        if not np.array_equal(idx + 1, locks[:, 0]) or len(set(idx)) != len(idx) or np.any(idx < 0) or np.any(idx >= len(last.param_values)):
            raise ValueError('Invalid fixed-parameter indices.')
        for path in files:
            model = OccamStartup.read(path) if path.stem == 'ITER00' else OccamIter.read(path)
            params = np.asarray(model.param_values)
            if not np.isfinite(params).all() or np.any(params < 0) or np.any(params > 5):
                raise ValueError('Invalid model parameter values.')
            if not np.allclose(params[idx], locks[:, 1], atol=1e-6, rtol=0):
                raise ValueError('Fixed values changed during inversion.')
        for name in names:
            if digest(work/name) != hashes[name]:
                raise ValueError('Frozen input changed during execution.')
        accepted = True
    except Exception as exc:
        warnings.append(str(exc))
    return InversionResult(method='mt', dimension='2d', backend='occam2d',
        status='converged' if accepted else 'needs_review', rms=rms, n_iter=iteration,
        workdir=str(work), warnings=warnings,
        files={p.name: str(p) for p in work.iterdir() if p.is_file()},
        metadata={'executed': code is not None or timed_out, 'timed_out': timed_out, 'exit_code': code, 'converged': accepted,
                  'numerical_acceptance': accepted, 'scientific_acceptance': None,
                  'binary': str(binary), 'binary_sha256': digest(binary),
                  'input_sha256': hashes, 'fixed_parameters': 'FixedParameters.dat'})
