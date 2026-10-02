# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Opt-in, disposable Docker fixture checks. Never run generated Python locally.

The caller supplies a trusted, locally installed Linux image by immutable ID
and a small JSON output contract. Inputs travel over the client's stdin, so
no host directory, project file or credential is mounted; the daemon may
be local, in WSL, or remote. ``PYCSAMT_DOCKER`` sets the client command
(for example ``wsl.exe -e docker``); by default ``docker`` on PATH is used. Container output is evidence from a generated program, not a proof
of scientific correctness or resistance to a malicious program's deception.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
import threading
import time
import uuid
from pathlib import Path, PurePosixPath

# This supervisor runs INSIDE Docker. A subprocess alone is not isolation.
_RUNNER = r'''
import json, stat, subprocess, sys
from pathlib import Path
bundle = json.loads(sys.stdin.read())
namespace = {}
exec(compile(bundle['values'], 'values.py', 'exec'), namespace)
validate_values = namespace['validate_values']
contract = bundle['contract']
# The script is held in memory and passed to a fresh interpreter's stdin.
run = subprocess.run(['python', '-I', '-'], input=bundle['script'].encode('utf-8'),
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
checks = {}
for name, rule in contract.items():
    try:
        path = Path('/work') / name
        # Reject links and special files before reading generated output.
        for part in [path, *path.parents]:
            if part == Path('/work'):
                break
            if part.is_symlink():
                raise ValueError('symlink output')
        info = path.stat()
        if not stat.S_ISREG(info.st_mode) or not 0 < info.st_size <= 1048576:
            raise ValueError('output must be a nonempty regular file <= 1 MiB')
        if rule['kind'] != 'file':
            value = json.loads(path.read_text())
            validate_values(value, rule)
        checks[name] = {'state': 'passed', 'reason': 'Fixture output meets its declared contract'}
    except Exception as exc:
        checks[name] = {'state': 'failed', 'reason': str(exc)[:300]}
print(json.dumps({'returncode': run.returncode, 'outputs': checks}))
'''


_BOOTSTRAP = 'import sys;exec(sys.stdin.readline())'


def _contract(outputs):
    if not isinstance(outputs, dict) or len(outputs) > 16:
        raise ValueError('Expected at most 16 output contracts')
    for name, rule in outputs.items():
        path = PurePosixPath(name)
        if (not name or len(name) > 200 or path.is_absolute()
                or '..' in path.parts or '\\' in name or ':' in name
                or str(path) != name or name == '.'):
            raise ValueError('Output names must be relative paths inside /work')
        if not isinstance(rule, dict) or rule.get('kind') not in {
            'file', 'positive_factors', 'array'
        }:
            raise ValueError('Unknown output contract')
        if rule['kind'] == 'array':
            shape = rule.get('shape')
            if (not isinstance(shape, list) or not 1 <= len(shape) <= 4
                    or any(type(n) is not int or not 1 <= n <= 1000 for n in shape)):
                raise ValueError('Array shape requires 1-4 bounded positive dimensions')
    return outputs


def docker_command():
    """Docker client command prefix, or ``None`` when no client is available."""
    configured = os.environ.get('PYCSAMT_DOCKER', '').strip()
    if configured:
        parts = shlex.split(configured, posix=os.name != 'nt')
        return parts if parts and shutil.which(parts[0]) else None
    docker = shutil.which('docker')
    return [docker] if docker else None


def _bounded_run(command, timeout, cancelled, stdin=None):
    """Bound Docker client output in memory and respond to cancellation."""
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               stdin=subprocess.PIPE if stdin is not None else None)
    chunks, size, overflow = [], [0], threading.Event()

    def read():
        while True:
            chunk = process.stdout.read(4096)
            if not chunk:
                break
            size[0] += len(chunk)
            if size[0] > 65536:
                overflow.set()
                break
            chunks.append(chunk)

    reader = threading.Thread(target=read, daemon=True)
    reader.start()
    if stdin is not None:
        try:
            process.stdin.write(stdin)
        finally:
            process.stdin.close()
    deadline = time.monotonic() + timeout
    try:
        while process.poll() is None:
            if cancelled and cancelled():
                raise RuntimeError('Fixture execution cancelled')
            if overflow.is_set():
                raise RuntimeError('Fixture report exceeded 64 KiB')
            if time.monotonic() >= deadline:
                raise RuntimeError('Fixture execution timed out')
            time.sleep(0.05)
        reader.join(timeout=1)
        if reader.is_alive() or overflow.is_set():
            raise RuntimeError('Fixture report incomplete or too large')
        return process.returncode, b''.join(chunks).decode('utf-8', errors='replace')
    finally:
        if process.poll() is None:
            process.kill()
        process.wait(timeout=5)
        reader.join(timeout=1)
        process.stdout.close()


def execute_fixture(code, *, image, outputs=None, timeout=20, cancelled=None):
    """Check an isolated script against caller-owned tiny output contracts.

    ``image`` must be a local Linux image's sha256 ID, containing Python and
    required dependencies. Pulling/building images is deliberately separate.
    Scripts use /work for output; optional synthetic input is embedded in code.
    Nothing is mounted from, or published to, the host filesystem.
    """
    skipped = {'state': 'not_checked', 'reason': 'No completed fixture execution'}
    result = {'executed': False, 'execution_attempted': False, 'execution': dict(skipped),
              'artifacts': dict(skipped), 'runtime_scientific': dict(skipped)}
    if not isinstance(image, str) or not re.fullmatch(r'sha256:[0-9a-f]{64}', image):
        result['execution'] = {'state': 'unverifiable', 'reason': 'A trusted local image sha256 ID is required'}
        return result
    docker = docker_command()
    if not docker:
        result['execution'] = {'state': 'unverifiable', 'reason': 'Docker unavailable; script was not executed'}
        return result
    outputs = _contract(outputs or {})
    if not 0 < timeout <= 60 or len(code.encode('utf-8')) > 262144:
        raise ValueError('Fixture timeout must be <= 60 s and code <= 256 KiB')
    if cancelled and cancelled():
        result['execution'] = {'state': 'unverifiable', 'reason': 'Cancelled before execution'}
        return result
    name = 'pycsamt-fixture-' + uuid.uuid4().hex
    bundle = json.dumps({
        'script': code,
        'values': Path(__file__).with_name('_fixture_values.py').read_text(encoding='utf-8'),
        'contract': outputs,
    })
    # Line 1 carries the supervisor, so the only argument is a quote-free
    # bootstrap that survives any client wrapper (e.g. wsl.exe).
    stdin = ('exec(' + json.dumps(_RUNNER) + ')\n' + bundle).encode('utf-8')
    command = [
        *docker, 'run', '--name', name, '--rm', '--pull=never', '--interactive',
        '--network=none', '--read-only', '--cap-drop=ALL',
        '--security-opt=no-new-privileges', '--user=65534:65534',
        '--memory=512m', '--memory-swap=512m', '--cpus=1', '--pids-limit=32',
        '--ulimit=cpu=30:30', '--ulimit=fsize=1048576:1048576',
        '--log-driver=none', '--workdir=/work',
        '--tmpfs=/work:rw,noexec,nosuid,nodev,size=16m,mode=1777',
        '--tmpfs=/tmp:rw,noexec,nosuid,nodev,size=16m,mode=1777',
        '--env=MPLBACKEND=Agg', '--env=MPLCONFIGDIR=/tmp/matplotlib',
        '--env=OPENBLAS_NUM_THREADS=1', '--entrypoint=python',
        image, '-I', '-c', _BOOTSTRAP,
    ]
    try:
        result['execution_attempted'] = True
        status, output = _bounded_run(command, timeout, cancelled, stdin=stdin)
        if status:
            raise RuntimeError(f'Container failed ({status}): {output[-1000:]}')
        data = json.loads(output)
        if (type(data.get('returncode')) is not int
                or set(data.get('outputs', {})) != set(outputs)
                or any(v.get('state') not in {'passed', 'failed'}
                       for v in data['outputs'].values())):
            raise ValueError('Invalid fixture report')
        result['executed'] = True
        result['execution'] = {
            'state': 'passed' if data['returncode'] == 0 else 'failed',
            'reason': f"Isolated fixture exit code {data['returncode']}",
            'image': image,
        }
        for key, names in (
            ('artifacts', list(outputs)),
            ('runtime_scientific', [n for n, r in outputs.items() if r['kind'] != 'file']),
        ):
            items = {n: data['outputs'][n] for n in names}
            result[key] = {
                'state': ('failed' if any(v['state'] == 'failed' for v in items.values())
                          else 'passed' if items else 'not_checked'),
                'reason': 'Declared fixture contracts only; not full scientific correctness',
                'items': items,
            }
    except (OSError, ValueError, RuntimeError, TypeError, AttributeError) as exc:
        result['execution'] = {'state': 'unverifiable', 'reason': str(exc)[:1200]}
    finally:
        # Killing the CLI alone does not stop the container.
        try:
            cleanup = subprocess.run([*docker, 'rm', '-f', name],
                                     stdout=subprocess.DEVNULL,
                                     stderr=subprocess.DEVNULL, timeout=5)
            if cleanup.returncode:
                result['cleanup'] = 'Container absent or removal failed; check Docker status'
        except (OSError, subprocess.TimeoutExpired):
            result['cleanup'] = 'Container cleanup could not be confirmed'
    return result
