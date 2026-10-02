# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Executor contracts use fake Docker transport; no local generated execution."""
import json
import os
from io import BytesIO
from types import SimpleNamespace

import pytest

from pycsamt.assistant.tools import fixture_execution as fixture
from pycsamt.assistant.tools._fixture_values import validate_values
from pycsamt.assistant.tools.validation_tools import validate_generated_code, validation_summary

IMAGE = 'sha256:' + 'a' * 64
LIVE_DOCKER = pytest.mark.skipif(
    not (os.environ.get('PYCSAMT_FIXTURE_IMAGE') and fixture.docker_command()),
    reason='Requires a configured local Docker image and client')


@pytest.mark.parametrize('value', [{}, {'s': 0}, {'s': -1}, {'s': float('nan')},
                                  {'s': float('inf')}, {'s': True}, {'s': '1'}])
def test_invalid_runtime_factors(value):
    with pytest.raises(ValueError):
        validate_values(value, {'kind': 'positive_factors'})


@pytest.mark.parametrize('value', [[[1], [2]], [[1, 2]], [[1, 2], [3, None]],
                                  [[1, 2], [3, float('nan')]]])
def test_invalid_runtime_arrays(value):
    with pytest.raises(ValueError):
        validate_values(value, {'kind': 'array', 'shape': [2, 2]})


def test_valid_runtime_values():
    validate_values({'s': 1.1}, {'kind': 'positive_factors'})
    validate_values([[1, 2], [3, 4]], {'kind': 'array', 'shape': [2, 2]})


@pytest.fixture
def transport(monkeypatch):
    commands = []
    monkeypatch.delenv('PYCSAMT_DOCKER', raising=False)
    monkeypatch.setattr(fixture.shutil, 'which', lambda name: '/docker')
    monkeypatch.setattr(fixture.subprocess, 'run',
                        lambda command, **kw: commands.append(command) or SimpleNamespace(returncode=0))
    return commands


def test_no_docker_never_executes(monkeypatch):
    monkeypatch.setattr(fixture.shutil, 'which', lambda name: None)
    result = fixture.execute_fixture('raise RuntimeError()', image=IMAGE)
    assert not result['executed']
    assert result['execution']['state'] == 'unverifiable'


@pytest.mark.parametrize('image', ['latest', '--privileged', None])
def test_requires_immutable_image(image):
    assert not fixture.execute_fixture('pass', image=image)['executed']


@pytest.mark.parametrize('name', ['../secret', '/etc/passwd', 'a\\b', '.', 'C:/file'])
def test_rejects_output_escape(name):
    with pytest.raises(ValueError):
        fixture._contract({name: {'kind': 'file'}})


def test_limits_and_failure_are_reported(monkeypatch, transport):
    def run(command, timeout, cancelled, stdin=None):
        for option in ('--network=none', '--read-only', '--cap-drop=ALL',
                       '--security-opt=no-new-privileges', '--user=65534:65534',
                       '--memory=512m', '--memory-swap=512m', '--pids-limit=32',
                       '--pull=never', '--cpus=1'):
            assert option in command
        # Inputs travel on stdin: no host path is mounted into the container.
        assert '--interactive' in command
        assert not [a for a in command if a.startswith(('--mount', '-v', '--volume'))]
        assert command[-3:] == ['-I', '-c', fixture._BOOTSTRAP]
        supervisor, bundle = stdin.decode('utf-8').split('\n', 1)
        assert supervisor == 'exec(' + json.dumps(fixture._RUNNER) + ')'
        assert json.loads(bundle)['script'] == 'pass'
        assert json.loads(bundle)['contract'] == {'factors.json': {'kind': 'positive_factors'}}
        return 0, json.dumps({'returncode': 1, 'outputs': {
            'factors.json': {'state': 'failed', 'reason': 'not positive'}}})
    monkeypatch.setattr(fixture, '_bounded_run', run)
    report = validate_generated_code('pass', execute_fixture=True, fixture_image=IMAGE,
                                    fixture_outputs={'factors.json': {'kind': 'positive_factors'}})
    assert report['executed'] and not report['ok']
    assert report['checks']['runtime_scientific']['state'] == 'failed'
    assert 'Script not executed' not in validation_summary(report)
    assert transport[-1][1:3] == ['rm', '-f']


def test_timeout_cleans_container(monkeypatch, transport):
    def timeout(*args, **kwargs):
        raise RuntimeError('Fixture execution timed out')
    monkeypatch.setattr(fixture, '_bounded_run', timeout)
    report = fixture.execute_fixture('pass', image=IMAGE)
    assert report['execution']['state'] == 'unverifiable'
    assert 'timed out' in report['execution']['reason']
    assert transport[-1][1:3] == ['rm', '-f']
    assert report['execution_attempted']


def test_client_cancellation_kills_client(monkeypatch):
    class Process:
        stdout = BytesIO(b'')
        killed = False

        def poll(self):
            return -1 if self.killed else None

        def kill(self):
            self.killed = True

        def wait(self, timeout):
            return -1

    process = Process()
    monkeypatch.setattr(fixture.subprocess, 'Popen', lambda *a, **kw: process)
    with pytest.raises(RuntimeError, match='cancelled'):
        fixture._bounded_run(['docker'], 20, lambda: True)
    assert process.killed


def test_client_rejects_oversized_report(monkeypatch):
    process = SimpleNamespace(stdout=BytesIO(b'x' * 70000), returncode=0,
                              poll=lambda: 0, wait=lambda **kw: 0)
    monkeypatch.setattr(fixture.subprocess, 'Popen', lambda *a, **kw: process)
    with pytest.raises(RuntimeError, match='too large'):
        fixture._bounded_run(['docker'], 20, None)


def test_container_supervisor_compiles():
    compile(fixture._RUNNER, '<fixture supervisor>', 'exec')
    compile(fixture._BOOTSTRAP, '<bootstrap>', 'exec')
    assert '"' not in fixture._BOOTSTRAP and "'" not in fixture._BOOTSTRAP


def test_configured_client_prefix(monkeypatch):
    monkeypatch.setenv('PYCSAMT_DOCKER', 'wsl.exe -e docker')
    monkeypatch.setattr(fixture.shutil, 'which', lambda name: 'C:/wsl.exe' if name == 'wsl.exe' else None)
    assert fixture.docker_command() == ['wsl.exe', '-e', 'docker']
    monkeypatch.setattr(fixture.shutil, 'which', lambda name: None)
    assert fixture.docker_command() is None
    monkeypatch.delenv('PYCSAMT_DOCKER')
    monkeypatch.setattr(fixture.shutil, 'which', lambda name: '/usr/bin/docker')
    assert fixture.docker_command() == ['/usr/bin/docker']


def test_client_receives_stdin(monkeypatch):
    received = {}

    class Pipe(BytesIO):
        def close(self):
            received['stdin'] = self.getvalue()

    class Process:
        stdout = BytesIO(b'{}')
        stdin = Pipe()
        returncode = 0

        def poll(self):
            return 0

        def wait(self, timeout):
            return 0

    monkeypatch.setattr(fixture.subprocess, 'Popen', lambda *a, **kw: Process())
    assert fixture._bounded_run(['docker'], 20, None, stdin=b'payload') == (0, '{}')
    assert received['stdin'] == b'payload'


def test_cancel_before_start(monkeypatch, transport):
    report = fixture.execute_fixture('pass', image=IMAGE, cancelled=lambda: True)
    assert not report['executed']
    assert not transport


def test_static_failure_does_not_execute(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('Static failure must prevent fixture execution')
    monkeypatch.setattr(fixture, 'execute_fixture', forbidden)
    report = validate_generated_code('def broken(:', execute_fixture=True, fixture_image=IMAGE)
    assert not report['ok'] and not report['executed']


def test_success_does_not_claim_unrequested_science(monkeypatch, transport):
    monkeypatch.setattr(fixture, '_bounded_run', lambda *args, **kwargs: (0, json.dumps({
        'returncode': 0, 'outputs': {}})))
    result = fixture.execute_fixture('pass', image=IMAGE)
    assert result['execution']['state'] == 'passed'
    assert result['runtime_scientific']['state'] == 'not_checked'


@LIVE_DOCKER
@pytest.mark.parametrize('code,contract,state', [
    ("import json; open('f.json', 'w').write(json.dumps({'s': 1.2}))",
     {'f.json': {'kind': 'positive_factors'}}, 'passed'),
    ("open('f.json', 'w').write('{\"s\": 0}')",
     {'f.json': {'kind': 'positive_factors'}}, 'failed'),
    ("open('a.json', 'w').write('[[1,2],[3,4]]')",
     {'a.json': {'kind': 'array', 'shape': [2, 2]}}, 'passed'),
    ("open('a.json', 'w').write('[[1],[2]]')",
     {'a.json': {'kind': 'array', 'shape': [2, 2]}}, 'failed'),
    ("open('a.json', 'w').write('[NaN]')",
     {'a.json': {'kind': 'array', 'shape': [1]}}, 'failed'),
    ('pass', {'missing.png': {'kind': 'file'}}, 'failed'),
    ("import os; os.symlink('/etc/passwd', 'link')",
     {'link': {'kind': 'file'}}, 'failed'),
])
def test_live_fixture_contracts(code, contract, state):
    result = fixture.execute_fixture(code, image=os.environ['PYCSAMT_FIXTURE_IMAGE'],
                                     outputs=contract)
    assert result['executed'], result
    assert result['execution']['state'] == 'passed', result
    assert result['artifacts']['state'] == state, result


@LIVE_DOCKER
@pytest.mark.parametrize('code', [
    "raise RuntimeError('fixture failure')",
    "open('/input/script.py', 'w').write('overwrite')",
    "open('/etc/fixture-output', 'w').write('escape')",
    "import socket; socket.create_connection(('1.1.1.1', 443), timeout=1)",
])
def test_live_execution_and_boundary_failures(code):
    result = fixture.execute_fixture(code, image=os.environ['PYCSAMT_FIXTURE_IMAGE'])
    assert result['executed'], result
    assert result['execution']['state'] == 'failed', result
