# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Provider-independent cancellation and stage boundaries."""
import threading

import pytest

from pycsamt.agents._base import AgentResult, BaseAgent
from pycsamt.agents._request import RequestCancelled, checkpoint, request_scope


class FakeAgent(BaseAgent):
    def execute(self, data):
        data['called'].append(self.name)
        data['cancel'].set()
        return AgentResult('success', 'finished')


def test_agent_discards_result_when_cancelled_in_operation():
    cancel, called = threading.Event(), []
    with request_scope(cancel.is_set), pytest.raises(RequestCancelled):
        FakeAgent('science').execute({'called': called, 'cancel': cancel})
    assert called == ['science']
    checkpoint()  # request scope is reset, not leaked into subsequent work


def test_cancelled_request_does_not_start_an_agent():
    with pytest.raises(RequestCancelled), request_scope(lambda: True):
        pytest.fail('Cancelled scope must not enter')


@pytest.mark.parametrize('provider,method', [
    ('claude', '_query_claude'), ('openai', '_query_openai'),
    ('gemini', '_query_gemini'), ('deepseek', '_query_deepseek'),
    ('minimax', '_query_minimax'),
])
def test_cloud_response_is_discarded_without_retry(monkeypatch, provider, method):
    cancel, calls = threading.Event(), []
    agent = FakeAgent('cloud', llm_provider=provider, api_key='fake-test-key')
    monkeypatch.setattr(type(agent), 'llm_available', property(lambda self: True))
    def transport(*args):
        calls.append(1)
        cancel.set()
        return 'late answer', 0
    monkeypatch.setattr(agent, method, transport)
    with request_scope(cancel.is_set), pytest.raises(RequestCancelled):
        agent.query_llm('test')
    assert calls == [1]


def test_coordinator_stops_before_next_stage_and_checkpoint(tmp_path):
    from pycsamt.agents.coordinator import AgentCoordinator
    cancel, called = threading.Event(), []
    coord = AgentCoordinator('test', checkpoint_dir=tmp_path)
    coord.add_step('first', FakeAgent('first'))
    coord.add_step('second', FakeAgent('second'))
    with request_scope(cancel.is_set), pytest.raises(RequestCancelled):
        coord.execute({'called': called, 'cancel': cancel})
    assert called == ['first']
    assert not list(tmp_path.glob('*.pkl'))
