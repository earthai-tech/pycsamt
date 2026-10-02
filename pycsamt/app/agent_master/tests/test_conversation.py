# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Multi-turn ownership, reset, artifact honesty and bounded context."""
import threading

from pycsamt.agents._generation import GenerationInput
from pycsamt.app.agent_master import _conversation as conversation
from pycsamt.app.agent_master.callbacks import chat
from pycsamt.assistant.memory import SessionState


def history_for(path, mid, **facts):
    state = SessionState(project_id=conversation.project_id({}), edi_path=path,
                         line=mid, facts=facts)
    return [{'role': 'user', 'content': 'load', 'mid': mid},
            {'role': 'assistant', 'content': 'loaded', 'memory': state.to_dict()}]


def test_parallel_conversations_have_independent_worker_memory(monkeypatch):
    barrier = threading.Barrier(2)
    def dispatch(jid, *args):
        state = chat._session()
        barrier.wait(timeout=5)
        chat._update_job(jid, status='done', result=state.edi_path, kind='answer')
    monkeypatch.setattr(chat, '_run_agent_impl', dispatch)
    jobs = [chat._new_job(), chat._new_job()]
    threads = [threading.Thread(target=chat._run_agent, args=(
        jid, 'explain', {}, {'provider': 'offline'}, {}, history_for(path, mid)))
        for jid, path, mid in zip(jobs, ['/a', '/b'], ['a', 'b'])]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=10)
        assert not thread.is_alive()
    assert [chat._get_job(j)['result'] for j in jobs] == ['/a', '/b']
    assert [chat._get_job(j)['memory']['line'] for j in jobs] == ['a', 'b']


def test_cancelled_late_result_does_not_publish_memory(monkeypatch):
    entered, release = threading.Event(), threading.Event()
    def dispatch(jid, *args):
        entered.set()
        release.wait(timeout=5)
        chat._session().set_data(edi_path='/late')
        chat._update_job(jid, status='done', result='late', kind='workflow')
    monkeypatch.setattr(chat, '_run_agent_impl', dispatch)
    jid = chat._new_job()
    thread = threading.Thread(target=chat._run_agent,
                              args=(jid, 'run', {}, {}, {}, history_for('/a', 'a')))
    thread.start()
    assert entered.wait(timeout=5)
    chat._stop_job_response([], [], {'jid': jid})
    release.set()
    thread.join(timeout=5)
    job = chat._get_job(jid)
    assert job['status'] == 'cancelled'
    assert job['result'] is None and 'memory' not in job


def test_new_chat_does_not_inherit_standalone_context(monkeypatch):
    chat._session().set_data(edi_path='/private')
    captured = []
    monkeypatch.setattr(chat, '_run_agent_impl', lambda *args: captured.append(chat._session().edi_path))
    chat._run_agent(chat._new_job(), 'hello', {}, {}, {}, [])
    assert captured == [None]
    chat._reset_session()


def test_project_switch_discards_prior_code_and_memory():
    history = history_for('/old', 'old')
    history[-1]['code'] = 'secret = 1'
    history.append({'role': 'user', 'content': 'change the script'})
    state, scoped = conversation.restore_memory(history, {'project_root': '/another'}, {})
    assert state.edi_path is None
    assert scoped == history[-1:]
    assert GenerationInput.from_chat('change the script', scoped).previous_code == ''
    assert conversation.scoped_data(history, {'project_root': '/another'}, {'path': '/old'}) == {}


def test_project_boundary_keeps_subsequent_turns():
    history = history_for('/old', 'old')
    project = conversation.project_id({'project_root': '/another'})
    state = SessionState(project_id=project)
    assert conversation.project_changed(history, project)
    history += [{'role': 'assistant', 'content': 'new project', 'context_reset': True,
                 'memory': state.to_dict()},
                {'role': 'user', 'content': 'keep this preference'},
                {'role': 'assistant', 'content': 'remembered', 'memory': state.to_dict()},
                {'role': 'user', 'content': 'continue'}]
    assert not conversation.project_changed(history, project)
    assert len(conversation.scoped_history(history, project)) == 4


def test_snapshot_preserves_workflow_evidence_across_questions():
    state = SessionState()
    conversation.remember_result(state, {'kind': 'workflow', 'workflow': 'qc',
                                         'result': '2 stations processed', 'figs': {'figure': {}}})
    snapshot = conversation.remember_result(state, {'kind': 'answer', 'result': 'Explanation'})
    restored = SessionState.from_dict(snapshot)
    assert restored.facts['last_workflow_result']['figure_count'] == 1
    assert '2 stations processed' in restored.context_summary()
    assert len(restored.context_summary()) <= 2000
    assert '2 stations processed' in conversation.recall_result('show the previous result', restored)
    assert conversation.recall_result('write code to summarize the previous result', restored) is None
    assert 'No completed workflow' in conversation.recall_result('show the previous result', SessionState())


def test_recent_runs_only_uses_the_current_conversation(agent_app):
    from .test_sidebar import _find
    fn = _find(agent_app, 'am-chat-window', 'am-sidebar-runs.children')
    history = history_for('/a', 'a', runs=[{'workflow': 'qc', 'summary': 'only a',
                                        'status': 'success', 'n_figures': 0}])
    assert 'only a' in str(fn([], history, {}))
    assert 'only a' not in str(fn([], [], {}))


def test_cancelled_progress_stays_cancelled():
    jid = chat._new_job()
    chat._update_job(jid, steps=[{'label': 'working', 'status': 'running'}])
    chat._update_job(jid, status='cancelled')
    chat._update_job(jid, status='done', result='late')
    assert chat._get_job(jid)['steps'][0]['status'] == 'cancelled'


def test_real_dispatch_recall_uses_only_supplied_snapshot():
    history = history_for('/a', 'a', last_workflow_result={
        'workflow': 'qc', 'summary': '2 stations passed', 'figure_count': 0})
    jid = chat._new_job()
    chat._run_agent(jid, 'show my previous result', {}, {'provider': 'offline'}, {}, history)
    job = chat._get_job(jid)
    assert job['status'] == 'done' and '2 stations passed' in job['result']
    assert 'no new computation' in job['execution']
    fresh = chat._new_job()
    chat._run_agent(fresh, 'show my previous result', {}, {'provider': 'offline'}, {}, [])
    assert '2 stations passed' not in chat._get_job(fresh)['result']
    assert 'No completed workflow' in chat._get_job(fresh)['result']


def test_answer_dispatch_suppresses_model_execution_claim(monkeypatch):
    from pycsamt.agents._base import AgentResult
    from pycsamt.agents.package_qa import PackageQAAgent
    monkeypatch.setattr(PackageQAAgent, 'execute', lambda *args: AgentResult(
        'success', 'answer', {'answer': 'I saved the plot to made-up.png.', 'source': 'llm'}))
    jid = chat._new_job()
    chat._dispatch_question(jid, 'how do I plot', llm_prov='claude', api_key=None,
                            sel_model=None, offline=True, history=[], step=lambda *args: None)
    assert 'made-up.png' not in chat._get_job(jid)['result']
    assert 'withheld' in chat._get_job(jid)['result']


def test_same_preview_does_not_replace_another_conversation(agent_app):
    from .test_sidebar import _find
    fn = _find(agent_app, 'am-store-messages', 'am-store-history.data')
    first = [{'role': 'user', 'content': 'hello', 'mid': 'one'},
             {'role': 'assistant', 'content': 'first'}]
    second = [{'role': 'user', 'content': 'hello', 'mid': 'two'},
              {'role': 'assistant', 'content': 'second'}]
    stored = fn(first, [])
    assert len(fn(second, stored)) == 2


def test_poll_does_not_insert_another_conversations_result(agent_app):
    from .test_sidebar import _find
    fn = _find(agent_app, 'am-interval-poll', 'am-chat-window.children')
    jid = chat._new_job()
    chat._update_job(jid, status='done', kind='answer', result='private', conversation_id='one')
    from dash import no_update
    result = fn(1, {'jid': jid}, [], {}, [{'role': 'user', 'mid': 'two'}])
    assert result[0] is no_update and result[3] is no_update
    assert result[1] is True


def test_old_script_remains_available_after_long_conversation():
    code = "output_dir = 'out'\nprint('keep me')"
    history = [{'role': 'assistant', 'code': code, 'generation': {'output_dir': 'out'}}]
    history += [{'role': 'user', 'content': 'explain ' + str(i)} for i in range(30)]
    generation = GenerationInput.from_chat('change only the output directory to new', history)
    assert generation.previous_code == code
    assert generation.previous_output_dir == 'out'
    assert generation.omitted_history_turns > 0
    assert len(generation.history_summary) <= 1000
    assert 'abbreviated' in generation.history_summary


def test_guidance_does_not_advertise_unproduced_artifacts(tmp_path):
    text = conversation.present_result({
        'result': 'Use this script.', 'kind': 'code',
        'script_path': str(tmp_path / 'missing.py'),
        'execution': conversation.result_evidence({'kind': 'code'}),
    })
    assert 'Saved script:' not in text
    assert 'no workflow computation performed' in text


def test_model_cannot_claim_a_saved_figure_in_an_explanation():
    answer = conversation.grounded_answer('I saved your plot to result.png.')
    assert 'withheld' in answer and 'result.png' not in answer
    instructions = 'Run this script to save a plot:\n```python\nfig.savefig("result.png")\n```'
    assert conversation.grounded_answer(instructions) == instructions


def test_execution_questions_are_answered_from_records():
    from types import SimpleNamespace

    empty = SimpleNamespace(facts={})
    script = [{'role': 'assistant', 'code': 'print(1)', 'content': 'Here is the script.'}]
    answer = conversation.execution_status(
        'You just wrote a script. Did it actually run and create the CSV?', empty, script)
    assert answer.startswith('No. The script in this conversation was generated, not executed')
    ran = SimpleNamespace(facts={'last_workflow_result': {'workflow': 'qc', 'figure_count': 2}})
    assert 'last recorded workflow run was qc' in conversation.execution_status(
        'Did the workflow run?', ran, [])
    assert conversation.execution_status('Did the workflow run?', empty, []).startswith(
        'No. Nothing has been executed')
    # Requests and explanations are not execution questions.
    for text in ('Run the QC workflow.', 'How do I save a figure?',
                 'Write a script that creates a CSV.'):
        assert conversation.execution_status(text, empty, script) is None


def test_passive_execution_claims_are_withheld():
    claim = 'Yes, the script was executed, and it created the CSV file qc.csv.'
    assert 'withheld' in conversation.grounded_answer(claim)
    honest = 'No. The script was not executed, so no CSV exists yet.'
    assert conversation.grounded_answer(honest) == honest
    notice = 'No workflow was run for this answer; here is how to run it.'
    assert conversation.grounded_answer(notice) == notice
    hedged = 'No, but the script was executed earlier.'
    assert 'withheld' in conversation.grounded_answer(hedged)


def test_result_uses_actual_sources_and_saved_script(tmp_path):
    script = tmp_path / 'script.py'
    script.write_text('pass')
    text = conversation.present_result({
        'result': 'Review this script.', 'script_path': str(script),
        'validation': {'api_evidence': [{'path': 'pycsamt/emtools/ss.py'}]},
    })
    assert 'pycsamt/emtools/ss.py' in text
    assert 'Saved script:' in text and 'python' in text


def test_new_chat_cancels_job_and_clears_pending_choices(agent_app):
    from .test_sidebar import _find
    fn = _find(agent_app, 'am-btn-new-chat', 'am-chat-window.children')
    jid = chat._new_job()
    result = fn(1, None, {'jid': jid})
    assert chat._get_job(jid)['status'] == 'cancelled'
    assert result[5:] == (True, {}, {}, {})
