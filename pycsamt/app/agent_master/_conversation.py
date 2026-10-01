# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Conversation-owned snapshots and factual result presentation."""
from __future__ import annotations

import hashlib
from pathlib import Path

from pycsamt.assistant.memory import SessionState


def project_id(settings):
    """An explicit project root/registry is the boundary, not the selected line."""
    settings = settings or {}
    root = settings.get('project_root') or Path.cwd()
    registry = settings.get('line_registry')
    if registry and ('\n' in registry or ': ' in registry):
        registry = 'inline-registry:' + hashlib.sha256(registry.encode('utf-8')).hexdigest()[:16]
    elif registry:
        registry = str(Path(registry).expanduser().resolve())
    return str(Path(root).expanduser().resolve()) + (
        '|' + registry if registry else '')


def scoped_history(history, project):
    """Discard earlier project context, keeping any current uncompleted turn."""
    messages = list(history or [])
    for index in range(len(messages) - 1, -1, -1):
        saved = messages[index].get('memory')
        if isinstance(saved, dict):
            if saved.get('project_id', '') != project:
                return messages[index + 1:]
            # A project-switch response becomes the new history boundary.
            if messages[index].get('context_reset'):
                start = index - 1 if index and messages[index - 1].get('role') == 'user' else index
                return messages[start:]
    return messages


def conversation_id(history):
    return next((m.get('mid') for m in history or [] if m.get('mid')), None)


def project_changed(history, project):
    for message in reversed(history or []):
        saved = message.get('memory')
        if isinstance(saved, dict):
            return saved.get('project_id', '') != project
    return False


def restore_memory(history, settings, edi_store):
    project = project_id(settings)
    history = scoped_history(history, project)
    state = SessionState.from_history(history, project_id=project)
    path = (edi_store or {}).get('path')
    if path and path != state.edi_path:
        state.set_data(edi_path=path)
        state.line = None
        state.last_summary = None
        state.facts.pop('execution', None)
    # The settings folder is a default, not a choice made in this conversation:
    # a folder named in a request must be able to override it.
    if settings.get('output_dir'):
        state.facts['default_output_dir'] = settings['output_dir']
    return state, history


def scoped_data(history, settings, edi_store):
    """Do not carry a previous project's still-loaded dataset into a new project."""
    project = project_id(settings)
    for message in reversed(history or []):
        saved = message.get('memory')
        if isinstance(saved, dict):
            if saved.get('project_id', '') != project:
                if (edi_store or {}).get('path') == saved.get('edi_path'):
                    return {}
            break
    return dict(edi_store or {})


def remember_result(state, job):
    """Only record facts supplied by the completed dispatch, not model prose."""
    state.last_summary = str(job.get('result') or job.get('error') or '')[:1200]
    state.facts['execution'] = job.get('execution', 'No computation evidence recorded')
    if job.get('script_path'):
        state.facts['script_path'] = job['script_path']
    if job.get('kind') == 'workflow':
        state.facts['last_workflow_result'] = {
            **state.facts.get('last_workflow_result', {}),
            'workflow': job.get('workflow'),
            'summary': state.last_summary,
            'figure_count': len(job.get('figs') or {}),
        }
    state.turns = []  # text already lives in STORE_MESSAGES; avoid recursive growth
    return state.to_dict()


def result_evidence(job):
    """Separate generated guidance from actual workflow dispatch results."""
    if job.get('execution'):
        return job['execution']
    kind = job.get('kind')
    if kind == 'code':
        return 'Script generated; no workflow computation performed'
    if kind == 'workflow':
        return 'Registered workflow returned; inspect its status and warnings'
    if kind == 'error' or job.get('status') == 'error':
        return 'Request failed; successful computation not established'
    return 'Explanation or suggestion only; no computation performed'


def recall_result(text, state):
    """Answer explicit result recall from recorded facts without another model run."""
    import re

    if re.search(r'\b(code|script|write|generate|rerun|execute|run again)\b', text, re.I):
        return None
    if not re.search(r'\b(last|previous|earlier)\s+(?:workflow|run|result|computation)s?\b', text, re.I):
        return None
    if not re.search(r'\b(what|show|summarize|recap|recall)\b', text, re.I):
        return None
    result = state.facts.get('last_workflow_result')
    if not result:
        return 'No completed workflow result is recorded in this conversation.'
    return (f"Previous {result.get('workflow') or 'workflow'} result:\n"
            + result.get('summary', '')
            + f"\nFigures returned in that run: {result.get('figure_count', 0)}. "
            'This recalls the recorded result; no new computation was run.')


def execution_status(text, state, history=None):
    """Answer "did it run / create the file?" from recorded evidence only.

    Whether something executed is a fact about this conversation, not
    something a model may assert. Generated scripts are never executed by
    Agent Master; only registered workflow runs are recorded as computation.
    """
    import re

    asks = re.search(
        r"(?:^|[.?!]\s+)(?:did|has|have|was|were|is)\b[^.?!]{0,60}?\b(?:run|ran|execut\w*|"
        r"creat\w*|produc\w*|sav\w*|writ\w*|wrote|complet\w*)\b",
        text, re.I)
    subject = re.search(r"\b(?:it|script|code|workflow|csv|file|figure|plot|output|run)\b",
                        text, re.I)
    if not (asks and subject):
        return None
    scripts = [m for m in (history or []) if m.get('role') == 'assistant' and m.get('code')]
    result = (state.facts.get('last_workflow_result') if state is not None else None)
    lines = []
    if scripts:
        lines.append('No. The script in this conversation was generated, not executed: Agent '
                     'Master does not run generated scripts, so it created no CSV, figure or '
                     'other output. Run it yourself after checking its paths and assumptions.')
    if result:
        lines.append(f"The last recorded workflow run was {result.get('workflow') or 'a workflow'}, "
                     f"with {result.get('figure_count', 0)} figure(s) returned; its summary is the "
                     'evidence for what it computed.')
    if not lines:
        lines.append('No. Nothing has been executed in this conversation, so no output file '
                     'exists from it.')
    return ' '.join(lines) + ' This answer comes from the recorded conversation; no new computation was run.'


def present_result(job):
    """Render authoritative status/sources around the agent's direct response."""
    text = job.get('result') or job.get('error') or 'No result was returned.'
    paths = []
    evidence = job.get('citations') or (job.get('validation') or {}).get('api_evidence', [])
    for item in evidence:
        path = (item.get('source_path') or item.get('path')) if isinstance(item, dict) else None
        if path and path not in paths and path not in text:
            paths.append(path)
    if paths:
        text += '\n\nSources: ' + ', '.join(f'`{p}`' for p in paths[:6]) + '.'
    path = job.get('script_path')
    if path and Path(path).is_file():
        text += f'\n\nSaved script: `{path}`. Run with `python "{path}"` after checking its paths and assumptions.'
    if job.get('execution'):
        text += '\n\n' + job['execution'] + '.'
    return text


def grounded_answer(answer):
    """Suppress definite first-person execution claims on an explanation-only path.

    This deliberately narrow guard is not a general truthfulness classifier.
    Instructions, code examples and cited statements are left for model evaluation.
    """
    import re

    prose = re.sub(r'```.*?```', '', answer, flags=re.S)
    claims = re.finditer(
        r"\b(?:I|we)(?:\s+have|'ve)?\s+(?:successfully\s+)?"
        r"(?:ran|executed|computed|calculated|plotted|saved)\b|"
        r"\b(?:I|we)(?:\s+have|'ve)?\s+(?:created|generated|produced)\s+"
        r"(?:the\s+|a\s+|your\s+)?(?:figure|plot|file|inversion|report)\b|"
        r"\b(?:figure|plot|file|inversion)\s+(?:has been|was)\s+"
        r"(?:saved|created|completed|generated)\b|"
        r"\b(?:script|code|workflow|analysis|inversion)\s+(?:has been|was)\s+"
        r"(?:successfully\s+)?(?:executed|run|completed)\b|"
        r"\b(?:script|code|workflow)\s+(?:ran|executed)\s+successfully\b|"
        r"\b(?:it|this|that)\s+(?:created|produced|wrote|saved)\s+"
        r"(?:the\s+|a\s+|an\s+|your\s+)?(?:csv|file|figure|plot|table|report)\b",
        prose, re.I,
    )
    # "No workflow was run" / "the script was not executed" are not claims.
    negated = re.compile(r"\b(?:no|not|never)\s+(?:\w+\s+){0,2}$", re.I)
    performed = any(not negated.search(prose[max(0, m.start() - 30):m.start()])
                    for m in claims)
    if performed:
        return ('The model response claimed work without execution evidence, so it was withheld. '
                'No workflow was run for this answer. Request a workflow to compute results, '
                'or ask for code and instructions to run it yourself.')
    return answer
