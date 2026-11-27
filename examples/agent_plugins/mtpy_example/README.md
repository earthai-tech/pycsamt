# MTpy example agent plugin

This shows how a package other than pyCSAMT can plug its own workflows
into the pyCSAMT agent loop. MTpy is **not** needed to import,
register, classify or validate; only the real step functions use it.

## The three pieces

A plugin is a `pycsamt.agents.plugins.DomainPlugin` bundling:

1. **`WorkflowRegistry`**: `descriptions`, `valid_workflows`,
   `dispatchable`, `classify(text)`, `normalise(wf)`, `family(wf)`
   and `steps(wf)`, which returns the `(step, agent, input, description)`
   tuples the orchestrator runs.
2. **`PlanValidator`**: `errors(plan)` (any error blocks the run),
   `needs_review(plan)` (`True` asks for approval of the plan hash),
   `cross_check(text, intent, wf)` and `route_flags(text, wf)` (a
   deterministic second opinion on the LLM's route).
3. **`record_factory(payload) -> RecordWriter`**: an object with
   `emit`, `snapshot`, `persist` and `corpus_manifest`. Subclassing
   `pycsamt.agents._provenance.RunRecord` gives you the first three.

The protocols are structural, so plain classes work; `DomainPlugin`
checks them when you build it.

## Using it

```python
import sys; sys.path.insert(0, "examples/agent_plugins")
from mtpy_example import make_plugin
from pycsamt.agents.plugins import register_plugin, get_plugin

register_plugin(make_plugin())
p = get_plugin("mtpy_example")
p.registry.classify("estimate strike for site 12")   # 'strike'
```

To ship a plugin from your own distribution, add an entry point:

```toml
[project.entry-points."pycsamt.agent_plugins"]
mtpy_example = "mtpy_example:make_plugin"
```

`get_plugin("mtpy_example")` then finds it automatically.
`get_plugin()` with no name always returns the built-in `pycsamt`
plugin, so existing behaviour does not change.
