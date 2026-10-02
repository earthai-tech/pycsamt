# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Static script checks with explicit, optional isolated fixture execution."""

from __future__ import annotations

import ast
import math
from pathlib import Path

from ._static_api import StaticAPI, signature

__all__ = [
    "validate_generated_code",
    "validate_script_file",
    "validation_summary",
]


def _check(state, reason, **details):
    return {"state": state, "reason": reason, **details}


def validate_generated_code(code, *, root=None, execute_fixture=False,
                            fixture_image=None, fixture_outputs=None,
                            fixture_timeout=20, cancelled=None):
    """``ok`` means no detected errors in performed checks, not full correctness.

    Arbitrary fixture execution fails closed without an isolated executor.
    """
    checks = {
        name: _check("not_checked", "Not reached")
        for name in (
            "syntax",
            "imports",
            "arguments",
            "scientific",
            "execution",
            "artifacts",
        )
    }
    checks["execution"] = _check(
        "unverifiable" if execute_fixture else "not_checked",
        "No isolated generated-code executor is configured; script was not executed.",
    )
    checks["artifacts"] = _check(
        "not_checked",
        "No isolated execution; output files and runtime arrays were not inspected.",
    )
    checks["runtime_imports"] = _check(
        "not_checked",
        "Modules were not imported; optional dependency availability and import side effects were not tested.",
    )
    report = {
        "ok": False,
        "status": "failed",
        "syntax_ok": False,
        "errors": [],
        "warnings": [],
        "checked": [],
        "checks": checks,
        "executed": False,
        "api_evidence": [],
    }
    try:
        if not (code or "").strip():
            raise SyntaxError("Empty code")
        tree = ast.parse(code)
        compile(tree, "<generated>", "exec")
    except (SyntaxError, ValueError, TypeError) as exc:
        report["errors"].append(f"SyntaxError: {exc}")
        checks["syntax"] = _check("failed", str(exc))
        return report
    report["syntax_ok"] = True
    checks["syntax"] = _check(
        "passed", "Parsed and compiled without execution"
    )
    api = StaticAPI(root)
    aliases, instances, imports, calls, scientific = {}, {}, [], [], []

    def resolve_name(node):
        if isinstance(node, ast.Name):
            return aliases.get(node.id) or instances.get(node.id)
        if isinstance(node, ast.Attribute):
            prefix = resolve_name(node.value)
            return prefix + "." + node.attr if prefix else None
        return None

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "pycsamt" or alias.name.startswith(
                    "pycsamt."
                ):
                    aliases[alias.asname or alias.name.split(".")[0]] = (
                        alias.name if alias.asname else "pycsamt"
                    )
                    state, _, evidence = api.resolve(alias.name)
                    imports.append(
                        _check(state, alias.name, evidence=evidence)
                    )
                else:
                    imports.append(
                        _check(
                            "not_checked",
                            f"External import {alias.name}: dependency availability was not tested",
                        )
                    )
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            for alias in node.names:
                name = module + "." + alias.name
                if node.level or alias.name == "*":
                    imports.append(
                        _check(
                            "unverifiable", f"Relative/wildcard import: {name}"
                        )
                    )
                elif module == "pycsamt" or module.startswith("pycsamt."):
                    state, _, evidence = api.resolve(name)
                    aliases[alias.asname or alias.name] = name
                    imports.append(_check(state, name, evidence=evidence))
                else:
                    imports.append(
                        _check(
                            "not_checked",
                            f"External import {name}: dependency availability was not tested",
                        )
                    )
    # Rebinding/shadowing makes global alias inference unsafe. Report it,
    # rather than applying the imported callable's signature to another object.
    shadowed = set()
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Name)
            and isinstance(node.ctx, ast.Store)
            and node.id in aliases
        ):
            shadowed.add(node.id)
        if (
            isinstance(
                node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
            )
            and node.name in aliases
        ):
            shadowed.add(node.name)
        if isinstance(node, ast.arg) and node.arg in aliases:
            shadowed.add(node.arg)
    for name in shadowed:
        aliases.pop(name, None)
        calls.append(
            _check(
                "unverifiable",
                f"Imported binding is reassigned or shadowed: {name}",
            )
        )
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            name = resolve_name(node.value.func)
            if name:
                state, target, _ = api.resolve(name)
                if state == "passed" and isinstance(target, ast.ClassDef):
                    for lhs in node.targets:
                        if isinstance(lhs, ast.Name):
                            instances[lhs.id] = name
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = resolve_name(node.func)
        if not name:
            continue
        state, target, evidence = api.resolve(name)
        bound = (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id in instances
        )
        if isinstance(target, ast.FunctionDef) and any(
            isinstance(d, ast.Name) and d.id == "classmethod"
            for d in target.decorator_list
        ):
            bound = True
        sig = signature(target, bound=bound) if target is not None else None
        if state == "failed":
            calls.append(
                _check("failed", f"Unknown API: {name}", line=node.lineno)
            )
        elif (
            sig is None
            or any(isinstance(a, ast.Starred) for a in node.args)
            or any(k.arg is None for k in node.keywords)
        ):
            calls.append(
                _check(
                    "unverifiable",
                    f"Dynamic or unresolved call: {name}",
                    line=node.lineno,
                )
            )
        else:
            try:
                sig.bind(
                    *[object() for _ in node.args],
                    **{k.arg: object() for k in node.keywords},
                )
                calls.append(_check("passed", name, line=node.lineno))
            except TypeError as exc:
                calls.append(
                    _check("failed", f"{name}: {exc}", line=node.lineno)
                )
        if evidence and sig is not None:
            report["api_evidence"].append(
                {
                    "symbol": name,
                    "signature": ast.unparse(target.args)
                    if isinstance(
                        target, (ast.FunctionDef, ast.AsyncFunctionDef)
                    )
                    else "constructor fields: " + ", ".join(sig.parameters),
                    **evidence,
                }
            )
        if name == "pycsamt.forward.synthetic.LayeredModel":
            values = {k.arg: k.value for k in node.keywords}
            values.update(
                {
                    key: value
                    for key, value in zip(
                        ("resistivity", "thickness"), node.args
                    )
                }
            )
            try:
                rho = ast.literal_eval(values["resistivity"])
                thickness = ast.literal_eval(values["thickness"])
                valid = (
                    isinstance(rho, (list, tuple))
                    and isinstance(thickness, (list, tuple))
                    and len(rho) >= 1
                    and len(thickness) == len(rho) - 1
                    and all(
                        not isinstance(x, bool) and math.isfinite(x) and x > 0
                        for x in [*rho, *thickness]
                    )
                )
                scientific.append(
                    _check(
                        "passed" if valid else "failed",
                        "Literal LayeredModel values must be finite positive vectors with n-1 thicknesses for n resistivities",
                    )
                )
            except (KeyError, ValueError, TypeError, OverflowError):
                scientific.append(
                    _check(
                        "unverifiable",
                        "LayeredModel arrays are computed at runtime",
                    )
                )
        if name == "pycsamt.emtools.ss.apply_ss_factors":
            value = next(
                (k.value for k in node.keywords if k.arg == "factors"),
                node.args[1] if len(node.args) > 1 else None,
            )
            try:
                factors = ast.literal_eval(value)
                valid = (
                    isinstance(factors, dict)
                    and bool(factors)
                    and all(
                        not isinstance(v, bool) and math.isfinite(v) and v > 0
                        for v in factors.values()
                    )
                )
                scientific.append(
                    _check(
                        "passed" if valid else "failed",
                        "Literal correction factors must be a nonempty mapping of finite positive values",
                    )
                )
            except (ValueError, TypeError, OverflowError):
                scientific.append(
                    _check(
                        "unverifiable",
                        "Correction factors are computed at runtime; positivity was not checked",
                    )
                )
    for name, entries in (
        ("imports", imports),
        ("arguments", calls),
        ("scientific", scientific),
    ):
        states = {e["state"] for e in entries}
        state = (
            "failed"
            if "failed" in states
            else "unverifiable"
            if states & {"unverifiable", "not_checked"}
            else "passed"
            if entries
            else "not_checked"
        )
        checks[name] = _check(
            state,
            "Static checkout inspection only; runtime behavior is not verified",
            items=entries,
        )
        for entry in entries:
            if entry["state"] == "failed":
                report["errors"].append(entry["reason"])
            elif entry["state"] in {"unverifiable", "not_checked"}:
                report["warnings"].append(entry["reason"])
            elif name == "imports":
                report["checked"].append(entry["reason"])
    report["ok"] = not report["errors"]
    report["status"] = (
        "failed"
        if report["errors"]
        else "unverifiable"
        if report["warnings"]
        else "passed"
    )
    report["scope"] = (
        "Static checks only. Execution, artifacts and runtime scientific invariants are not certified."
    )
    if execute_fixture and fixture_image and report["ok"]:
        from .fixture_execution import execute_fixture as run_fixture

        runtime = run_fixture(code, image=fixture_image, outputs=fixture_outputs,
                              timeout=fixture_timeout, cancelled=cancelled)
        report["executed"] = runtime["executed"]
        report["execution_attempted"] = runtime.get("execution_attempted", False)
        for name in ("execution", "artifacts", "runtime_scientific"):
            checks[name] = runtime[name]
            if runtime[name]["state"] == "failed":
                report["errors"].append(f"{name}: {runtime[name]['reason']}")
                report["errors"].extend(
                    f"{path}: {item['reason']}"
                    for path, item in runtime[name].get("items", {}).items()
                    if item["state"] == "failed"
                )
            elif runtime[name]["state"] == "unverifiable":
                report["warnings"].append(runtime[name]["reason"])
        if runtime.get("cleanup"):
            report["warnings"].append(runtime["cleanup"])
        report["ok"] = not report["errors"]
        report["status"] = ("failed" if report["errors"] else
                            "unverifiable" if report["warnings"] else "passed")
        report["scope"] = "Static checks and optional isolated fixture contracts; not full scientific correctness."
    return report


def validation_summary(report):
    checks = report.get("checks", {})
    if not checks:
        if not report.get("syntax_ok", False):
            return "Validation: syntax error — " + "; ".join(
                report.get("errors", [])
            )
        return (
            "Validation: some symbols could not be verified — "
            + "; ".join(report.get("errors", []))
            if not report.get("ok")
            else "Validation: legacy report; detailed checks unavailable. Script not executed."
        )
    return (
        "Validation: "
        + "; ".join(
            f"{name}: {check['state']}" for name, check in checks.items()
        )
        + ".\n"
        + "\n".join(report.get("errors", [])[:5])
        + (
            "\nRepair: "
            + str(len(report["repair"].get("attempts", [])))
            + " attempt(s); "
            + report["repair"].get("stop_reason", "unknown")
            + "."
            if report.get("repair")
            else ""
        )
        + "\nChecks do not establish full scientific correctness. "
        + ("Script executed only in an isolated fixture; outputs were not published."
           if report.get("executed") else
           "Fixture execution was attempted but completion could not be verified."
           if report.get("execution_attempted") else "Script not executed.")
    )


def validate_script_file(path):
    try:
        return validate_generated_code(Path(path).read_text(encoding="utf-8"))
    except (OSError, UnicodeError) as exc:
        result = validate_generated_code("")
        result["errors"] = [f"Cannot read script: {exc}"]
        return result
