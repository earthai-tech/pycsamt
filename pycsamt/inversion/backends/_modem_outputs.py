"""Explicit, fresh 3-D ModEM output validation; no directory discovery."""
from pathlib import Path
import hashlib
import re
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_outputs(raw, workdir, inputs):
    """Reject reused/aliased output paths before a process is launched."""
    if not isinstance(raw, dict) or set(raw) != {"model", "response", "log"}:
        raise ValueError("outputs must declare model, response and log paths")
    root = Path(workdir).resolve()
    paths = {}
    for role, name in raw.items():
        p = Path(name)
        p = (root / p).resolve() if not p.is_absolute() else p.resolve()
        if not p.is_relative_to(root) or p == root or p.exists():
            raise ValueError(f"Output must be absent and inside workdir: {role}")
        paths[role] = p
    if len(set(paths.values())) != 3:
        raise ValueError("Output paths must be distinct")
    protected = {Path(p).resolve(): digest(p) for p in inputs}
    if set(paths.values()) & set(protected):
        raise ValueError("Output aliases an input")
    return paths, protected


def _rows(path):
    rows = {}
    for line in Path(path).read_text().splitlines():
        t = line.split()
        if not t or t[0].startswith(("#", ">")) or line.strip().startswith("ModEM impedance responses for"):
            continue
        if len(t) != 11:
            raise ValueError("Malformed ModEM response row")
        period = float(t[0])
        xyz = np.array([float(v) for v in t[4:7]])
        value = complex(float(t[8]), float(t[9]))
        error = float(t[10])
        if period <= 0 or not np.isfinite([period, *xyz, value.real, value.imag, error]).all():
            raise ValueError("Nonfinite or invalid response row")
        key = period, t[1], t[7]
        if key in rows:
            raise ValueError("Duplicate response key")
        rows[key] = value, error, xyz
    if not rows:
        raise ValueError("Empty response data")
    return rows


def _headers(path):
    """Require complete block metadata instead of accepting parser defaults."""
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith(("#", "ModEM impedance responses for"))]
    blocks = {}
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            raise ValueError("Response rows require a complete block header")
        meta = lines[i:i+6]
        if len(meta) != 6 or not all(v.startswith(">") for v in meta):
            raise ValueError("Incomplete ModEM block header")
        meta = [v[1:].strip() for v in meta]
        family, sign, units = meta[:3]
        if family in blocks or not sign or not units:
            raise ValueError("Duplicate or empty ModEM block metadata")
        rotation = float(meta[3]); origin = tuple(map(float, meta[4].split()))
        counts = tuple(map(int, meta[5].split()))
        if len(origin) != 2 or len(counts) != 2 or min(counts) <= 0 or not np.isfinite([rotation,*origin]).all():
            raise ValueError("Invalid ModEM block geometry/counts")
        i += 6; rows = []
        while i < len(lines) and not lines[i].startswith(">"):
            t = lines[i].split()
            if len(t) != 11:
                raise ValueError("Malformed ModEM block row")
            rows.append(t); i += 1
        if (len({float(t[0]) for t in rows}), len({t[1] for t in rows})) != counts:
            raise ValueError("Header counts disagree with rows")
        blocks[family] = (sign, units, rotation, origin, counts,
                          frozenset((float(t[0]),t[1],t[7]) for t in rows))
    if not blocks:
        raise ValueError("Missing ModEM block headers")
    return blocks


def compatible_headers(observation, prediction):
    a, b = _headers(observation), _headers(prediction)
    if a.keys() != b.keys():
        return False
    for family, x in a.items():
        y = b[family]
        # Native output prints origin metadata to three decimal degrees.
        if x[:3] != y[:3] or x[4:] != y[4:] or not np.allclose(x[3],y[3],atol=5e-4,rtol=0):
            return False
    return True


def native_control(path):
    """Validate the supported eight-record a36/g15.7/i4 native layout."""
    from ...models.modem.control import _KEYS
    lines = Path(path).read_text().splitlines()
    if len(lines) != len(_KEYS):
        raise ValueError("Native control requires exactly eight ordered records")
    values = {}
    for index, (line, (attr, label)) in enumerate(zip(lines, _KEYS)):
        if line[:36].strip() != label + ":":
            raise ValueError("Native control label/order/column mismatch")
        field = line[36:]
        if index == 0:
            if not field.strip() or len(field) > 80:
                raise ValueError("Invalid native output stem")
            values[attr] = field.strip()
        elif index == 7:
            if len(field) > 4 or not field.strip().isdigit() or int(field) < 1:
                raise ValueError("Native iteration limit must fit i4")
            values[attr] = int(field)
        else:
            if len(field) > 15 or "." not in field:
                raise ValueError("Native real field requires decimal point within g15.7")
            value = float(field.replace("D", "E").replace("d", "e"))
            if not np.isfinite(value) or value <= 0:
                raise ValueError("Native real controls must be finite and positive")
            values[attr] = value
    return values


def classify_stop(log_path, control_path, iteration, rms):
    """Describe source-supported stopping evidence without asserting convergence."""
    controls = native_control(control_path)
    target = controls["target_rms"]
    limit_raw = controls["max_iterations"]
    text = Path(log_path).read_text()
    pattern = r"NLCG iterations:\s*(\d+)\s+function evaluations:\s*(\d+)"
    endings = list(re.finditer(pattern, text))
    local = "Unable to get out of a local minimum" in text
    info = {"target_rms": target, "max_iterations": int(limit_raw)}
    if not np.isfinite(rms) or rms < 0 or iteration < 0 or iteration > limit_raw:
        raise ValueError("Invalid final iteration or RMS")
    if len(endings) > 1:
        raise ValueError("Multiple native completion footers")
    if local:
        if endings:
            raise ValueError("Contradictory local-minimum return and completion footer")
        return {**info, "reason": "local_minimum"}
    if not endings:
        return {**info, "reason": "completion_unverified"}
    footer = endings[0]
    if int(footer.group(1)) != iteration or int(footer.group(2)) < iteration:
        raise ValueError("Native completion counters disagree")
    if text[footer.end():].strip():
        raise ValueError("Unexpected content after native completion footer")
    # Use the same conservative uncertainty budget as printed-response/log
    # RMS consistency. Do not turn sub-print-precision differences into success.
    margin = 1e-4 + 1e-3 * max(rms, target)
    info["rms_threshold_margin"] = margin
    if abs(rms-target) <= margin:
        reason = "threshold_ambiguous"
    elif rms + margin < target:
        reason = "target_rms_reached"
    elif iteration == limit_raw:
        reason = "iteration_limit"
    else:
        reason = "other_stop"
    return {**info, "reason": reason}


def validate_outputs(paths, protected, observation, control=None, initial_model=None):
    """Validate fresh triplet and RMS; target/stopping qualification is separate."""
    from ...models.modem.model3d import ModEmModel3D
    from ...models.modem.log import ModEmLog
    for p, expected in protected.items():
        if digest(p) != expected:
            raise ValueError("ModEM input changed during execution")
    for p in paths.values():
        if not p.is_file() or p.stat().st_size == 0 or p.is_symlink():
            raise ValueError("Missing, empty or linked declared output")
    def iteration(p):
        match = re.search(r"_(\d{3,})$", p.stem)
        if not match:
            raise ValueError("Output must have explicit numbered iteration")
        return int(match.group(1))
    n = iteration(paths["model"])
    if paths["model"].stem != paths["response"].stem or iteration(paths["response"]) != n:
        raise ValueError("Model/response iteration mismatch")
    model = ModEmModel3D.read(paths["model"])
    if model.rho_loge.size == 0 or not np.isfinite(model.rho_loge).all():
        raise ValueError("Nonfinite model")
    if initial_model is not None:
        initial = ModEmModel3D.read(initial_model)
        for attr in ["x_widths", "y_widths", "z_widths", "origin"]:
            a,b = getattr(initial,attr),getattr(model,attr)
            if a.shape != b.shape or not np.allclose(a,b,rtol=1e-6,atol=1e-3):
                raise ValueError("Output model geometry differs from input")
        if initial.rho_loge.shape != model.rho_loge.shape or not np.isclose(initial.rotation,model.rotation):
            raise ValueError("Output model shape/rotation differs from input")
    log = ModEmLog.read(paths["log"])
    if not log.n_iter or int(log.iterations[-1]) != n or not np.isfinite(log.rms).all():
        raise ValueError("Log/model iteration mismatch or invalid RMS")
    if any(not np.isfinite(v).all() for v in
           [log.objective, log.model_norm, log.lagrange, log.alpha]) or (np.diff(log.iterations) <= 0).any():
        raise ValueError("Invalid iteration history")
    if not compatible_headers(observation, paths["response"]):
        raise ValueError("Observed/predicted block metadata differ")
    obs, pred = _rows(observation), _rows(paths["response"])
    if obs.keys() != pred.keys():
        raise ValueError("Observed/predicted keys differ")
    residuals = []
    for key, (value, error, xyz) in obs.items():
        pv, _, pxyz = pred[key]
        if not np.allclose(xyz, pxyz, atol=1e-3, rtol=0):
            raise ValueError("Response coordinates differ")
        if error <= 0:
            raise ValueError("Nonpositive observational error")
        if error < 1e12:
            residuals.append((value-pv)/error)
    if not residuals:
        raise ValueError("No active observations")
    rms = float(np.sqrt(np.mean(np.abs(residuals)**2)/2))
    if not np.isclose(rms, log.final_rms, rtol=1e-3, atol=1e-4):
        raise ValueError("Recomputed RMS disagrees with log")
    return {"model": model, "log": log, "rms": rms, "iteration": n,
            "n_real_scalars": 2*len(residuals),
            "sha256": {role: digest(p) for role,p in paths.items()},
            "stopping": (classify_stop(paths["log"], control, n, rms)
                         if control is not None else {"reason": "control_unverified"}),
            "converged": False}
