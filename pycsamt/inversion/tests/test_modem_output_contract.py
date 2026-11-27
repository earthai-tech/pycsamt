"""Invocation-bound output checks, without native process execution."""
import numpy as np
import pytest
from pycsamt.inversion.backends._modem_outputs import prepare_outputs, validate_outputs
from pycsamt.models.modem.control import ModEmControl
from pycsamt.models.modem.model3d import ModEmModel3D

HEADER = "> Full_Impedance\n> exp(+i omega t)\n> [mV/km]/[nT]\n> 0\n> 0 0\n> 1 1\n"


@pytest.fixture
def contract(tmp_path):
    obs = tmp_path / "obs.dat"
    obs.write_text(HEADER + "1 S 0 0 0 0 0 ZXY 2 2 1\n")
    paths, protected = prepare_outputs(
        {"model": "run_001.rho", "response": "run_001.dat", "log": "run.log"},
        tmp_path, [obs])
    model = ModEmModel3D()
    model.x_widths = model.y_widths = model.z_widths = np.array([100.])
    model.rho_loge = np.log(np.full((1,1,1), 100.))
    model.write(paths["model"])
    paths["response"].write_text(HEADER + "1 S 0 0 0 0 0 ZXY 1 1 1e15\n")
    paths["log"].write_text("Completed NLCG iteration 1\nwith: f= 1 m2= 1 rms= 1 lambda= 1 alpha= 1\n")
    return paths, protected, obs


def test_fresh_triplet_and_real_scalar_rms(contract):
    paths, protected, obs = contract
    result = validate_outputs(paths, protected, obs)
    assert result["rms"] == pytest.approx(1.)
    assert result["n_real_scalars"] == 2
    assert result["iteration"] == 1
    assert result["converged"] is False
    assert set(result["sha256"]) == {"model", "response", "log"}


@pytest.mark.parametrize("damage", ["missing", "empty", "input", "iteration", "rms", "keys", "coords", "duplicate", "nan"])
def test_invalid_triplet_rejected(contract, damage):
    paths, protected, obs = contract
    if damage == "missing": paths["model"].unlink()
    if damage == "empty": paths["log"].write_text("")
    if damage == "input": obs.write_text("changed")
    if damage == "iteration": paths["log"].write_text("START: f= 1 m2= 1 rms= 1 lambda= 1 alpha= 1")
    if damage == "rms": paths["log"].write_text("Completed NLCG iteration 1\nwith: f= 1 m2= 1 rms= 0.1 lambda= 1 alpha= 1")
    if damage == "keys": paths["response"].write_text(HEADER + "1 OTHER 0 0 0 0 0 ZXY 1 1 1e15")
    if damage == "coords": paths["response"].write_text(HEADER + "1 S 0 0 100 0 0 ZXY 1 1 1e15")
    if damage == "duplicate": paths["response"].write_text(paths["response"].read_text()*2)
    if damage == "nan": paths["response"].write_text(HEADER + "1 S 0 0 0 0 0 ZXY nan 1 1e15")
    with pytest.raises((ValueError, FileNotFoundError)):
        validate_outputs(paths, protected, obs)


@pytest.mark.parametrize("kind", ["existing", "outside", "duplicate"])
def test_invalid_output_paths_rejected_before_launch(tmp_path, kind):
    obs = tmp_path / "obs"; obs.write_text("input")
    paths = {"model": "m_001.rho", "response": "m_001.dat", "log": "run.log"}
    if kind == "existing": (tmp_path / paths["model"]).write_text("stale")
    if kind == "outside": paths["model"] = "../outside.rho"
    if kind == "duplicate": paths["log"] = paths["response"]
    with pytest.raises(ValueError): prepare_outputs(paths, tmp_path, [obs])


@pytest.mark.parametrize("write_outputs", [True, False])
def test_adapter_validates_declared_outputs_after_execution(tmp_path, monkeypatch, write_outputs):
    import pycsamt.models.modem as modem
    from pycsamt.inversion.backends.modem import ModEMBackend
    from pycsamt.inversion.config import InversionConfig
    obs = tmp_path / "obs.dat"
    obs.write_text(HEADER + "1 S 0 0 0 0 0 ZXY 2 2 1\n")
    for name in ["input.rho", "control.inv", "cov.cov"]:
        (tmp_path / name).write_text("Exit search when rms is less than: 1.05\nMaximum number of iterations: 10\n")
    def run(self, *args, **kwargs):
        assert kwargs["load_result"] is False
        if write_outputs:
            model = ModEmModel3D()
            model.x_widths = model.y_widths = model.z_widths = np.array([100.])
            model.rho_loge = np.log(np.full((1,1,1),100.))
            model.write(tmp_path / "run_001.rho")
            (tmp_path / "run_001.dat").write_text(HEADER + "1 S 0 0 0 0 0 ZXY 1 1 1e15")
            (tmp_path / "run.log").write_text("Completed NLCG iteration 1\nwith: f= 1 m2= 1 rms= 1 lambda= 1 alpha= 1")
    monkeypatch.setattr(modem.ModEmRunner, "run", run)
    cfg = InversionConfig(method="mt", dimension="3d", backend="modem", workdir=str(tmp_path),
        run_external=True, backend_options={
        "files": {"model": "input.rho", "data": "obs.dat", "control": "control.inv", "covariance": "cov.cov"},
        "runner": {"outputs": {"model": "run_001.rho", "response": "run_001.dat", "log": "run.log"}}})
    ctrl = ModEmControl(); ctrl.target_rms=1.05; ctrl.max_iterations=10
    ctrl.write(tmp_path / "control.inv")
    initial = ModEmModel3D()
    initial.x_widths=initial.y_widths=initial.z_widths=np.array([100.])
    initial.rho_loge=np.log(np.full((1,1,1),100.))
    initial.write(tmp_path / "input.rho")
    result = ModEMBackend(cfg).run()
    assert result.metadata["executed"] is True
    assert result.metadata["converged"] is False
    if write_outputs:
        assert result.status == "loaded"
        assert result.metadata["output_provenance"] == "validated_triplet"
        assert result.rms == pytest.approx(1.)
    else:
        assert result.status == "needs_review"
        assert np.isnan(result.rms)


@pytest.mark.parametrize("change", ["units", "sign", "rotation", "origin", "counts", "missing"])
def test_header_mismatches_rejected(contract, change):
    paths, protected, obs = contract
    text = paths["response"].read_text()
    replacements = {"units": ("[mV/km]/[nT]", "ohm"), "sign": ("+i", "-i"),
                    "rotation": ("> 0\n", "> 30\n"), "origin": ("> 0 0", "> 1 0"),
                    "counts": ("> 1 1", "> 2 1"), "missing": (HEADER, "")}
    a,b = replacements[change]
    paths["response"].write_text(text.replace(a,b))
    with pytest.raises(ValueError): validate_outputs(paths, protected, obs)


@pytest.mark.parametrize("rms,iteration,end,local,expected", [
    (1.,1,True,False,"target_rms_reached"),
    (1.05,10,True,False,"threshold_ambiguous"),
    (2.,10,True,False,"iteration_limit"),
    (2.,2,False,True,"local_minimum"),
    (1.,1,False,False,"completion_unverified"),
    (2.,2,True,False,"other_stop")])
def test_native_stop_classification(tmp_path,rms,iteration,end,local,expected):
    from pycsamt.inversion.backends._modem_outputs import classify_stop
    ctrl=tmp_path/"control.inv"
    c=ModEmControl(); c.target_rms=1.05; c.max_iterations=10; c.write(ctrl)
    log=tmp_path/"run.log"
    log.write_text((f"NLCG iterations: {iteration} function evaluations: 20" if end else "")+
                   ("\nUnable to get out of a local minimum. Exiting..." if local else ""))
    assert classify_stop(log,ctrl,iteration,rms)["reason"] == expected


def test_native_title_continuation_and_origin_print_precision(contract):
    from pycsamt.inversion.backends._modem_outputs import compatible_headers
    paths, protected, obs = contract
    obs.write_text(obs.read_text().replace("> 0 0", "> -31.95556 41.53481"))
    paths["response"].write_text("# ModEM impedance responses for\n ModEM impedance responses for title\n"+
        paths["response"].read_text().replace("> 0 0", "> -31.956 41.535"))
    assert compatible_headers(obs, paths["response"])


@pytest.mark.parametrize("damage", ["order", "width", "decimal", "limit"])
def test_native_control_layout_rejected(tmp_path, damage):
    from pycsamt.inversion.backends._modem_outputs import native_control
    path=tmp_path/"control.inv"; ModEmControl().write(path)
    lines=path.read_text().splitlines()
    if damage == "order": lines[1],lines[2]=lines[2],lines[1]
    if damage == "width": lines[1]="    "+lines[1]
    if damage == "decimal": lines[1]=lines[1][:36]+"10"
    if damage == "limit": lines[-1]=lines[-1][:36]+"10000"
    path.write_text("\n".join(lines))
    with pytest.raises(ValueError): native_control(path)


@pytest.mark.parametrize("changed", [True, False])
def test_output_model_geometry_matches_initial(contract,tmp_path,changed):
    paths,protected,obs=contract
    initial=ModEmModel3D.read(paths["model"])
    if changed: initial.x_widths *= 2
    initial.write(tmp_path/"initial.ws")
    if changed:
        with pytest.raises(ValueError,match="geometry"):
            validate_outputs(paths,protected,obs,initial_model=tmp_path/"initial.ws")
    else:
        assert validate_outputs(paths,protected,obs,initial_model=tmp_path/"initial.ws")["iteration"] == 1


@pytest.mark.parametrize("status", ["loaded", "prepared", "success", "executed"])
def test_explicit_nonconvergence_overrides_legacy_status(status):
    from pycsamt.inversion.results import InversionResult
    result=InversionResult("mt","3d","modem",status=status,metadata={"converged":False})
    assert result.converged is False


@pytest.mark.parametrize("text,iteration,rms", [
    ("NLCG iterations: 1 function evaluations: 2\nNLCG iterations: 1 function evaluations: 2",1,1.),
    ("NLCG iterations: 1 function evaluations: 2\nStarting line search...",1,1.),
    ("NLCG iterations: 2 function evaluations: 2",1,1.),
    ("NLCG iterations: 1 function evaluations: 0",1,1.),
    ("Unable to get out of a local minimum\nNLCG iterations: 1 function evaluations: 2",1,1.),
    ("NLCG iterations: 11 function evaluations: 20",11,1.),
    ("NLCG iterations: 1 function evaluations: 2",1,float("nan"))])
def test_contradictory_terminal_evidence_rejected(tmp_path,text,iteration,rms):
    from pycsamt.inversion.backends._modem_outputs import classify_stop
    ctrl=tmp_path/"control.inv"; c=ModEmControl();c.max_iterations=10;c.write(ctrl)
    log=tmp_path/"run.log";log.write_text(text)
    with pytest.raises(ValueError):classify_stop(log,ctrl,iteration,rms)


@pytest.mark.parametrize("rms", [1.0499999,1.05,1.0500001])
def test_print_threshold_neighbours_not_promoted(tmp_path,rms):
    from pycsamt.inversion.backends._modem_outputs import classify_stop
    ctrl=tmp_path/"control.inv";c=ModEmControl();c.target_rms=1.05;c.write(ctrl)
    log=tmp_path/"run.log";log.write_text("NLCG iterations: 1 function evaluations: 2")
    assert classify_stop(log,ctrl,1,rms)["reason"]=="threshold_ambiguous"
