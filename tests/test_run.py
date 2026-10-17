"""Step G: chaining, resume from checkpoints, gates and the Slurm script, against a fake gmx."""

import json
from pathlib import Path

import pytest

from gromatrix import checks
from gromatrix import mdp as mdpmod
from gromatrix import run as runmod
from gromatrix.forcefields import acquire_lock, resolve

FAKE_GMX = str(Path(__file__).with_name("fake_gmx.py"))


def _system(b, sid="S1", sha="aaa", smoke=None):
    wd = b / "06-systems" / sid
    wd.mkdir(parents=True)
    stamp = {"id": sid, "status": "ok", "inputs_sha": sha, "role": "complex"}
    if smoke is not None:
        stamp["smoke"] = {"ok": smoke}
    (wd / "system.json").write_text(json.dumps(stamp))
    mdpmod.write_all(resolve("amber"), wd / "mdp", 300.0, 0.02)  # 10 000 production steps
    for f in ("solv_ions.gro", "topol.top", "index.ndx"):
        (wd / f).write_text(f"{f}\n")
    return wd


def _run(b, *extra, env=None, monkeypatch=None):
    for k, v in (env or {}).items():
        monkeypatch.setenv(k, v)
    rc = runmod.main(["--batch-dir", str(b), "--gmx", FAKE_GMX, "--poll", "0.01", *extra])
    for k in env or {}:
        monkeypatch.delenv(k)
    return rc


def _calls(wd):
    return [json.loads(x) for x in (wd / "gmx_calls.jsonl").read_text().splitlines()]


@pytest.fixture
def batch(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    acquire_lock(tmp_path, "amber")
    return tmp_path


def test_full_protocol_chains_coordinates_and_checkpoints(batch, capsys):
    wd = _system(batch)
    assert _run(batch) == 0
    lines = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    assert lines[-1]["ok"] and lines[-1]["done"] == ["S1"]
    assert [x["phase"] for x in lines if x.get("status") == "done"] == ["em", "nvt", "npt", "md"]
    assert lines[-2]["overall_pct"] == 100.0

    grompp = {c[c.index("-o") + 1]: c for c in _calls(wd) if c[0] == "grompp"}
    assert ["-c", "em.gro"] == grompp["nvt.tpr"][3:5] and "-r" in grompp["nvt.tpr"]
    npt = grompp["npt.tpr"]
    assert npt[npt.index("-t") + 1] == "nvt.cpt" and npt[npt.index("-r") + 1] == "nvt.gro"
    md = grompp["md.tpr"]
    assert md[md.index("-t") + 1] == "npt.cpt" and "-r" not in md  # production is unrestrained

    st = json.loads((wd / "run.json").read_text())
    assert st["status"] == "done" and st["inputs_sha"] == "aaa" and st["trajectory"] == "md.xtc"
    assert set(st["phases"]) == {"em", "nvt", "npt", "md"}

    # A second run has nothing to do: every phase is complete.
    capsys.readouterr()
    assert _run(batch) == 0
    lines = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    skipped = [x for x in lines if x.get("status") == "skipped"]
    assert [x["phase"] for x in skipped] == ["em", "nvt", "npt", "md"]


def test_crash_resumes_from_checkpoint(batch, monkeypatch):
    wd = _system(batch)
    assert _run(batch, env={"FAKE_GMX_KILL": "npt"}, monkeypatch=monkeypatch) == 1
    st = json.loads((wd / "run.json").read_text())
    assert st["status"] == "failed" and st["phases"]["npt"]["status"] == "failed"
    assert (wd / "npt.cpt").exists() and not (wd / "npt.gro").exists()
    assert runmod.phase_state(wd, runmod.PHASES[2], "npt") == "partial"

    (wd / "gmx_calls.jsonl").unlink()
    assert _run(batch) == 0
    calls = _calls(wd)
    # em/nvt untouched; npt continued with -cpi and NOT re-gromped; md fresh.
    assert not any(c[0] == "grompp" and c[-3] in ("em.tpr", "nvt.tpr", "npt.tpr") for c in calls)
    npt_run = next(c for c in calls if c[0] == "mdrun" and "npt" in c)
    assert npt_run[npt_run.index("-cpi") + 1] == "npt.cpt"
    assert any(c[0] == "grompp" and "md.tpr" in c for c in calls)
    assert json.loads((wd / "run.json").read_text())["status"] == "done"


def test_rerunning_a_phase_restarts_everything_after_it(batch, monkeypatch):
    wd = _system(batch)
    assert _run(batch) == 0
    (wd / "nvt.gro").unlink()  # NVT output lost: NPT and production depend on it
    (wd / "nvt.cpt").unlink()
    (wd / "gmx_calls.jsonl").unlink()
    assert _run(batch) == 0
    runs = [c[c.index("-deffnm") + 1] for c in _calls(wd) if c[0] == "mdrun"]
    assert runs == ["nvt", "npt", "md"] and not any("-cpi" in c for c in _calls(wd))


def test_gates_em_convergence_and_lincs(batch, monkeypatch, capsys):
    wd = _system(batch, "S1")
    _system(batch, "S2")
    assert _run(batch, env={"FAKE_GMX_EM_FMAX": "5000"}, monkeypatch=monkeypatch) == 1
    out = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert {f["id"] for f in out["failed"]} == {"S1", "S2"}  # one failure does not stop the batch
    assert "EM did not converge" in json.loads((wd / "run.json").read_text())["reason"]
    assert runmod.phase_state(wd, runmod.PHASES[0], "em") == "new" and (wd / "em.failed.gro").exists()

    for s in ("S1", "S2"):
        for f in (batch / "06-systems" / s).glob("em.*"):
            f.unlink()
    assert _run(batch, "--only", "S1", env={"FAKE_GMX_LINCS": "nvt"}, monkeypatch=monkeypatch) == 1
    assert "LINCS" in json.loads((wd / "run.json").read_text())["reason"]
    # LINCS warnings in production are reported, not fatal.
    for f in wd.glob("nvt.*"):
        f.unlink()
    assert _run(batch, "--only", "S1", env={"FAKE_GMX_LINCS": "md"}, monkeypatch=monkeypatch) == 0
    assert json.loads((wd / "run.json").read_text())["warnings"] == ["md: 1 LINCS warnings"]


def test_custom_deffnm_and_mdrun_flags(batch):
    wd = _system(batch)
    assert _run(batch, "--deffnm", "prod", "--threads", "4", "--gpu-id", "0", "--maxh", "1.5") == 0
    md = next(c for c in _calls(wd) if c[0] == "mdrun" and "prod" in c)
    assert md[md.index("-nt") + 1] == "4" and md[md.index("-gpu_id") + 1] == "0" and "-maxh" in md
    assert (wd / "prod.xtc").exists() and (wd / "nvt.gro").exists()


def test_selection_skips_failed_smoke_and_unknown_ids(batch, capsys):
    _system(batch, "OK1", smoke=True)
    _system(batch, "BAD", smoke=False)
    assert _run(batch, "--only", "OK1", "BAD", "NOPE") == 0  # skipped systems are reported, not run
    out = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert out["done"] == ["OK1"]
    assert {s["id"]: s["reason"][:5] for s in out["skipped"]} == {"BAD": "smoke", "NOPE": "no su"}


def test_slurm_array_script(batch, capsys):
    _system(batch, "A1")
    _system(batch, "A2")
    assert (
        runmod.main(
            ["--batch-dir", str(batch), "--slurm", "--time", "1-00:00:00", "--threads", "16", "--partition", "gpu",
             "--gpus-per-task", "1", "--prelude", "module load gromacs", "--max-parallel", "4"]
        )
        == 0
    )  # fmt: skip
    out = json.loads(capsys.readouterr().out)
    script = Path(out["script"]).read_text()
    assert out["systems"] == ["A1", "A2"] and out["next"].startswith("sbatch ")
    assert "#SBATCH --array=0-1%4" in script and "#SBATCH --cpus-per-task=16" in script
    assert "#SBATCH --partition=gpu" in script and "#SBATCH --gpus-per-task=1" in script
    assert "module load gromacs" in script and "SYSTEMS=(A1 A2)" in script
    assert '--only "$ID"' in script and "--maxh 22.633" in script  # 24 h * 0.95 - 10 min
    assert "--signal=B:TERM@300" in script and "--force" not in script
    runmod.main(["--batch-dir", str(batch), "--slurm", "--force"])
    assert "--force" in Path(json.loads(capsys.readouterr().out)["script"]).read_text()


def test_log_helpers(tmp_path):
    log = tmp_path / "md.log"
    log.write_text(
        "Log file opened\n           Step           Time\n          100       0.2\n\nReceived the TERM signal\n"
        "Finished mdrun\nLog file opened\n           Step           Time\n          400       0.8\n\nFinished mdrun\n"
    )
    seg = checks.last_segment(log.read_text())
    assert runmod.last_step(seg) == 400 and "TERM" not in seg  # only the newest segment counts
    md = runmod.PHASES[3]
    (tmp_path / "md.gro").write_text("x")
    assert runmod.phase_state(tmp_path, md, "md") == "done"
    log.write_text(log.read_text() + "Log file opened\nStep 9: Run time exceeded 0.990 hours\nFinished mdrun\n")
    (tmp_path / "md.cpt").write_text("x")
    assert runmod.phase_state(tmp_path, md, "md") == "partial"  # -maxh stop is not completion
    # A phase whose checkpoint the next phase reads is not done without it.
    (tmp_path / "nvt.log").write_text("Log file opened\nFinished mdrun\n")
    (tmp_path / "nvt.gro").write_text("x")
    assert runmod.phase_state(tmp_path, runmod.PHASES[1], "nvt") == "new"


@pytest.mark.parametrize(
    "walltime,hours",
    [("90", 1.5), ("30:00", 0.5), ("24:00:00", 24), ("1-12", 36), ("1-12:30", 36.5), ("2-01:30:00", 49.5)],
)
def test_slurm_walltime_formats(walltime, hours):
    assert runmod._hours(walltime) == pytest.approx(hours)


def test_failed_gate_is_rerun_fresh_not_resumed(batch, monkeypatch):
    wd = _system(batch)
    assert _run(batch, env={"FAKE_GMX_LINCS": "nvt"}, monkeypatch=monkeypatch) == 1
    assert (wd / "nvt.failed.log").exists() and not (wd / "nvt.cpt").exists()
    (wd / "gmx_calls.jsonl").unlink()
    assert _run(batch) == 0
    calls = _calls(wd)
    nvt = next(c for c in calls if c[0] == "mdrun" and "nvt" in c)
    assert "-cpi" not in nvt and any(c[0] == "grompp" and "nvt.tpr" in c for c in calls)


def test_em_left_by_smoke_test_is_gated(batch, monkeypatch):
    wd = _system(batch)
    # step F's smoke test ran EM in the system directory and it did not converge
    (wd / "em.gro").write_text("x")
    (wd / "em.log").write_text(
        "Log file opened\nPotential Energy  = -1.0e+05\nMaximum force = 5.0e+04\nFinished mdrun\n"
    )
    assert _run(batch) == 1
    st = json.loads((wd / "run.json").read_text())
    assert "EM did not converge" in st["reason"] and not (wd / "nvt.tpr").exists()


def test_crash_is_failed_even_after_an_earlier_signal_stop(batch, monkeypatch):
    wd = _system(batch)
    assert _run(batch, env={"FAKE_GMX_MAXH": "md"}, monkeypatch=monkeypatch) == 1
    st = json.loads((wd / "run.json").read_text())
    assert st["status"] == "interrupted" and runmod.phase_state(wd, runmod.PHASES[3], "md") == "partial"
    # Resume fails before mdrun writes anything; the old "Run time exceeded" must not mask it.
    assert _run(batch, env={"FAKE_GMX_FAIL": "md"}, monkeypatch=monkeypatch) == 1
    assert json.loads((wd / "run.json").read_text())["status"] == "failed"
    assert _run(batch) == 0  # and a clean resume finishes it


def test_maxh_is_one_budget_for_the_invocation(batch, monkeypatch, capsys):
    _system(batch, "S1")
    _system(batch, "S2")
    clock = iter(range(0, 10**6, 600))  # every mdrun "takes" 10 minutes
    monkeypatch.setattr(runmod.time, "monotonic", lambda: float(next(clock)))
    assert _run(batch, "--maxh", "0.5") == 1  # em, nvt, npt fit; production does not
    out = json.loads(capsys.readouterr().out.splitlines()[-1])
    maxh = [float(c[c.index("-maxh") + 1]) for c in _calls(batch / "06-systems" / "S1") if c[0] == "mdrun"]
    assert maxh == sorted(maxh, reverse=True) and maxh[0] < 0.5  # what is left, not the full hour
    assert out["interrupted"] == ["S1"] and out["not_started"] == ["S2"]  # budget ends the batch


def test_signal_between_phases_stops_the_batch(batch, monkeypatch):
    wd = _system(batch)
    real = runmod.Runner.grompp

    def grompp_then_term(self, phase, prev):
        real(self, phase, prev)
        if phase.name == "nvt":
            self.request_stop()  # SIGTERM arrives while no mdrun is running

    monkeypatch.setattr(runmod.Runner, "grompp", grompp_then_term)
    assert _run(batch) == 1
    st = json.loads((wd / "run.json").read_text())
    assert st["status"] == "interrupted" and st["phases"]["em"]["status"] == "done"
    assert not any(c[0] == "mdrun" and "nvt" in c for c in _calls(wd))


def test_corrupt_run_json_does_not_block(batch):
    wd = _system(batch)
    (wd / "run.json").write_text('{"status": "runn')  # killed mid-write
    assert _run(batch) == 0


def test_abbreviated_flags_are_rejected():
    with pytest.raises(SystemExit):
        runmod.main(["--batch-dir", "x", "--backg"])  # would otherwise re-spawn itself in background
