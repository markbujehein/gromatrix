"""Step G: chaining, resume from checkpoints, gates and the Slurm script, against a fake gmx."""

import json
from pathlib import Path

import pytest

from gromatrix import mdp as mdpmod
from gromatrix import run as runmod
from gromatrix.forcefields import acquire_lock, resolve

FAKE_GMX = str(Path(__file__).with_name("fake_gmx.py"))


def _system(b, sid="S1", sha="aaa", smoke=None):
    wd = b / "systems" / sid
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
    assert runmod.phase_state(wd, "npt") == "partial"

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


def test_rebuilt_system_restarts_from_em(batch):
    wd = _system(batch)
    assert _run(batch) == 0
    stamp = json.loads((wd / "system.json").read_text())
    (wd / "system.json").write_text(json.dumps({**stamp, "inputs_sha": "bbb"}))
    (wd / "gmx_calls.jsonl").unlink()
    assert _run(batch) == 0
    assert [c[c.index("-deffnm") + 1] for c in _calls(wd) if c[0] == "mdrun"] == ["em", "nvt", "npt", "md"]
    assert json.loads((wd / "run.json").read_text())["inputs_sha"] == "bbb"


def test_gates_em_convergence_and_lincs(batch, monkeypatch, capsys):
    wd = _system(batch, "S1")
    _system(batch, "S2")
    assert _run(batch, env={"FAKE_GMX_EM_FMAX": "5000"}, monkeypatch=monkeypatch) == 1
    out = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert {f["id"] for f in out["failed"]} == {"S1", "S2"}  # one failure does not stop the batch
    assert "EM did not converge" in json.loads((wd / "run.json").read_text())["reason"]
    assert runmod.phase_state(wd, "em") != "done" and (wd / "em.failed.gro").exists()

    for s in ("S1", "S2"):
        for f in (batch / "systems" / s).glob("em.*"):
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
    assert "--signal=B:TERM@300" in script


def test_log_helpers(tmp_path):
    log = tmp_path / "md.log"
    log.write_text(
        "Log file opened\n           Step           Time\n          100       0.2\n\nReceived the TERM signal\n"
        "Finished mdrun\nLog file opened\n           Step           Time\n          400       0.8\n\nFinished mdrun\n"
    )
    seg = runmod.last_segment(log)
    assert runmod.last_step(seg) == 400 and "TERM" not in seg  # only the newest segment counts
    (tmp_path / "md.gro").write_text("x")
    assert runmod.phase_state(tmp_path, "md") == "done"
    log.write_text(log.read_text() + "Log file opened\nStep 9: Run time exceeded 0.990 hours\nFinished mdrun\n")
    (tmp_path / "md.cpt").write_text("x")
    assert runmod.phase_state(tmp_path, "md") == "partial"  # -maxh stop is not completion
    assert runmod._hours("2-01:30:00") == pytest.approx(49.5) and runmod._hours("90") == 1.5


def test_progress_lines_while_mdrun_runs(batch, monkeypatch, capsys):
    _system(batch)
    assert _run(batch, env={"FAKE_GMX_SLOW": "md"}, monkeypatch=monkeypatch) == 0
    lines = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    running = [x for x in lines if x.get("status") == "running" and x["phase"] == "md"]
    r = next(x for x in running if x["progress_pct"] == 50.0)  # polled while mdrun sat at nsteps/2
    assert r["step"] == "run" and r["system"] == "S1"
    # overall: nvt + npt (50 000 steps each) done, half of the 10 000-step production
    assert r["overall_pct"] == pytest.approx(100 * 105000 / 110000, abs=0.1)
