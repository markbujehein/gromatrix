"""Where and how the MD of a batch runs: the one switch between local and remote compute.

    gromatrix-compute --batch-dir B detect
    gromatrix-compute --batch-dir B set --mode local --threads 8 [--gpu-id 0] [--maxh 20]
    gromatrix-compute --batch-dir B set --mode slurm --time 24:00:00 --partition gpu --gpus-per-task 1
    gromatrix-compute --batch-dir B set --mode remote --time 48:00:00
    gromatrix-compute --batch-dir B show
    gromatrix-compute --batch-dir B run-args        # the exact flags for `gromatrix-run`

Everything up to step F (structure prep, validation, GROMACS builds, smoke test) needs only CPU
and GROMACS. Only the production MD (G) is expensive, and it can go three ways:

  local    detached on this machine (`gromatrix-run --background`); CPU-only is fine for short runs
  slurm    write the Slurm array script here and submit it from this machine (a login node)
  remote   write the Slurm array script, then copy the batch to another machine and submit there;
           copy back 06-systems/ afterwards and run step H here

The choice is stored in 00-config/compute.json so the MD agents and `gromatrix-status` read it
instead of asking again. Prints one JSON line {"ok": ..., ...}.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

from . import layout, result

MODES = ("local", "slurm", "remote")
DEFAULTS = {
    "mode": "local",
    "threads": None,
    "gpu_id": None,
    "maxh": None,
    "time": "24:00:00",
    "partition": None,
    "account": None,
    "gpus_per_task": None,
    "max_parallel": None,
}


def load(batch_dir: Path) -> dict | None:
    try:
        return json.loads((Path(batch_dir) / layout.COMPUTE).read_text())
    except (OSError, ValueError):
        return None


def detect() -> dict:
    """What this machine offers; never fails (a missing tool is just false)."""
    return {
        "cpus": os.cpu_count(),
        "gpu": shutil.which("nvidia-smi") is not None,
        "slurm": shutil.which("sbatch") is not None,
        "gmx": shutil.which("gmx") is not None or shutil.which("gmx_mpi") is not None,
    }


def run_args(cfg: dict) -> list[str]:
    """Flags for `gromatrix-run` that realise `cfg` (CPU-only when gpu_id is unset)."""
    a: list[str] = []
    if cfg["mode"] == "local":
        a += ["--threads", str(cfg.get("threads") or os.cpu_count() or 1)]
        if cfg.get("gpu_id") is not None:
            a += ["--gpu-id", str(cfg["gpu_id"])]
        if cfg.get("maxh"):
            a += ["--maxh", str(cfg["maxh"])]
        return a + ["--background"]
    a += ["--slurm", "--time", cfg.get("time") or DEFAULTS["time"]]
    for flag, key in (
        ("--partition", "partition"),
        ("--account", "account"),
        ("--gpus-per-task", "gpus_per_task"),
        ("--max-parallel", "max_parallel"),
    ):
        if cfg.get(key) is not None:
            a += [flag, str(cfg[key])]
    return a


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gromatrix-compute", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("detect", "show", "run-args"):
        sub.add_parser(name)
    s = sub.add_parser("set")
    s.add_argument("--mode", choices=MODES, required=True)
    s.add_argument("--threads", type=int)
    s.add_argument("--gpu-id", type=int, help="omit for CPU-only")
    s.add_argument("--maxh", type=float, help="stop cleanly with a checkpoint after this many hours")
    s.add_argument("--time", help="Slurm walltime, HH:MM:SS")
    s.add_argument("--partition")
    s.add_argument("--account")
    s.add_argument("--gpus-per-task", type=int)
    s.add_argument("--max-parallel", type=int)
    a = ap.parse_args(argv)

    if a.cmd == "detect":
        return result.finish({"ok": True, **detect()})
    cfg = load(a.batch_dir)
    if a.cmd == "set":
        cfg = {**DEFAULTS, **{k: v for k, v in vars(a).items() if k in DEFAULTS and v is not None}}
        path = a.batch_dir / layout.COMPUTE
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cfg, indent=2) + "\n")
    if cfg is None:
        return result.finish(
            result.fail("compute_unset", f"no {layout.COMPUTE}; run `gromatrix-compute set --mode ...`")
        )
    out = {"ok": True, "compute": cfg}
    if a.cmd == "run-args":
        out["args"] = run_args(cfg)
        out["submit"] = cfg["mode"] == "slurm"  # remote: copy the batch and submit on the other machine
    return result.finish(out)


if __name__ == "__main__":
    sys.exit(main())
