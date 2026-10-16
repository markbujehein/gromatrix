#!/usr/bin/env python3
"""A stand-in for `gmx` in tests: writes the files real GROMACS writes (tpr,
log with Step/Time blocks, gro, cpt, xtc) so step G's chaining, resume and gates
run without GROMACS. Every call is appended to gmx_calls.jsonl in the cwd.

Fault injection (environment):
  FAKE_GMX_KILL=<name>   mdrun of <name> writes half its steps + a checkpoint, then exits 1
  FAKE_GMX_EM_FMAX=<x>   EM ends with this Fmax (default 500)
  FAKE_GMX_LINCS=<name>  mdrun of <name> logs a LINCS warning
  FAKE_GMX_SLOW=<name>   mdrun of <name> logs step nsteps/2, pauses 0.5 s, then finishes
  FAKE_GMX_FAIL=<name>   mdrun of <name> exits 1 before writing anything (bad flag, GPU error)
  FAKE_GMX_MAXH=<name>   mdrun of <name> stops at nsteps/2 like -maxh: checkpoint, exit 0
"""

import json
import os
import sys
import time
from pathlib import Path


def opt(args, flag, default=None):
    return args[args.index(flag) + 1] if flag in args else default


def main():
    args = [a for a in sys.argv[1:] if a not in ("-quiet", "-nobackup")]
    with open("gmx_calls.jsonl", "a") as fh:
        fh.write(json.dumps(args) + "\n")
    sub, args = args[0], args[1:]
    if sub == "grompp":
        mdp = {}
        for line in Path(opt(args, "-f")).read_text().splitlines():
            if "=" in line and not line.startswith(";"):
                k, v = line.split("=", 1)
                mdp[k.strip()] = v.strip()
        for f in (opt(args, "-c"), opt(args, "-r"), opt(args, "-t")):
            if f and not Path(f).exists():
                print(f"Fatal error: File '{f}' does not exist", file=sys.stderr)
                return 1
        tpr = {"nsteps": int(float(mdp.get("nsteps", 0))), "c": opt(args, "-c"), "em": mdp.get("integrator") == "steep"}
        Path(opt(args, "-o")).write_text(json.dumps(tpr))
        return 0
    if sub == "mdrun":
        name = opt(args, "-deffnm")
        if os.environ.get("FAKE_GMX_FAIL") == name:
            print("Fatal error: GPU not available", file=sys.stderr)
            return 1
        tpr = json.loads(Path(f"{name}.tpr").read_text())
        n = tpr["nsteps"]
        start = 0
        if "-cpi" in args and Path(f"{name}.cpt").exists():
            start = json.loads(Path(f"{name}.cpt").read_text())["step"]
        kill = os.environ.get("FAKE_GMX_KILL") == name
        maxh = os.environ.get("FAKE_GMX_MAXH") == name
        stop = n // 2 if kill or maxh else n
        lines = ["Log file opened on today", f"Started mdrun from step {start}"]
        for s in sorted({start, (start + stop) // 2, stop}):
            lines += ["           Step           Time", f"          {s}       {s * 0.002:.5f}", ""]
        if os.environ.get("FAKE_GMX_LINCS") == name:
            lines.append("Step 10, time 0.02 (ps)  LINCS WARNING")
        Path(f"{name}.cpt").write_text(json.dumps({"step": stop}))
        if os.environ.get("FAKE_GMX_SLOW") == name:
            block = f"           Step           Time\n          {n // 2}    0.1\n"
            Path(f"{name}.log").write_text("Log file opened\n" + block)
            time.sleep(0.5)
        if kill:
            with open(f"{name}.log", "a") as fh:
                fh.write("\n".join(lines) + "\n")
            return 1
        if tpr["em"]:
            fmax = float(os.environ.get("FAKE_GMX_EM_FMAX", 500))
            lines += ["Potential Energy  = -1.0e+05", f"Maximum force     =  {fmax:.5e} on atom 12"]
        if maxh:
            lines.append(f"Step {stop}: Run time exceeded 0.990 hours, will terminate the run")
        lines.append("Finished mdrun on rank 0")
        with open(f"{name}.log", "a") as fh:
            fh.write("\n".join(lines) + "\n")
        Path(f"{name}.gro").write_text(Path(tpr["c"]).read_text())  # mdrun writes confout on -maxh too
        Path(f"{name}.xtc").write_text("frames")
        return 0
    print(f"fake gmx: unsupported subcommand {sub}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
