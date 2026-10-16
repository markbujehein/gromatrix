"""MDP generation. Non-bonded settings come only from the locked ForceField so
every system in a matrix uses identical, forcefield-correct cut-off schemes."""

from __future__ import annotations

from pathlib import Path

from .forcefields import ForceField


def _fmt(d: dict) -> str:
    return "\n".join(f"{k:<24s} = {v}" for k, v in d.items()) + "\n"


def protocols(ff: ForceField, temp: float = 300.0, md_ns: float = 100.0) -> dict[str, dict]:
    nb = dict(ff.nonbonded, nstlist=10, pbc="xyz")
    em = {"integrator": "steep", "emtol": 1000.0, "emstep": 0.01, "nsteps": 50000, **nb}
    dyn = {
        "integrator": "md",
        "dt": 0.002,
        "constraints": "h-bonds",
        "constraint-algorithm": "lincs",
        "tcoupl": "V-rescale",
        "tc-grps": "System",
        "tau_t": 0.1,
        "ref_t": temp,
        "nstxout-compressed": 5000,
        "nstenergy": 5000,
        "nstlog": 5000,
        **nb,
    }
    nvt = {
        "define": "-DPOSRES",
        "nsteps": 50000,
        "gen_vel": "yes",
        "gen_temp": temp,
        "gen_seed": -1,
        "pcoupl": "no",
        **dyn,
    }
    npt = {
        "define": "-DPOSRES",
        "nsteps": 50000,
        "continuation": "yes",
        "gen_vel": "no",
        "pcoupl": "C-rescale",
        "pcoupltype": "isotropic",
        "tau_p": 2.0,
        "ref_p": 1.0,
        "compressibility": 4.5e-5,
        "refcoord-scaling": "com",
        **dyn,
    }
    md = {
        **{k: v for k, v in npt.items() if k not in ("define", "refcoord-scaling")},
        "nsteps": int(md_ns * 1000 / 0.002),
        "tau_p": 5.0,
    }
    ions = {**em, "nsteps": 0}
    return {"ions": ions, "em": em, "nvt": nvt, "npt": npt, "md": md}


def write_all(ff: ForceField, outdir: Path, temp: float, md_ns: float) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for name, params in protocols(ff, temp, md_ns).items():
        (outdir / f"{name}.mdp").write_text(f"; {name} - forcefield {ff.key}\n" + _fmt(params))


def emtol(mdp_dir: Path) -> float:
    return float(read_mdp(Path(mdp_dir) / "em.mdp").get("emtol", 1000))


def read_mdp(path: Path) -> dict[str, str]:
    out = {}
    for line in Path(path).read_text().splitlines():
        line = line.split(";")[0]
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out
