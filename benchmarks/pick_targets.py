"""Batch driver for the PXR/CAR panel: choose a target ligand per entry, run step C.

The target is the largest non-ion, non-water, non-additive HET (by heavy atoms)
in the fixed structure. For PXR/CAR deposits that is the drug/agonist.
"""

import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

from gromatrix.fix import ADDITIVES
from gromatrix.pdbio import ION_RES, WATER_RES, read_pdb
from gromatrix.structure import POLYMER_RES

B = Path(sys.argv[1])
targets, info = [], {}
for pdb in sorted((B / "fixed").glob("*.pdb")):
    pid = pdb.stem
    sizes = Counter()
    seen = set()
    for a in read_pdb(pdb):
        if a.resname in POLYMER_RES or a.resname in WATER_RES or a.resname in ION_RES or a.resname in ADDITIVES:
            continue
        sizes[a.resname] += 1
        seen.add((a.chain, a.resname, a.resseq))
    copies = Counter(r for _, r, _ in seen)
    per_copy = {r: sizes[r] / copies[r] for r in sizes}
    meta = json.loads(pdb.with_suffix(".meta.json").read_text())
    ents = [e["description"] for e in meta["entities"].values() if e["type"] == "polymer"]
    receptor = "CAR" if any("NR1I3" in d.upper() or "ANDROSTANE" in d.upper() for d in ents) else "PXR"
    if any("NR1I2" in d.upper() or "PREGNANE" in d.upper() for d in ents) and receptor == "CAR":
        receptor = "PXR+CAR"
    t = max(per_copy, key=per_copy.get) if per_copy else None
    info[pid] = {"receptor": receptor, "het": {r: int(per_copy[r]) for r in per_copy}, "target": t, "entities": ents}
    if t:
        targets.append(f"{t}@{pid}")
(B / "panel_targets.json").write_text(json.dumps(info, indent=1))
print(
    json.dumps(
        {"entries": len(info), "with_target": len(targets), "receptors": Counter(v["receptor"] for v in info.values())}
    )
)
r = subprocess.run(
    [sys.executable, "-m", "gromatrix.extract", "--batch-dir", str(B), "--targets", ",".join(targets)],
    capture_output=True,
    text=True,
)
(B / "extract.json").write_text(r.stdout)
print("extract exit", r.returncode)
