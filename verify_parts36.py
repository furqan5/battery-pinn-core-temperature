"""Re-execute Parts 3 and 6 in THIS repository, to upgrade their numbers from
'transcribed from another machine' to 'verified here'.

The notebooks discover their dataset by recursive search from '.', so they are
copied to the repo root and executed with cwd = repo root, where
Data_Sets\\Li-ion Battery Dataset from NASA PCoE\\...\\B0005.mat is reachable.

allow_errors=True deliberately: a cell that fails is information (it tells us
which numbers are re-verifiable and which are not), and stopping at the first
failure would hide the rest.
"""
import os
import shutil
import sys
import time

import nbformat
from nbclient import NotebookClient

REPO = r"C:\Users\Nouman\Desktop\Furqan's Docs\battery_pinn"
SRC = os.path.join(REPO, "CoreField_Battery_Package_v1.0",
                   "CoreField_Battery_Package_v1.0", "notebooks")
OUT = os.path.join(REPO, "verify")
os.makedirs(OUT, exist_ok=True)

targets = sys.argv[1:] or ["Part3_CRLB_Classical_Inverse.ipynb",
                           "Part6_Real_Data_PINN.ipynb"]

for name in targets:
    src = os.path.join(SRC, name)
    tmp = os.path.join(REPO, "_run_" + name)
    shutil.copy2(src, tmp)
    nb = nbformat.read(tmp, as_version=4)
    ncode = sum(1 for c in nb.cells if c.cell_type == "code")
    print(f"\n{'='*90}\nEXECUTING {name}: {len(nb.cells)} cells ({ncode} code)\n{'='*90}",
          flush=True)

    client = NotebookClient(nb, timeout=3000, kernel_name="python3",
                            allow_errors=True,
                            resources={"metadata": {"path": REPO}})
    t0 = time.time()
    try:
        client.execute()
    except Exception as e:
        print(f"  client-level failure: {type(e).__name__}: {e}", flush=True)

    dst = os.path.join(OUT, name.replace(".ipynb", "_executed.ipynb"))
    nbformat.write(nb, dst)
    os.remove(tmp)

    errs, ran = [], 0
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code":
            continue
        outs = c.get("outputs", [])
        if outs:
            ran += 1
        for o in outs:
            if o.get("output_type") == "error":
                errs.append((i, o.get("ename"), (o.get("evalue") or "")[:110]))
    print(f"\n  {time.time()-t0:.0f} s | {ran}/{ncode} code cells produced output | "
          f"{len(errs)} cells errored", flush=True)
    for i, en, ev in errs:
        print(f"    cell {i:3d}  {en}: {ev}", flush=True)
    print(f"  wrote {dst}", flush=True)
