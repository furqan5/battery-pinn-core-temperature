"""Part 6 seed sweep: how stable is the fixed-residual fit actually?

Motivation. The notebook quotes R0 = 144.78 mOhm at w = 20 and scores prediction
G2 ("within +/-3 % of the cycle-0 classical 144.35") as a HIT at +0.30 %.
Re-executing the same seeded cell in this repository returned 152.83 mOhm, which
is +5.87 % and a MISS. Both runs are n = 1. Neither is wrong; the quantity was
simply never characterised.

Note the notebook DOES set torch.manual_seed(seed), so the divergence is across
environments (torch build, thread count, reduction order), not across seeds.
That is itself worth reporting: a seeded notebook that is not reproducible
across machines cannot support a value quoted to five significant figures.

Method: execute the notebook's setup cells to build the state, extract the
training function, then call it across seeds at both data weights.
"""
import json
import os
import re
import sys
import time

import nbformat
from nbclient import NotebookClient

REPO = r"C:\Users\Nouman\Desktop\Furqan's Docs\battery_pinn"
SRC = os.path.join(REPO, "verify", "Part6_Real_Data_PINN_executed.ipynb")
N_SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
WEIGHTS = [float(w) for w in (sys.argv[2].split(",") if len(sys.argv) > 2
                              else ["20", "200"])]

nb = nbformat.read(SRC, as_version=4)

# Cell 14 defines train_real_I and then runs it twice. Keep only the definition
# so the sweep does not pay for two throwaway trainings.
src14 = "".join(nb.cells[14]["source"]) if isinstance(nb.cells[14]["source"], list) \
    else nb.cells[14]["source"]
cut = src14.find('print(f"\\n{\'variant\'')
if cut < 0:
    cut = src14.find("print(f\"\\n{'variant'")
if cut < 0:
    raise SystemExit("could not locate the end of the training-function block")
defs_only = src14[:cut]
if "def train_real_I" not in defs_only:
    raise SystemExit("train_real_I definition not captured")

sweep = f'''
import numpy as _np, json as _json, time as _t
_res = []
for _w in {WEIGHTS!r}:
    for _s in range({N_SEEDS}):
        _t0 = _t.time()
        _r = train_real_I(P, seed=_s, wd2=_w)
        _r["w_data"] = _w; _r["seed"] = _s; _r["wall"] = _t.time() - _t0
        _res.append(_r)
        print(f"w={{_w:6.0f}} seed={{_s}}  R0={{_r['R0_mOhm']:8.3f}} mOhm  "
              f"c1={{_r['c1']:+7.3f}}  c2={{_r['c2']:+7.3f}}  "
              f"RMSE={{_r['rmse']:.4f}}  Lr={{_r['Lr']:.3e}}  "
              f"[{{_r['wall']:.0f}} s]", flush=True)
open(r"{os.path.join(REPO, 'results', 'p6_seeds.json')}", "w").write(_json.dumps(_res))

_target = c0f["R0_mOhm"]
print()
print("=" * 92)
print(f"CYCLE-0 CLASSICAL TARGET: R0 = {{_target:.3f}} mOhm, c1 = {{c0f['c1']:.3f}}, "
      f"c2 = {{c0f['c2']:.3f}}, RMSE = {{c0f['rmse']:.4f}} K")
print("=" * 92)
for _w in {WEIGHTS!r}:
    _g = [r for r in _res if r["w_data"] == _w]
    _R = _np.array([r["R0_mOhm"] for r in _g])
    _c1 = _np.array([r["c1"] for r in _g])
    _rm = _np.array([r["rmse"] for r in _g])
    _dev = 100 * (_R - _target) / _target
    print(f"\\n  w_data = {{_w:.0f}}   n = {{len(_g)}}")
    print(f"    R0    mean {{_R.mean():8.3f}}  sd {{_R.std(ddof=1):6.3f}}  "
          f"min {{_R.min():8.3f}}  max {{_R.max():8.3f}}  spread {{_R.ptp():6.3f}} mOhm")
    print(f"    vs classical: mean {{_dev.mean():+6.2f}} %  range "
          f"[{{_dev.min():+6.2f}}, {{_dev.max():+6.2f}}] %")
    print(f"    c1    mean {{_c1.mean():+7.3f}}  sd {{_c1.std(ddof=1):6.3f}}  "
          f"(all negative: {{bool((_c1 < 0).all())}})")
    print(f"    RMSE  mean {{_rm.mean():7.4f}}  sd {{_rm.std(ddof=1):6.4f}}")
    _hit = _np.abs(_dev) <= 3.0
    print(f"    G2 (|dev| <= 3 %%): {{int(_hit.sum())}}/{{len(_g)}} seeds would score HIT")
'''

nb.cells = nb.cells[:14] + [
    nbformat.v4.new_code_cell(defs_only),
    nbformat.v4.new_code_cell(sweep),
]

print(f"executing setup + {N_SEEDS} seeds x {len(WEIGHTS)} weights", flush=True)
t0 = time.time()
client = NotebookClient(nb, timeout=20000, kernel_name="python3",
                        allow_errors=False, resources={"metadata": {"path": REPO}})
client.execute()
nbformat.write(nb, os.path.join(REPO, "verify", "Part6_seed_sweep.ipynb"))
print(f"\ndone in {(time.time()-t0)/60:.1f} min", flush=True)
