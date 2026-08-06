"""Execute a notebook top to bottom and write outputs back in.

    python run_notebook.py [notebook.ipynb]

This is the check that the deliverable actually runs as claimed, rather than
being a document assembled from scripts that were each run separately.
Any cell that raises stops the run and reports which one.
"""
import sys
import time
import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

NB = sys.argv[1] if len(sys.argv) > 1 else "Part7_CoreValidation.ipynb"

nb = nbformat.read(NB, as_version=4)
print(f"executing {NB}: {len(nb.cells)} cells "
      f"({sum(1 for c in nb.cells if c.cell_type=='code')} code)")

client = NotebookClient(nb, timeout=5000, kernel_name="python3",
                        allow_errors=False, resources={"metadata": {"path": "."}})

t0 = time.time()
try:
    client.execute()
except CellExecutionError as e:
    # locate the offending cell
    for i, c in enumerate(nb.cells):
        if c.cell_type == "code" and any(
                o.get("output_type") == "error" for o in c.get("outputs", [])):
            print(f"\nFAILED in cell {i}:")
            print("  source:", c.source.strip().split("\n")[0][:100])
            for o in c["outputs"]:
                if o.get("output_type") == "error":
                    print("  ", o.get("ename"), ":", o.get("evalue"))
                    print("\n".join("   " + l for l in o.get("traceback", [])[-12:]))
            break
    nbformat.write(nb, NB)
    print(f"\npartial outputs written after {(time.time()-t0)/60:.1f} min")
    sys.exit(1)

nbformat.write(nb, NB)
n_err = sum(1 for c in nb.cells if c.cell_type == "code"
            and any(o.get("output_type") == "error" for o in c.get("outputs", [])))
print(f"\nOK: executed cleanly in {(time.time()-t0)/60:.1f} min, {n_err} cells with errors")
