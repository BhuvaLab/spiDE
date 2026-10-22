"""Fit spidepy once on the toy data and pickle (adata, fit, res) to _cache/toy_spide.pkl, so the
executed examples of the Python reference share one fit instead of refitting on every page.
Run by build_reference.sh from site/.
"""
import pickle
import time
from pathlib import Path

import spidepy as sp

t = time.time()
adata = sp.datasets.toy_spide()
sp.build_niches(adata, sigma=30)
res = sp.spide(adata, condition="condition", sigma=30, verbose=False)
Path("_cache").mkdir(exist_ok=True)
with open("_cache/toy_spide.pkl", "wb") as f:
    pickle.dump((adata, res.fit, res), f)
print(f"make_python_cache.py: toy fit cached in {time.time() - t:.0f} s")
