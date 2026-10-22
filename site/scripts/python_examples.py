"""Append an executed Examples section to every generated spidepy reference page.

spidepy's docstrings carry no Examples yet, so the example code lives here, keyed by the
page name quartodoc gives (``build_niches``, ``pl.plot_niche_map``, ...). When a docstring gains
runnable Examples, delete its entry here (quartodoc then renders the docstring's own).

Every example starts from the shared toy objects built once per build by make_python_cache.py:
``adata`` (toy_spide() with the sigma = 30 niches), ``fit`` and ``res`` (the slopes engine on all
index types, with the condition). Functions that mutate their input get a copy first.

Usage: python scripts/python_examples.py        (run by build_reference.sh, from site/)
Entries in SKIP are shown unexecuted; there are none today.
"""
import re
import sys
from pathlib import Path

SKIP: dict[str, str] = {}

SETUP = '''```{python}
#| echo: false
import os, pickle, matplotlib
import matplotlib.pyplot as plt
import spidepy as sp
with open(os.path.join(os.environ.get("QUARTO_PROJECT_DIR", "../.."), "_cache", "toy_spide.pkl"), "rb") as f:
    adata, fit, res = pickle.load(f)
```'''

EX = {
    "build_niches": '''adata0 = sp.datasets.toy_spide()
adata0 = sp.build_niches(adata0, sigma=30)
adata0.obsm["Niche30"].head()''',
    "merge_niches": '''a = sp.merge_niches(adata, groups={"AC": ["A", "C"]}, sigma=30, copy=True)
a.obsm["Niche30"].columns.tolist()''',
    "fit_spide": '''fit = sp.fit_spide(adata, index="A", condition="condition", sigma=30, verbose=False)
fit''',
    "test_spide": '''res_a = sp.test_spide(fit)
res_a''',
    "spide": '''a = adata.copy()
res_a = sp.spide(a, condition="condition", sigma=30, index="A", verbose=False)
sp.results(res_a, test="pooled").head()''',
    "results": '''sp.results(res, test="pooled").head()''',
    "patient_slopes": '''sp.patient_slopes(fit, gene="G1", index="A", niche="B").head()''',
    "patient_intercepts": '''sp.patient_intercepts(fit, gene="G1", index="A").head()''',
    "SpiDEFit": '''fit''',
    "SpiDEResults": '''res''',
    "spi_gsea": '''sets = {"first": [f"G{i}" for i in range(1, 6)], "second": [f"G{i}" for i in range(6, 13)]}
sp.spi_gsea(res, sets, min_size=3)''',
    "test_niche_abundance": '''na = sp.test_niche_abundance(adata, condition="condition", sigma=30)
na.sort_values("p").head()''',
    "toy_spide": '''a = sp.datasets.toy_spide()
a''',
    "datasets.toy_spide": '''a = sp.datasets.toy_spide()
a''',
    "datasets.simulate_toy": '''a = sp.datasets.simulate_toy(n_samples=4, n_per=60, seed=1)
a''',
    "datasets.simulate_niche": '''a = sp.datasets.simulate_niche(seed=1)
a''',
    "datasets.simulate_clustered": '''a = sp.datasets.simulate_clustered(seed=1)
a''',
    "pl.plot_niche_map": '''sp.pl.plot_niche_map(adata, niche="B", sample="S1", index="A")''',
    "pl.plot_niche_response": '''sp.pl.plot_niche_response(adata, res, gene="G1", index="A", niche="B")''',
    "pl.plot_patient_slopes": '''sp.pl.plot_patient_slopes(res, gene="G1", index="A", niche="B")''',
    "pl.plot_patient_effects": '''sp.pl.plot_patient_effects(res, index="A")''',
    "pl.plot_call_map": '''sp.pl.plot_call_map(res)''',
    "pl.plot_top_triplets": '''sp.pl.plot_top_triplets(res, test="pooled", n=15, adata=adata)''',
    "pl.plot_spillover": '''sp.pl.plot_spillover(res, adata)''',
    "pl.plot_triplet_heatmap": '''sp.pl.plot_triplet_heatmap(res, test="pooled", n=10)''',
    "pl.plot_p_values": '''sp.pl.plot_p_values(res)''',
    "pl.spide_colours": '''sp.pl.spide_colours("condition")''',
    "pl.style": '''with sp.pl.style():
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot([0, 1], [0, 1])
fig''',
}


def clean(txt: str) -> str:
    """Docstrings are written for Sphinx; show their roles and double backticks as Markdown."""
    txt = re.sub(r":(?:func|class|meth|mod|data|attr):`~?([\w.]+)`",
                 lambda m: "`" + m.group(1).split(".")[-1] + "`", txt)
    return re.sub(r"``([^`\n]+)``", r"`\1`", txt)


def main(ref: Path) -> int:
    missing = []
    done = 0
    for name, code in EX.items():
        page = ref / f"{name}.qmd"
        if not page.exists():
            continue
        if name in SKIP:
            page.write_text(clean(page.read_text(encoding="utf-8")), encoding="utf-8")
            continue
        txt = clean(page.read_text(encoding="utf-8"))
        # an Examples section from the docstring itself wins
        if "## Examples" in txt:
            continue
        txt += f"\n\n## Examples\n\n{SETUP}\n\n```{{python}}\n{code}\n```\n"
        page.write_text(txt, encoding="utf-8")
        done += 1
    pages = {p.stem for p in ref.glob("*.qmd") if p.stem not in ("index",)}
    missing = sorted(pages - set(EX) - set(SKIP))
    if missing:
        print("python_examples.py: pages without an example entry:", ", ".join(missing), file=sys.stderr)
        return 1
    print(f"python_examples.py: examples appended to {done} pages")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path("reference/python")))
