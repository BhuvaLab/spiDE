"""Build the Python reference pages with quartodoc, then make every docstring Examples section
executable. Run by build_reference.sh from site/.

1. quartodoc 0.11 cannot render a numpy ``Yields`` section (``pl.style`` has one); it is rendered
   here like ``Returns``.
2. quartodoc prints a docstring's ``>>>`` Examples as one unexecuted ``python`` block. Each is
   rewritten as executed Quarto chunks: the statements run on spidepy's own toy_spide() (the
   docstrings are self-contained and doctest-checked in spidepy's CI), and the expected output
   lines of the doctest are dropped, since Quarto shows the real output after the chunk that
   produces it.

Usage: python scripts/build_python_reference.py        (from site/; config in _quarto.yml)
"""
import re
import sys
from pathlib import Path

from quartodoc.__main__ import build as qd_build
from quartodoc._griffe_compat import docstrings as ds
from quartodoc.renderers.md_renderer import MdRenderer

REF = Path("reference/python")


def _render_yields(self, el: ds.DocstringSectionYields):
    rows = [self.render(ds.DocstringReturn(y.name, annotation=y.annotation, description=y.description))
            for y in el.value]
    return self._render_table(rows, ["Name", "Type", "Description"], "returns")


MdRenderer.render.register(_render_yields)

# resolution of the example plots (spidepy sets each figure's own size)
SETUP = '''```{python}
#| echo: false
import matplotlib
matplotlib.rcParams["figure.dpi"] = 90
```'''

BLOCK = re.compile(r"(## Examples[^\n]*\n\n)```python\n(.*?)\n```", re.S)


def show_figures(doctest: str) -> str:
    """The docstrings check a plot with ``type(fig).__name__`` and close it (right for doctest);
    on the page show the figure instead."""
    doctest = re.sub(r"^>>> type\((\w+)\)\.__name__\n'Figure'\n", "", doctest, flags=re.M)
    return re.sub(r"^>>> (?:plt|matplotlib\.pyplot)\.close\((\w+)\)$", r">>> \1", doctest, flags=re.M)


def chunks(doctest: str):
    """Split a doctest into executable statement groups, each ending where output is expected."""
    groups, cur, had_out = [], [], False
    for line in doctest.splitlines():
        if line.startswith((">>> ", "... ")) or line in (">>>", "..."):
            if had_out:  # output of the previous statement is over: close the group
                groups.append(cur)
                cur, had_out = [], False
            cur.append(line[4:])
        elif line.strip():
            had_out = True  # an expected-output line: the group above prints something
        # blank lines inside a doctest only separate; drop
    if cur:
        groups.append(cur)
    return ["```{python}\n" + "\n".join(g) + "\n```" for g in groups]


def main() -> int:
    qd_build.main(["--config", "_quarto.yml"], standalone_mode=False)
    n, bare = 0, []
    for page in sorted(REF.glob("*.qmd")):
        if page.stem == "index":
            continue
        txt = page.read_text(encoding="utf-8")
        m = BLOCK.search(txt)
        if not m or ">>>" not in m.group(2):
            bare.append(page.stem)
            continue
        body = "\n\n".join(chunks(show_figures(m.group(2))))
        txt = txt[:m.start()] + m.group(1) + SETUP + "\n\n" + body + txt[m.end():]
        page.write_text(txt, encoding="utf-8")
        n += 1
    print(f"build_python_reference.py: executable examples on {n} pages")
    if bare:
        print("pages without a docstring Examples section:", ", ".join(bare), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
