This website documents two packages that implement the same method: the R package **spiDE** and the Python package **spidepy**. This page walks through the R package. The same walk-through for Python, with its outputs, is the [Python tutorial](python-tutorial.qmd), and [R vs Python](r-vs-python.qmd) lists how the two correspond and where they differ. Both packages can be installed as follows:

::: {.panel-tabset group="language"}

## R

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BhuvaLab/spiDE")
```

## Python

spidepy is not on PyPI yet; install it from GitHub. The extra index gives a CPU-only build of PyTorch, which it needs through spanormpy:

```bash
pip install --extra-index-url https://download.pytorch.org/whl/cpu "spidepy @ git+https://github.com/BhuvaLab/spidepy.git"
```

:::
