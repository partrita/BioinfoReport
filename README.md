# Therapeutic protein development analysis
[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://partrita.github.io)

The reasons for starting this project are as follows.

> Therapeutic proteins are proteins that are engineered in the laboratory for pharmaceutical use, and it's very nasty to make.

Bioinformatics tools are widely used to develop protein pharmaceuticals. With so many tools being used, I've found that they are often applied based on one's intuition and experience without a proper analysis pipeline.

Then, I saw a demo report on an antibody analysis service in a company called [amber biolog](https://www.amberbiology.com/), and thought it would be good to write reproducible code using Jupyter notebook.

# How to use

```bash
# build html
$jupyter-book build demo-report/Aflibercept --all
# build pdf by latex
$jupyter-book build demo-report/Aflibercept --builder pdflatex

```

# Tools I used

- [PDM](https://pdm.fming.dev/): Python package manager with PEP 582 support.
- [Jupyer book](https://jupyterbook.org/en/stable/intro.html)
- Biopython
- [Biotite](https://www.biotite-python.org/tutorial/target/index.html)
- netMHCpan 4.1
- netMHCIIPan 4.0
