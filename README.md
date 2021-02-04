This repository contains supplementary data for 
"Histone dynamics mediate DNA unwrapping and sliding in nucleosomes: insights from multi-microsecond molecular dynamics simulations" by Armeev et al.

Please, use the interactive web site rendered through GitHub pages from this repository at http://intbio.github.io/Armeev_et_al_2021/

## Trajectories availability
Filtered trajectories are available [here](trj)

Full trajectories are available per request.

## Analysis scripts
[Trajectory analysis pipeline example](https://nbviewer.jupyter.org/github/intbio/Armeev_et_al_2021/blob/main/analysis_scripts_examples/pynucl_analysis.ipynb)

[Fiber analysis pipeline example](https://nbviewer.jupyter.org/github/intbio/Armeev_et_al_2021/blob/main/analysis_scripts_examples/fiber_analysis.ipynb)

## Installation
Our analysis pipelines require python version 3.7 and Jupyter notebooks
#### Required libraries (tested versions):
* matplotlib (3.1.3)
* seaborn (0.10.0)
* plotnine (0.6.0)
* numpy (1.17.3)
* scipy (1.3.2)
* pandas (1.0.1)
* tqdm (4.42.1)
* requests (2.23.0)
* nglview (2.7.1)
* Bio (0.3.0)
* MDAnalysis (1.0.0)
* [seqplot](https://github.com/intbio/seqplot)
* [DNAtools](https://github.com/intbio/DNAtools)
* [pymolint](https://github.com/intbio/pymolint)
* [pytexshade](https://github.com/intbio/pytexshade)

#### Required software
* [3DNA](https://x3dna.org/) should be in PATH for the python interpreter.