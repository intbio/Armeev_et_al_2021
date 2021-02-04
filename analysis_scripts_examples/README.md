## Analysis scripts
### Nucleosome MD trajectory analysis code examples

Custom analysis scripts and pipelines were written in Python 3 using MDAnalysis (coordinate manipulation, 3D alignment), and 3DNA (determination of DNA base pair centers, calculation of base pair and base pair step parameters)
The following Jupyter notebook file provides an example of key analysis procedures used in our study. Please refer to the installation instructions below to set up a working python environment.

[trajectory_analysis.ipynb](https://nbviewer.jupyter.org/github/intbio/Armeev_et_al_2021/blob/main/analysis_scripts_examples/trajectory_analysis.ipynb)

#### Installation
Jupyter notebook interactive environment with appropriate python version 3.7 kernel is required. 
Please refer to installation instructions here https://jupyter.org/install (we suggest using conda).
The list of required python libraries is provided below. 
To use the DNA geometry and unwrapping analysis functionality [3DNA](https://x3dna.org/) should also be installed and callable from the python interpreter.
##### Required libraries (tested versions):
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

### [NCP fiber analysis pipeline example](https://nbviewer.jupyter.org/github/intbio/Armeev_et_al_2021/blob/main/analysis_scripts_examples/fiber_analysis.ipynb)