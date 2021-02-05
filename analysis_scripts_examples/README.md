# Analysis scripts
## Nucleosome MD trajectory analysis code examples

Custom analysis scripts and pipelines were written in Python 3 using MDAnalysis (coordinate manipulation, 3D alignment), and 3DNA (determination of DNA base pair centers, calculation of base pair and base pair step parameters)
The following Jupyter notebook file provides an example of key analysis procedures used in our study. Please refer to the installation instructions below to set up a working python environment.

To perform sample analysis please open and run the following jupyter notebook (expected run time - 10 minutes), expected output is already present in the notebook.

[trajectory_analysis.ipynb](html/trajectory_analysis.html)

## Nucleosome fiber generation and analysis algorithms
This python script use snapshots from MD trajectories to connect nucleosomes with straight B-DNA linkers.

To perform sample analysis please open and run the following jupyter notebook (expected run time - 5 minutes), expected output is already present in the notebook.

[fiber_generation_with_analysis.ipynb](html/fiber_generation_with_analysis.html)

## Installation
Software environment was tested on Ubuntu server 18.04.
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
* [pynucl](https://github.com/intbio/pynucl)
* [seqplot](https://github.com/intbio/seqplot)
* [DNAtools](https://github.com/intbio/DNAtools)
* [pymolint](https://github.com/intbio/pymolint)
* [pytexshade](https://github.com/intbio/pytexshade)
* [pynamod](https://github.com/intbio/pynamod) (required for fiber generation and analysis)

##### Install time
Around 30 minutes
