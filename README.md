# HBV GWS project

## Automatically run notebooks

1. `cd` into the root directory of the project
1. run the bash script `src/run.sh`
1. look at the results in the `results/` directory

## Project structure

### Data processing & analysis

The order matters: some notebooks store processed data, 

1. Clinical data notebook: process the clinical data from the `csv` file. Stores a DataFrame object.
1. Viral data notebook
1. Joint viral and clinical data: combine the two datasets. PCA colored with genotypes.

### Documentation & tutorials

* Statistics notebook: put altogether all relevant information about statistical theoretical background
* Plink introduction notebook: basic procedures and commands of plink. Mainly follows the official tutorial.
* Plink and Python notebook: how to import plink files into python.
* Scitas tutorial
* tutorial/HapMap notebook: processing example data (official tutorial).

## TO DO / OPTIMIZE:

* **Find a way to manage paths **
* Assemble scripts together ?
* Make a kind of setup file (python, bash, both ???)
* Manage figure storage (some kind of option at beggining of the notebook, or make it separately again ??)

If more time:
* Look at makefiles
