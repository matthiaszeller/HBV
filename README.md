# HBV GWS project

## Requirements

* plink version 1.9
* plink version 2
* have the data stored in a structure that is consistent with the paths defined in the `src/setup.py` file.
 
## Automatically run notebooks

Method to execute all notebooks in order (!) and output the result as notebooks converted to html. Make sure the data is consistent with the paths defined in the setup file. Note that running all notebooks can take a while.

1. `cd` into the root directory of the project
1. run the bash script `src/run.sh`
1. look at the results in the `results/` directory

## Project structure

### Data processing & analysis

The order matters: some notebooks store processed data, 

1. Clinical data notebook: process the clinical data from the `csv` file. Stores a DataFrame object.
1. Viral data notebook
1. Joint viral and clinical data notebook: combine the two datasets. PCA colored with genotypes.
1. Host genotype data preparation notebook: quality control, application of 

### Documentation & tutorials

* Statistics notebook: put altogether all relevant information about statistical theoretical background
* Plink introduction notebook: basic procedures and commands of plink. Mainly follows the official tutorial.
* Plink and Python notebook: how to import plink files into python.
* Scitas tutorial
* tutorial/HapMap notebook: processing example data (official tutorial).
