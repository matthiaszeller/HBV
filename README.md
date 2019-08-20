# HBV GWS project

## About

This repository is the result of my bachelor project as a three-week long internship in *Fellay Lab, EPFL* under the supervision of Sina Rüeger. This project aims to perform a genome-wide association study (GWAS) of HBV-infected individuals from different populations, i.e. find associations between host single nucleotide polymorphisms (SNPs) and viral amino acid variants. The project eventually focuses on a genetically-related asian subpopulation.

## Requirements

* plink version 1.9
* plink version 2
* [Pyhattan](https://github.com/Pudkip/Pyhattan)
* stanard data science packages in Python (pandas, numpy, matplotlib, seaborn)
* have the data stored in a structure that is consistent with the paths defined in the `src/setup.py` file.

## Automatically run notebooks

Method to execute all notebooks in order (!) and output the result as notebooks converted to html. Make sure the data is consistent with the paths defined in the setup file. Note that running all notebooks can take a while. Moreover, you will probably not be able to run it on a computer that has only a few cores and a small amount of RAM, especially for the last notebook.

1. `cd` into the root directory of the project
1. run the bash script `src/run.sh`
1. look at the results in the `results/` directory

## Project structure

The analyses is performed directly inside notebooks, and some of them store processed data. Thus the order of the notebooks (see below) matters. One can simply automatically run notebooks (see above) once the setup file is consistently defined. 

### Data processing & analysis

1. Clinical data notebook: process the clinical data from the `csv` file. Stores a DataFrame object.
1. Viral data notebook
1. Joint viral and clinical data notebook: combine the two datasets. PCA colored with genotypes.
1. Host genotype data preparation notebook: quality control, application of filters
1. Host genotype data analysis notebook
1. GWAS of asian subpopulation

### Documentation & tutorials

* tutorial/Statistics notebook: put altogether relevant information about statistical theoretical background
* tutorial/Plink introduction notebook: basic procedures and commands of plink. Mainly follows the official tutorial.
* tutorial/Plink and Python notebook: how to import plink files into python.
* tutorial/Scitas tutorial
* tutorial/HapMap notebook: processing example data (official tutorial).
