# HBV GWS project

## About

This repository is the result of my bachelor project as a three-week long internship in *Fellay Lab, EPFL* under the supervision of Sina Rüeger. This project aims to perform a genome-to-genome  (G2G) study of HBV-infected individuals from different populations, i.e. find associations between host single nucleotide polymorphisms (SNPs) and viral amino acid variants. The project eventually focuses on a genetically-related asian subpopulation.

## Requirements

* [plink version 1.9](https://www.cog-genomics.org/plink2/)
* [plink version 2](https://www.cog-genomics.org/plink/2.0/)
* [Pyhattan](https://github.com/Pudkip/Pyhattan)
* [assocplots](https://github.com/khramts/assocplots)
* Standard data science packages in Python (pandas, numpy, scipy, matplotlib, seaborn)
* Consistent data storage accoding to the paths defined in [`src/setup.py`](src/setup.py)

## Project structure

The analyses is performed directly inside notebooks, and some of them store processed data. Thus the order of the notebooks (see below) matters. One can convert the notebooks to PDF with `nbconvert`, optionally with `--execute` to re-run the notebooks.

### Data desription & format

* Clinical data: `.csv` file
* Viral sequencing data: `.csv` file
* Host exome sequencing data: `.ped` and `.map` files

### Data processing & analysis

All computations are performed in notebooks, which one has to run in the following order:

1. **Clinical data notebook**: process the clinical data from the `csv` file. Stores a DataFrame binary object.
1. **Viral data notebook**: process viral data from a `csv`. Stores a processed DataFrame in a binary file. 
1. **Joint viral and clinical data notebook**: combine the two datasets. PCA colored with genotypes
1. **Host genotype data preparation notebook**: quality control, application of filters
1. **Host genotype data analysis notebook**: PCA, association analyses, clustering
1. **G2G of asian subpopulation**: prepare new dataset, try monovariate models, implement multivariate models
1. **G2G computer**: multivariate models computation, analysis of results
1. **Interpretation of results**: extract and analyse significant associations

### Documentation & tutorials

Useful resources and references:

* tutorial/Statistics notebook: put altogether relevant information about statistical theoretical background
* tutorial/Plink introduction notebook: basic procedures and commands of plink. Mainly follows the official tutorial.
* tutorial/Plink and Python notebook: how to import plink files into python
* tutorial/Scitas tutorial: run jobs and launch Jupyter on remote server
* tutorial/HapMap notebook: processing example data (official tutorial)
