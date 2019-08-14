#!/bin/bash

jupyter nbconvert --execute "Clinical data.ipynb" --to html --output-dir=results
jupyter nbconvert --execute "Viral data.ipynb" --to html --output-dir=results
jupyter nbconvert --execute "Joint viral and clinical data.ipynb" --to html --output-dir=results
jupyter nbconvert --execute "Host genotypes - Data preparation.ipynb" --to html --output-dir=results
jupyter nbconvert --execute "Host genotypes - Analysis.ipynb" --to html --output-dir=results
