#!/bin/sh

#SBATCH --time 6:00:00

module load gcc
module load python
module list

echo "LAUNCHING JUPYTER NBCONVERT\n\n"

jupyter nbconvert --execute "/home/mazeller/HBV/G2G computer.nbconvert.ipynb" --to notebook --ExecutePreprocessor.timeout=-1

echo "\n\n\nJUPYTER NBCONVERT HAS ENDED"

