#!/bin/bash

### Job name
#SBATCH -J build_site 
#SBATCH -e ./build_site.txt
#SBATCH -o ./build_site.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 120:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=360G -c 64 -x big-node001

source ~/.bashrc
conda activate r-4.1

R -e "pkgdown::build_site(preview = FALSE)"
