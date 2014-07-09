#!/bin/sh
#
#PBS -N arj32
#PBS -l nodes=1,walltime=06:00:00


module add R/3.0.1
cd /share/scratch/arj32/cancerdata
Rscript convertRawData/BRCA/methylPromoters.R
