#!/bin/sh
#
#PBS -N arj32
#PBS -l nodes=1,walltime=03:00:00


module add R/3.0.1
cd /share/scratch/arj32/cancerdata
R CMD BATCH correlationAnalysis/GBM/GBMmethylRnaseqLm.R
