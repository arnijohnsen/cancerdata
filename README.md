# Introduction
This repository contains R scripts which analyse data from
The Cancer Genome Atlas, currently:

- Breast invasive carcinoma (BRCA)
    - 450K Methylation
    - RNASeq
- Glioblastoma (GBM)
    - 450K Methylation
    - RNASeqV2
- Kidney renal clear cell carcinoma (KIRC)
    - 450K Methylation
    - RNASeq
- Ovarian serous cystadenocarcinoma (OV)
    - 27K Methylation
    - RNASeqV2
- Prostate adenocarcinoma (PRAD)
    - 450K Methylation
    - RNASeqV2
- Liver hepatocellular carcinoma (LIHC)
    - 450K Methylation
    - RNASeqV2

Data is available for download at https://tcga-data.nci.nih.gov/tcga/

Author: Arni Johnsen, arni.johnsen@gmail.com

# Installation 

Installation: 

    git clone https://arnijohnsen@bitbucket.org/arnijohnsen/cancerdata.git
    mkdir rawdata
    mkdir Rdata

The Rscripts rely on the following packages:

    WGCNA

# Data structure

Each cancer type (e.g. BRCA) there are several data directories:

    .
    |-- rawdata
    |   `-- BRCA
    |       |-- Methylation
    |       `-- RNASeq
    `-- Rdata
        |-- BRCA
	    |-- calc
	    |-- data
	    |-- info
            `-- tmp

tar files from TCGA should be untarred to one of

- `rawdata/BRCA/Methylation`
- `rawdata/BRCA/RNASeq`

File reading scripts will save .Rdata files to subdirectories of Rdata, 
which contain

- `Rdata/BRCA/calc` result of analysis, such as linear models
- `Rdata/BRCA/data` data.frames containing data read from rawdata
- `Rdata/BRCA/info` information data, such as sample names
- `Rdata/BRCA/tmp ` temporary files, users should never use these

