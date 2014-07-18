# Introduction
This repository contains R scripts which analyse data from
The Cancer Genome Atlas, currently:

- Breast invasive carcinoma (BRCA)
    - Uses 450K Methylation and RNASeq
- Colon adenocarcinoma (COAD)
    - Uses 450K Methylation and RNASeqV2
- Glioblastoma (GBM)
    - Uses 450K Methylation and RNASeqV2
- Kidney renal clear cell carcinoma (KIRC)
    - Uses 450K Methylation and RNASeq
- Liver hepatocellular carcinoma (LIHC)
    - Uses 450K Methylation and RNASeqV2
- Lung adenocarinoma (LUAD)
    - Uses 450K Methylation and RNASeqV2
- Lung squamous cell carcinoma(LUSC)
    - Uses 450K Methylation and RNASeqV2
- Ovarian serous cystadenocarcinoma (OV)
    - Uses 27K Methylation and RNASeqV2
- Prostate adenocarcinoma (PRAD)
    - Uses 450K Methylation and RNASeqV2
- Rectal adenocarcinoma (READ)
    - Uses 450K Methylation and RNASeqV2

Data is available for download at https://tcga-data.nci.nih.gov/tcga/

Author: Arni Johnsen, arni.johnsen@gmail.com

# Installation 

Installation: 

    git clone https://arnijohnsen@bitbucket.org/arnijohnsen/cancerdata.git
    mkdir rawdata
    mkdir Rdata

The Rscripts rely on the following packages:

    WGCNA

# Structure

## Directory structure

Each cancer type (e.g. BRCA) there are several data directories:

    .
    |-- rawdata
    |   `-- BRCA
    |       |-- Methylation
    |       `-- RNASeq
    `-- Rdata
        `-- BRCA
            |-- calc
            |-- data
            |-- info
            `-- tmp

tar files from TCGA should be untarred to one of

- `rawdata/BRCA/Methylation`
- `rawdata/BRCA/RNASeq`

File reading scripts will save .Rdata files to subdirectories of Rdata, 
which contain

- `Rdata/BRCA/calc` Result of analysis, such as linear models
- `Rdata/BRCA/data` Data.frames containing data read from rawdata
- `Rdata/BRCA/info` Information data, such as sample names
- `Rdata/BRCA/tmp ` Temporary files, users should never use these

## Data frame structure

Data contained in `Rdata/XXXX/data` has a consistant nomenclature and format: 

- Data files are name `XXXX-YYY.Rdata`, where 
    - `XXXX` is the cancer type (e.g. BRCA)
    - `YYY` is one of the following
        - `CMP` Methylation of Promoters in Cancer samples
        - `NMP` Methylation of Promoters in Normal samples
        - `CEA` Expression of All gene in Cancer samples
        - `NEA` Expression of All gene in Normal samples
- Each data file contains one `data.frame`, named `XXXX.YYY`, following the same nomenclature
- Each row of a data frame represents one sample and is named by its TCGA barcode (first 14 letters)
- Each column of a data frame represents either
    - one probe methylation, and is named by its cg number (e.g. `cg04582861`)
    - one gene expression, and is named by its HGNC name (e.g. `BRCA1`)

A data frame linking probes to genes can be found in `Rdata/XXXX/info/XXXX-linked-probes-genes.Rdata`. 
It contains two columns, `probe` and `gene`, each row corresponds to a probe in a promoter region for a gene.

