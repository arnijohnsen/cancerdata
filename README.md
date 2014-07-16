# Introduction
This repository contains R scripts which analyse data from
The Cancer Genome Atlas, currently:

- BRCA (Breast invasive carcinoma):        450K / RNASeq
- GBM  (Glioblastoma)                      450K / RNASeqV2
- KIRC (Kidney renal clear cell carcinoma) 450K / RNASeq
- OV   (Ovarian serous cystadenocarcinoma)  27K / RNASeqV2
- PRAD (Prostate adenocarcinoma)           450K / RNASeqV2
- LICH (Liver hepatocellular carcinoma)    450K / RNASeqV2

Data is available for download at https://tcga-data.nci.nih.gov/tcga/

Installation: 
    git clone https://arnijohnsen@bitbucket.org/arnijohnsen/cancerdata.git
    mkdir rawdata
    mkdir Rdata

The Rscripts rely on the following packages:
    WGCNA

The data sets used are Methylation and RNA squencing, 
which should be placed in a directory
    ../rawdata/
relative to where the source code was cloned. 
This directory should also contain a subdirectory with probe
annotation files, ../rawdata/AnnotationFiles/

Scripts found in convertRawData/ convert this data to .Rdata files, 
which are saved in to a directory called
    ../Rdata/

Data files from TCGA must be extracted using tar, prior to execution. 

A sample file structure looks as follows: 

arj32
|-- cancerdata
|   |-- convertRawData
|   |   `-- BRCA
|   |       |-- generateProbeGenesList.R
|   |       |-- methylPromoters.R
|   |       `-- rnaseqAllGenes.R
|   |-- correlationAnalysis
|   |   `-- methylRnaseqCorrelation.R
|   |-- plotUtils
|   |   `-- plotMethylationExpression.R
|   `-- README
|-- rawdata
|   |-- AnnotationFiles
|   |   `-- GenomeStudioProbeAnnotations.txt
|   `-- BRCA
|       |-- Methylation
|       |   |-- DNA_Methylation
|       |   |   `-- JHU_USC__HumanMethylation450
|       |   |       `-- Level_3
|       |   |           |-- jhu-usc.edu_BRCA.HumanMethylation450.10.lvl-3.TCGA-A2-A1FV-01A-11D-A13K-05.txt
|       |   |           |-- jhu-usc.edu_BRCA.HumanMethylation450.10.lvl-3.TCGA-A2-A1FW-01A-11D-A13K-05.txt
|       |   |           |-- jhu-usc.edu_BRCA.HumanMethylation450.10.lvl-3.TCGA-A2-A1FX-01A-11D-A13K-05.txt
|       |   |               ... REST OF BRCA DNA METHYLATION FILES ...
|       |   |-- file_manifest.txt
|       |   |-- FILE_SAMPLE_MAP.txt
|       |   |-- METADATA
|       |   |   `-- JHU_USC__HumanMethylation450
|       |   |       |-- jhu-usc.edu_BRCA.HumanMethylation450.1.8.0.idf.txt
|       |   |       `-- jhu-usc.edu_BRCA.HumanMethylation450.1.8.0.sdrf.txt
|       |   `-- README_DCC.txt
|       `-- RNASeq
|           |-- file_manifest.txt
|           |-- FILE_SAMPLE_MAP.txt
|           |-- METADATA
|           |   `-- UNC__IlluminaHiSeq_RNASeq
|           |       |-- unc.edu_BRCA.IlluminaHiSeq_RNASeq.2.idf.txt
|           |       `-- unc.edu_BRCA.IlluminaHiSeq_RNASeq.2.sdrf.txt
|           |-- README_DCC.txt
|           `-- RNASeq
|               `-- UNC__IlluminaHiSeq_RNASeq
|                   `-- Level_3
|                       |-- UNCID_1109810.TCGA-GM-A2D9-01A-11R-A18M-07.111101_UNC13-SN749_0133_BC04U8ABXX.4_2.trimmed.annotated.translated_to_genomic.spljxn.quantification.txt
|                       |-- UNCID_1109879.TCGA-GM-A2DC-01A-11R-A18M-07.111101_UNC13-SN749_0133_BC04U8ABXX.7_2.trimmed.annotated.gene.quantification.txt
|                       |-- UNCID_1110044.TCGA-GM-A2D9-01A-11R-A18M-07.111101_UNC13-SN749_0133_BC04U8ABXX.4_2.trimmed.annotated.gene.quantification.txt
|                       |-- UNCID_1110094.TCGA-GM-A2DB-01A-31R-A18M-07.111101_UNC13-SN749_0133_BC04U8ABXX.6_2.trimmed.annotated.translated_to_genomic.spljxn.quantification.txt
|                       |-- UNCID_1110118.TCGA-GM-A2DB-01A-31R-A18M-07.111101_UNC13-SN749_0133_BC04U8ABXX.6_2.trimmed.annotated.gene.quantification.txt
|                       |-- UNCID_1110143.TCGA-GM-A2D9-01A-11R-A18M-07.111101_UNC13-SN749_0133_BC04U8ABXX.4_2.trimmed.annotated.translated_to_genomic.exon.quantification.txt
|                           ... REST OF BRCA RNA SEQUENCING FILES ...
`-- Rdata
    `-- BRCA
        |-- cancerRnaseqAllGenes.Rdata
        |-- genesList.Rdata
        |-- normalRnaseqAllGenes.Rdata
        `-- rnaseqSampList.Rdata

Author: Arni Johnsen, arni.johnsen@gmail.com
