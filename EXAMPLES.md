# Examples
This file includes examples on how to use data created by the scripts in this repository

## Start R on Jotunn
Request an interactive session (not necessary, but it's favourable to running R on login node): 

    qsub -I -X

Enter the correct directory, on my system: 

    cd /share/scratch/arj32/cancerdata

Add `R` module and run `R`

    module add R
    R

I recommend using Xforwarding to view plots created in interactive sessions. On Windows, this can be achived by running Xming on your computer and enabling Xforwarding in Putty.

## Create a simple methylation plot
Start by loading expression and methylation data for cancer samples from the *BRCA* data set

    > load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
    > load("../Rdata/BRCA/data/BRCA-CEA.Rdata")

Check if data loaded correctly

    > dim(BRCA.CMP)
    [1]    732 107599
    > dim(BRCA.CEA)
    [1]  1059 19901

Note that the methylation set `BRCA.CMP` has 732 samples (rows) but the expression set `BRCA.CEA` only has 1059 samples (rows), so we must (for this analysis) use only samples which exist in both sets:

    > cancer.samples <- intersect(rownames(BRCA.CMP), rownames(BRCA.CEA))
    > length(cancer.samples)
    [1] 731

So 731 samples are included in both data sets. 

Now we need to choose a probe or gene, for this we need to load a linked list of probe and gene names

    > load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")
    > head(BRCA.linked.probes.genes)
    probes  genes
    1 cg00000622  NIPA2
    2 cg00000734   CNBP
    3 cg00000769  DDX55
    4 cg00001245 MRPS25
    5 cg00001349   MAEL
    6 cg00001582  ZMIZ1

This plot will include data for the *BRCA1* gene, so I search for `BRCA1` in the gene list

    > grep("BRCA1", BRCA.linked.probes.genes$genes)
     [1]  17885  19966  20238  30347  36048  38510  40265  45944  53914  57092
    [11]  63895  68680  69963  71184  77779  79370  81943  84024  85813  99076
    [21] 107102

There are 21 probes for the BRCA1 gene, so one must be selected. I'll select the first result, number 17885

    > i <- 17885
    > BRCA.linked.probes.genes$probes[i]
    [1] cg04110421
    105182 Levels: cg00000622 cg00000734 cg00000769 cg00001245 ... ch.X.881546R
    > BRCA.linked.probes.genes$genes[i]
    [1] BRCA1
    17692 Levels: A1BG A1CF A2BP1 A2LD1 A2ML1 A4GALT A4GNT AAA1 AAAS AACS ... psiTPTE22

Now we select the data

    > x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[i])]
    > y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[i])]
    > str(x)
     num [1:731] 0.0333 0.0529 0.0274 0.0289 0.0274 ...
    > str(y)
     num [1:731] 606 329 981 426 675 ...

Finally, we can create a plot

    > plot(x,y, pch=20, xlab="Methylation", ylab="Expression")
