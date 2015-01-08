data.file.dir = "../rawdata/BRCA/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
output.dir    = "../Rdata/BRCA-new/methylation/parsed-data/"

load("../Rdata/BRCA-new/methylation/file-lists/matched-normal-file-names.Rdata")
load("../Rdata/BRCA-new/methylation/file-lists/matched-cancer-file-names.Rdata")
load("../Rdata/BRCA-new/methylation/file-lists/unmatched-cancer-file-names.Rdata")

load("../Rdata/BRCA-new/methylation/probe-annotations/body-island-probes.Rdata")
load("../Rdata/BRCA-new/methylation/probe-annotations/body-shore-probes.Rdata")
load("../Rdata/BRCA-new/methylation/probe-annotations/body-none-probes.Rdata")
load("../Rdata/BRCA-new/methylation/probe-annotations/promoter-probes.Rdata")
load("../Rdata/BRCA-new/methylation/probe-annotations/enhancer-probes.Rdata")

chromosomes <- c(1:22, "X", "Y")

##########################
# Matched normal samples #
##########################

for(j in 1:length(chromosomes)){
  write(paste(body.island.probes[[j]], collapse="\t"), file=paste(output.dir,"matched-normal/body-island/chr-",chromosomes[j],".txt",sep=""))
  write(paste(body.shore.probes[[j]],  collapse="\t"), file=paste(output.dir,"matched-normal/body-shore/chr-",chromosomes[j],".txt",sep=""))
  write(paste(body.none.probes[[j]],   collapse="\t"), file=paste(output.dir,"matched-normal/body-none/chr-",chromosomes[j],".txt",sep=""))
  write(paste(promoter.probes[[j]],    collapse="\t"), file=paste(output.dir,"matched-normal/promoter/chr-",chromosomes[j],".txt",sep=""))
  write(paste(enhancer.probes[[j]],    collapse="\t"), file=paste(output.dir,"matched-normal/enhancer/chr-",chromosomes[j],".txt",sep=""))
}

for(i in 1:length(matched.normal.file.names)){
  cat("Reading matched normal file", i, "of", length(matched.normal.file.names), "\n")
  raw.file.data <- read.table(paste(data.file.dir, matched.normal.file.names[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1, stringsAsFactors=F)
  sample.name <- sub(".txt","", sub(".*TCGA","TCGA", matched.normal.file.names[i]))

  for(j in 1:length(chromosomes)){
    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.island.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-normal/body-island/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.shore.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-normal/body-shore/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.none.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-normal/body-none/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% promoter.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-normal/promoter/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% enhancer.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-normal/enhancer/chr-",chromosomes[j],".txt",sep=""), append=T)
  }
}


##########################
# Matched cancer samples #
##########################

for(j in 1:length(chromosomes)){
  write(paste(body.island.probes[[j]], collapse="\t"), file=paste(output.dir,"matched-cancer/body-island/chr-",chromosomes[j],".txt",sep=""))
  write(paste(body.shore.probes[[j]],  collapse="\t"), file=paste(output.dir,"matched-cancer/body-shore/chr-",chromosomes[j],".txt",sep=""))
  write(paste(body.none.probes[[j]],   collapse="\t"), file=paste(output.dir,"matched-cancer/body-none/chr-",chromosomes[j],".txt",sep=""))
  write(paste(promoter.probes[[j]],    collapse="\t"), file=paste(output.dir,"matched-cancer/promoter/chr-",chromosomes[j],".txt",sep=""))
  write(paste(enhancer.probes[[j]],    collapse="\t"), file=paste(output.dir,"matched-cancer/enhancer/chr-",chromosomes[j],".txt",sep=""))
}

for(i in 1:length(matched.cancer.file.names)){
  cat("Reading matched cancer file", i, "of", length(matched.cancer.file.names), "\n")
  raw.file.data <- read.table(paste(data.file.dir, matched.cancer.file.names[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1, stringsAsFactors=F)
  sample.name <- sub(".txt","", sub(".*TCGA","TCGA", matched.cancer.file.names[i]))

  for(j in 1:length(chromosomes)){
    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.island.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-cancer/body-island/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.shore.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-cancer/body-shore/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.none.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-cancer/body-none/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% promoter.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-cancer/promoter/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% enhancer.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"matched-cancer/enhancer/chr-",chromosomes[j],".txt",sep=""), append=T)
  }
}


############################
# Unmatched cancer samples #
############################

for(j in 1:length(chromosomes)){
  write(paste(body.island.probes[[j]], collapse="\t"), file=paste(output.dir,"unmatched-cancer/body-island/chr-",chromosomes[j],".txt",sep=""))
  write(paste(body.shore.probes[[j]],  collapse="\t"), file=paste(output.dir,"unmatched-cancer/body-shore/chr-",chromosomes[j],".txt",sep=""))
  write(paste(body.none.probes[[j]],   collapse="\t"), file=paste(output.dir,"unmatched-cancer/body-none/chr-",chromosomes[j],".txt",sep=""))
  write(paste(promoter.probes[[j]],    collapse="\t"), file=paste(output.dir,"unmatched-cancer/promoter/chr-",chromosomes[j],".txt",sep=""))
  write(paste(enhancer.probes[[j]],    collapse="\t"), file=paste(output.dir,"unmatched-cancer/enhancer/chr-",chromosomes[j],".txt",sep=""))
}

for(i in 1:length(unmatched.cancer.file.names)){
  cat("Reading unmatched cancer file", i, "of", length(unmatched.cancer.file.names), "\n")
  raw.file.data <- read.table(paste(data.file.dir, unmatched.cancer.file.names[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1, stringsAsFactors=F)
  sample.name <- sub(".txt","", sub(".*TCGA","TCGA", unmatched.cancer.file.names[i]))

  for(j in 1:length(chromosomes)){
    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.island.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"unmatched-cancer/body-island/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.shore.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"unmatched-cancer/body-shore/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% body.none.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"unmatched-cancer/body-none/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% promoter.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"unmatched-cancer/promoter/chr-",chromosomes[j],".txt",sep=""), append=T)

    tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% enhancer.probes[[j]] ]
    write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir,"unmatched-cancer/enhancer/chr-",chromosomes[j],".txt",sep=""), append=T)
  }
}
quit(save="no")
