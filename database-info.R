data.sets.list <- c("BRCA", "COAD", "GBM", "KIRC", "LIHC", "LUAD", "LUSC", "OV", "PRAD", "READ")
n <- length(data.sets.list)

CMP.no.samples <- list()
NMP.no.samples <- list()
CEA.no.samples <- list()
NEA.no.samples <- list()
no.probes      <- list()
no.genes       <- list()
C.co.samples   <- list()
N.co.samples   <- list()
for (i in 1:n){
  set <- data.sets.list[i]
  if (file.exists(paste("../Rdata/", set, "/info/", set, "-CMP-samples.Rdata", sep=""))){
    load(paste("../Rdata/", set, "/info/", set, "-CMP-samples.Rdata", sep=""))
    CMP.no.samples[[set]] <- length(get(paste(set,".CMP.samples", sep="")))
  }else{
    CMP.no.samples[[set]] <- 0
  }
  if (file.exists(paste("../Rdata/", set, "/info/", set, "-NMP-samples.Rdata", sep=""))){
    load(paste("../Rdata/", set, "/info/", set, "-NMP-samples.Rdata", sep=""))
    NMP.no.samples[[set]] <- length(get(paste(set,".NMP.samples", sep="")))
  }else{
    NMP.no.samples[[set]] <- 0
  }
  if (file.exists(paste("../Rdata/", set, "/info/", set, "-CEA-samples.Rdata", sep=""))){
    load(paste("../Rdata/", set, "/info/", set, "-CEA-samples.Rdata", sep=""))
    CEA.no.samples[[set]] <- length(get(paste(set,".CEA.samples", sep="")))
  }else{
    CEA.no.samples[[set]] <- 0
  }
  if (file.exists(paste("../Rdata/", set, "/info/", set, "-NEA-samples.Rdata", sep=""))){
    load(paste("../Rdata/", set, "/info/", set, "-NEA-samples.Rdata", sep=""))
    NEA.no.samples[[set]] <- length(get(paste(set,".NEA.samples", sep="")))
  }else{
    NEA.no.samples[[set]] <- 0
  }
  if (file.exists(paste("../Rdata/", set, "/info/", set, "-probes.Rdata", sep=""))){
    load(paste("../Rdata/", set, "/info/", set, "-probes.Rdata", sep=""))
    no.probes[[set]] <- length(get(paste(set,".probes", sep="")))
  }else{
    no.probes[[set]] <- 0
  }
  if (file.exists(paste("../Rdata/", set, "/info/", set, "-genes.Rdata", sep=""))){
    load(paste("../Rdata/", set, "/info/", set, "-genes.Rdata", sep=""))
    no.genes[[set]] <- length(get(paste(set,".genes", sep="")))
  }else{
    no.genes[[set]] <- 0
  }
  C.co.samples[[set]] <- length(intersect(get(paste(set, ".CMP.samples", sep="")),
                                          get(paste(set, ".CEA.samples", sep=""))))
  if (exists(paste(set,".NMP.samples", sep="")) && exists(paste(set,".NEA.samples", sep=""))){
    N.co.samples[[set]] <- length(intersect(get(paste(set, ".NMP.samples", sep="")),
                                            get(paste(set, ".NEA.samples", sep=""))))
  }else{
    N.co.samples[[set]] <- 0
  }
}
all.data.stats <- data.frame(CMP.samples  = unlist(CMP.no.samples),
                             CEA.samples  = unlist(CEA.no.samples),
                             C.co.samples = unlist(C.co.samples),
                             NMP.samples  = unlist(NMP.no.samples),
                             NEA.samples  = unlist(NEA.no.samples),
                             N.co.samples = unlist(N.co.samples),
                             num.probes   = unlist(no.probes),
                             num.genes    = unlist(no.genes))

sink("database-info.txt")
options("width"=120)
print(t(all.data.stats))
sink()
quit(save="no")
