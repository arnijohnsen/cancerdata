# Ask user which data sets he wants to use
cat("Enter data sets you want to use, seperated by spaces\n")
cat("(Available are BRCA GBM KIRC OV)\n")
data.sets.input <- readline()
data.sets.list <- unlist(strsplit(data.sets.input, " "))

for (i in 1:length(data.sets.list)){
  if(!exists(paste(data.sets.list[i], ".linked.probes.genes", sep=""))){
    cat("Loading", data.sets.list[i], "probe annotations\n")
    load(paste("../Rdata/", data.sets.list[i], "/info/", data.sets.list[i], "-linked-probes-genes.Rdata", sep=""))
  }
  if(!exists(paste(data.sets.list[i], ".CMP", sep=""))){
    cat("Loading", data.sets.list[i], "cancer methylation\n")
    load(paste("../Rdata/", data.sets.list[i], "/data/", data.sets.list[i], "-CMP.Rdata", sep=""))
  }
  if(!exists(paste(data.sets.list[i], ".CEA", sep=""))){
    cat("Loading", data.sets.list[i], "cancer expression\n")
    load(paste("../Rdata/", data.sets.list[i], "/data/", data.sets.list[i], "-CEA.Rdata", sep=""))
  }
}

search <- readline("Enter gene or probe name: ")
search <- paste("^", search, "$", sep="")

n.list <- list()


for(i in 1:length(data.sets.list)){
  n.list[[data.sets.list[i]]] <- union(
         grep(search, get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))$probes), 
         grep(search, get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))$genes))
}

all.probes <- character(0)
all.genes  <- character(0)
for(i in 1:length(data.sets.list)){
  all.probes <- c(all.probes, 
    as.character(get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))[n.list[[data.sets.list[i]]], ]$probes))
  all.genes <- c(all.genes, 
    as.character(get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))[n.list[[data.sets.list[i]]], ]$genes))
}

unique.probes.genes <- strsplit(unique(paste(all.probes, all.genes, sep="-")), "-")

layout(matrix(seq(1,2*length(data.sets.list)), ncol=length(data.sets.list)))

for(i in 1:length(unique.probes.genes)){
  probe <- unique.probes.genes[[i]][1]
  gene  <- unique.probes.genes[[i]][2]
  cat(probe, gene, "\n")
  for(i in 1:length(data.sets.list)){
    cancer.samples <- intersect(rownames(get(paste(data.sets.list[i], ".CMP", sep=""))), 
                                rownames(get(paste(data.sets.list[i], ".CEA", sep=""))))
    x.cancer <- get(paste(data.sets.list[i], ".CMP", sep=""))[cancer.samples, probe]
    y.cancer <- get(paste(data.sets.list[i], ".CEA", sep=""))[cancer.samples, gene]
    if(is.null(x.cancer) || is.null(y.cancer) || length(x.cancer) == 0 || length(y.cancer) == 0){
      plot.new()
      plot.new()
    }else{
      plot(x.cancer, y.cancer, main=data.sets.list[i], xlim=c(0,1))
      plot(log(x.cancer), log(y.cancer), main=data.sets.list[i])
    }
  }
  # Ask user for input to plot next plot
  cat ("Press [enter] for next plot")
  line <- readline()
}

