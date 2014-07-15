# Ask user which data sets he wants to use
cat("Enter data sets you want to use, seperated by spaces\n")
cat("(Available are BRCA GBM KIRC OV PRAD)\n")
data.sets.input <- readline()
data.sets.list <- unlist(strsplit(data.sets.input, " "))

# Remove unavailable data sets
available.sets <- c("BRCA", "GBM", "KIRC", "OV", "PRAD")
cat("Removing invalid sets:", data.sets.list[!(data.sets.list %in% available.sets)], "\n")
data.sets.list <- data.sets.list[data.sets.list %in% available.sets]
normal.sets <- c("BRCA")

# Load missing data sets
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
  if(data.sets.list[i] %in% normal.sets){
    if(!exists(paste(data.sets.list[i], ".NMP", sep=""))){
      cat("Loading", data.sets.list[i], "normal methylation\n")
      load(paste("../Rdata/", data.sets.list[i], "/data/", data.sets.list[i], "-NMP.Rdata", sep=""))
    }
    if(!exists(paste(data.sets.list[i], ".NEA", sep=""))){
      cat("Loading", data.sets.list[i], "normal expression\n")
      load(paste("../Rdata/", data.sets.list[i], "/data/", data.sets.list[i], "-NEA.Rdata", sep=""))
    }
  }
}

# Ask user for probe or gene to plot
search <- readline("Enter gene or probe name: ")
search <- paste("^", search, "$", sep="")

# Create list which includes all numbers from XXXX.linked.probes.genes
# which are to be used
n.list <- list()

for(i in 1:length(data.sets.list)){
  n.list[[data.sets.list[i]]] <- union(
         grep(search, get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))$probes), 
         grep(search, get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))$genes))
}

# Create list of all probes and genes which are to be used (duplicated allowed)
all.probes <- character(0)
all.genes  <- character(0)
for(i in 1:length(data.sets.list)){
  all.probes <- c(all.probes, 
    as.character(get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))[n.list[[data.sets.list[i]]], ]$probes))
  all.genes <- c(all.genes, 
    as.character(get(paste(data.sets.list[i], ".linked.probes.genes", sep=""))[n.list[[data.sets.list[i]]], ]$genes))
}

# Filter duplicates from all.probes/all.genes and truncate to one list
unique.probes.genes <- strsplit(unique(paste(all.probes, all.genes, sep="-")), "-")

# Set plot layout
layout(matrix(seq(1,2*length(data.sets.list)), ncol=length(data.sets.list)))

# Loop over all probes/genes that will be plotted
for(i in 1:length(unique.probes.genes)){
  # Name of probe and gene that will be used
  probe <- unique.probes.genes[[i]][1]
  gene  <- unique.probes.genes[[i]][2]
  cat(probe, gene, "\n")
  for(i in 1:length(data.sets.list)){
    cancer.samples <- intersect(rownames(get(paste(data.sets.list[i], ".CMP", sep=""))), 
                                rownames(get(paste(data.sets.list[i], ".CEA", sep=""))))
    x.cancer <- get(paste(data.sets.list[i], ".CMP", sep=""))[cancer.samples, probe]
    y.cancer <- get(paste(data.sets.list[i], ".CEA", sep=""))[cancer.samples, gene]
    if(data.sets.list[i] %in% normal.sets){
      normal.samples <- intersect(rownames(get(paste(data.sets.list[i], ".NMP", sep=""))), 
                                rownames(get(paste(data.sets.list[i], ".NEA", sep=""))))
      x.normal <- get(paste(data.sets.list[i], ".NMP", sep=""))[normal.samples, probe]
      y.normal <- get(paste(data.sets.list[i], ".NEA", sep=""))[normal.samples, gene]

    }
    if(is.null(x.cancer) || is.null(y.cancer) || length(x.cancer) == 0 || length(y.cancer) == 0){
      plot.new()
      plot.new()
    }else{
      plot(x.cancer, y.cancer, main=data.sets.list[i], xlim=c(0,1),
           xlab="Methylation", ylab="Expression", 
	   col="red", pch=20)
      if(data.sets.list[i] %in% normal.sets){
        points(x.normal, y.normal, col="blue", pch=20)
      }
      plot(log(x.cancer), log(y.cancer), main=data.sets.list[i], 
           xlab="log Methylation", ylab="log Expression", 
	   col="red", pch=20)
      if(data.sets.list[i] %in% normal.sets){
        points(log(x.normal), log(y.normal), col="blue", pch=20)
      }
    }
  }
  # Ask user for input to plot next plot
  if(i != length(unique.probes.genes)){
    cat("Press [enter] for next plot")
    line <- readline()
  }else{
    cat("Last plot reached\n")
  }
}

