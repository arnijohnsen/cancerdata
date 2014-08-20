# Ask user if he wants abline to be drawn
cat("Do you want abline to be drawn? (y/n)\n")
abline.answer <- readline()
if( abline.answer %in% c("y", "Y", "yes", "Yes", "YES", "true", "TRUE", "t", "T")){
  draw.abline = TRUE
} else {
  draw.abline = FALSE
}

# Ask user which data sets he wants to use
cat("Enter data sets you want to use, seperated by spaces\n")
cat("(Available are BRCA COAD GBM KIRC LIHC LUAD LUSC OV PRAD READ)\n")
data.sets.input <- readline()
data.sets.list <- unlist(strsplit(data.sets.input, " "))

# Remove invalid data sets
available.sets <- c("BRCA", "COAD", "GBM", "KIRC", "LIHC", "LUAD", "LUSC", "OV", "PRAD", "READ")
if(!all(data.sets.list %in% available.sets)){
  cat("Removing invalid sets:", data.sets.list[!(data.sets.list %in% available.sets)], "\n")
}
data.sets.list <- data.sets.list[data.sets.list %in% available.sets]
normal.sets <- c("BRCA", "LIHC", "LUAD", "LUSC", "PRAD", "READ")

if(length(data.sets.list) == 0){
  stop("No valid data set entered\n")
}

# Load missing data sets, previously loaded sets are not reloaded
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

# Create variable to add to 0 values (to avoid log(0) errors)
eps <- 1e-6

# Loop over all probes/genes that will be plotted
if(length(unique.probes.genes) == 0){
  cat("No hits\n")
}else{
  # Set plot layout, need 2 for each selected data set
  layout(matrix(seq(1,2*length(data.sets.list)), ncol=length(data.sets.list)))

  # Loop over list of probes and/or genes
  for(i in 1:length(unique.probes.genes)){
    # Name of probe and gene that will be used
    probe <- unique.probes.genes[[i]][1]
    gene  <- unique.probes.genes[[i]][2]
    cat(probe, gene, "\n")
    for(i in 1:length(data.sets.list)){
      # Select which samples to use and get x (methylation) and y (expression) data 
      # for those samples
      cancer.samples <- intersect(rownames(get(paste(data.sets.list[i], ".CMP", sep=""))), 
                                  rownames(get(paste(data.sets.list[i], ".CEA", sep=""))))
      x.cancer <- get(paste(data.sets.list[i], ".CMP", sep=""))[cancer.samples, probe]
      y.cancer <- get(paste(data.sets.list[i], ".CEA", sep=""))[cancer.samples, gene]
      if(data.sets.list[i] %in% normal.sets){
        # If the data set has normal values, then repeat the process for normal data sets
        normal.samples <- intersect(rownames(get(paste(data.sets.list[i], ".NMP", sep=""))), 
                                  rownames(get(paste(data.sets.list[i], ".NEA", sep=""))))
        x.normal <- get(paste(data.sets.list[i], ".NMP", sep=""))[normal.samples, probe]
        y.normal <- get(paste(data.sets.list[i], ".NEA", sep=""))[normal.samples, gene]
      }
      # If no data or NULL data is found for cancer sets, then an empty plot is created
      # This happens if data for a certain probe is in some (but not all) selected data sets
      if(is.null(x.cancer) || is.null(y.cancer) || length(x.cancer) == 0 || length(y.cancer) == 0){
        plot.new()
        plot.new()
      }else{
        # If the data exists, start by doing a linear regression on available data
        # (either just cancer or normal and cancer)
        if(data.sets.list[i] %in% normal.sets){
          lin <- lm(c(y.normal,y.cancer) ~ c(x.normal, x.cancer))
          sum <- summary(lin)
        }else{
          lin <- lm(y.cancer ~ x.cancer)
          sum <- summary(lin)
        }
        # Create plot for cancer data, with red color
        plot(x.cancer, y.cancer, main=data.sets.list[i], xlim=c(0,1),
             xlab="Methylation", ylab="Expression", 
             sub = paste("R^2=", format(sum$r.squared, digits=3),
                         ", r=", format(sqrt(sum$r.squared)*sign(sum$coefficients[2]), digits=3), 
                         ", p=", format(sum$coefficients[8], digits=3), sep=""),
             col="red", pch=20)
        if(data.sets.list[i] %in% normal.sets){
          # Add normal data, if it exists
          points(x.normal, y.normal, col="blue", pch=20)
        }
        if(draw.abline){
          # Draw abline, if it was requested
          abline(a=sum$coefficients[1], b=sum$coefficients[2])
        }

        # Repeat the process for log transformed data
        if(data.sets.list[i] %in% normal.sets){
          eps <- min(c(y.normal,y.cancer)[c(y.normal,y.cancer) != 0], na.rm=T)
          lin <- lm(log(c(y.normal+eps,y.cancer+eps)) ~ log(c(x.normal+eps, x.cancer+eps)))
          sum <- summary(lin)
        }else{
          eps <- min(y.cancer[y.cancer != 0], na.rm=T)
          lin <- lm(log(y.cancer+eps) ~ log(x.cancer+eps))
          sum <- summary(lin)
        }
        plot(log(x.cancer+eps), log(y.cancer+eps), main=data.sets.list[i], 
             xlab="log Methylation", ylab="log Expression", 
             sub = paste("R^2=", format(sum$r.squared, digits=3),
                         ", r=", format(sqrt(sum$r.squared)*sign(sum$coefficients[2]), digits=3), 
                         ", p=", format(sum$coefficients[8], digits=3), sep=""),
                   col="red", pch=20)
        if(data.sets.list[i] %in% normal.sets){
          points(log(x.normal+eps), log(y.normal+eps), col="blue", pch=20)
        }
        if(draw.abline){
          abline(a=sum$coefficients[1], b=sum$coefficients[2])
        }
      }
    }
    # Add annotation to plot with gene and probe name
    mtext(paste(probe, gene, sep=" "), side=3, line=-1.5, outer=TRUE, font=2)
    # Ask user for input to plot next plot
    cat("Press [enter] for next plot")
    line <- readline()
  }
}
