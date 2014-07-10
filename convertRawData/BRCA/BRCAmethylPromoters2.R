library(WGCNA)
# Load data from previous script
load(file="../Rdata/BRCA/tmp/normalBetaValues.Rdata")
load(file="../Rdata/BRCA/tmp/cancerBetaValues.Rdata")

# Transpose and filter out bad data
ggNormal <- goodGenes(t(normalBetaValues), verbose=3)
ggCancer <- goodGenes(t(cancerBetaValues), verbose=3)

normalMethyl <- as.data.frame(t(normalBetaValues)[, ggNormal & ggCancer])
cancerMethyl <- as.data.frame(t(cancerBetaValues)[, ggNormal & ggCancer])
probeList <- colnames(normalMethyl)
normalMethylSampList <- rownames(normalMethyl)
cancerMethylSampList <- rownames(cancerMethyl)

# Save to files
save(normalMethyl, file="../Rdata/BRCA/normalMethylPromoters.Rdata")
save(cancerMethyl, file="../Rdata/BRCA/cancerMethylPromoters.Rdata")
save(probeList,    file="../Rdata/BRCA/probeList.Rdata")
save(normalMethylSampList, cancerMethylSampList, file="../Rdata/BRCA/methylSampList.Rdata")
quit(save="no")
