load("../Rdata/BRCA/calc/BRCA-statistics.Rdata")

names <- rownames(BRCA.statistics)[BRCA.statistics$slope < 100 &
                                   BRCA.statistics$slope > -650 &
                                   BRCA.statistics$r < -0.25 &
                                   BRCA.statistics$p.adj < 0.01 &
                                   BRCA.statistics$NMP.medians < 0.1 &
                                   BRCA.statistics$CMP.99quant > 0.4]

BRCA.filtered.list <- data.frame(probe = gsub("\\..*", "", names), gene = gsub(".*\\.", "", names))
write.table(BRCA.filtered.list, file="../Rdata/BRCA/calc/BRCA-filtered-list.txt", row.names=F, quote=F)
