load("../Rdata/OV/calc/OV-statistics.Rdata")

names <- rownames(OV.statistics)[OV.statistics$slope < 100 &
                                   OV.statistics$slope > -650 &
                                   OV.statistics$r < -0.25 &
                                   OV.statistics$p.adj < 0.01 &
                                   OV.statistics$NMP.medians < 0.1&
                                   OV.statistics$CMP.99quant > 0.4]

OV.filtered.list <- data.frame(probe = gsub("\\..*", "", names), gene = gsub(".*\\.", "", names))
write.table(OV.filtered.list, file="../Rdata/OV/calc/OV-filtered-list.txt", row.names=F, quote=F)
