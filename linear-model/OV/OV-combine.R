# Load files
load("../Rdata/OV/calc/OV-linear-ME.Rdata")
#load("../Rdata/OV/calc/OV-nm-expr.Rdata")
load("../Rdata/OV/calc/OV-quantiles.Rdata")
#load("../Rdata/OV/calc/OV-diff-fold.Rdata")
#load("../Rdata/OV/calc/OV-new-vars.Rdata")
load("../Rdata/OV/calc/OV-expr-medians.Rdata")

OV.statistics <- cbind(OV.linear.ME, OV.quantiles, OV.expr.medians)
save(OV.statistics, file="../Rdata/OV/calc/OV-statistics.Rdata")
quit(save="no")
