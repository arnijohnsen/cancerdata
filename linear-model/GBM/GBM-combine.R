# Load files
load("../Rdata/GBM/calc/GBM-linear-ME.Rdata")
#load("../Rdata/GBM/calc/GBM-nm-expr.Rdata")
load("../Rdata/GBM/calc/GBM-quantiles.Rdata")
#load("../Rdata/GBM/calc/GBM-diff-fold.Rdata")
#load("../Rdata/GBM/calc/GBM-new-vars.Rdata")
load("../Rdata/GBM/calc/GBM-expr-medians.Rdata")

GBM.statistics <- cbind(GBM.linear.ME, GBM.quantiles, GBM.expr.medians)
save(GBM.statistics, file="../Rdata/GBM/calc/GBM-statistics.Rdata")
quit(save="no")
