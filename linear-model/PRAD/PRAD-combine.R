# Load files
load("../Rdata/PRAD/calc/PRAD-linear-ME.Rdata")
#load("../Rdata/PRAD/calc/PRAD-nm-expr.Rdata")
load("../Rdata/PRAD/calc/PRAD-quantiles.Rdata")
#load("../Rdata/PRAD/calc/PRAD-diff-fold.Rdata")
#load("../Rdata/PRAD/calc/PRAD-new-vars.Rdata")
load("../Rdata/PRAD/calc/PRAD-expr-medians.Rdata")

PRAD.statistics <- cbind(PRAD.linear.ME, PRAD.quantiles, PRAD.expr.medians)
save(PRAD.statistics, file="../Rdata/PRAD/calc/PRAD-statistics.Rdata")
quit(save="no")
