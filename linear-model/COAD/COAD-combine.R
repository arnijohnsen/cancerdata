# Load files
load("../Rdata/COAD/calc/COAD-linear-ME.Rdata")
#load("../Rdata/COAD/calc/COAD-nm-expr.Rdata")
load("../Rdata/COAD/calc/COAD-quantiles.Rdata")
#load("../Rdata/COAD/calc/COAD-diff-fold.Rdata")
#load("../Rdata/COAD/calc/COAD-new-vars.Rdata")
load("../Rdata/COAD/calc/COAD-expr-medians.Rdata")

COAD.statistics <- cbind(COAD.linear.ME, COAD.quantiles, COAD.expr.medians)
save(COAD.statistics, file="../Rdata/COAD/calc/COAD-statistics.Rdata")
quit(save="no")
