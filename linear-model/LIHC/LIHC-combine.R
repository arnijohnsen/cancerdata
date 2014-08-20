# Load files
load("../Rdata/LIHC/calc/LIHC-linear-ME.Rdata")
#load("../Rdata/LIHC/calc/LIHC-nm-expr.Rdata")
load("../Rdata/LIHC/calc/LIHC-quantiles.Rdata")
#load("../Rdata/LIHC/calc/LIHC-diff-fold.Rdata")
#load("../Rdata/LIHC/calc/LIHC-new-vars.Rdata")
load("../Rdata/LIHC/calc/LIHC-expr-medians.Rdata")

LIHC.statistics <- cbind(LIHC.linear.ME, LIHC.quantiles, LIHC.expr.medians)
save(LIHC.statistics, file="../Rdata/LIHC/calc/LIHC-statistics.Rdata")
quit(save="no")
