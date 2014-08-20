# Load files
load("../Rdata/LUSC/calc/LUSC-linear-ME.Rdata")
#load("../Rdata/LUSC/calc/LUSC-nm-expr.Rdata")
load("../Rdata/LUSC/calc/LUSC-quantiles.Rdata")
#load("../Rdata/LUSC/calc/LUSC-diff-fold.Rdata")
#load("../Rdata/LUSC/calc/LUSC-new-vars.Rdata")
load("../Rdata/LUSC/calc/LUSC-expr-medians.Rdata")

LUSC.statistics <- cbind(LUSC.linear.ME, LUSC.quantiles, LUSC.expr.medians)
save(LUSC.statistics, file="../Rdata/LUSC/calc/LUSC-statistics.Rdata")
quit(save="no")
