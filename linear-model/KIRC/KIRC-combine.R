# Load files
load("../Rdata/KIRC/calc/KIRC-linear-ME.Rdata")
#load("../Rdata/KIRC/calc/KIRC-nm-expr.Rdata")
load("../Rdata/KIRC/calc/KIRC-quantiles.Rdata")
#load("../Rdata/KIRC/calc/KIRC-diff-fold.Rdata")
#load("../Rdata/KIRC/calc/KIRC-new-vars.Rdata")
load("../Rdata/KIRC/calc/KIRC-expr-medians.Rdata")

KIRC.statistics <- cbind(KIRC.linear.ME, KIRC.quantiles, KIRC.expr.medians)
save(KIRC.statistics, file="../Rdata/KIRC/calc/KIRC-statistics.Rdata")
quit(save="no")
