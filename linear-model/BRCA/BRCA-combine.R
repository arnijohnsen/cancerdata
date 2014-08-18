# Load files
load("../Rdata/BRCA/calc/BRCA-linear-ME.Rdata")
#load("../Rdata/BRCA/calc/BRCA-nm-expr.Rdata")
load("../Rdata/BRCA/calc/BRCA-quantiles.Rdata")
#load("../Rdata/BRCA/calc/BRCA-diff-fold.Rdata")
#load("../Rdata/BRCA/calc/BRCA-new-vars.Rdata")
load("../Rdata/BRCA/calc/BRCA-expr-medians.Rdata")

BRCA.statistics <- cbind(BRCA.linear.ME, BRCA.quantiles, BRCA.expr.medians)
save(BRCA.statistics, file="../Rdata/BRCA/calc/BRCA-statistics.Rdata")
quit(save="no")
