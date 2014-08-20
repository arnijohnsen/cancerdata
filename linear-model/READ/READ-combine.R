# Load files
load("../Rdata/READ/calc/READ-linear-ME.Rdata")
#load("../Rdata/READ/calc/READ-nm-expr.Rdata")
load("../Rdata/READ/calc/READ-quantiles.Rdata")
#load("../Rdata/READ/calc/READ-diff-fold.Rdata")
#load("../Rdata/READ/calc/READ-new-vars.Rdata")
load("../Rdata/READ/calc/READ-expr-medians.Rdata")

READ.statistics <- cbind(READ.linear.ME, READ.quantiles, READ.expr.medians)
save(READ.statistics, file="../Rdata/READ/calc/READ-statistics.Rdata")
quit(save="no")
