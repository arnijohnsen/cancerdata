# Load files
load("../Rdata/LUAD/calc/LUAD-linear-ME.Rdata")
#load("../Rdata/LUAD/calc/LUAD-nm-expr.Rdata")
load("../Rdata/LUAD/calc/LUAD-quantiles.Rdata")
#load("../Rdata/LUAD/calc/LUAD-diff-fold.Rdata")
#load("../Rdata/LUAD/calc/LUAD-new-vars.Rdata")
load("../Rdata/LUAD/calc/LUAD-expr-medians.Rdata")

LUAD.statistics <- cbind(LUAD.linear.ME, LUAD.quantiles, LUAD.expr.medians)
save(LUAD.statistics, file="../Rdata/LUAD/calc/LUAD-statistics.Rdata")
quit(save="no")
