library(data.table)

frax <- fread("../results/table12_GWAS_Fx.csv")
bmd <- fread("../results/table12_GWAS_BMD.csv")
frax[P.I <= 0.0001, .N]
bmd[P.NI <= 0.0001, .N]
frax[, .N]
bmd[, .N]

frax[P.I <= 0.001, .N]/frax[, .N]
bmd[P.NI <= 0.001, .N]/bmd[, .N]

