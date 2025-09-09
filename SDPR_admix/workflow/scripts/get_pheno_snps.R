library(data.table)
library(optparse)

option_list <- list(
  make_option("--beta", type = "character"),
  make_option("--rep", type = "integer"),
  make_option("--out", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

dt <- fread(opt$beta)
dt <- dt[, .SD, .SDcols = c(1, (2 * opt$rep):(2 * opt$rep + 1))]
dt <- dt[!is.na(dt[[2]]), ]

writeLines(dt$snp, con = opt$out)
