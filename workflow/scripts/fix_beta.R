library(data.table)
library(optparse)

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--output", type = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

dt_beta <- fread(opt$input)
ind_0 <- lapply(1:10, function(k) {
  which(rowSums(dt_beta[, (2 * k):(2 * k + 1)] == 0) == 2)
})

ind_0_sample <- lapply(1:10, function(k) {
  sort(sample(ind_0[[k]], size = nrow(dt_beta) - 1000))
})

for (k in 1:10) {
  dt_beta[ind_0_sample[[k]], (2 * k):(2 * k + 1)] <- NA
}

fwrite(dt_beta, file=opt$output)
