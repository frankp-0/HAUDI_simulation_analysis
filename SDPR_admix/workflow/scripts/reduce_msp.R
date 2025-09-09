library(data.table)
library(optparse)

option_list <- list(
  make_option("--msp_in", type = "character"),
  make_option("--msp_out", type = "character"),
  make_option("--samples", type = "character"),
  make_option("--snps_in", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

reduce_msp <- function(msp_in, msp_out, snps_in, samples) {
  snps <- readLines(snps_in)
  samples <- readLines(samples)
  pos_keep <- sapply(strsplit(snps, ":"), function(x) x[[2]]) |>
    as.integer() |>
    sort()

  hdr1 <- readLines(msp_in, n = 1)
  hdr2 <- c("#chm", "spos", "epos", "sgpos", "egpos", "n snps")
  hdr2 <- paste(c(hdr2, samples), collapse = "\t")
  hdr <- c(hdr1, hdr2)

  msp <- fread(msp_in, skip = 1, sep = "\t")
  colkeep <- c(colnames(msp)[1:6], paste0(samples, rep(c(".0", ".1"), each = length(samples))))
  msp <- msp[, ..colkeep]

  vcf_start <- pos_keep[1]
  vcf_stop <- pos_keep[length(pos_keep)]
  if (!(vcf_start %in% msp$spos)) {
    pos_keep <- c(max(msp$spos[msp$spos < vcf_start]), pos_keep)
  }
  if (!(vcf_stop %in% msp$epos)) {
    pos_keep <- c(pos_keep, min(msp$epos[msp$epos > vcf_stop]))
  }
  idx_s <- which(msp$spos %in% pos_keep)
  idx_e <- which(msp$epos %in% pos_keep)
  if (max(idx_e) != nrow(msp)) {
    idx_s <- idx_s[-length(idx_s)]
  }
  if (min(idx_s) != 1) {
    idx_e <- idx_e[-1]
  }
  epos <- msp[idx_e, ]$epos
  egpos <- msp[idx_e, ]$egpos
  msp <- msp[idx_s, ]
  msp$epos <- epos
  msp$egpos <- egpos
  msp$spos[1] <- vcf_start
  msp$epos[length(msp$epos)] <- vcf_stop

  writeLines(hdr, msp_out)
  fwrite(msp, file = msp_out, sep = "\t", append = TRUE, col.names = F)
}

reduce_msp(opt$msp_in, opt$msp_out, opt$snps_in, opt$samples)
