library(VariantAnnotation)
library(optparse)

option_list <- list(
  make_option("--vcf_file",
    type = "character",
    help = "path to VCF file"
  ),
  make_option("--lanc_file",
    type = "character",
    help = "path to lanc file"
  ),
  make_option("--out_file",
    type = "character",
    help = "file path to output file (should be .vcf)"
  )
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

#### HELPER FUNCTIONS ####

# gets GRanges object for switch points
get_switch_ranges <- function(lanc_file, vcf_file) {
  # get the index of variants where lanc "switches"
  # for at least one subject
  ind_where_switch <- readLines(lanc_file, warn = FALSE)[-1] |>
    strsplit(" ") |>
    unlist() |>
    sub(pattern = "\\:.*", replacement = "") |>
    as.numeric() |>
    unique() |>
    sort()

  # get ranges for all variants in target VCF
  ranges_param <- VariantAnnotation::ScanVcfParam(geno = NA)
  ranges_all <- VariantAnnotation::readVcf(vcf_file, param = ranges_param) |>
    MatrixGenerics::rowRanges()

  # subset to ranges where switch occurs in >= 1 subject
  ranges_where_switch <- ranges_all[ind_where_switch]
  return(ranges_where_switch)
}


# gets ancestry information for switch points
get_anc_switch <- function(lanc_file) {
  # list containing (for each subject) a vector of switch points
  lanc <- readLines(lanc_file, warn = FALSE)[-1] |> strsplit(" ")

  # vector of sample index for each switch
  idx_sample <- rep.int(seq_along(lanc), times = sapply(lanc, length))

  # vector of variant index for each switch
  idx_var <- lanc |>
    unlist() |>
    sub(pattern = "\\:.*", replacement = "") |>
    as.numeric()

  # ancestries for each switch
  an <- lanc |>
    unlist() |>
    sub(pattern = ".*\\:", replacement = "")
  an1 <- substr(an, start = 1, stop = 1) |> as.integer()
  an2 <- substr(an, start = 2, stop = 2) |> as.integer()

  # list of switch point information
  anc_switch <- list(
    idx_sample = idx_sample,
    idx_var = idx_var,
    an1 = an1,
    an2 = an2
  )
  return(anc_switch)
}

upper_fill_ancestry <- function(x) {
  idx_na <- which(is.na(x))
  idx_no_na <- which(!is.na(x))

  upper <- findInterval(idx_na, idx_no_na) + 1

  x[idx_na] <- x[idx_no_na[upper]]
  return(x)
}

make_lanc_vcf <- function(vcf_switch, anc_switch) {
  ## Make AN1, AN2
  an1 <- matrix(NA_integer_,
    nrow = nrow(vcf_switch),
    ncol = ncol(vcf_switch)
  )

  an2 <- matrix(NA_integer_,
    nrow = nrow(vcf_switch),
    ncol = ncol(vcf_switch)
  )

  uniq_idx_var <- anc_switch$idx_var |>
    unique() |>
    sort()

  an1[cbind(
    match(anc_switch$idx_var, uniq_idx_var),
    anc_switch$idx_sample
  )] <- anc_switch$an1

  an2[cbind(
    match(anc_switch$idx_var, uniq_idx_var),
    anc_switch$idx_sample
  )] <- anc_switch$an2

  ## Interpolate ancestry
  an1 <- sapply(seq_len(ncol(an1)), function(j) {
    upper_fill_ancestry(x = an1[, j])
  })

  an2 <- sapply(seq_len(ncol(an2)), function(j) {
    upper_fill_ancestry(x = an2[, j])
  })

  geno(vcf_switch)$AN1 <- an1
  geno(vcf_switch)$AN2 <- an2

  hdr <- VariantAnnotation::header(vcf_switch)
  hdr_an <- data.frame(
    Number = rep("1", 2), Type = rep("Integer", 2),
    Description = c(
      "Ancestry of first haplotype",
      "Ancestry of second haplotype"
    ),
    row.names = c("AN1", "AN2")
  )

  geno(hdr) <- rbind(geno(hdr), hdr_an)
  header(vcf_switch) <- hdr

  return(vcf_switch)
}


#### ANALYSIS

# get ranges for switch points
cat("Getting ranges\n")
ranges_where_switch <- get_switch_ranges(
  lanc_file = opt$lanc_file,
  vcf_file = opt$vcf_file
)

# read in VCF for switch points
switch_param <- VariantAnnotation::ScanVcfParam(which = ranges_where_switch)
cat(sprintf("%s switch points\n", length(ranges_where_switch)))
cat("Reading switch VCF\n")
vcf_switch <- VariantAnnotation::readVcf(
  file = opt$vcf_file,
  param = switch_param
)

# get ancestry information for switch points
cat("Getting ancestry information")
anc_switch <- get_anc_switch(opt$lanc_file)

cat("Making switch VCF\n")
vcf_switch <- make_lanc_vcf(
  vcf_switch = vcf_switch, anc_switch = anc_switch
)


# write
cat("Writing switch VCF\n")
VariantAnnotation::writeVcf(
  obj = vcf_switch,
  filename = opt$out_file,
  index = FALSE
)
