library(data.table)
library(optparse)

option_list <- list(
  make_option("--plink_prefix", type = "character"),
  make_option("--msp_out", type = "character"),
  make_option("--map_file", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

read_lanc <- function(lanc_file) {
  lanc <- readLines(lanc_file)[-1]
  dims <- sapply(strsplit(readLines(lanc_file)[1], " "), as.numeric)[, 1]
  anc <- array(dim = c(dims[2], dims[1], 2))
  for (i in 1:length(lanc)) {
    switches <- strsplit(lanc[i], " ")[[1]]
    start <- 1
    for (switch in switches) {
      stop <- as.integer(strsplit(switch, ":")[[1]][[1]])
      lanc1i <- as.integer(substr(strsplit(switch, ":")[[1]][[2]], 1, 1))
      lanc2i <- as.integer(substr(strsplit(switch, ":")[[1]][[2]], 2, 2))
      anc[i, start:stop, 1] <- lanc1i
      anc[i, start:stop, 2] <- lanc2i
      start <- stop + 1
    }
  }
  anc
}

get_cm <- function(pvar_file, map_file) {
  dt_pvar <- fread(pvar_file)
  pos <- dt_pvar$POS
  map <- fread(map_file)
  map <- map[
    gsub(pattern = "chr", replacement = "", map[[1]]) %in% unique(dt_pvar$`#CHROM`),
  ]

  cm <- map[[3]][match(pos, map[[4]])]
  for (i in which(is.na(cm))) {
    posi <- pos[i]
    idx1 <- max(which(posi >= map[[4]]))
    idx2 <- min(which(posi <= map[[4]]))
    cmi <- map$V3[idx1] + (posi - map[[4]][idx1]) / (map[[4]][idx2] - map[[4]][idx1]) * (map[[3]][idx2] - map[[3]][idx1])
    cm[i] <- cmi
  }
  cm
}

convert_lanc_to_msp <- function(plink_prefix, map_file, msp_out) {
  lines <- readLines(sprintf("%s.lanc", plink_prefix))
  header <- scan(text = lines[1], what = integer(), quiet = TRUE)
  n_indiv <- header[2]

  samples <- fread(sprintf("%s.psam", plink_prefix))[[1]]

  hdr_codes <- "#Subpopulation order/codes: AFR=0\tEUR=1"
  write(file = msp_out, hdr_codes)
  hdr_cols <- paste0(
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\t",
    paste(paste0(rep(samples, each = 2), c(".0", ".1")), collapse = "\t")
  )
  write(file = msp_out, hdr_cols, append = TRUE)
  pvar <- fread(sprintf("%s.pvar", plink_prefix))
  pvar$`#CHROM` <- as.character(pvar$`#CHROM`)
  map <- fread(map_file)
  cm <- get_cm(sprintf("%s.pvar", plink_prefix), map_file)
  anc <- read_lanc(sprintf("%s.lanc", plink_prefix))
  is_switch <- c(FALSE, apply(anc[, -1, ] != anc[, -dim(anc)[2], ], 2, function(a) sum(a) != 0))
  is_switch[length(is_switch)] <- TRUE
  idx_switch <- which(is_switch)

  lines <- vector(length = length(idx_switch))
  i <- 1
  is_first <- TRUE
  for (eidx in idx_switch) {
    if (is_first) {
      sidx <- 1
      is_first <- FALSE
    }
    info_string <- c(
      pvar$`#CHROM`[sidx],
      pvar$POS[sidx],
      pvar$POS[eidx],
      round(cm[sidx], 2),
      round(cm[eidx], 2),
      eidx - sidx
    ) |> paste(collapse = "\t")
    anc_string <- paste(t(anc[, sidx, ]), collapse = "\t")
    line <- paste(info_string, anc_string, sep = "\t")
    lines[i] <- line
    i <- i + 1
    sidx <- eidx
  }
  out_con <- file(msp_out, "a")
  writeLines(lines, out_con)
  close(out_con)
}

convert_lanc_to_msp(opt$plink_prefix, opt$map_file, opt$msp_out)
