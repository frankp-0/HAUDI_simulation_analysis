library(HAUDI)
library(data.table)

vcf_file <- "./data/CEU-YRI.lanc.vcf.gz"
fbm_pref <- "./data/CEU-YRI"
chunk_size <- 400
anc_names <- c("CEU", "YRI")

fbm_list <- make_fbm(
  vcf_file = vcf_file,
  fbm_pref = fbm_pref,
  chunk_size = chunk_size,
  anc_names = anc_names,
  geno_format = "GT",
  min_ac = 50
)

saveRDS(fbm_list$FBM, file = "./data/CEU-YRI.rds")
fwrite(fbm_list$info, file = "./data/CEU-YRI.fbm_info.txt", sep = "\t")
