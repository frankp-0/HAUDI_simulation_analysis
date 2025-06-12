library(optparse)
library(HAUDI)
library(data.table)

option_list <- list(
  make_option("--pheno_file", type = "character"),
  make_option("--beta_file", type = "character"),
  make_option("--results_file", type = "character"),
  make_option("--model_file", type = "character"),
  make_option("--rep", type = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

fbm <- readRDS("./data/CEU-YRI.rds")
fbm_info <- data.table::fread("./data/CEU-YRI.fbm_info.txt")
y <- data.table::fread(opt$pheno_file)[[opt$rep + 1]]

dt_beta <- data.table::fread(opt$beta_file)[,
  .SD,
  .SDcols = c(1, (2 * opt$rep):(2 * opt$rep + 1))
]

dt_beta <- na.omit(dt_beta)
snps <- dt_beta$snp
snps <- intersect(snps, fbm_info$rsid)

ind_test_all <- split(1:10000, cut(1:10000, 10, labels = FALSE))
ind_test <- ind_test_all[[opt$rep]]
ind_train <- (1:10000)[-ind_test]

time_run <- system.time({
  model <- HAUDI::lasso(
    fbm_obj = fbm, fbm_info = fbm_info,
    y = y, ind_train = ind_train, K = 5,
    family = "gaussian",
    snps = snps
  )
})

y_predicted <- predict(model, fbm)
r2_test <- cor(y_predicted[ind_test], y[ind_test])^2
mse_test <- mean((y[ind_test] - y_predicted[ind_test])^2)
loss_validation <- -summary(model)$validation_loss


param_vec <- strsplit(opt$pheno_file, split = "_")[[1]][2:4]
param_vec[3] <- strsplit(param_vec[3], ".pheno")[[1]]

df_result <- data.frame(
  heritability = param_vec[1],
  gen_cor = param_vec[2],
  n_causal = param_vec[3],
  rep = opt$rep,
  r2_test = r2_test,
  loss_val = loss_validation,
  mse_test = mse_test,
  gamma = NA,
  time_run = time_run[3]
)

write.table(
  x = df_result, file = opt$results_file, sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = TRUE
)

saveRDS(model, file = opt$model_file)
