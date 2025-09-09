library(data.table)

ind_test_all <- split(1:10000, cut(1:10000, 10, labels = FALSE))
for (hsq in c(0.2, 0.6)) {
  for (gcor in c(0.5, 0.75, 1)) {
    for (ncaus in c(50, 200)) {
      pheno_file <- sprintf("../data/pheno/pheno_%s_%s_%s.pheno", hsq, gcor, ncaus)
      dt_pheno <- fread(pheno_file)
      for (rep in 1:10) {
        ind_test <- ind_test_all[[rep]]
        ind_train <- (1:10000)[-ind_test]
        file_train <- sprintf("data/pheno/train/train_%s_%s_%s_rep%s.pheno", hsq, gcor, ncaus, rep)
        file_test <- sprintf("data/pheno/test/test_%s_%s_%s_rep%s.pheno", hsq, gcor, ncaus, rep)
        dt_pheno[ind_train, .SD, .SDcols = c(1, 1, rep + 1)] |>
          fwrite(file = file_train, sep = " ", quote = FALSE, col.names = FALSE)
        dt_pheno[ind_test, .SD, .SDcols = c(1, 1, rep + 1)] |>
          fwrite(file = file_test, sep = " ", quote = FALSE, col.names = FALSE)
      }
    }
  }
}
