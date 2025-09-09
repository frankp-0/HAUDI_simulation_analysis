samples <- paste0("Sample_", 1:10000)
ind_test_all <- split(1:10000, cut(1:10000, 10, labels = FALSE))
for (rep in 1:10) {
  ind_test <- ind_test_all[[rep]]
  ind_train <- (1:10000)[-ind_test]
  train_file <- sprintf("data/samples/samples_train_rep%s.txt", rep)
  test_file <- sprintf("data/samples/samples_test_rep%s.txt", rep)
  writeLines(samples[ind_train], train_file)
  writeLines(samples[ind_test], test_file)
}
