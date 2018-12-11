#Integrated non negative matrix factorization 
library(IntNMF)


args <- commandArgs(trailingOnly = TRUE)

K <- args[1]

X_skcm_matrix_no_zero <- readRDS('~/X_skcm_matrix_no_zero.rds')
X_skcm_gbm_matrix_no_zero <- readRDS('~/X_skcm_gbm_matrix_no_zero.rds')

test_NMF_K <- nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=K, n.ini=5)

saveRDS(test_NMF_K, paste('~/test_NMF_',K,'.rds', collapse='', sep=''))
