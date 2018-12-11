#Integrated non negative matrix factorization 
library(RcppCNPy)
library(tidyverse)
library(curatedTCGAData)
library(HGNChelper)
library(plyr)
library(IntNMF)
library(edgeR)
library(enrichR)
library(gplots)
library(heatmap3)

X_skcm <- read_csv('~/Python/6_867_project/X_sckm_counts.csv', col_names = FALSE)

X_skcm_to_gbm <- read_csv('~/Python/6_867_project/X_sckm_to_gbm_counts.csv', col_names = FALSE)

X_skcm_to_gbm_25 <- read_csv('~/Python/6_867_project/X_sckm_to_gbm_25_counts.csv', col_names = FALSE)
X_skcm_to_gbm_50 <- read_csv('~/Python/6_867_project/X_sckm_to_gbm_50_counts.csv', col_names = FALSE)
X_skcm_to_gbm_75 <- read_csv('~/Python/6_867_project/X_sckm_to_gbm_75_counts.csv', col_names = FALSE)

X_skcm_to_gbm[X_skcm_to_gbm<0] <- 0

X_skcm_to_gbm_25[X_skcm_to_gbm_25<0] <- 0
X_skcm_to_gbm_50[X_skcm_to_gbm_50<0] <- 0
X_skcm_to_gbm_75[X_skcm_to_gbm_75<0] <- 0

#raw data 
skcm <- curatedTCGAData('SKCM', 'RNASeq2GeneNorm',F)
gene_names <- skcm@ExperimentList@listData$`SKCM_RNASeq2GeneNorm-20160128`@NAMES
gene_names <- tibble(Gene=gene_names)

#write gene names out to get HGNC online
write.csv(gene_names, '~/R/gene_names.csv', sep=',', row.names = F, col.names = F)

#read in table I made with https://www.genenames.org/tools/multi-symbol-checker/
gene_info <- read_csv('~/R/hgnc-symbol-check.csv')
names(gene_info) <- c('Gene', 'Status', 'Approved', 'ApprovedName', 'HGNC', 'Location')

look_up <- function(gene){
  index <- which(gene_info$Gene==gene)[1]
  return(gene_info$HGNC[index])
}

hgnc_names <- unlist(lapply(gene_names$Gene, function(gene) look_up(gene)))

names(X_skcm) <- hgnc_names
names(X_skcm_to_gbm) <- hgnc_names
names(X_skcm_to_gbm_25) <- hgnc_names
names(X_skcm_to_gbm_50) <- hgnc_names
names(X_skcm_to_gbm_75) <- hgnc_names

##### edgeR #####
X <- rbind(X_skcm, X_skcm_to_gbm)
X <- t(X)
colnames(X) <- c(unlist(lapply(seq(1:95), 
                function(i) paste(i,'SKCM', collapse='', sep=''))), 
                unlist(lapply(seq(1:95), 
                function(i) paste(i,'GBM', collapse='', sep=''))))
rownames(X) <- gene_names$Gene

dgList <- DGEList(counts=X, genes=rownames(X))

countsPerMillion <- cpm(dgList)
summary(countsPerMillion)

countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
dgList <- calcNormFactors(dgList, method="TMM")

#hgnc update
hgnc_names <- hgnc_names[keep]
gene_names <- gene_names$Gene[keep]

sampleType <- rep(c("SKCM", "GBM"), each=dim(X)[2]/2)

designMat <- model.matrix(~sampleType)

dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=2)

edgeR_result <- topTags(lrt)

deGenes <- as.logical(decideTestsDGE(lrt, p=.05/20501))
deGenesNames <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenesNames, xlab='Average log CountsPerMillion',
          ylab='Average Log Fold Change', main='Differential Expression SCKM to GBM')
abline(h=c(-1, 1), col=3)


##### NNMF #####

X_skcm_matrix <- as.matrix(X_skcm[,keep]) #removed non-keepers
X_skcm_gbm_matrix <- as.matrix(X_skcm_to_gbm[,keep])

X_skcm_gbm_25_matrix <- as.matrix(X_skcm_to_gbm_25[,keep])
X_skcm_gbm_50_matrix <- as.matrix(X_skcm_to_gbm_50[,keep])
X_skcm_gbm_75_matrix <- as.matrix(X_skcm_to_gbm_75[,keep])

#HGNC update 
hgnc_names <- hgnc_names[deGenes]
gene_names <- gene_names[deGenes]

#remove non DE genes 
X_skcm_matrix_DE <- X_skcm_matrix[,deGenes]
X_skcm_gbm_matrix_DE <- X_skcm_gbm_matrix[,deGenes]

X_skcm_gbm_25_matrix_DE <- X_skcm_gbm_25_matrix[,deGenes]
X_skcm_gbm_50_matrix_DE <- X_skcm_gbm_50_matrix[,deGenes]
X_skcm_gbm_75_matrix_DE <- X_skcm_gbm_75_matrix[,deGenes]

#remove zero columns 
zero_columns <- c(which(colSums(X_skcm_matrix_DE)==0), which(colSums(X_skcm_gbm_matrix_DE)==0),  
                  which(colSums(X_skcm_gbm_25_matrix_DE)==0),  which(colSums(X_skcm_gbm_50_matrix_DE)==0), 
                  which(colSums(X_skcm_gbm_75_matrix_DE)==0))


X_skcm_matrix_no_zero <- X_skcm_matrix_DE[,-zero_columns]
X_skcm_gbm_matrix_no_zero <- X_skcm_gbm_matrix_DE[,-zero_columns]

X_skcm_gbm_25_matrix_no_zero <- X_skcm_gbm_25_matrix_DE[,-zero_columns]
X_skcm_gbm_50_matrix_no_zero <- X_skcm_gbm_50_matrix_DE[,-zero_columns]
X_skcm_gbm_75_matrix_no_zero <- X_skcm_gbm_75_matrix_DE[,-zero_columns]

saveRDS(X_skcm_gbm_25_matrix_no_zero, 'X_skcm_gbm_25_matrix_no_zero.rds')
saveRDS(X_skcm_gbm_50_matrix_no_zero, 'X_skcm_gbm_50_matrix_no_zero.rds')
saveRDS(X_skcm_gbm_75_matrix_no_zero, 'X_skcm_gbm_75_matrix_no_zero.rds')
saveRDS(X_skcm_matrix_no_zero, 'X_skcm_matrix_no_zero.rds')
saveRDS(X_skcm_gbm_matrix_no_zero, 'X_skcm_gbm_matrix_no_zero.rds')

hgnc_names <- hgnc_names[-zero_columns]
gene_names <- gene_names[-zero_columns]


X_skcm_matrix_no_zero <- t(X_skcm_matrix_no_zero)
X_skcm_gbm_matrix_no_zero <- t(X_skcm_gbm_matrix_no_zero)


print(dim(X_skcm_matrix_no_zero))
# system.time(nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=2, n.ini=1))
# 
#test_NMF2 <- nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=2, n.ini=1)
# test_NMF3 <- nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=3, n.ini=10)
# test_NMF4 <- nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=4, n.ini=10)
# test_NMF5 <- nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=5, n.ini=10)
# test_NMF6 <- nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=6, n.ini=10)
# test_NMF10 <- nmf.mnnals(dat=list(X_skcm_matrix_no_zero, X_skcm_gbm_matrix_no_zero), k=10, n.ini=10)


test_NMF2 <- readRDS('~/Downloads/5fracs_test_NMF_2.rds')
test_NMF3 <- readRDS('~/Downloads/5fracs_test_NMF_3.rds')
test_NMF6 <- readRDS('~/Downloads/5fracs_test_NMF_6.rds')
test_NMF7 <- readRDS('~/Downloads/5fracs_test_NMF_7.rds')
test_NMF8 <- readRDS('~/Downloads/5fracs_test_NMF_8.rds')
test_NMF9 <- readRDS('~/Downloads/5fracs_test_NMF_9.rds')
test_NMF10 <- readRDS('~/Downloads/5fracs_test_NMF_10.rds')

test_NMF20 <- readRDS('~/Downloads/5fracs_test_NMF_20.rds')
test_NMF30 <- readRDS('~/Downloads/5fracs_test_NMF_30.rds')
test_NMF40 <- readRDS('~/Downloads/5fracs_test_NMF_40.rds')
test_NMF50 <- readRDS('~/Downloads/5fracs_test_NMF_50.rds')
test_NMF60 <- readRDS('~/Downloads/5fracs_test_NMF_60.rds')


min.f.WH <- list(test_NMF2$min.f.WH, test_NMF3$min.f.WH, 
              test_NMF6$min.f.WH, test_NMF7$min.f.WH,test_NMF8$min.f.WH,
              test_NMF9$min.f.WH,test_NMF10$min.f.WH,  test_NMF20$min.f.WH,
              test_NMF30$min.f.WH, test_NMF40$min.f.WH,
              test_NMF50$min.f.WH, test_NMF60$min.f.WH)


WH_lengths <- unlist(lapply(min.f.WH, function(m) length(m)))


K <- rep(c(2,3,6,7,8,9,10,20, 30, 40,50, 60), WH_lengths)
WH_df <- data.frame(WH=unlist(min.f.WH), k=K)
ggplot(data=WH_df, aes(x=K, y=WH, color=factor(K)))+geom_point()+
  stat_summary(fun.y="mean", geom="line", aes(group=1), color='black')+theme_bw()+
  ylab('Reconstruction Error')+
  xlab('Number of Metagenes')+
  guides(color=FALSE)+
  ggtitle('Reconstruction Error for SKCM to GBM')

#make heatmap of W 

NMF <- test_NMF60
W <- NMF$W
W_normalized <- prop.table(W,2)

heatmapW <- heatmap3(W, Rowv = NA, Colv=NA, verbose=T, useRaster = T, main='Metagenes Heatmap', xlab='Metagene', ylab='Gene')

heatmapW_normalized <- heatmap3(W_normalized, scale='row', Rowv = NA, Colv=NA, verbose=T, useRaster = T, main='Metagenes Heatmap', xlab='Metagene', ylab='Gene')

H <- NMF$H

heatmapH1 <- heatmap3(t(H[[1]]), scale='row', Rowv = NA, Colv=NA, verbose=T, useRaster = T, main='SKCM Heatmap', xlab='Metagene', ylab='Sample')



##### GO terms ##### 

dbs <- listEnrichrDbs()

dbs <- c("GO_Biological_Process_2015")

#remove no name matches 
NA_names <- is.na(rownames(test_NMF2$W))


iNNMF_df <- data.frame(genes=rownames(test_NMF2$W)[-NA_names], scores=test_NMF2$W[-NA_names,1])


enriched <- enrichr(iNNMF_df, dbs)


