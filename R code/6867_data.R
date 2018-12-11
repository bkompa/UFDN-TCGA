library(curatedTCGAData)
brcatest <- curatedTCGAData('DLBC', 'RNASeq2GeneNorm', F)
brcadata <- brcatest@ExperimentList@listData[[1]]@assays[[1]]
write.csv(brcadata, '~/Downloads/dlbc.csv')

cancers <- c('ACC' ,'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA' ,'GBM', 'HNSC', 'KICH',
             'KIRC', 'KIRP', 'LAML', 'LGG' ,'LIHC', 'LUAD', 'LUSC' , 'MESO', 'OV', 'PAAD' ,'PCPG',
             'PRAD', 'READ' ,'SARC' , 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC' , 'UCS', 'UVM')

for(subtype in cancers){
    for(assay in c('RNASeq2GeneNorm')){
      #get the name of the assay 
      available_assays <- names(curatedTCGAData(subtype))
      rnaseq_name_index <- which(grepl(assay, available_assays))
      rnaseq_name <- available_assays[rnaseq_name_index]
      print(paste('Downloading .... ', rnaseq_name))
      data <- curatedTCGAData(subtype, assay, F)
      #assay_data <- data@ExperimentList@listData[[1]]@assays[[1]]
      path <- write.csv(data, paste0(c('~','Downloads','TCGA',subtype, 
                                      paste(subtype,'_',assay,'.csv', sep='')), sep='/'))
    }
}

RNA_seq_normalized_cancers <- c()
for(subtype in cancers){
  assays <- names(curatedTCGAData(subtype))
  if(any(grepl('RNASeq2GeneNorm', assays))){
    RNA_seq_normalized_cancers <- c(RNA_seq_normalized_cancers, subtype)
  }
}