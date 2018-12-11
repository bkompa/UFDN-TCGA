


library(topGO)

W <- readRDS('~/6_867_project/20_W.rds')
geneMap <- readRDS('~/6_867_project/gene2GO.rds')
genes <- readRDS('~/6_867_project/gene_names.rds')

geneList <- W[,1]
names(geneList) <-genes


GOdata <- new("topGOdata", 
              ontology='MF', 
              allGenes=geneList, 
              annot=annFUN.gene2GO, 
              gene2GO=geneMap, 
              geneSel=function(allScore){TRUE}, 
              nodeSize=5)


test.stat <- new("classicScore", testStatistic=GOKSTest, name="KS tests")
resultKS <- getSigGroups(GOdata, test.stat)
allRes <- GenTable(GOdata, classicKS=resultKS, topNodes=10)