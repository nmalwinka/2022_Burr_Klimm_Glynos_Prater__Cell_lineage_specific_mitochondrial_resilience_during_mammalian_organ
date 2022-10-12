#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()

library("Seurat") 
library("anndata")
library("SCENIC")
library("R2HTML")
library("doMC")
library("BiocParallel")


Project        <- "MBU_spb54_005__SCENIC"
baseDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005"
setwd(baseDir)


message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description', "hgnc_symbol", "hgnc_id"), mart = ensembl)          


message("+-------------------------------------------------------------------------------")
message("+                    extract counts from h5ad files                             ")
message("+-------------------------------------------------------------------------------")

Counts_raw <- anndata::read_h5ad("Input/10xMutatorMiceE85_raw.h5ad")
cell_info <- Counts_raw$obs
Counts <- Counts_raw$X
dim(Counts)
Counts[1:20,1:50]
Counts_t <- WGCNA::transposeBigData(as.matrix(Counts), blocksize = 10000)
Counts_t <- Matrix::Matrix(Counts_t, sparse = TRUE)


exprMat <- GetAssayData(Seurat_obj, slot = "data") # used for SCENIC


message("+-------------------------------------------------------------------------------")
message("+                     create Seurat objectand QC                                ")
message("+-------------------------------------------------------------------------------")

Seurat_obj <- CreateSeuratObject(Counts_t, project = Project, meta.data = cell_info)
head(Seurat_obj@meta.data)

Seurat_obj@meta.data$celltype <- cellInfo[match( rownames(Seurat_obj@meta.data), cellInfo$X),]$celltype
table(Seurat_obj@meta.data$celltype)

# store mitochondrial percentage in object meta data
Seurat_obj <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-", col.name = "percent.mt")

Idents(Seurat_obj) <- Seurat_obj@meta.data$batch
VlnPlot(Seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)




message("--------------------------------------------------------------------------------")
message("+                           cell cycle scoring                                  ")
message("+-------------------------------------------------------------------------------")

cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes <- tolower(s.genes)
g2m.genes <- tolower(g2m.genes)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

s.genes <- firstup(s.genes)
g2m.genes <- firstup(g2m.genes)

Seurat_obj <- CellCycleScoring(Seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



message("--------------------------------------------------------------------------------")
message("+                        process with Seurat                                    ")
message("+-------------------------------------------------------------------------------")

celltype_cols <- c( "amnion"="plum1", "mesodermProgenitors"="darkolivegreen1", "neuralTube"="plum4","mixedMesoderm"="gold2", "neuralCrest"="purple4", "midHindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomiticMesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngealMesoderm"="royalblue","foregut"="violetred3", "extraembryonicMesoderm"="steelblue3", "cardiac"="firebrick" , "somiticMesoderm"= "darkgreen",  "endothelial"="orange1", "midHindgut"="violetred2", "blood"="black")

celltype_cols <- c( "amnion"="plum1", "mesoderm progenitors"="darkolivegreen1", "neural tube"="plum4","mixed mesoderm"="gold2", "neural crest"="purple4", "mid hindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomitic mesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngeal mesoderm"="royalblue","foregut"="violetred3", "extraembryonic mesoderm"="steelblue3", "cardiac"="firebrick" , "somitic mesoderm"= "darkgreen",  "endothelial"="orange1", "mid hindgut"="violetred2", "blood"="black")

mouse_cols <- c( "m.5024C>T"= "firebrick", "WT"="darkolivegreen4", "m.5019A>G"="dodgerblue4")
Seurat_obj@meta.data$mouse <- gsub(  "A5019G", "m.5019A>G", Seurat_obj@meta.data$mouse )
Seurat_obj@meta.data$mouse <- gsub(  "C5024T", "m.5024C>T", Seurat_obj@meta.data$mouse )


Seurat_obj  <- SCTransform(Seurat_obj,  verbose = TRUE, return.only.var.genes = FALSE, vars.to.regress = c("batch") )
Seurat_obj  <- RunPCA(Seurat_obj, npcs = 30, verbose = TRUE) 
Seurat_obj  <- RunUMAP(Seurat_obj, dims = 1:30)
Seurat_obj  <- FindNeighbors(Seurat_obj, dims = 1:30, k.param = 30)
Seurat_obj  <- FindClusters(Seurat_obj, resolution = 1, random.seed = 1) 


message("+-------------------------------------------------------------------------------")
message("+                                  HARMONY                                       ")
message("+-------------------------------------------------------------------------------")

harmony_matrix.su <- RunHarmony(Seurat_obj, group.by.vars = "batch", reduction = "pca", assay.use="SCT")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:20)
harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony", dims = 1:10) 
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
head(harmony_matrix.su@meta.data,2)






message("+-------------------------------------------------------------------------------")
message("+                          initialize SCENIC                                    ")
message("+-------------------------------------------------------------------------------")

scenicOptions <- initializeScenic(datasetTitle= paste0("SCENIC_SCT_batch_regr", "__", Project), org="mgi", dbDir=dbDir, dbs=c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"),  nCores=2) 
scenicOptions@inputDatasetInfo$cellInfo <- cell_info
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 2
scenicOptions@settings$seed <- 123

#saveRDS(scenicOptions, file=paste0("int/scenicOptions_SCT_batch_regr", Project ,".Rds") )

exprMat  <- readRDS( "Input/MBU_spb54_005_Flo_SCENIC__SCT_counts.Rds")




genesKept <- geneFiltering(as.matrix(exprMat), scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(as.matrix(exprMat_filtered), scenicOptions)
exprMat_log <- log2(exprMat_filtered+1)

# genie3 takes forever so do arboreto again!
exportsForArboreto(as.matrix(exprMat_log), scenicOptions, dir = "int")




message("--------------------------------------------------------------------------------")
message("+     RUN ARBORETO on HPC and return with RGN output                            ")
message("+-------------------------------------------------------------------------------")

# dask no longer working so run arboreto with multiprocessing instead.

RGN_linklist <- read.table("int/network_arboreto_numpy_output_SCT_NOT_REGRESSED.tsv", sep = "\t")
colnames(RGN_linklist) <- c("TF", "Target", "weight")
length(unique(RGN_linklist$TF))     #  1499
length(unique(RGN_linklist$Target)) #  14496

corrMat <- readRDS("int/1.2_corrMat.Rds")
dim(corrMat)

RGN_linklist2 <- RGN_linklist[RGN_linklist$Target %in% rownames(corrMat),]
dim(RGN_linklist)
dim(RGN_linklist2)

#saveRDS(RGN_linklist2, 'int/1.4_GENIE3_linkList.Rds')

runSCENIC_1_coexNetwork2modules(scenicOptions)




message("--------------------------------------------------------------------------------")
message("+               Post-Arboreto analysis- step2                                   ")
message("+-------------------------------------------------------------------------------")

nCores <- 2
register(MulticoreParam(workers=nCores), default = TRUE)
register(SnowParam(workers=nCores), default = TRUE)

scenicOptions@settings$db_mcVersion
scenicOptions@settings$dbs
motifAnnot <- getDbAnnotations(scenicOptions)

scenicOptions@settings$nCores <- nCores
runSCENIC_2_createRegulons(scenicOptions) 




message("--------------------------------------------------------------------------------")
message("+                step  3   -->   score cells                                    ")
message("+-------------------------------------------------------------------------------")

scenicOptions@settings$nCores <- 1 # not higher as otherwise will give error in score cells...
runSCENIC_3_scoreCells(scenicOptions, exprMat)

runSCENIC_4_aucell_binarize(scenicOptions)














