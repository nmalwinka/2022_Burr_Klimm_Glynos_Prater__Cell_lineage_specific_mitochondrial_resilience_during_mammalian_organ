#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()


##############################################################################################################
#                                                                                                            #
#   Project: MBU_spb54_005_MouseEmbryo                                                                       #
#   Malwina Prater (mn367@cam.ac.uk), 2022                                                                   #
#   MRC MBU, University of Cambridge                                                                         #
#   Script: scRNA-seq mouse dataset - integration with Tabula Muris                                          # 
#                                                                                                            #
##############################################################################################################


message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Seurat")
library("dplyr")
library("biomaRt")
library("ComplexHeatmap")
})
theme_set(theme_cowplot())


Project        <- "MBU_spb54_005__CELL_REV_adult_mouse_tabula_muris_"
baseDir        <- "/Users/xxx/Documents/xxx/xxx/xxx" # replace with your path
resDir        <- "/Users/xxx/Documents/xxx/xxx/xxx/GSE132042_tabula_muris" # replace with your path
setwd(resDir)


# get gene annotation from ensembl:
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', "transcript_length"), mart = ensembl)



message("+-------------------------------------------------------------------------------")
message("+                         Load in mitocarta                                     ")
message("+-------------------------------------------------------------------------------")

# download Mitocarta v3 first from: https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways

MitoCarta3 <- read.csv("/Users/mn367/Documents/MBU-Projects/Databases/Mouse.MitoCarta3.0_summarised.csv")
MitoCarta3_pathways <- read.csv("/Users/mn367/Documents/MBU-Projects/Databases/Mouse.MitoCarta3.0_pathways.csv")
MitoCarta3_genes_mouse <- unique(MitoCarta3$Symbol)

MitoCarta_oxPhos <- MitoCarta3[grep("OXPHOS", MitoCarta3$MitoCarta3.0_MitoPathways), c("Symbol","EnsemblGeneID","MitoCarta3.0_MitoPathways","Synonyms","Description","MitoCarta3.0_SubMitoLocalization","Tissues")]

MitoCarta_oxPhos$ComplexOxphos <- gsub( "OXPHOS > " , "", MitoCarta_oxPhos$MitoCarta3.0_MitoPathways)
MitoCarta_oxPhos$ComplexOxphos <- gsub( " \\| OXPHOS subunits.*" , "", MitoCarta_oxPhos$ComplexOxphos)
MitoCarta_oxPhos$ComplexOxphos <- gsub( " \\| Metabolism.*" , "", MitoCarta_oxPhos$ComplexOxphos)
MitoCarta_oxPhos$ComplexOxphos <- gsub( " > .*" , "", MitoCarta_oxPhos$ComplexOxphos)
oxPhos_genes <- unique(MitoCarta_oxPhos$Symbol)




message("+-------------------------------------------------------------------------------")
message("+                              ok lets try bulk:                                ")
message("+-------------------------------------------------------------------------------")

meta <- read.csv("Input/GSE132040_MACA_Bulk_metadata.csv", header = 1)
colnames(meta) <- c("Sample_name", "title" ,  "sample_classification" ,"organism" , "age", "dev_stage","sex" ,  "molecule", "desc" ,"processed.data.file" , "raw.file", "BioSample","Instrument.Model" )
meta$tissue <- gsub( "_[0-9]{2}", "", meta$sample_classification)
meta$tissue <- gsub( "_[0-9]", "", meta$tissue)
meta$tissue <- gsub( "[0-9]", "", meta$tissue)
table(meta$age, meta$tissue)
idx <- grep("NA", meta$age)
meta <- meta[-idx,]

meta_final <- meta[meta$age %in% c("9" , "1"),]#
table(meta_final$sex, meta_final$tissue)

meta_final <- meta_final[,c("Sample_name", "age","dev_stage","sex" ,"sample_classification", "tissue")]
unique(meta_final$tissue)
meta_final$tissue_Lung <- ifelse(meta_final$tissue == "Lung", "Lung" , "other")
meta_final$tissue_Heart <- ifelse(meta_final$tissue == "Heart", "Heart" , "other")
meta_final$tissue_Kidney <- ifelse(meta_final$tissue == "Kidney", "Kidney" , "other")
meta_final$tissue_Brain <- ifelse(meta_final$tissue == "Brain", "Brain" , "other")
meta_final$tissue_Bone <- ifelse(meta_final$tissue == "Bone", "Bone" , "other")
meta_final$tissue_Limb_Muscle <- ifelse(meta_final$tissue == "Limb_Muscle", "Limb_Muscle" , "other")
meta_final$tissue_Spleen <- ifelse(meta_final$tissue == "Spleen", "Spleen" , "other")
meta_final$tissue_Pancreas <- ifelse(meta_final$tissue == "Pancreas", "Pancreas" , "other")
meta_final$tissue_Liver <- ifelse(meta_final$tissue == "Liver", "Liver" , "other")

Counts <- read.csv("bulk/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv", header = 1, row.names = 1)
idx <- grep("^_", rownames(Counts))
Counts <- Counts[-idx,]

colnames(Counts) <- gsub(".gencode.vM19" , "", colnames(Counts))
colnames(Counts) %in% meta$Sample_name

Counts <- Counts[,colnames(Counts) %in% meta$Sample_name]
Counts_final <- Counts[,colnames(Counts) %in% meta_final$Sample_name]

all_genes <- data.frame(genes=rownames(Counts_final))



meta2 <- meta[meta$age %in% c("1","9"),]
Counts2 <- Counts[rowMeans(Counts) > 10,]
Counts2 <- Counts2[,colnames(Counts2) %in% meta2$Sample_name]


message("+-------------------------------------------------------------------------------")
message("+ Create ddsHTSeq object")
message("+-------------------------------------------------------------------------------")

ddsHTSeq <- DESeqDataSetFromMatrix(Counts2, meta2, design = ~sex + tissue  )

colData(ddsHTSeq)
dds <- DESeq(ddsHTSeq)
resultsNames(dds)
colData(dds)

rld <- rlogTransformation(dds)


message("+-------------------------------------------------------------------------------")
message("+                   Heatmap with isoforms                  ")
message("+-------------------------------------------------------------------------------")

# load in expression matrix from Mouse Embryo!!
Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))
Seurat_obj_WT <- subset(Seurat_obj, mouse == "WT")
rm(Seurat_obj) # we are using only WT data so full dataset not needed.

celltype_cols <- c( "amnion"="plum1", "mesoderm progenitors"="darkolivegreen1", "neural tube"="plum4","mixed mesoderm"="gold2", "neural crest"="purple4", "mid hindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomitic mesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngeal mesoderm"="royalblue","foregut"="violetred3", "extraembryonic mesoderm"="steelblue3", "cardiac"="firebrick" , "somitic mesoderm"= "darkgreen",  "endothelial"="orange1", "mid hindgut"="violetred2", "blood"="black")


Seurat_obj_WT@meta.data$celltype2 <- Seurat_obj_WT@meta.data$celltype
Seurat_obj_WT@meta.data$celltype2 <- gsub( "_", " ", Seurat_obj_WT@meta.data$celltype2)
Seurat_obj_WT@meta.data$celltype2 <- gsub( "Hind", " hind", Seurat_obj_WT@meta.data$celltype2)
Seurat_obj_WT@meta.data$celltype2 <- gsub( "Crest", " crest", Seurat_obj_WT@meta.data$celltype2)
Seurat_obj_WT@meta.data$celltype2 <- gsub( "Progenitors", " progenitors", Seurat_obj_WT@meta.data$celltype2)
Seurat_obj_WT@meta.data$celltype2 <- gsub( "Tube", " tube", Seurat_obj_WT@meta.data$celltype2)
Seurat_obj_WT@meta.data$celltype2 <- gsub( "Mesoderm", " mesoderm", Seurat_obj_WT@meta.data$celltype2)
unique(Seurat_obj_WT@meta.data$celltype2)


Idents(Seurat_obj_WT) <- Seurat_obj_WT@meta.data$celltype2
expr_mat <- AverageExpression(Seurat_obj_WT, assay = "SCT", slot = "scale.data")[[1]]
expr_mat2 <- expr_mat
expr_mat2[1:5,1:5]
summary(expr_mat2)


dds <- readRDS("dds_bulk_1m_9m_all_tissues_Tabula_Muris.Rds")
rld <- vst(dds)
#rld <- readRDS("rld_9m_all_tissues.Rds")



# select genes for plotting on the heatmaps. Select from 2 options:

selected_genes <-  c("Cox6a2", "Cox7a1","Cox6a1", "Cox7a2","Cox8a","Cox8b","Cox8c", "Cox4i1", "Cox4i2", "Ndufa4","Ndufa4l2") # just selected isoforms

selected_genes <- oxPhos_genes # all oxphos genes





mat <- assay(rld)[rownames(assay(rld)) %in% selected_genes ,]

mat_embryo <- expr_mat2[rownames(expr_mat2) %in% selected_genes,]

dim(mat_embryo)
mat_embryo <- mat_embryo[match(rownames(mat), rownames(mat_embryo)),]
rownames(mat) == rownames(mat_embryo)
rownames(mat_embryo) <- rownames(mat)
rownames(mat) ==  rownames(mat_embryo)

mat_anno <- data.frame(sample= colnames(mat)) 
mat_anno$tissue <- meta[match( mat_anno$sample, meta$Sample_name),]$tissue
mat_anno$age <- meta[match( mat_anno$sample, meta$Sample_name),]$age

anno_embryo <- data.frame(celltype= colnames(mat_embryo)) 
anno_embryo$celltype <- gsub( "Hind", " hind", anno_embryo$celltype)
anno_embryo$celltype <- gsub( "Mesoderm", " mesoderm", anno_embryo$celltype)
anno_embryo$celltype <- gsub( "Crest", " crest", anno_embryo$celltype)
anno_embryo$celltype <- gsub( "Progenitors", " progenitors", anno_embryo$celltype)
anno_embryo$celltype <- gsub( "Tube", " tube", anno_embryo$celltype)
col_split <- anno_embryo$celltype

ha = HeatmapAnnotation(CellType = col_split, col = list(CellType = celltype_cols))
ha1 = HeatmapAnnotation(Age = mat_anno$age, col = list(Age = c("1"= "green4", "9"="hotpink1"))) 
ht_adult <- ComplexHeatmap::Heatmap(matrix = mat, name= "Expression in adult", show_row_names = TRUE, show_column_names = FALSE, 
                                    cluster_rows = TRUE, cluster_columns = TRUE, column_dend_side = "bottom",
                                    column_split = mat_anno$tissue, top_annotation = ha1,  row_title_rot = 180, column_title_rot = 90,
                                    height = nrow(mat)*unit(5, "mm"))  



f1 = circlize::colorRamp2( c(-1.5, 0, 1.5), c("blue", "white","red"), space = "RGB") 
ht_embryo <- ComplexHeatmap::Heatmap(matrix = mat_embryo, name= "Expression in embryo", column_dend_side = "bottom", 
                                     show_row_names = TRUE, show_column_names = TRUE, 
                                     cluster_rows = FALSE, cluster_columns = TRUE, 
                                     column_names_side = "top", col = f1,
                                     row_title_rot = 90, column_title_rot = 90, top_annotation = ha,
                                     height = nrow(mat)*unit(6, "mm"), width = unit(6, "cm"),  row_names_gp = gpar(fontface = "italic")) 

ht_list <- ht_adult + ht_embryo


#pdf(paste("ComplexHeatmap",  "GSE132040", "TabulaMuris_adult_tissues_ALLisoform_expression_clust_vs_EMBRYO", ".pdf", sep="_"), onefile=FALSE, width=15, height=40)
pdf(paste("ComplexHeatmap",  "GSE132040", "TabulaMuris_adult_tissues_DIFFisoform_expression_clust_vs_EMBRYO", "v2.pdf", sep="_"), onefile=FALSE, width=7, height=13)
par(bg=NA)
draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()










