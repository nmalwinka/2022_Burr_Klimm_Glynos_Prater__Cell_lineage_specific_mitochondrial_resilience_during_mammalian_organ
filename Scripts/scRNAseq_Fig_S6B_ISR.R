#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()


library("ggplot2")
library("ggrepel")
library("cowplot")
library("Seurat")
library("dplyr")
theme_set(theme_cowplot())
library("biomaRt")
library(circlize)
library(ComplexHeatmap)
library(SeuratDisk)
library(enrichR)
library(GO.db)
library(org.Mm.eg.db)



Project        <- "MBU_spb54_005__CELL_REV_IntegratedStressResponse_"
baseDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005/"
resDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005/CELL_REVISION_mouse_organogenesis"
setwd(resDir)




go_id = GOID( GOTERM[ Term(GOTERM) == "integrated stress response signaling"])
allegs = get(go_id, org.Mm.egGO2ALLEGS)
genes = unlist(mget(allegs,org.Mm.egSYMBOL))
genes

ISR_genes <- unique(genes)

ISR_genes <- c("Atf4","Eif2s1","Eif2ak3","Eif2ak1","Impact","Nck1","Nck2","Nfe2l2","Ptpn1","Ptpn2","Agr2","Eif2ak4","Abca7","Bok","Tmed2","Dele1","Oma1","Tmem33","Qrich1","Chop","Ddit3","Asns","Chac1","Pck2","Psph" )



message("+-------------------------------------------------------------------------------")
message("+                          load MARKERS tables                                  ")
message("+-------------------------------------------------------------------------------")


Markers_C5024T <- read.csv("../Input/MBU_spb54_005_Markers_C5024T_Seurat.csv", row.names = 1)
Markers_A5019G <- read.csv("../Input/MBU_spb54_005_Markers_A5019G_Seurat.csv", row.names = 1)
Markers_A5019G_C5024T <- read.csv("../Input/MBU_spb54_005_Markers_C5024T_A5019G_vs_WT_Seurat.csv", row.names = 1)
Markers_C5024T[rownames(Markers_C5024T) %in% ISR_genes,]
Markers_A5019G[rownames(Markers_A5019G) %in% ISR_genes,]
Markers_A5019G_C5024T[rownames(Markers_A5019G_C5024T) %in% ISR_genes,]



Markers_amnion_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/amnion_A5019G.csv"), row.names = 1)
Markers_amnion_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/amnion_C5024T.csv"), row.names = 1)
Markers_cardiac_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/cardiac_A5019G.csv"), row.names = 1)
Markers_cardiac_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/cardiac_C5024T.csv"), row.names = 1)
Markers_endothelial_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/endothelial_A5019G.csv"), row.names = 1)
Markers_endothelial_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/endothelial_C5024T.csv"), row.names = 1)
Markers_extraembryonicMesoderm_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/extraembryonicMesoderm_A5019G.csv"), row.names = 1)
Markers_extraembryonicMesoderm_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/extraembryonicMesoderm_C5024T.csv"), row.names = 1)
Markers_forebrain_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/forebrain_A5019G.csv"), row.names = 1)
Markers_forebrain_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/forebrain_C5024T.csv"), row.names = 1)
Markers_foregut_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/foregut_A5019G.csv"), row.names = 1)
Markers_foregut_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/foregut_C5024T.csv"), row.names = 1)
Markers_mesodermProgenitors_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/mesodermProgenitors_A5019G.csv"), row.names = 1)
Markers_mesodermProgenitors_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/mesodermProgenitors_C5024T.csv"), row.names = 1)
Markers_midHindbrain_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/midHindbrain_A5019G.csv"), row.names = 1)
Markers_midHindbrain_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/midHindbrain_C5024T.csv"), row.names = 1)
Markers_midHindgut_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/midHindgut_A5019G.csv"), row.names = 1)
Markers_midHindgut_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/midHindgut_C5024T.csv"), row.names = 1)
Markers_mixedMesoderm_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/mixedMesoderm_A5019G.csv"), row.names = 1)
Markers_mixedMesoderm_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/mixedMesoderm_C5024T.csv"), row.names = 1)
Markers_neuralCrest_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/neuralCrest_A5019G.csv"), row.names = 1)
Markers_neuralCrest_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/neuralCrest_C5024T.csv"), row.names = 1)
Markers_neuralTube_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/neuralTube_A5019G.csv"), row.names = 1)
Markers_neuralTube_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/neuralTube_C5024T.csv"), row.names = 1)
Markers_notochord_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/notochord_A5019G.csv"), row.names = 1)
Markers_notochord_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/notochord_C5024T.csv"), row.names = 1)
Markers_pharyngealMesoderm_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/pharyngealMesoderm_A5019G.csv"), row.names = 1)
Markers_pharyngealMesoderm_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/pharyngealMesoderm_C5024T.csv"), row.names = 1)
Markers_placodes_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/placodes_A5019G.csv"), row.names = 1)
Markers_placodes_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/placodes_C5024T.csv"), row.names = 1)
Markers_presomiticMesoderm_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/presomiticMesoderm_A5019G.csv"), row.names = 1)
Markers_presomiticMesoderm_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/presomiticMesoderm_C5024T.csv"), row.names = 1)
Markers_somiticMesoderm_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/somiticMesoderm_A5019G.csv"), row.names = 1)
Markers_somiticMesoderm_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/somiticMesoderm_C5024T.csv"), row.names = 1)


Markers_celltype_genotype_list <- list(Markers_amnion_A5019G=Markers_amnion_A5019G,   
                                       Markers_amnion_C5024T=Markers_amnion_C5024T,
                                       Markers_cardiac_A5019G=Markers_cardiac_A5019G, 
                                       Markers_cardiac_C5024T=Markers_cardiac_C5024T,
                                       Markers_endothelial_A5019G=Markers_endothelial_A5019G, 
                                       Markers_endothelial_C5024T=Markers_endothelial_C5024T,
                                       Markers_extraembryonicMesoderm_A5019G=Markers_extraembryonicMesoderm_A5019G, 
                                       Markers_extraembryonicMesoderm_C5024T=Markers_extraembryonicMesoderm_C5024T,
                                       Markers_forebrain_A5019G=Markers_forebrain_A5019G,
                                       Markers_forebrain_C5024T=Markers_forebrain_C5024T,
                                       Markers_foregut_A5019G=Markers_foregut_A5019G,
                                       Markers_foregut_C5024T=Markers_foregut_C5024T,
                                       Markers_mesodermProgenitors_A5019G=Markers_mesodermProgenitors_A5019G,
                                       Markers_mesodermProgenitors_C5024T=Markers_mesodermProgenitors_C5024T,
                                       Markers_midHindbrain_A5019G=Markers_midHindbrain_A5019G,
                                       Markers_midHindbrain_C5024T=Markers_midHindbrain_C5024T,
                                       Markers_midHindgut_A5019G=Markers_midHindgut_A5019G,
                                       Markers_midHindgut_C5024T=Markers_midHindgut_C5024T,
                                       Markers_mixedMesoderm_A5019G=Markers_mixedMesoderm_A5019G,
                                       Markers_mixedMesoderm_C5024T=Markers_mixedMesoderm_C5024T,
                                       Markers_neuralCrest_A5019G=Markers_neuralCrest_A5019G,
                                       Markers_neuralCrest_C5024T=Markers_neuralCrest_C5024T,
                                       Markers_neuralTube_A5019G=Markers_neuralTube_A5019G,
                                       Markers_neuralTube_C5024T=Markers_neuralTube_C5024T,
                                       Markers_notochord_A5019G=Markers_notochord_A5019G,
                                       Markers_notochord_C5024T=Markers_notochord_C5024T,
                                       Markers_pharyngealMesoderm_A5019G=Markers_pharyngealMesoderm_A5019G,
                                       Markers_pharyngealMesoderm_C5024T=Markers_pharyngealMesoderm_C5024T,
                                       Markers_placodes_A5019G=Markers_placodes_A5019G,
                                       Markers_placodes_C5024T=Markers_placodes_C5024T,
                                       Markers_presomiticMesoderm_A5019G=Markers_presomiticMesoderm_A5019G,
                                       Markers_presomiticMesoderm_C5024T=Markers_presomiticMesoderm_C5024T,
                                       Markers_somiticMesoderm_A5019G=Markers_somiticMesoderm_A5019G,
                                       Markers_somiticMesoderm_C5024T=Markers_somiticMesoderm_C5024T)



names(Markers_celltype_genotype_list)


for(i in names(Markers_celltype_genotype_list)){
  Markers_celltype_genotype_list[[i]]$gene <- rownames(Markers_celltype_genotype_list[[i]])
}

View(Markers_celltype_genotype_list[["endothelial C5024T"]])




message("+-------------------------------------------------------------------------------")
message("+                             ISR enrichment                                    ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 0.25


library(clusterProfiler)
db <- c("GO_Biological_Process_2018") # "WikiPathways_2019_Mouse"


identify_ISR_enrichment <- function(markers_table=NULL, l2fc_cutoff=0.25, terms_of_interest="integrated stress response"){
  
  genes_for_enrichment <- rownames(markers_table[ abs(markers_table$avg_log2FC) > l2fc_cutoff, ])
  
  if( length(genes_for_enrichment) > 0){
    
    ego2 <- enrichGO(gene = genes_for_enrichment, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont= "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.1)
    ego2_df <- as.data.frame(ego2)
    for(i in terms_of_interest){
      matches1 <- unique (grep(paste(terms_of_interest,collapse="|"), 
                               ego2_df$Description, value=TRUE))
      
      DF_selected_terms <- ego2_df[ego2_df$Description %in% matches1,]
      return(DF_selected_terms)
    } 
  } else { return(NA)}
}


identify_ISR_genes <- function(markers_table=NULL, l2fc_cutoff=0.25, ISR_genes=ISR_genes){
  markers <- rownames(markers_table[ abs(markers_table$avg_log2FC) > l2fc_cutoff, ])
  markers_ISR <- markers[markers %in% ISR_genes]
  return(markers_ISR)
}

ISR_results_list <- list()
ISR_genes_list <- list()

for (i in seq_along(names(Markers_celltype_genotype_list))){
  name <- names(Markers_celltype_genotype_list)[[i]]
  print(name)
  ISR_results_list[[i]] <- identify_ISR_enrichment(markers_table=Markers_celltype_genotype_list[[i]], l2fc_cutoff=0.25)
  
  markers_table <- Markers_celltype_genotype_list[[i]]
  markers <- rownames(markers_table[ abs(markers_table$avg_log2FC) > l2fc_cutoff, ])
  markers_ISR <- markers[markers %in% ISR_genes]
  ISR_genes_list[[i]] <- markers_ISR
}

ISR_results_list
ISR_genes_list


names(ISR_genes_list) <- names(Markers_celltype_genotype_list)


names(Markers_celltype_genotype_list)[[9]]
ISR_genes_list[[9]]

names(Markers_celltype_genotype_list)[[19]]
ISR_genes_list[[19]]






Atf4_reg_genes <- c( "Nrf2", "Clock","Tfeb", "Tfe3", "Nupr1","Cebpg","Cep290","Cenpf","Fos","Ep300","Trib3","Egln3","Egln1")
ISR_genes2 <- c(ISR_genes, Atf4_reg_genes )



message("+-------------------------------------------------------------------------------")
message("+                     heatmaps                            ")
message("+-------------------------------------------------------------------------------")

Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))

Seurat_obj@meta.data$CellType_Genotype <- paste(Seurat_obj@meta.data$celltype, Seurat_obj@meta.data$mouse, sep = "_")

Seurat_obj@meta.data$CellType_Genotype2 <- gsub( "_", " ", Seurat_obj@meta.data$CellType_Genotype)
Seurat_obj@meta.data$CellType_Genotype2 <- gsub( "Hind", " hind", Seurat_obj@meta.data$CellType_Genotype2)
Seurat_obj@meta.data$CellType_Genotype2 <- gsub( "Crest", " crest", Seurat_obj@meta.data$CellType_Genotype2)
Seurat_obj@meta.data$CellType_Genotype2 <- gsub( "Progenitors", " progenitors", Seurat_obj@meta.data$CellType_Genotype2)
Seurat_obj@meta.data$CellType_Genotype2 <- gsub( "Tube", " tube", Seurat_obj@meta.data$CellType_Genotype2)
Seurat_obj@meta.data$CellType_Genotype2 <- gsub( "Meso", " meso", Seurat_obj@meta.data$CellType_Genotype2)
unique(Seurat_obj@meta.data$CellType_Genotype2)
unique(Seurat_obj@meta.data$CellType_Genotype)


celltype_cols <- c( "amnion"="plum1", "mesoderm progenitors"="darkolivegreen1", "neural tube"="plum4","mixed mesoderm"="gold2", "neural crest"="purple4", "mid hindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomitic mesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngeal mesoderm"="royalblue","foregut"="violetred3", "extraembryonic mesoderm"="steelblue3", "cardiac"="firebrick" , "somitic mesoderm"= "darkgreen",  "endothelial"="orange1", "mid hindgut"="violetred2", "blood"="black", "PGCs"="grey10")


Idents(Seurat_obj) <- Seurat_obj@meta.data$CellType_Genotype2
pseudobulk_avg_expr <- AverageExpression(Seurat_obj, assays = "SCT", slot = "data")[[1]]
pseudobulk_avg_expr2 <- t(scale(t(pseudobulk_avg_expr)))
head(pseudobulk_avg_expr2)




mat <- pseudobulk_avg_expr2[rownames(pseudobulk_avg_expr2) %in% ISR_genes,]
idx <- which(colnames(mat) == "blood A5019G" )
mat <- mat[,-idx]

mat_anno <- data.frame(group=colnames(mat))
mat_anno$tissue <- gsub( " A5019G", "", mat_anno$group)
mat_anno$tissue <- gsub( " C5024T", "", mat_anno$tissue)
mat_anno$tissue <- gsub( " WT", "", mat_anno$tissue)


mat_anno$mouse <- gsub( ".* ", "", mat_anno$group)
mat_anno$mouse <- gsub( "C5024T", "m.5024C>T", mat_anno$mouse)
mat_anno$mouse <- gsub( "A5019G", "m.5019A>G", mat_anno$mouse)

mat_anno$mouse <- factor(mat_anno$mouse, levels = c( "WT", "m.5024C>T", "m.5019A>G" ))
mat_anno <- mat_anno[order(mat_anno$mouse),]
mat_anno <- mat_anno[order(mat_anno$tissue),]


names(Markers_celltype_genotype_list) <- gsub( "Markers_", "", names(Markers_celltype_genotype_list))
names(Markers_celltype_genotype_list) <- gsub( "_", " ", names(Markers_celltype_genotype_list))
names(Markers_celltype_genotype_list) <- gsub( "Meso", " meso", names(Markers_celltype_genotype_list))
names(Markers_celltype_genotype_list) <- gsub( "Prog", " prog", names(Markers_celltype_genotype_list))
names(Markers_celltype_genotype_list) <- gsub( "Hind", " hind", names(Markers_celltype_genotype_list))
names(Markers_celltype_genotype_list) <- gsub( "Tube", " tube", names(Markers_celltype_genotype_list))
names(Markers_celltype_genotype_list) <- gsub( "Crest", " crest", names(Markers_celltype_genotype_list))

mat_fdr <- mat
mat_fdr <- mat_fdr*0
colnames(mat_fdr)[colnames(mat_fdr) %in% names(Markers_celltype_genotype_list)]

for(i in names(Markers_celltype_genotype_list)){
  
  for(gene in rownames(mat_fdr)){
    mat_fdr[gene, i] <- ifelse( gene %in% rownames(Markers_celltype_genotype_list[[i]] ) , TRUE, FALSE)
  }
  
}

mat <- mat[, match( mat_anno$group, colnames(mat))]
mat_fdr <- mat_fdr[, match( mat_anno$group, colnames(mat_fdr))]

ha = HeatmapAnnotation(CellType = mat_anno$tissue, col = list(CellType = celltype_cols))

mat_anno$tissue <- gsub( " mesoderm", "\nmesoderm", mat_anno$tissue)
mat_anno$tissue <- gsub( " progenitors", "\nprogenitors", mat_anno$tissue)
colnames(mat) <- gsub( ".* ", "", colnames(mat))


cell_fun2 = function(j, i, x, y, w, h, fill) {
  if(mat_fdr[i, j] == TRUE) {
    grid.text("*", x, y)
  } 
}



row_split <- data.frame(gene =rownames(mat))
row_split$split <- ifelse(row_split$gene %in% ISR_genes,   "ISR", "Atf4 regulator")


colnames(mat) <- gsub( "C5024T", "m.5024C>T", colnames(mat))
colnames(mat) <- gsub( "A5019G", "m.5019A>G", colnames(mat))

f1 = circlize::colorRamp2( c(-3, 0,3), c("green3", "white", "magenta4"), space = "RGB") #  colorRampPalette(c("blue","white","red"))(100)
ht1 <- ComplexHeatmap::Heatmap(mat, name="Expression", col = f1,  cluster_columns = FALSE, show_row_names = TRUE, column_split = mat_anno$tissue , top_annotation = ha, cluster_column_slices=TRUE, column_title_rot = 90, height  = unit(10, "cm"), width = unit(20, "cm"),  cell_fun=cell_fun2, row_names_gp = gpar(fontface = "italic"), split = row_split$split) 


ht1

pdf(paste( Project, "_ComplexHeatmap_", "ISR_genes_",  "_scaledExpr",".pdf", sep="_"), width=15,height=10) # "_celltype_regulators",
par(bg=NA)
ht1
dev.off()


Markers_A5019G[rownames(Markers_A5019G) %in% ISR_genes,]
Markers_C5024T[rownames(Markers_C5024T) %in% ISR_genes,]








