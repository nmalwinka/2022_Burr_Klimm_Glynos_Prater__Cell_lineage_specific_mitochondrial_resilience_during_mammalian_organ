#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()

library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Seurat") 
library("biomaRt")

Project        <- "MBU_spb54_005__fig__Buffering"
baseDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005"
setwd(baseDir)

Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))



message("+-------------------------------------------------------------------------------")
message("+                       add extra labeling                                      ")
message("+-------------------------------------------------------------------------------")


celltype_cols <- c( "amnion"="plum1", "mesoderm progenitors"="darkolivegreen1", "neural tube"="plum4","mixed mesoderm"="gold2", "neural crest"="purple4", "mid hindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomitic mesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngeal mesoderm"="royalblue","foregut"="violetred3", "extraembryonic mesoderm"="steelblue3", "cardiac"="firebrick" , "somitic mesoderm"= "darkgreen",  "endothelial"="orange1", "mid hindgut"="violetred2", "blood"="black")

Seurat_obj@meta.data$celltype2 <- Seurat_obj@meta.data$celltype
Seurat_obj@meta.data$celltype2 <- gsub( "Hind", " hind", Seurat_obj@meta.data$celltype2)
Seurat_obj@meta.data$celltype2 <- gsub( "Mesoderm", " mesoderm", Seurat_obj@meta.data$celltype2)
Seurat_obj@meta.data$celltype2 <- gsub( "Crest", " crest", Seurat_obj@meta.data$celltype2)
Seurat_obj@meta.data$celltype2 <- gsub( "Progenitors", " progenitors", Seurat_obj@meta.data$celltype2)
Seurat_obj@meta.data$celltype2 <- gsub( "Tube", " tube", Seurat_obj@meta.data$celltype2)
Seurat_obj@meta.data$celltype2 <- factor(Seurat_obj@meta.data$celltype2, levels = c(
  "amnion", "mesoderm progenitors","neural tube","mixed mesoderm", "neural crest",  "mid hindbrain", "notochord", "placodes","presomitic mesoderm",  "forebrain",  "pharyngeal mesoderm",  "foregut", "extraembryonic mesoderm",   "cardiac",  "somitic mesoderm" , "endothelial", "mid hindgut", "blood"))


Seurat_obj@meta.data$celltype_genotype <- paste(Seurat_obj@meta.data$celltype2, Seurat_obj@meta.data$mouse, sep = "_")
Seurat_obj@meta.data$celltype_genotype <- gsub( "_", " ", Seurat_obj@meta.data$celltype_genotype)
Seurat_obj@meta.data$celltype_genotype <- gsub( "A5019G", "m.5019A>G", Seurat_obj@meta.data$celltype_genotype)
Seurat_obj@meta.data$celltype_genotype <- gsub( "C5024T", "m.5024C>T", Seurat_obj@meta.data$celltype_genotype)

unique(Seurat_obj@meta.data$celltype_genotype)
Seurat_obj@meta.data$celltype_genotype <- factor(Seurat_obj@meta.data$celltype_genotype, levels = c(
  "blood m.5019A>G",
  "amnion WT", "amnion m.5024C>T" ,"amnion m.5019A>G" ,
  "mesoderm progenitors WT" ,"mesoderm progenitors m.5024C>T","mesoderm progenitors m.5019A>G",
  "neural tube WT","neural tube m.5024C>T" , "neural tube m.5019A>G",
  "mixed mesoderm WT",  "mixed mesoderm m.5024C>T" ,"mixed mesoderm m.5019A>G" ,
  "neural crest WT", "neural crest m.5024C>T"  ,"neural crest m.5019A>G", 
  "mid hindbrain WT", "mid hindbrain m.5024C>T", "mid hindbrain m.5019A>G", 
  "notochord WT","notochord m.5024C>T", "notochord m.5019A>G",          
  "placodes WT", "placodes m.5024C>T"   ,   "placodes m.5019A>G",     
  "presomitic mesoderm WT",  "presomitic mesoderm m.5024C>T", "presomitic mesoderm m.5019A>G",  
  "forebrain WT","forebrain m.5024C>T" ,"forebrain m.5019A>G" ,         
  "pharyngeal mesoderm WT", "pharyngeal mesoderm m.5024C>T", "pharyngeal mesoderm m.5019A>G", 
  "foregut WT", "foregut m.5024C>T", "foregut m.5019A>G" ,  
  "extraembryonic mesoderm WT", "extraembryonic mesoderm m.5024C>T", "extraembryonic mesoderm m.5019A>G",
  "cardiac WT", "cardiac m.5024C>T", "cardiac m.5019A>G" ,       
  "somitic mesoderm WT" , "somitic mesoderm m.5024C>T" ,"somitic mesoderm m.5019A>G",    
  "endothelial WT", "endothelial m.5024C>T", "endothelial m.5019A>G",    
  "mid hindgut WT","mid hindgut m.5024C>T", "mid hindgut m.5019A>G" ))
Idents(Seurat_obj) <- Seurat_obj@meta.data$celltype_genotype

Seurat_obj@meta.data$mouse <- gsub( "A5019G", "m.5019A>G", Seurat_obj@meta.data$mouse)
Seurat_obj@meta.data$mouse <- gsub( "C5024T", "m.5024C>T", Seurat_obj@meta.data$mouse)

Seurat_obj@meta.data$celltype2 <- factor(Seurat_obj@meta.data$celltype2, levels = c("blood", "amnion","cardiac","endothelial","extraembryonic mesoderm","forebrain" ,"foregut" , "mesoderm progenitors", "mid hindbrain", "mid hindgut","mixed mesoderm",  "neural crest" ,"neural tube" , "notochord", "pharyngeal mesoderm", "placodes", "presomitic mesoderm","somitic mesoderm"    ))

Seurat_obj <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-", col.name = "percent.mt")



message("+-------------------------------------------------------------------------------")
message("+                     load in mitocarta genes                                   ")
message("+-------------------------------------------------------------------------------")

MitoCarta3 <- read.csv("/Users/mn367/Documents/MBU-Projects/Databases/Mouse.MitoCarta3.0_summarised.csv")
MitoCarta3_pathways <- read.csv("/Users/mn367/Documents/MBU-Projects/Databases/Mouse.MitoCarta3.0_pathways.csv")
MitoCarta3_genes_mouse <- MitoCarta3$Symbol



message("+-------------------------------------------------------------------------------")
message("+                            load in Mootha genes                               ")
message("+-------------------------------------------------------------------------------")

genes_mootha <- read.csv("/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005/Input/Mootha_Gene_List.csv")

homolog_genes <- human2mouse(genes_mootha$gene, db = homologene::homologeneData)
genes_mootha$mouse_gene <- homolog_genes[match(genes_mootha$gene , homolog_genes$humanGene),]$mouseGene

genes_mootha$buff <- ifelse(genes_mootha$category == "Mootha Buffering/supressor genes", TRUE, FALSE)
genes_mootha$leth <- ifelse(genes_mootha$category == "Mootha Synthetic Sick/Lethal", TRUE, FALSE)

mootha_mouse_genes <- unique(genes_mootha$mouse_gene)
genes_mootha_buff <- unique(genes_mootha[genes_mootha$buff == TRUE,]$mouse_gene)
genes_mootha_leth <- unique(genes_mootha[genes_mootha$leth == TRUE,]$mouse_gene)



message("+-------------------------------------------------------------------------------")
message("+                            load in Larsson genes                              ")
message("+-------------------------------------------------------------------------------")

Larsson_IMT_vs_ctrl <- read.csv("/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_003_001/Larsson_2021_IMT1/embr202153054-sup-0003-datasetev2_CRISPR-Cas9_Screen.csv")
Larsson_IMT_vs_ctrl$significant <- ifelse(Larsson_IMT_vs_ctrl$logP > 1.3, TRUE, FALSE) # threshold 1.3 for padj 0.05 and 2 for padj 0.01
Larsson_IMT_vs_ctrl <- Larsson_IMT_vs_ctrl[Larsson_IMT_vs_ctrl$significant == TRUE,]

Larsson_IMT_vs_ctrl$binary <- ifelse(Larsson_IMT_vs_ctrl$LFC > 0, "POSITIVE_HIT_resistance_to_MTI1", "NEGATIVE_HIT_susceptibility_to_MTI1")
table(Larsson_IMT_vs_ctrl[abs(Larsson_IMT_vs_ctrl$LFC) > 0.6,]$binary)

homolog_genes2 <- human2mouse(Larsson_IMT_vs_ctrl$gene_name, db = homologene::homologeneData)
Larsson_IMT_vs_ctrl$mouse_gene <- homolog_genes2[match(Larsson_IMT_vs_ctrl$gene_name , homolog_genes2$humanGene),]$mouseGene

idx <- grep("Control_.*", Larsson_IMT_vs_ctrl$gene_name)
Larsson_IMT_vs_ctrl <- Larsson_IMT_vs_ctrl[-idx,]
Larsson_IMT_vs_ctrl_l2fc0.6 <- Larsson_IMT_vs_ctrl[abs(Larsson_IMT_vs_ctrl$LFC) >= 0.6,] 

Larsson_IMT_vs_ctrl_BUFF <- Larsson_IMT_vs_ctrl[Larsson_IMT_vs_ctrl$binary == "POSITIVE_HIT_resistance_to_MTI1",] 
genes_larsson_res <- Larsson_IMT_vs_ctrl_l2fc0.6[Larsson_IMT_vs_ctrl_l2fc0.6$binary == "POSITIVE_HIT_resistance_to_MTI1",]$mouse_gene
genes_larsson_sus <- Larsson_IMT_vs_ctrl_l2fc0.6[Larsson_IMT_vs_ctrl_l2fc0.6$binary == "NEGATIVE_HIT_susceptibility_to_MTI1",]$mouse_gene



message("+-------------------------------------------------------------------------------")
message("+                         Load in markers                                       ")
message("+-------------------------------------------------------------------------------")

Markers_C5024T <- read.csv("Input/MBU_spb54_005_Markers_C5024T_Seurat.csv", row.names = 1)
Markers_A5019G <- read.csv("Input/MBU_spb54_005_Markers_A5019G_Seurat.csv", row.names = 1)

Markers_genotype_list <- list(Markers_C5024T=Markers_C5024T,   
                              Markers_A5019G=Markers_A5019G)

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

for (i in seq_along(Markers_celltype_genotype_list) ){
  Markers_celltype_genotype_list[[i]]$gene <- rownames(Markers_celltype_genotype_list[[i]])
  Markers_celltype_genotype_list[[i]]$celltype_genotype <- names(Markers_celltype_genotype_list)[[i]]
  
}

Markers_A5019G_celltype <- rbind(Markers_celltype_genotype_list[[1]], 
                                   Markers_celltype_genotype_list[[3]],  
                                   Markers_celltype_genotype_list[[5]], 
                                   Markers_celltype_genotype_list[[7]],
                                   Markers_celltype_genotype_list[[9]],
                                   Markers_celltype_genotype_list[[11]],
                                   Markers_celltype_genotype_list[[13]],
                                   Markers_celltype_genotype_list[[15]],
                                   Markers_celltype_genotype_list[[17]],
                                   Markers_celltype_genotype_list[[19]],
                                   Markers_celltype_genotype_list[[21]],
                                   Markers_celltype_genotype_list[[23]],
                                   Markers_celltype_genotype_list[[25]],
                                   Markers_celltype_genotype_list[[27]],
                                   Markers_celltype_genotype_list[[29]],
                                   Markers_celltype_genotype_list[[31]],
                                   Markers_celltype_genotype_list[[33]])


Markers_C5024T_celltype <- rbind(Markers_celltype_genotype_list[[2]], 
                                   Markers_celltype_genotype_list[[4]],  
                                   Markers_celltype_genotype_list[[6]], 
                                   Markers_celltype_genotype_list[[8]],
                                   Markers_celltype_genotype_list[[10]],
                                   Markers_celltype_genotype_list[[12]],
                                   Markers_celltype_genotype_list[[14]],
                                   Markers_celltype_genotype_list[[16]],
                                   Markers_celltype_genotype_list[[18]],
                                   Markers_celltype_genotype_list[[20]],
                                   Markers_celltype_genotype_list[[22]],
                                   Markers_celltype_genotype_list[[24]],
                                   Markers_celltype_genotype_list[[26]],
                                   Markers_celltype_genotype_list[[28]],
                                   Markers_celltype_genotype_list[[30]],
                                   Markers_celltype_genotype_list[[32]],
                                   Markers_celltype_genotype_list[[34]])






message("+-------------------------------------------------------------------------------")
message("+           percent Larsson and Mootha genes in cell type markers               ")
message("+-------------------------------------------------------------------------------")

total_genes <- length(unique( rownames(GetAssayData(Seurat_obj))))

CellType_name_list <- list()
CellType_genes_list <- list()
CellType_genes_count <- list()
CellType_Mootha_buffering <- list()
CellType_Mootha_lethal <- list()
CellType_Larsson_susceptibility <- list()
CellType_Larsson_resistance <- list()

names(Markers_celltype_genotype_list)

for ( i in seq_along(Markers_celltype_genotype_list)){
  
  CellType_name <- names(Markers_celltype_genotype_list)[i]
  CellType_name <- gsub( "Markers_", "", CellType_name)
  CellType_markers <- Markers_celltype_genotype_list[[i]]

  CellType_markers <- rownames(CellType_markers[CellType_markers$avg_log2FC < 0, ])
  
  CellType_Mootha_buffering[[i]] <- length(CellType_markers[CellType_markers %in% unique(as.character(genes_mootha_buff)) ]  )
  CellType_Mootha_lethal[[i]] <-  length(CellType_markers[CellType_markers %in% unique(as.character(genes_mootha_leth)) ]  )
  CellType_Larsson_susceptibility[[i]] <- length(CellType_markers[CellType_markers %in% Larsson_IMT_vs_ctrl[Larsson_IMT_vs_ctrl$binary == "NEGATIVE_HIT_susceptibility_to_MTI1" ,]$mouse_gene ] )
  CellType_Larsson_resistance[[i]] <- length(CellType_markers[CellType_markers %in% Larsson_IMT_vs_ctrl[Larsson_IMT_vs_ctrl$binary == "POSITIVE_HIT_resistance_to_MTI1" ,]$mouse_gene] )
  CellType_name_list[[i]] <- CellType_name
  CellType_genes_list[[i]] <- CellType_markers
  CellType_genes_count[[i]] <- length(unique(CellType_markers))
}


annotation_Mootha_Larsson <- data.frame(Cell_type=unlist(CellType_name_list), 
                                        Markers_size=unlist(CellType_genes_count),
                                        Mootha_buffering=unlist(CellType_Mootha_buffering),
                                        Mootha_lethal=unlist(CellType_Mootha_lethal),
                                        Larsson_susceptibility=unlist(CellType_Larsson_susceptibility),
                                        Larsson_resistance=unlist(CellType_Larsson_resistance))
annotation_Mootha_Larsson$Cell_type <- gsub( "_C5024T", " m.5024C>T", annotation_Mootha_Larsson$Cell_type)
annotation_Mootha_Larsson$Cell_type <- gsub( "_A5019G", " m.5019A>G", annotation_Mootha_Larsson$Cell_type)
rownames(annotation_Mootha_Larsson) <- annotation_Mootha_Larsson$Cell_type


annotation_Mootha_Larsson_DOWN_MARKERS <- annotation_Mootha_Larsson



message("--------------------------------------------------------------------------------")
message("+                        HYPERGEOMETRIC TEST                                    ")
message("+-------------------------------------------------------------------------------")

total_genes <- length(unique( rownames(GetAssayData(Seurat_obj))))


Gene_Proportions_in_Patterns <- annotation_Mootha_Larsson_DOWN_MARKERS[,c("Markers_size","Mootha_buffering" ,"Mootha_lethal","Larsson_susceptibility", "Larsson_resistance" )]


# How to use phyper in R: HYPERGEOMETRIC TEST
####### https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/

GEN_SET_NAME <- "Mootha_buffering"
GEN_SET_GENES <- genes_mootha_buff

#GEN_SET_NAME <- "Larsson_susceptibility"
#GEN_SET_GENES <- genes_larsson_sus



hypergeometric_tests_list <- list()
fisher_exact_tests_list <- list()
count <- 1

for(CELL_TYPE in rownames(Gene_Proportions_in_Patterns)){
  GEN_SET_GENES <- GEN_SET_GENES[GEN_SET_GENES %in% unique( rownames(GetAssayData(Seurat_obj)))]
  group1 = length(unique(GEN_SET_GENES))
  group2 = Gene_Proportions_in_Patterns[rownames(Gene_Proportions_in_Patterns) == CELL_TYPE,]$Markers_size
  Overlap = Gene_Proportions_in_Patterns[CELL_TYPE, GEN_SET_NAME]
  Total= total_genes
  # Test for over-representation (enrichment)
  phyper(Overlap-1, group2, Total-group2, group1, lower.tail= FALSE)
  #phyper(q=Overlap -1, m=group1, n=Total-group1, k=group2, lower.tail=FALSE)
  hypergeometric_tests_list[[count]] <- phyper(Overlap-1, group2, Total-group2, group1, lower.tail= FALSE)
  # FIsher exact test
  contingency.table <- data.frame(matrix(nrow=2, ncol=2))
  rownames(contingency.table) <- c("predicted.target", "non.predicted")
  colnames(contingency.table) <- c("class.member", "non.member")
  contingency.table["predicted.target", "class.member"] <- Overlap ## Number of marked genes in the selection
  contingency.table["predicted.target", "non.member"] <- group2 - Overlap ## Number of non-marked genes in the selection
  contingency.table["non.predicted", "class.member"] <- group1 - Overlap ## Number of marked genes outside of the selection
  contingency.table["non.predicted", "non.member"] <- Total - (group2 - Overlap) ## Number of non-marked genes in the selection
  (contingency.row.sum <- apply(contingency.table, 1, sum))
  (contingency.col.sum <- apply(contingency.table, 2, sum))
  contingency.table.margins <- cbind(contingency.table, contingency.row.sum)
  contingency.table.margins <- rbind(contingency.table.margins, apply(contingency.table.margins, 2, sum))
  names(contingency.table.margins) <- c(names(contingency.table), "total")
  rownames(contingency.table.margins) <- c(rownames(contingency.table), "total")
  print(contingency.table.margins)
  print(sum(contingency.table)) ## The value shoudl equal N, since every
  count <- count + 1
}

names(hypergeometric_tests_list) <- rownames(Gene_Proportions_in_Patterns)



message("+-------------------------------------------------------------------------------")
message("+                       Figure 4 G                              ")
message("+-------------------------------------------------------------------------------")

hypergeom_testing_df <- data.frame( pval = unlist(hypergeometric_tests_list))
hypergeom_testing_df$padj <- p.adjust(hypergeom_testing_df$pval, method = "BH" )
hypergeom_testing_df$star <- ifelse(hypergeom_testing_df$padj < 0.05, "*", "ns")

mat_heatmap_hyper <- data.frame(Downregulated=hypergeom_testing_df$padj)
rownames(mat_heatmap_hyper) <- rownames(hypergeom_testing_df)
rownames(mat_heatmap_hyper) <- gsub( "Hind", " hind", rownames(mat_heatmap_hyper))
rownames(mat_heatmap_hyper) <- gsub( "Mesoderm", " mesoderm", rownames(mat_heatmap_hyper))
rownames(mat_heatmap_hyper) <- gsub( "Crest", " crest", rownames(mat_heatmap_hyper))
rownames(mat_heatmap_hyper) <- gsub( "Progenitors", " progenitors", rownames(mat_heatmap_hyper))
rownames(mat_heatmap_hyper) <- gsub( "Tube", " tube", rownames(mat_heatmap_hyper))
rownames(mat_heatmap_hyper)

mat_heatmap_hyper$cell_type <- gsub(" m\\..*", "", rownames(mat_heatmap_hyper)) 
mat_heatmap_hyper$genotype <- gsub(".* m.", "m.", rownames(mat_heatmap_hyper)) 
mat_heatmap_hyper_molten <- reshape2::melt(mat_heatmap_hyper)

mat_heatmap_hyper_cast <- reshape2::dcast(mat_heatmap_hyper_molten[,-c(3)], formula= cell_type ~ genotype)
rownames(mat_heatmap_hyper_cast) <- mat_heatmap_hyper_cast$cell_type

mat_heatmap_hyper_cast <- mat_heatmap_hyper_cast[,c(2,3)]
head(mat_heatmap_hyper_cast)


ht_cols = circlize::colorRamp2( c(0,  0.051,  1 ), c( "green4",  "grey95", "grey95"), space = "RGB") #  v"#edf8fb", "#8c6bb1", 
#ht_cols = circlize::colorRamp2( c(0,  0.051,  1 ), c( "#6e016b",  "grey95", "grey95"), space = "RGB") #  v"#edf8fb", "#8c6bb1", 


ht_celltype_DOWN = ComplexHeatmap:: Heatmap(mat_heatmap_hyper_cast,  col = ht_cols, name = "Downregulated",  row_title = "", column_title = "Downregulated", show_row_names = TRUE, show_column_names = TRUE, heatmap_legend_param = list(title = paste0("P.value (adj)"), legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", column_names_rot = 0, width = unit(5, "cm"),column_names_side = "bottom", column_title_rot = 0, column_names_centered = TRUE) 


pdf(paste( Project, "ComplexHeatmap", "Mootha_gene_enrichment", "DOWN_regulated_gene_markers_","METHOD", "hypergeometric_test","refined_Mootha_genes_v3.pdf", sep="_"), width=5,height=6) 
par(bg=NA)
ht_celltype_DOWN
dev.off()





message("--------------------------------------------------------------------------------")
message("+                                Fig 4 C VENN                                   ")
message("+-------------------------------------------------------------------------------")

library("venneuler")
library("VennDiagram")

Markers_A5019G_celltype_genes <- unique(c(rownames(Markers_amnion_A5019G), rownames(Markers_cardiac_A5019G),  rownames(Markers_endothelial_A5019G), rownames(Markers_extraembryonicMesoderm_A5019G), rownames(Markers_forebrain_A5019G), rownames(Markers_foregut_A5019G), rownames(Markers_mesodermProgenitors_A5019G), rownames(Markers_midHindbrain_A5019G), rownames(Markers_midHindgut_A5019G), rownames(Markers_mixedMesoderm_A5019G), rownames(Markers_neuralCrest_A5019G), rownames(Markers_neuralTube_A5019G), rownames(Markers_notochord_A5019G), rownames(Markers_pharyngealMesoderm_A5019G), rownames(Markers_placodes_A5019G), rownames(Markers_presomiticMesoderm_A5019G), rownames(Markers_somiticMesoderm_A5019G) )) 

Markers_C5024T_celltype_genes <- unique(c(rownames(Markers_amnion_C5024T), rownames(Markers_cardiac_C5024T),  rownames(Markers_endothelial_C5024T), rownames(Markers_extraembryonicMesoderm_C5024T), rownames(Markers_forebrain_C5024T), rownames(Markers_foregut_C5024T), rownames(Markers_mesodermProgenitors_C5024T), rownames(Markers_midHindbrain_C5024T), rownames(Markers_midHindgut_C5024T), rownames(Markers_mixedMesoderm_C5024T), rownames(Markers_neuralCrest_C5024T), rownames(Markers_neuralTube_C5024T), rownames(Markers_notochord_C5024T), rownames(Markers_pharyngealMesoderm_C5024T), rownames(Markers_placodes_C5024T), rownames(Markers_presomiticMesoderm_C5024T), rownames(Markers_somiticMesoderm_C5024T) )) 




# for mootha:::
gen_modif <- unique(c(genes_mootha_buff,genes_mootha_leth))

# for larsson:::
gen_modif <- unique(c(genes_larsson_res,genes_larsson_sus))




gen_modif <- gen_modif[!is.na(gen_modif)]

venn_df <- data.frame(m5019AG=Markers_A5019G_celltype_genes)
#venn_df <- data.frame(m5019AG=rownames(Markers_A5019G))
rownames(venn_df) <- venn_df$m5019AG
tmp_venn <- data.frame(m5024CT=Markers_C5024T_celltype_genes)
#tmp_venn <- data.frame(m5024CT=rownames(Markers_C5024T))
rownames(tmp_venn) <- tmp_venn$m5024CT
venn_df <- merge(venn_df, tmp_venn, by = "row.names", all= TRUE)
tmp_venn <- data.frame(gen_mod=gen_modif)
rownames(tmp_venn) <- tmp_venn$gen_mod
venn_df <- merge(venn_df, tmp_venn, by.x = "Row.names", by.y = "row.names", all= TRUE)
rownames(venn_df) <- venn_df$Row.names
venn_df <- venn_df[,-1]

foo<- function(x) { sum(!is.na(x)) }
purrr::map_dbl(venn_df, foo) # gives number of genes in venn for each group.
# prepare venn table now:
for (i in 1:length(venn_df)) {
  venn_df[[i]] <- sapply(venn_df[[i]], foo, simplify=T)  # 1=yes, 0=NA
}

venn_counts <- limma::vennCounts(venn_df)

#VennDiagram::venn.diagram(venn_counts)

venn_counts



grid.newpage()                    # Create new plotting page
draw.triple.venn(area1 = sum(venn_df$m5019AG), 
                 area2 = sum(venn_df$m5024CT), 
                 area3 = sum(venn_df$gen_mod), 
                 n12 = venn_counts[7,4] + venn_counts[8,4], #(1547+ 22), 
                 n23 = venn_counts[4,4] + venn_counts[8,4], #(41 + 22),
                 n13 = venn_counts[6,4] + venn_counts[8,4], #(18 + 22),
                 n123 = venn_counts[8,4] , # 22, 
                 col = "red",
                 fill = c("pink", "blue", "grey"))



library(venneuler)

venn_plot <- plot(venneuler(c( m5019AG=sum(venn_df$m5019AG),
                               m5024CT=sum(venn_df$m5024CT),
                               gen_mod= sum(venn_df$gen_mod),
                               "m5019AG&m5024CT"=venn_counts[7,4] + venn_counts[8,4],
                               "m5019AG&gen_mod"= venn_counts[4,4] + venn_counts[8,4],
                               "gen_mod&m5024CT"= venn_counts[6,4] + venn_counts[8,4],
                               "gen_mod&m5024CT&m5019AG"= venn_counts[8,4]
)))



venn_plot <- plot(venneuler(c( m5019AG=sum(venn_df$m5019AG),
                               m5024CT=sum(venn_df$m5024CT),
                               gen_mod= sum(venn_df$gen_mod),
                               "m5019AG&m5024CT"=(1547+ 22),
                               "m5019AG&gen_mod"= (18 + 22),
                               "gen_mod&m5024CT"= (41 + 22),
                               "gen_mod&m5024CT&m5019AG"= 22
)))


pdf(paste(Project, "VENNeuler", "_compare_marker_genes_with_Mootha_buff_modifiers", "venn.pdf", sep="_"), width=6, height=5, onefile=FALSE)
par(bg=NA)
plot(venneuler(c( m5019AG=sum(venn_df$m5019AG),
                  m5024CT=sum(venn_df$m5024CT),
                  gen_mod= sum(venn_df$gen_mod),
                  "m5019AG&m5024CT"=venn_counts[7,4] + venn_counts[8,4],
                  "m5019AG&gen_mod"= venn_counts[4,4] + venn_counts[8,4],
                  "gen_mod&m5024CT"= venn_counts[6,4] + venn_counts[8,4],
                  "gen_mod&m5024CT&m5019AG"= venn_counts[8,4]
)))

dev.off()




message("--------------------------------------------------------------------------------")
message("+                 Suppl figures                                    ")
message("+-------------------------------------------------------------------------------")

tbl_S3_genes <- rownames(venn_df[venn_df$m5019AG == 1 & venn_df$m5024CT == 1 & venn_df$gen_mod == 1,])
tbl_S3_genes <- venn_df[(venn_df$m5019AG == 1 | venn_df$m5024CT == 1) & venn_df$gen_mod == 1,]
tbl_S3 <- as.data.frame(tbl_S3_genes)
colnames(tbl_S3) <- "gene"

#write.csv(tbl_S3, "MBU_spb54_005_tbl_S3_overlap_genes_mootha_markers_v2.csv")
#write.csv(tbl_S3, "MBU_spb54_005_tbl_S4_overlap_genes_larsson_markers_v2.csv")


Tbl_S2_5019 <- read.csv("MBU_spb54_005_Markers_A5019G_pseudobulk_Seurat.csv", row.names = 1)
Tbl_S2_5024 <- read.csv("MBU_spb54_005_Markers_C5024T_pseudobulk_Seurat.csv", row.names = 1)

Tbl_S2_5019_mito <- Tbl_S2_5019[Tbl_S2_5019$gene %in% Mito_genes,]
Tbl_S2_5024_mito <- Tbl_S2_5024[Tbl_S2_5024$gene %in% Mito_genes,]

write.csv(Tbl_S2_5019_mito, "MBU_spb54_005_Markers_A5019G_pseudobulk_Seurat_mito.carta3.csv")
write.csv(Tbl_S2_5024_mito, "MBU_spb54_005_Markers_C5024T_pseudobulk_Seurat_mito.carta3.csv")







message("+-------------------------------------------------------------------------------")
message("+           Figure S6 A heatmap mootha buff all                                 ")
message("+-------------------------------------------------------------------------------")

exprMat_avg <- AverageExpression(Seurat_obj, slot = "scale.data", assay= "SCT")
exprMat_avg <- exprMat_avg[[1]]
exprMat_avg[1:10,1:10]
exprMat_avg <- as.data.frame(exprMat_avg)

selected_genes <- unique(c( genes_mootha_buff  ))
Heatmap_name <- "Mootha_buffering_genes"



exprMat_ht <- exprMat_avg[rownames(exprMat_avg) %in% selected_genes,-1]

col_split <- gsub( " WT", "", colnames(exprMat_ht))
col_split <- gsub( " m.5024C>T", "", col_split)
col_split <- gsub( " m.5019A>G", "", col_split)

exprMat_ht[1:10,1:10]
#write.csv(exprMat_ht, "MBU_spb54_005_Tbl_S6_for_Fig_S6A_mat.csv")

exprMat_ht2 <- exprMat_ht



colnames(exprMat_ht) <- gsub(  ".* WT", "WT", colnames(exprMat_ht))
colnames(exprMat_ht) <- gsub(  ".* m.5019A>G", "m.5019A>G", colnames(exprMat_ht))
colnames(exprMat_ht) <- gsub(  " m.5024C>T", "  m.5024C>T", colnames(exprMat_ht))

dend1 = cluster_between_groups(exprMat_ht, col_split)
ha = HeatmapAnnotation(CellType = col_split, col = list(CellType = celltype_cols))

          

ht_colss = circlize::colorRamp2( c(-1.2, 0,  1.2), c( "blue",  "white", "red"), space = "sRGB") #  v"#edf8fb", "#8c6bb1", 

ht_genes = ComplexHeatmap:: Heatmap(exprMat_ht,  col = ht_colss, name = "Mootha buffering genes",  
                                    row_title = "", column_title = "",
                                    show_row_names = TRUE, show_column_names = TRUE, 
                                    cluster_columns = dend1, cluster_rows = TRUE,
                                    heatmap_legend_param = list(title = paste0("Scaled expression"), legend_height = unit(3, "cm"),
                                                                title_position = "topleft"),  
                                    row_title_rot = 0, column_names_rot = 90 , row_names_rot = 0, column_title_rot = 90 , 
                                    row_names_side ="right", row_dend_side="left",column_names_side = "bottom",column_title_side ="top", 
                                    row_names_max_width = max_text_width(colnames(exprMat_ht)), bottom_annotation = ha,
                                    width = unit(18, "cm"),  height = unit(50, "cm"), row_names_gp = gpar(fontface="italic") ) 

ht_genes



pdf(paste( Project, "ComplexHeatmap", Heatmap_name, "celltype_genotype_slot_scale.data_v2.pdf", sep="_"), width=12, height=25) # "_celltype_regulators",
par(bg=NA)
ht_genes
dev.off()





message("+-------------------------------------------------------------------------------")
message("+           Figure 5 A heatmap mootha buff most variable                        ")
message("+-------------------------------------------------------------------------------")


exprMat_for_split <- as.data.frame(t(exprMat_ht2))
exprMat_for_split
rownames(exprMat_for_split)



exprMat_split <- split(exprMat_for_split, col_split)
names(exprMat_split)
exprMat_split[[1]][,1:5]

KEEP_GENES <- list()
DIFF <- -0.2

for (i in seq_along(exprMat_split)){
  
  df <- exprMat_split[[i]]
  df <- as.data.frame(t(df))
  colnames(df) <- gsub( ".* ", "", colnames(df))
  
  df$`m.5024C>T_vs_WT` <- df$`m.5024C>T` - df$WT
  df$`m.5019A>G_vs_WT` <- df$`m.5019A>G` - df$WT
  KEEP_GENES[[i]] <- rownames(df[df$`m.5024C>T_vs_WT` < DIFF | df$`m.5019A>G_vs_WT` < DIFF,])
  
}


KEEP_GENES <- unique(unlist(KEEP_GENES))




MIN_MAX <- 0.5
exprMat_ht3 <- exprMat_ht[rownames(exprMat_ht) %in% KEEP_GENES,]
exprMat_ht3 <- exprMat_ht3[apply(exprMat_ht3, 1, function(x) max(x)) > MIN_MAX,]


ht_colss = circlize::colorRamp2( c(-1, 0,  1), c( "blue",  "white", "red"), space = "sRGB") #  v"#edf8fb", "#8c6bb1", 



colnames(exprMat_ht3) <- gsub(  ".* WT", "WT", colnames(exprMat_ht3))
colnames(exprMat_ht3) <- gsub(  ".* m.5019A>G", "m.5019A>G", colnames(exprMat_ht3))
colnames(exprMat_ht3) <- gsub(  " m.5024C>T", "   m.5024C>T", colnames(exprMat_ht3))



ht_genes = ComplexHeatmap:: Heatmap(exprMat_ht3,  col = ht_colss, name = "Mootha buffering genes",  row_title = "", column_title = "",
                                    show_row_names = TRUE, show_column_names = TRUE, 
                                    cluster_columns = dend1, cluster_rows = TRUE,
                                    heatmap_legend_param = list(title = paste0("Scaled expression"), legend_height = unit(3, "cm"),
                                                                title_position = "topleft"),  
                                    row_title_rot = 0, column_names_rot = 90 , row_names_rot = 0, column_title_rot = 90 , 
                                    row_names_side ="right", row_dend_side="left",column_names_side = "bottom",column_title_side ="top", 
                                    row_names_max_width = max_text_width(colnames(exprMat_ht)), bottom_annotation = ha,
                                    width = unit(23, "cm"),  height = unit(17, "cm"), row_names_gp = gpar(fontface="italic") ) 
# width = unit(6, "cm"),  column_names_max_height=max_text_height(rownames(exprMat_ht)), cluster_column_slices = TRUE,
ht_genes



pdf(paste( "MBU_spb54_005__fig6", "ComplexHeatmap", Heatmap_name, "celltype_genotype_slot_scale.data", "FILT_diff", DIFF, "max", MIN_MAX , "v2.pdf", sep="_"), width=15, height=15) # "_celltype_regulators",
par(bg=NA)
ht_genes
dev.off()







message("+-------------------------------------------------------------------------------")
message("+                                enrich R                                       ")
message("+-------------------------------------------------------------------------------")

db <- c("WikiPathways_2019_Mouse", "GO_Biological_Process_2018")


l2fc_cutoff <- 0.25
names(Markers_celltype_genotype_list)
Enrichment_celltype_genotype_list <- list()

for (i in names(Markers_celltype_genotype_list)){
  
  markers_table <- Markers_celltype_genotype_list[[i]]
  celltype_name <- i
  #colnames(markers_table)[ colnames(markers_table) == "avg_log2FC"] <- "logfoldchanges"
  #markers_table$Direction <- ifelse(markers_table$logfoldchanges > 0, "Overexpressed in mutator", "Underexpressed in mutator")
  #genes_for_enrichment <- unique(markers_table[abs(markers_table$logfoldchanges) > l2fc_cutoff & markers_table$Direction =="Overexpressed in mutator",]$gene.name)
  genes_for_enrichment <- unique( rownames(markers_table[(markers_table$avg_log2FC) > l2fc_cutoff ,]) )
  
  print(length(genes_for_enrichment))
  enrich_RES <- enrichr(genes_for_enrichment, databases = db )
  enrich_RES_WP <- enrich_RES[[1]][enrich_RES[[1]]$Adjusted.P.value < 0.05,]
  if(nrow(enrich_RES_WP)>0){
    enrich_RES_WP$db <- "WP"
    write.csv(enrich_RES_WP, paste0(Project, "_enrich_res_Seurat", "WP_", celltype_name, "_vs_WT_l2fc0.25", ".csv"))
  }
  enrich_RES_GOBP <- enrich_RES[[2]][enrich_RES[[2]]$Adjusted.P.value < 0.05,]
  if(nrow(enrich_RES_GOBP)>0){
    enrich_RES_GOBP$db <- "GOBP"
    write.csv(enrich_RES_GOBP, paste0(Project, "_enrich_res_Seurat", "GOBP_", celltype_name, "_vs_WT_l2fc0.25", ".csv"))
  }
  enrich_RES_merged <- rbind(enrich_RES_WP, enrich_RES_GOBP)
  Enrichment_celltype_genotype_list[[i]] <- enrich_RES_merged
}

names(Enrichment_celltype_genotype_list)


GO_matrix <- data.frame(c("xx", "xxx"))
colnames(GO_matrix) <- "Term"

for (i in names(Enrichment_celltype_genotype_list)){
  comparison_name <- gsub( "Markers_", "", i)
  enrich_RES_merged <- Enrichment_celltype_genotype_list[[i]]
  GO_matrix_to_add <- enrich_RES_merged[,c("Term", "Adjusted.P.value")]
  colnames(GO_matrix_to_add)[2] <- comparison_name
  GO_matrix <- unique(merge(GO_matrix, GO_matrix_to_add, by.x = "Term", by.y ="Term", all.x = TRUE, all.y = TRUE ))  
}

dim(GO_matrix) #
GO_matrix[,1] <- gsub(",",";", GO_matrix[,1])
GO_matrix$db <- GO_matrix$Term
GO_matrix$db <- gsub( ".* WP.*", "WP", GO_matrix$db)
GO_matrix$db <- gsub( ".* \\(GO.*", "GOBP", GO_matrix$db)
GO_matrix <- GO_matrix[GO_matrix$Term != "xx",]
GO_matrix <- GO_matrix[GO_matrix$Term != "xxx",]
rownames(GO_matrix) <- GO_matrix$Term

colnames(GO_matrix) <- gsub( "_C5024T", " m.5024C>T", colnames(GO_matrix))
colnames(GO_matrix) <- gsub( "_A5019G", " m.5019A>G", colnames(GO_matrix))

WP_matrix <- GO_matrix[GO_matrix$db == "WP", -c(1,ncol(GO_matrix))]
GOBP_matrix <- GO_matrix[GO_matrix$db == "GOBP", -c(1,ncol(GO_matrix))]





message("+-------------------------------------------------------------------------------")
message("+                       Figure 4 F heatmap WP                                   ")
message("+-------------------------------------------------------------------------------")


message("+                  enrichment heatmap for WP                   ")

GO_matrix3  <- WP_matrix
GO_res_name <- "WikiPathways"

rownames(GO_matrix3) <- gsub( "regulation" ,  "reg." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "regulated" ,  "reg." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "dependent" ,  "dep." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "positive" ,  "+" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "negative" ,  "-" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "pathway" ,  "path." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "chemical" ,  "chem." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "polymerase" ,  "pol" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "cotranslational" ,  "cotransl." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "biosynthetic" ,  "biosynth." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "mitochondrial" ,  "MT" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "calcium" ,  "Ca" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "Homo sapiens.*" ,  "" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( " WP.*" ,  "" , rownames(GO_matrix3))

GO_matrix3[is.na(GO_matrix3)] <- 1

GO_matrix3 <- GO_matrix3[,order(colnames(GO_matrix3))]


colnames(GO_matrix3)  <- gsub( "Hind", " hind", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Crest", " crest", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Progenitors", " progenitors", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Tube", " tube", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Mesoderm", " mesoderm", colnames(GO_matrix3) )

split = gsub(".* ", "", colnames(GO_matrix3))
colnames(GO_matrix3)  <- gsub( " m.*", "", colnames(GO_matrix3) )

f1 = colorRamp2( c(0, 0.0001, 0.05, 0.051, 0.5, 1), c("#006d2c",  "#2ca25f", "#66c2a4", "white","white", "lightgrey"), space = "RGB") 
ht1 = Heatmap(as.matrix(GO_matrix3),  col = f1, name = GO_res_name,  row_title = "", column_title = "", show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE ,  row_dend_side = "right", row_names_side ="left", width = unit(ncol(GO_matrix3)/2.4, "cm"), height = unit(nrow(GO_matrix3)/2.4, "cm"), column_split = split, row_names_gp = gpar(fontsize = 12)) 


print(ht1)

pdf(paste( Project, "Figure_4F", GO_res_name, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix3)/2+5), height=nrow(GO_matrix3)/3) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()





message("+-------------------------------------------------------------------------------")
message("+                       Figure S5 B heatmap BP                                  ")
message("+-------------------------------------------------------------------------------")

message("+               rrvgo reduce terms                ")

library(rrvgo)
Threshold <- 0.7
GO_res_name <- "GOBP"

GOBP_ids <- (rownames(GOBP_matrix))
GOBP_ids <- gsub( "\\).*", "", GOBP_ids)
GOBP_ids <- gsub( ".*\\(", "", GOBP_ids)

simMatrix_BP <-   calculateSimMatrix(GOBP_ids,
                                     orgdb="org.Mm.eg.db",
                                     ont=c("BP"),
                                     method="Rel")

reducedTerms_BP <- reduceSimMatrix(simMatrix_BP, 
                                   scores = NULL,
                                   threshold=Threshold,
                                   orgdb="org.Mm.eg.db")



length(unique(reducedTerms_BP$parentTerm)) # 76
length(unique(reducedTerms_BP$term)) # 760



dim(GOBP_matrix)
GO_matrix3 <- GOBP_matrix[apply(GOBP_matrix, 1, function(y) !all(is.na(y))),]
dim(GO_matrix3)


GO_matrix3$ID <- gsub( ".*\\(","" ,rownames(GO_matrix3))
GO_matrix3$ID <- gsub( "\\).*","" ,GO_matrix3$ID)
GO_matrix3$parentTerm <- reducedTerms_BP[match(GO_matrix3$ID, reducedTerms_BP$go),]$parentTerm
GO_matrix3$na_count <- apply(GO_matrix3[,-c(ncol(GO_matrix3), ncol(GO_matrix3)-1)], 1, function(x) sum(is.na(x)))
GO_matrix3 <- GO_matrix3[order(GO_matrix3$na_count, decreasing = FALSE),]
GO_matrix3 <- GO_matrix3[!duplicated(GO_matrix3$parentTerm),]

GO_matrix3[is.na(GO_matrix3)] <- 1
GO_matrix3$parentID <- reducedTerms_BP[match(GO_matrix3$ID, reducedTerms_BP$go),]$parent
rownames(GO_matrix3) <- paste0(GO_matrix3$parentTerm, " ", GO_matrix3$parentID)
GO_matrix3 <- GO_matrix3[, -c(ncol(GO_matrix3), ncol(GO_matrix3)-1, ncol(GO_matrix3)-2, ncol(GO_matrix3)-3)]
GO_matrix3 <- GO_matrix3[rownames(GO_matrix3) != "1 NA",]


rownames(GO_matrix3) <- gsub( " \\(GO:.*" ,  "" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "regulation" ,  "reg." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "dependent" ,  "dep." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "positive" ,  "+" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "negative" ,  "-" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "pathway" ,  "path." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "chemical" ,  "chem." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "polymerase" ,  "pol" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "cotranslational" ,  "cotransl." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "biosynthetic" ,  "biosynth." , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "mitochondrial" ,  "MT" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "calcium" ,  "Ca" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "Homo sapiens.*" ,  "" , rownames(GO_matrix3))
rownames(GO_matrix3) <- gsub( "GO.*" ,  "" , rownames(GO_matrix3))

GO_matrix3 <- GO_matrix3[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33, 2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34)]

colnames(GO_matrix3)  <- gsub( "Hind", " hind", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Crest", " crest", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Progenitors", " progenitors", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Tube", " tube", colnames(GO_matrix3) )
colnames(GO_matrix3)  <- gsub( "Mesoderm", " mesoderm", colnames(GO_matrix3) )


split = gsub(".* ", "", colnames(GO_matrix3))
colnames(GO_matrix3)  <- gsub( " m.*", "", colnames(GO_matrix3) )

f1 = colorRamp2( c(0, 0.0001, 0.05, 0.051, 0.5, 1), c("#006d2c",  "#2ca25f", "#66c2a4", "white","white", "lightgrey"), space = "RGB") 

ht1 = Heatmap(as.matrix(GO_matrix3),  col = f1, name = GO_res_name,  row_title = "", column_title = "", show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE ,  row_dend_side = "right",column_split = split, row_names_side ="left", width = unit(ncol(GO_matrix3)/1.5, "cm"), height = unit(nrow(GO_matrix3)/1.5, "cm")) # width = unit(140, "cm"),
print(ht1)

pdf(paste( Project, "SEURAT_ComplexHeatmap_RRVGO_reduced", Threshold, GO_res_name, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, height=(nrow(GO_matrix3)/2+6), width=ncol(GO_matrix3)/1.5) # /2 for wikipathways
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


ht1t = Heatmap(as.matrix(t(GO_matrix3)),  col = f1, name = GO_res_name,  row_title = "", column_title = "", show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = FALSE ,  row_dend_side = "right", row_split = split, row_names_side ="right", height = unit(ncol(GO_matrix3)/2, "cm"),width  = unit(nrow(GO_matrix3)/2, "cm")) # width = unit(140, "cm"),
print(ht1t)


pdf(paste( Project, "SEURAT_ComplexHeatmap_RRVGO_reduced", Threshold, GO_res_name, "l2fc", l2fc_cutoff, "flipped.pdf", sep="_"), onefile=FALSE, width=(nrow(GO_matrix3)/2+6), height=ncol(GO_matrix3)/1.5) # /2 for wikipathways
par(bg=NA)
draw(ht1t, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()













