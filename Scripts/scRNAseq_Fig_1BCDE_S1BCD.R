#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()


##############################################################################################################
#                                                                                                            #
#   Project: MBU_spb54_005_MouseEmbryo                                                                       #
#   Malwina Prater (mn367@cam.ac.uk), 2022                                                                   #
#   MRC MBU, University of Cambridge                                                                         #
#   Script: scRNA-seq mouse dataset - Wildtype embryos                                                       # 
#                                                                                                            #
##############################################################################################################


message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
library("ggplot2")
library("SCENIC")
library("Seurat")  
library("ggchromatic")
library("ggtree")
library("patchwork") 
})


Project        <- "MBU_spb54_005_Fig1"
baseDir        <- "/Users/xxx/Documents/xxx/xxx/xxx" # replace with your path
setwd(baseDir)


# load in Seurat object with whole mouse dataset
Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))


# correct mt transcript
Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-")


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





message("--------------------------------------------------------------------------------")
message("+                         Fig 1 C - UMAP                                        ")
message("+-------------------------------------------------------------------------------")

#  UMAP With cell types - just WT samples  
Seurat_obj_WT <- subset(Seurat_obj, mouse == "WT")
rm(Seurat_obj)


pdf(paste0( Project, "_C__UMAP_WT_", "celltype", "_HARMONY_SCT_batch_regressed", "_no_labs", ".pdf"), onefile=FALSE, width=5, height=4) 
par(bg=NA)
DimPlot(Seurat_obj_WT, reduction = "umap", group.by = "celltype2", pt.size = 1, cols = celltype_cols, label= FALSE, label.size = 5, repel = TRUE)   + coord_fixed(ratio = 1) + scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank())
dev.off()



message("--------------------------------------------------------------------------------")
message("+                           Fig S4 E & D                                        ")
message("+-------------------------------------------------------------------------------")


pdf(paste0( Project, "_S4E__UMAP_ALL_", "celltype", "_HARMONY_SCT_batch_regressed", "_with_legend", ".pdf"), onefile=FALSE, width=7, height=4) 
par(bg=NA)
DimPlot(Seurat_obj, reduction = "umap", group.by = "celltype2", pt.size = 1, cols = celltype_cols, label= FALSE, label.size = 5, repel = TRUE)   + coord_fixed(ratio = 1) + scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL) + theme( plot.title = element_blank())
dev.off()


pdf(paste0( Project, "_S4D__UMAP_ALL_", "mouseBatch", "_HARMONY_SCT_batch_regressed", "_with_legend", ".pdf"), onefile=FALSE, width=7, height=4) 
par(bg=NA)
DimPlot(Seurat_obj, reduction = "umap", group.by = "mouseBatch", pt.size = 0.2, label= FALSE, label.size = 5, repel = TRUE)   + coord_fixed(ratio = 1) + scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL) + theme( plot.title = element_blank())
dev.off()




message("--------------------------------------------------------------------------------")
message("+                        Fig 1 C -  number of cells                             ")
message("+-------------------------------------------------------------------------------")

celltype_cols_WT <- celltype_cols[1:17]
celltype_cols_WT2 <- celltype_cols_WT[c(1,14,16,13,10,12,2,6,17,4,5,3,7,11,8,9,15)]


df_for_barplot <- as.data.frame(table(Seurat_obj_WT@meta.data$celltype2))
df_for_barplot <- df_for_barplot[df_for_barplot$Freq > 0,]
colnames(df_for_barplot) <- c("cell_type", "number_of_cells")
df_for_barplot$ordering <- order(match(names(celltype_cols_WT2), df_for_barplot$cell_type))
df_for_barplot <- df_for_barplot[order(df_for_barplot$ordering, decreasing = FALSE),]

plt_bar <- ggplot(data = df_for_barplot, aes(x = number_of_cells, y = reorder(cell_type, -ordering) )) +
  geom_bar(stat="identity", fill = (celltype_cols_WT2)) + labs(x ="number of cells", y = "") 


pdf(paste0( Project, "_C__barplot_WT_", "celltype", "_HARMONY_SCT_batch_regressed", "", ".pdf"), onefile=FALSE, width=5.5, height=4) 
par(bg=NA)
plt_bar
dev.off()





message("--------------------------------------------------------------------------------")
message("+                        Fig 1 D -  mtDNA transcripts                           ")
message("+-------------------------------------------------------------------------------")

Idents(Seurat_obj_WT) <- Seurat_obj_WT@meta.data$celltype2
Idents(Seurat_obj_WT) <- factor(Idents(Seurat_obj_WT), levels = rev(c("amnion", "cardiac" , "endothelial" ,            "extraembryonic mesoderm", "forebrain" , "foregut", "mesoderm progenitors", "mid hindbrain",  "mid hindgut", "mixed mesoderm"  ,"neural crest", "neural tube", "notochord","pharyngeal mesoderm" , "placodes" ,"presomitic mesoderm" ,"somitic mesoderm" )))

df_mtDNA <- Seurat_obj_WT@meta.data[,c("celltype2", "percent.mt")]

plt_vln <- VlnPlot(Seurat_obj_WT, features = "percent.mt", cols = celltype_cols_WT2, pt.size = 0) + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) + coord_flip() + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) +  labs(y ="mtDNA transcript levels", x = "cell type") +  theme( axis.text.y=element_blank(),axis.ticks.y=element_blank())


pdf(paste0( Project, "_D__violin_WT_", "celltype", "_HARMONY_SCT_batch_regressed", "xx", ".pdf"), onefile=FALSE, width=4, height=4) 
par(bg=NA)
plt_vln
dev.off()




message("--------------------------------------------------------------------------------")
message("+                        Fig 1 D -  ANOVA for mtDNA transcripts                 ")
message("+-------------------------------------------------------------------------------")

summary(df_mtDNA$percent.mt)
df_mtDNA$celltype2 <- factor(df_mtDNA$celltype2, levels = c("amnion", "cardiac" , "endothelial" , "extraembryonic mesoderm", "forebrain" , "foregut", "mesoderm progenitors", "mid hindbrain",  "mid hindgut", "mixed mesoderm"  ,"neural crest", "neural tube", "notochord","pharyngeal mesoderm" , "placodes" ,"presomitic mesoderm" ,"somitic mesoderm" ))


test_t_1celltype_vs_all <- function(dataset, variable, celltype){
  # Levene test for equal variance between the groups:
  dataset2 <- dataset
  dataset2$celltype2 <- ifelse(dataset2$celltype2 == celltype, celltype, "rest")
  Levene_test_res <- car::leveneTest(dataset2[[variable]] ~ dataset2[["celltype2"]], data = dataset2)
  # ANOVA if equal variance:
  if (Levene_test_res$`Pr(>F)`[1] > 0.05){
    print("Levene test Pval > 0.05: variance between groups is equal. Can proceed to ANOVA.")
    res_pairwise_t_test <- pairwise.t.test(dataset2[[variable]], dataset2[["celltype2"]], p.adjust.method = "BH", pool.sd = TRUE)
    return(res_pairwise_t_test$p.value[1])
  } else {
    print("Levene test Pval < 0.05: variance between groups is NOT equal. Can proceed to alternatives of ANOVA: Welch test, Kruskal-Wallis rank sum and pairwise t-test without assumption of equal variances. ")
    res_pairwise_t_test <- pairwise.t.test(dataset2[[variable]], dataset2[["celltype2"]], p.adjust.method = "BH", pool.sd = FALSE)
    return( res_pairwise_t_test$p.value[1])
  }
}


fig_1D_stats <- data.frame(celltype = c("amnion" ,"mesoderm progenitors" , "neural tube" , "mixed mesoderm" ,  "neural crest" , "mid hindbrain", "notochord" , "placodes", "presomitic mesoderm",     "forebrain" ,"pharyngeal mesoderm", "foregut", "extraembryonic mesoderm", "cardiac" ,"somitic mesoderm" , "endothelial" , "mid hindgut" ))
fig_1D_stats$t_test_p.value <- ""
fig_1D_stats[1,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "amnion")
fig_1D_stats[2,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "mesoderm progenitors")
fig_1D_stats[3,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "neural tube")
fig_1D_stats[4,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "mixed mesoderm")
fig_1D_stats[5,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "neural crest")
fig_1D_stats[6,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "mid hindbrain")
fig_1D_stats[7,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "notochord")
fig_1D_stats[8,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "placodes")
fig_1D_stats[9,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "presomitic mesoderm")
fig_1D_stats[10,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "forebrain")
fig_1D_stats[11,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "pharyngeal mesoderm")
fig_1D_stats[12,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "foregut")
fig_1D_stats[13,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "extraembryonic mesoderm")
fig_1D_stats[14,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "cardiac")
fig_1D_stats[15,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "somitic mesoderm")
fig_1D_stats[16,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "endothelial")
fig_1D_stats[17,2] <- test_t_1celltype_vs_all(dataset = df_mtDNA, variable = "percent.mt", celltype= "mid hindgut")

fig_1D_stats$t_test_p.adjust_BH <- p.adjust(fig_1D_stats$t_test_p.value, method = "BH")
fig_1D_stats$star <- ifelse(fig_1D_stats$t_test_p.adjust_BH < 0.05, "*", "NS")

write.csv(fig_1D_stats, "fig_1D_stats.csv")






message("--------------------------------------------------------------------------------")
message("+               related to Fig 1E & S1B -  stats - which celltype DEGs?         ")
message("+-------------------------------------------------------------------------------")


Idents(Seurat_obj_WT) <- Seurat_obj_WT@meta.data$celltype2

Markers_WT_celltype <- FindAllMarkers(Seurat_obj_WT)
Markers_WT_celltype2 <- Markers_WT_celltype[Markers_WT_celltype$p_val_adj < 0.05,]

Markers_WT_celltype2 <- Markers_WT_celltype2[order(Markers_WT_celltype2$p_val_adj),]
Markers_WT_celltype2 <- Markers_WT_celltype2[order(Markers_WT_celltype2$cluster),]


MitoCarta3_genes_mouse2 <- MitoCarta3_genes_mouse[MitoCarta3_genes_mouse %in% rownames(GetAssayData(Seurat_obj_WT))]
mito_genes_not_changing <- unique(MitoCarta3_genes_mouse2[!MitoCarta3_genes_mouse2 %in% Markers_WT_celltype2$gene])
mito_genes_not_changing <- mito_genes_not_changing[order(mito_genes_not_changing)]
mito_genes_changing     <- unique(MitoCarta3_genes_mouse2[MitoCarta3_genes_mouse2 %in% Markers_WT_celltype2$gene])
mito_genes_changing     <- mito_genes_changing[order(mito_genes_changing)]
length(mito_genes_changing)
length(mito_genes_not_changing)
mito_genes_changing <- as.data.frame(mito_genes_changing)
mito_genes_not_changing <- as.data.frame(mito_genes_not_changing)

Markers_WT_celltype2 <- Markers_WT_celltype2[order(Markers_WT_celltype2$gene),]
Markers_WT_celltype_selected <- Markers_WT_celltype2[Markers_WT_celltype2$gene %in% c("Timm50", "Ndufa8", "Cox8a", "Polg", "Gldc", "Bnip3"),]

write.csv(mito_genes_changing, "MBU_spb54_005__Seurat_obj_WT__Fig1E_mito_genes_changing_celltype_padj0.05.csv")
write.csv(mito_genes_not_changing, "MBU_spb54_005__Seurat_obj_WT__Fig1E_mito_genes_not_changing_celltype_padj0.05.csv")


Markers_WT_celltype_mitocarta <- Markers_WT_celltype2[Markers_WT_celltype2$gene %in% MitoCarta3_genes_mouse2,]
write.csv(Markers_WT_celltype_mitocarta, "MBU_spb54_005__Seurat_obj_WT__all_markers_MITOCARTA_celltype_padj0.05.csv")
table(Markers_WT_celltype_mitocarta$gene)
Markers_WT_celltype_mitocarta <- Markers_WT_celltype_mitocarta[order(Markers_WT_celltype_mitocarta$p_val_adj),]




message("--------------------------------------------------------------------------------")
message("+                Fig 1 E -  mtDNA related gene expression                       ")
message("+-------------------------------------------------------------------------------")

Markers <- read.csv( "Input/MBU_spb54_005__Seurat_obj_WT__all_markers_celltype_padj0.05.csv")

metadata_WT <- Seurat_obj_WT@meta.data
metadata_WT$cell <- rownames(metadata_WT)

selected_genes <- c("Timm44", "Cox10", "Ndufb5")
selected_genes <- c("Grsf1", "Slc25a4", "Chchd10")

col_genes <- c("Grsf1"="darkred", "Slc25a4"="darkred", "Chchd10"="darkred", "Timm44"="blue", "Ndufb5"="blue", "Cox10"="blue")

expr_mat_scaled <- GetAssayData(Seurat_obj_WT, assay = "SCT", slot = "scale.data")
expr_mat_scaled_Vln <- expr_mat_scaled[rownames(expr_mat_scaled) %in% selected_genes,]
expr_mat_scaled_Vln <- as.data.frame(t(expr_mat_scaled_Vln))
expr_mat_scaled_Vln[1:5,]
expr_mat_scaled_Vln$cell_type <- metadata_WT[match( rownames(expr_mat_scaled_Vln), rownames(metadata_WT)),]$celltype2

expr_mat_scaled_Vln_molten <- reshape2::melt(expr_mat_scaled_Vln)
colnames(expr_mat_scaled_Vln_molten)[2] <- "gene"
expr_mat_scaled_Vln_molten$cell_type <- factor(expr_mat_scaled_Vln_molten$cell_type, levels = (c("amnion", "cardiac" , "endothelial" ,            "extraembryonic mesoderm", "forebrain" , "foregut", "mesoderm progenitors", "mid hindbrain",  "mid hindgut", "mixed mesoderm"  ,"neural crest", "neural tube", "notochord","pharyngeal mesoderm" , "placodes" ,"presomitic mesoderm" ,"somitic mesoderm" )))


#expr_mat_scaled_Vln_molten$fill <- ifelse(expr_mat_scaled_Vln_molten$gene %in% c("Timm44","Ndufb5","Cox10"), "darkred", "blue")
expr_mat_scaled_Vln_molten$fill <- ifelse(expr_mat_scaled_Vln_molten$gene %in% c("Grsf1","Slc25a4","Chchd10"), "darkred", "blue")


min(expr_mat_scaled_Vln_molten$value)
max(expr_mat_scaled_Vln_molten$value)

vln1 <- ggplot(expr_mat_scaled_Vln_molten, aes(x=cell_type, y=value, fill = fill)) + 
  geom_violin( ) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs( x = "", y = "") + ylim(min(expr_mat_scaled_Vln_molten$value), 5) +
  theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) + 
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~gene, ncol = 1, strip.position = "right", scales = "fixed") 

vln1

#pdf(paste0( Project, "_E__violin_WT_", "celltype", "_HARMONY_SCT_batch_regressed", "_v3", ".pdf"), onefile=FALSE, width=4.5, height=5) 
#par(bg=NA)
vln1
#dev.off()


vln1 <- ggplot(expr_mat_scaled_Vln_molten, aes(x=cell_type, y=value)) + 
  geom_violin( fill = "blue") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs( x = "", y = "") + ylim(min(expr_mat_scaled_Vln_molten$value), 5) +
  theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) + 
  theme( axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~gene, ncol = 1, strip.position = "left", scales = "fixed") + 
  theme(strip.text.y.left = element_text(angle = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank() ) 

vln1

pdf(paste0( Project, "_E__violin_WT_", "celltype", "_HARMONY_SCT_batch_regressed", "Grsf1_Chchd10_Slc25a4", ".pdf"), onefile=FALSE, width=4.5, height=5) 
par(bg=NA)
vln1
dev.off()







message("--------------------------------------------------------------------------------")
message("+                Fig 1 B -  Pdgra Pou3f1 Foxa1 UMAP                             ")
message("+-------------------------------------------------------------------------------")

expr_mat_scaled <- GetAssayData(Seurat_obj_WT, slot = "scale.data", assay = "SCT")
expr_mat_scaled <- as.data.frame(t(expr_mat_scaled))
expr_mat_scaled[1:5,1:5]


matrix.umap <- as.data.frame(Embeddings(object=Seurat_obj_WT, reduction="umap"))
matrix.meta <- Seurat_obj_WT@meta.data
matrix.umap$mouse <- matrix.meta$mouse
matrix.umap$celltype <- matrix.meta$celltype2
head(matrix.umap)

sum(rownames(matrix.umap) == rownames(expr_mat_scaled))
matrix.umap$Pou3f1 <- expr_mat_scaled$Pou3f1
matrix.umap$Pdgfra <- expr_mat_scaled$Pdgfra
matrix.umap$Foxa1 <- expr_mat_scaled$Foxa1

summary(matrix.umap$Pou3f1)
summary(matrix.umap$Pdgfra)
summary(matrix.umap$Foxa1)

matrix.umap$Pou3f1 <- ifelse(matrix.umap$Pou3f1 > 12, 12, matrix.umap$Pou3f1)
matrix.umap$Pdgfra <- ifelse(matrix.umap$Pdgfra > 7, 7, matrix.umap$Pdgfra)
matrix.umap$Foxa1 <- ifelse(matrix.umap$Foxa1 > 13, 13, matrix.umap$Foxa1)



plt_rgb <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour = rgb_spec(Pou3f1, Pdgfra, Foxa1))) + 
  geom_point() + coord_fixed() +
  #scale_y_continuous(breaks=c(-10,0,10)) + scale_x_continuous(breaks=c(-10,0,10)) +
  theme(legend.position="none") + 
  annotate("text", x = -10.8, y = 7, label = "Pou3f1", colour = "red", size = 5) +
  annotate("text", x = 13, y = 2, label = "Pdgfra", colour = "green3", size = 5) +
  annotate("text", x = 2, y = -12.5, label = "Foxa1", colour = "blue", size = 5) +
  scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank())




plt_cmy <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour = cmy_spec(Pou3f1, Pdgfra, Foxa1))) + 
  geom_point(size = 1) + coord_fixed() +
  scale_y_continuous(breaks=c(-10,0,10)) + scale_x_continuous(breaks=c(-10,0,10)) +
  theme(legend.position="none") + 
  annotate("text", x = -10.8, y = 7, label = "Pou3f1", colour = "cyan3", size = 5) +
  annotate("text", x = 13, y = 2, label = "Pdgfra", colour = "orchid", size = 5) +
  annotate("text", x = 2, y = -12.5, label = "Foxa1", colour = "yellow3", size = 5) +
  scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank())




pdf(paste(Project, "UMAP_WT_", "Pou3f1_Pdgfra_Foxa1","cmy_spec_enh.pdf", sep="_"), width=5, height=4)
par(bg=NA)
plt_cmy
dev.off()

pdf(paste(Project, "UMAP_WT_", "Pou3f1_Pdgfra_Foxa1","rgb_spec_enh.pdf", sep="_"), width=5, height=4)
par(bg=NA)
plt_rgb
dev.off()







message("--------------------------------------------------------------------------------")
message("+                Fig 1 F - Cox6a1, Cox6a2 Dotplot                             ")
message("+-------------------------------------------------------------------------------")

selected_genes <- c( "Cox6a1", "Cox6a2")

expr_mat_scaled <- GetAssayData(Seurat_obj_WT, slot = "scale.data", assay = "SCT")
expr_mat_scaled[1:5,1:5]

AvgExpr_scaled <- AverageExpression(Seurat_obj_WT, slot = "scale.data", assay = "SCT")
AvgExpr_scaled <- AvgExpr_scaled[[1]]
AvgExpr_scaled_Dotplot <- AvgExpr_scaled[rownames(AvgExpr_scaled) %in% selected_genes,]


clust <- hclust(dist(t(AvgExpr_scaled_Dotplot) %>% as.matrix(), method = "canberra"), method="complete") # hclust with distance matrix

ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram, right = FALSE)
ggtree_plot
clust$labels[clust$order]

pdf(paste(Project, "DotPlot_WT", "Cox6a1_Cox6a2","with_dendro.pdf", sep="_"), width=6, height=5)
par(bg=NA)
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.2,2), align = 'h')
dev.off()


Idents(Seurat_obj_WT) <-  Seurat_obj_WT@meta.data$celltype2
dotplot <- DotPlot(Seurat_obj_WT, features = selected_genes, cols = c("white","darkred")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "")   + theme(axis.text.x=element_text( face = "italic"))

pdf(paste(Project, "DotPlot_WT", "Cox6a1_Cox6a2","no_dendro.pdf", sep="_"), width=5, height=5)
par(bg=NA)
dotplot
dev.off()


#We also have Cox7a1/Cox7a2/Cox7a2l and Cox8a/Cox8b/Cox8c dot plots in the supplement. I was wondering whether the Cox7a one might look better with just Cox7a1 and Cox7a2 (Cox7a2l is thought to have a slightly different function, so we can probably just compare the ubiquitous (7a2) and muscle/heart (7a1) isoforms on the dot plot.


dotplot <- DotPlot(Seurat_obj_WT, features = c("Cox7a1", "Cox7a2"), cols = c("white","darkred")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "")   + theme(axis.text.x=element_text( face = "italic")) 

pdf(paste(Project, "DotPlot_WT", "Cox7a1_Cox7a2_Cox7a2l","no_dendro.pdf", sep="_"), width=5, height=5)
par(bg=NA)
dotplot
dev.off()


dotplot <- DotPlot(Seurat_obj_WT, features = c("Cox8a"), cols = c("white","darkred")) + theme(axis.text.x=element_text( face = "italic")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "")    # , "Cox8b", "Cox8c"

pdf(paste(Project, "DotPlot_WT", "Cox8a","no_dendro.pdf", sep="_"), width=5, height=5)
par(bg=NA)
dotplot
dev.off()


dotplot <- DotPlot(Seurat_obj_WT, features = c("Cox4i1", "Cox4i2"), cols = c("white","darkred")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "")    + theme(axis.text.x=element_text( face = "italic"))

pdf(paste(Project, "DotPlot_WT", "Cox4i1_Cox4i2","no_dendro.pdf", sep="_"), width=5, height=5)
par(bg=NA)
dotplot
dev.off()





selected_genes2 <-  c("Cox6a1","Cox6a2", "Cox7a2","Cox7a1", "Cox8a", "Cox4i1", "Cox4i2", "Ndufa4","Ndufa4l2")# "Cox8b","Cox8c",

DotPlot(Seurat_obj_WT, features = selected_genes2, cols = c("white","darkred")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "")    + theme(axis.text.x=element_text( face = "italic"))

















