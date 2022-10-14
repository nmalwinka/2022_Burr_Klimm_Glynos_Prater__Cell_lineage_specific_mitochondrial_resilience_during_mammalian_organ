#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()

library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Seurat") 
library("biomaRt")



Project        <- "MBU_spb54_005__fig__PSEUDOBULK"
baseDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005"
setwd(baseDir)

Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))

mouse_cols <- c( "m.5024C>T"= "firebrick", "WT"="darkolivegreen4", "m.5019A>G"="dodgerblue4")




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
message("+                         Load in pseudobulk markers                            ")
message("+-------------------------------------------------------------------------------")

Idents(Seurat_obj) <- Seurat_obj@meta.data$mouse

#Markers_C5024T <- FindMarkers(Seurat_obj, ident.1 = "m.5024C>T", ident.2 = "WT")
#Markers_A5019G <- FindMarkers(Seurat_obj, ident.1 = "m.5019A>G", ident.2 = "WT")
#Markers_C5024T_vs_A5019G <- FindMarkers(Seurat_obj, ident.1 = "m.5024C>T", ident.2 = "m.5019A>G")

Markers_C5024T <- read.csv("Input/MBU_spb54_005_Markers_C5024T_Seurat.csv", row.names = 1)
Markers_A5019G <- read.csv("Input/MBU_spb54_005_Markers_A5019G_Seurat.csv", row.names = 1)
Markers_C5024T_vs_A5019G <- read.csv("Input/MBU_spb54_005_Markers_C5024T_vs_A5019G_Seurat.csv", row.names = 1)
Markers_C5024T$gene <- rownames(Markers_C5024T)
Markers_A5019G$gene <- rownames(Markers_A5019G)
Markers_C5024T_vs_A5019G$gene <- rownames(Markers_C5024T_vs_A5019G)

All_mouse_markers <- rbind(Markers_A5019G, Markers_C5024T)



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


# just downregulated marker genes:::
selected_genes <- genes_mootha_buff[genes_mootha_buff %in% All_mouse_markers[All_mouse_markers$avg_log2FC < 0,]$gene]




message("+-------------------------------------------------------------------------------")
message("+                          Fig 4C   dotplot                                     ")
message("+-------------------------------------------------------------------------------")

Idents(Seurat_obj) <- Seurat_obj@meta.data$mouse
#Idents(Seurat_obj) <- factor(Idents(Seurat_obj), levels = c("WT", "m.5024C>T", "m.5019A>G"))
Idents(Seurat_obj) <- factor(Idents(Seurat_obj), levels = c( "m.5019A>G", "m.5024C>T","WT"))


selected_genes2 <- selected_genes[selected_genes %in% Markers_C5024T$gene | selected_genes %in% Markers_A5019G$gene]
selected_genes_df <- data.frame(gene=selected_genes2)
selected_genes_df$col_genes <- ifelse(selected_genes_df$gene %in% c(MitoCarta3_genes_mouse, "Atp5b", "Atp5d", "Atp5a1"), "green4", "black")

plt_dot <- DotPlot(Seurat_obj, features = selected_genes_df$gene, cols = c("white", "darkred") ) + theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) + labs(x = "", y = "") + theme(axis.text.x=element_text(colour=selected_genes_df$col_genes , face = "italic", angle = 90, vjust = 1, hjust=1 ))# + coord_flip()


pdf(paste( Project, "DotPlot","fig_4C", "Mootha_buff_genes", "markers_down_shared.pdf", sep="_"), width=8,height=2)
par(bg=NA)
plt_dot
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                           Fig 3D                                              ")
message("+-------------------------------------------------------------------------------")

Seurat_obj_placodes <- subset(Seurat_obj, celltype2 == "placodes")
Seurat_obj_forebrain <- subset(Seurat_obj, celltype2 == "forebrain")
Seurat_obj_somitic_mesoderm <- subset(Seurat_obj, celltype2 == "somitic mesoderm")


PT_SIZE <- 0.0
SLOT    <-"data" 


vln1  <- VlnPlot(Seurat_obj_placodes, features = "Wnt6", cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none") +  labs(title= "placodes", y ="Wnt6", x = "") + geom_boxplot(width=0.1)

vln2  <- VlnPlot(Seurat_obj_forebrain, features = "Six3", cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none") +  labs(title= "forebrain", y ="Six3", x = "") + geom_boxplot(width=0.1)

vln3  <- VlnPlot(Seurat_obj_somitic_mesoderm, features = "Tcf15", cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none") +  labs(title= "somitic mesoderm", y ="Tcf15", x = "") + geom_boxplot(width=0.1)


pdf(paste( Project, "3D_VlnPlot_data",".pdf", sep="_"), width=7,height=3.5) # "_celltype_regulators",
par(bg=NA)
plot_grid(vln1,vln2,vln3 , ncol = 3)
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                           Fig 4E                                              ")
message("+-------------------------------------------------------------------------------")

Idents(Seurat_obj) <- Seurat_obj@meta.data$mouse
Idents(Seurat_obj) <- factor(Idents(Seurat_obj), levels = c("WT", "m.5024C>T", "m.5019A>G"))

selected_genes <- c("mt-Co1", "mt-Rnr1", "mt-Rnr2", "Sdha", "Uqcrc1", "Mtch2")

PT_SIZE <- 0.0
SLOT    <-"data" 


vln1  <- VlnPlot(Seurat_obj, features = selected_genes[1], cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) +  geom_boxplot(width=0.1) + ggtitle( selected_genes[1]) +  labs(y="", x = selected_genes[1])

vln2  <- VlnPlot(Seurat_obj, features = selected_genes[2], cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) +  geom_boxplot(width=0.1) + ggtitle( selected_genes[1]) +  labs(y="", x = selected_genes[2])

vln3  <- VlnPlot(Seurat_obj, features = selected_genes[3], cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) +  geom_boxplot(width=0.1) + ggtitle( selected_genes[1]) +  labs(y="", x = selected_genes[3])

vln4  <- VlnPlot(Seurat_obj, features = selected_genes[4], cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) +  geom_boxplot(width=0.1) + ggtitle( selected_genes[1]) +  labs(y="", x = selected_genes[4])

vln5  <- VlnPlot(Seurat_obj, features = selected_genes[5], cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) +  geom_boxplot(width=0.1) + ggtitle( selected_genes[1]) +  labs(y="", x = selected_genes[5])

vln6  <- VlnPlot(Seurat_obj, features = selected_genes[6], cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank()) +  geom_boxplot(width=0.1) + ggtitle( selected_genes[1]) +  labs(y="",  x = selected_genes[6])


pdf(paste( Project, "4E_VlnPlot_data","2.pdf", sep="_"), width=14,height=3.5) # "_celltype_regulators",
par(bg=NA)
plot_grid(vln1,vln2,vln3,vln4,vln5,vln6 , ncol = 6)
dev.off()

vln_grid  <- VlnPlot(Seurat_obj, features = selected_genes, cols = mouse_cols, slot = SLOT, assay = "SCT", pt.size= PT_SIZE, ncol = 6) +  geom_boxplot(width=0.1) #+ labs(x = "")   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.title = element_blank(), legend.position="none", plot.title = element_blank())  +

pdf(paste( Project, "4E_VlnPlot_data","2.pdf", sep="_"), width=14,height=4) # "_celltype_regulators",
par(bg=NA)
vln_grid
dev.off()


message("+-------------------------------------------------------------------------------")
message("+                           Fig 5A stats                                        ")
message("+-------------------------------------------------------------------------------")

sel_df <- All_mouse_markers[All_mouse_markers$gene %in%  c("mt-Co1", "mt-Rnr1", "mt-Rnr2", "Sdha", "Uqcrc1", "Mtch2"),]
sel_df <- sel_df[order(sel_df$gene),-7]

write.csv(sel_df , "Fig_5A_vln_stats.csv")









message("+-------------------------------------------------------------------------------")
message("+                           Fig 3B   mtDNA levels                               ")
message("+-------------------------------------------------------------------------------")

# store mitochondrial percentage in object meta data
Seurat_obj <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-", col.name = "percent.mt")

Seurat_obj@meta.data$mouse <- factor(Seurat_obj@meta.data$mouse, levels = c( "WT", "m.5024C>T", "m.5019A>G" ))
Seurat_obj@meta.data$mouse2 <- gsub( "WT", "aWT", Seurat_obj@meta.data$mouse)


Idents(Seurat_obj) <- Seurat_obj@meta.data$celltype_genotype
Seurat_obj_no_blood <- subset(Seurat_obj, celltype2 != "blood")
Seurat_obj_no_blood@meta.data$celltype2 <- factor(Seurat_obj_no_blood@meta.data$celltype2, levels = c( "amnion","cardiac","endothelial","extraembryonic mesoderm","forebrain" ,"foregut" , "mesoderm progenitors", "mid hindbrain", "mid hindgut","mixed mesoderm",  "neural crest" ,"neural tube" , "notochord", "pharyngeal mesoderm", "placodes", "presomitic mesoderm","somitic mesoderm"    ))
Seurat_obj_no_blood@meta.data$mouse <- factor(Seurat_obj_no_blood@meta.data$mouse, levels = c( "WT", "m.5024C>T", "m.5019A>G" ))

pdf(paste( Project, "3B_VlnPlot_data",".pdf", sep="_"), width=9,height=6) 
par(bg=NA)
VlnPlot(Seurat_obj_no_blood, features = "percent.mt",   slot = "data", assay = "RNA", pt.size= PT_SIZE,  group.by =  "celltype2", split.by = "mouse", cols = rep(c("darkolivegreen4",  "firebrick", "dodgerblue4"), 17)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  labs(title= "", y ="mtDNA transcript levels", x = "") 
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                     Fig 3B mtDNA levels STATS                                 ")
message("+-------------------------------------------------------------------------------")

# within cell type only!

df_mtDNA <- Seurat_obj_no_blood@meta.data[,c("celltype2", "percent.mt", "mouse")]
summary(df_mtDNA$percent.mt)

test_t_within_genotype <- function(dataset, variable, celltype){
  dataset <- dataset[dataset$celltype2 == celltype,]
  Levene_test_res <- car::leveneTest(dataset[[variable]] ~ dataset[["mouse"]], data = dataset)
  if (Levene_test_res$`Pr(>F)`[1] > 0.05){
    res_pairwise_t_test <- pairwise.t.test(dataset[[variable]], dataset[["mouse"]], p.adjust.method = "BH", pool.sd = TRUE)
    res_pairwise_t_test <- res_pairwise_t_test$p.value
    t_test_vs_WT <- as.data.frame(t(res_pairwise_t_test[,1]))
    rownames(t_test_vs_WT) <- celltype
    return(t_test_vs_WT)
  } else {
     res_pairwise_t_test <- pairwise.t.test(dataset[[variable]], dataset[["mouse"]], p.adjust.method = "BH", pool.sd = FALSE)
    res_pairwise_t_test <- res_pairwise_t_test$p.value
    t_test_vs_WT <- as.data.frame(t(res_pairwise_t_test[,1]))
    rownames(t_test_vs_WT) <- celltype
    return(t_test_vs_WT)
  }
}


fig_3B_stats <- test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "amnion")
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "mesoderm progenitors"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "neural tube"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "mixed mesoderm"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "neural crest"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "mid hindbrain"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "placodes"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "notochord"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "presomitic mesoderm"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "forebrain"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "pharyngeal mesoderm"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "foregut"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "extraembryonic mesoderm"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "cardiac"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "somitic mesoderm"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "endothelial"))
fig_3B_stats <- rbind(fig_3B_stats, test_t_within_genotype(dataset = df_mtDNA, variable = "percent.mt",  celltype= "mid hindgut"))

fig_3B_stats$star_5019_adj <- p.adjust(fig_3B_stats$`m.5019A>G`)
fig_3B_stats$star_5024_adj <- p.adjust(fig_3B_stats$`m.5024C>T`)
fig_3B_stats$star_5019 <- ifelse(fig_3B_stats$star_5019_adj < 0.05, "*", "ns")
fig_3B_stats$star_5024 <- ifelse(fig_3B_stats$star_5024_adj < 0.05, "*", "ns")
fig_3B_stats <- fig_3B_stats[order(rownames(fig_3B_stats)),]

#write.csv(fig_3B_stats, "fig_3B_stats.csv")




message("+-------------------------------------------------------------------------------")
message("+          Stacked bar plot for cell type proportions- fig 3C                   ")
message("+-------------------------------------------------------------------------------")

metadata <- Seurat_obj@meta.data[,c("orig.ident","batch","mouse","celltype")]
metadata$celltype_genotype <- paste0(metadata$celltype, "_",metadata$mouse)

names(celltype_cols)
CellTypesList <- c("placodes", "forebrain",  "mid hindbrain","neural tube","neural crest",  "mesoderm progenitors" ,"presomitic mesoderm" ,"mixed mesoderm" , "somitic mesoderm", "pharyngeal mesoderm", "extraembryonic mesoderm", "endothelial" ,"cardiac","amnion" ,"notochord","foregut", "mid hindgut"      )

CellTypeFILL <- celltype_cols[names(celltype_cols) != "blood"]
CellTypeFILL <- CellTypeFILL[match(CellTypesList, names(CellTypeFILL))]

celltype.means.melt <- as.data.frame(table(metadata[,c("celltype", "mouse")]))
colnames(celltype.means.melt) <- c("CellType", "Mouse", "Proportion_of_cells")
celltype.means.melt <- celltype.means.melt[celltype.means.melt$CellType != "blood",]
celltype.means.melt$CellType <- gsub( "Hind", " hind", celltype.means.melt$CellType)
celltype.means.melt$CellType <- gsub( "Mesoderm", " mesoderm", celltype.means.melt$CellType)
celltype.means.melt$CellType <- gsub( "Crest", " crest", celltype.means.melt$CellType)
celltype.means.melt$CellType <- gsub( "Progenitors", " progenitors", celltype.means.melt$CellType)
celltype.means.melt$CellType <- gsub( "Tube", " tube", celltype.means.melt$CellType)

celltype.means.melt$CellType <- factor(celltype.means.melt$CellType, levels =  CellTypesList )
celltype.means.melt <- celltype.means.melt[order(celltype.means.melt$CellType),] 
celltype.means.melt$Mouse <- factor(celltype.means.melt$Mouse, levels = c( "WT","m.5019A>G", "m.5024C>T"))
celltype.means.melt$CellType <- factor(celltype.means.melt$CellType, levels =  CellTypesList )

p4 <- ggplot() + geom_bar(aes(y = Proportion_of_cells, x = Mouse, fill = CellType), colour = "black", data = celltype.means.melt,
                          stat="identity", position="fill") +  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) + theme(legend.title = element_blank()) + labs(x="", y="Fraction of cells") + scale_fill_manual(values=(CellTypeFILL)) +
  theme(axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = 0.5, vjust = 0, face = "plain")) +
  theme(legend.text=element_text(size=16))
p4


pdf(paste(Project, "_fig_3C_CellType_stacked_barplot.pdf", sep=""), width=6.5,height=5, onefile=FALSE)
par(bg=NA)
p4
dev.off()








message("+-------------------------------------------------------------------------------")
message("+                           Fig S4F   heatmap                                   ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 0.4


Markers_mito <-  All_mouse_markers[abs(All_mouse_markers$avg_log2FC) > l2fc_cutoff & All_mouse_markers$gene %in% MitoCarta3_genes_mouse,]



selected_genes <- unique(Markers_mito$gene)
selected_genes <- unique(c(selected_genes, "mt-Rnr1", "mt-Rnr2"))


expr_mat <- GetAssayData(Seurat_obj, assay = "SCT", slot = "data")
expr_mat[1:5,1:5]


expr_mat <- as.matrix(expr_mat[rownames(expr_mat) %in% selected_genes,])
mat <- as.data.frame(t(scale(t(expr_mat))))
min(mat)
max(mat)
mat[1:5,1:5]

length(selected_genes) == nrow(expr_mat)
mat <- mat[match(selected_genes,  rownames(expr_mat)),]

mat <- mat[,order(colnames(mat))]

col_split <- as.data.frame(colnames(mat))
colnames(col_split) <- "cell"
col_split$mouse <- Seurat_obj@meta.data[match(col_split$cell, rownames(Seurat_obj@meta.data)),]$mouse
col_split$mouseBatch <- Seurat_obj@meta.data[match(col_split$cell, rownames(Seurat_obj@meta.data)),]$mouseBatch
unique(col_split$mouse)
col_split$mouse <- factor(col_split$mouse, levels = c("WT","m.5024C>T","m.5019A>G" ))

ha = HeatmapAnnotation(Mouse = col_split$mouse, col = list(Mouse = mouse_cols))

set.seed(14)
f1 = circlize::colorRamp2( c(-1.5, 0, 1.5), c("blue","white","red"), space = "LUV") 
ht1 <- ComplexHeatmap::Heatmap(mat,col = f1, name="Expression", cluster_rows = TRUE,  cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE, column_split = col_split$mouse, top_annotation = ha, cluster_column_slices = TRUE, row_title_side = "left", row_names_side = "left" , row_names_gp = gpar(fontface = "italic")) 

pdf(paste( Project, "SCTBatchRegr","fig_S4F", "___heatmap_",  "_SCT_data_scaled", "mito_genes_rnr", ".pdf", sep="_"), width=10,height=6) # "_celltype_regulators",
par(bg=NA)
ht1
dev.off()







message("+-------------------------------------------------------------------------------")
message("+                           Fig 3A   KEGG                                       ")
message("+-------------------------------------------------------------------------------")


ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)



CalculateKeggEnrichment <- function(RESULTS_TABLE, LOG2FOLDCHANGE=0.25 ){
  Kegg_genes <- na.omit(RESULTS_TABLE)
  Kegg_genes$Gene_ID <- ensEMBL2id[match( Kegg_genes$gene, ensEMBL2id$external_gene_name),]$entrezgene_id
  Kegg_genes <- na.omit(Kegg_genes)
  Kegg_genes <- Kegg_genes[abs(Kegg_genes$avg_log2FC) > LOG2FOLDCHANGE,]
  
  foldchanges = Kegg_genes$avg_log2FC
  names(foldchanges) = Kegg_genes$Gene_ID
  foldchanges <- sort(foldchanges, decreasing = T)
  head(foldchanges)
  
  kk_all <- enrichKEGG(names(foldchanges), organism="mmu", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg")
  #head(summary(kk_all))
  kk_all <- setReadable(kk_all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  
  kk_down <- enrichKEGG(names(foldchanges[foldchanges<0]), organism="mmu", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg") 
  #head(summary(kk_down))
  kk2_down <- setReadable(kk_down, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kk_res_down <- as.data.frame(kk2_down)
  kk_res_down$direction <- "down"
  
  kk_up <- enrichKEGG(names(foldchanges[foldchanges>0]), organism="mmu", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg", ) 
  #head(summary(kk_up))
  kk2_up <- setReadable(kk_up, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kk_res_up <- as.data.frame(kk2_up)
  kk_res_up$direction <- "up"
  
  
  kk_results <- rbind(kk_res_up, kk_res_down)
  
  return(as.data.frame(kk_all))
}

kk_all_Markers_A5019G_vs_WT <- CalculateKeggEnrichment(Markers_A5019G, LOG2FOLDCHANGE=0.3)
kk_all_Markers_C5024T_vs_WT <- CalculateKeggEnrichment(Markers_C5024T, LOG2FOLDCHANGE=0.3)
kk_all_Markers_C5024T_vs_A5019G <- CalculateKeggEnrichment(Markers_C5024T_vs_A5019G, LOG2FOLDCHANGE=0.3)


write.csv(CalculateKeggEnrichment(Markers_A5019G, LOG2FOLDCHANGE=0.3), paste0("MBU_spb54_005_", "_enrichKEGG_", "Markers_A5019G", "_l2fc_","0.3" ,".csv"))
write.csv(CalculateKeggEnrichment(Markers_C5024T, LOG2FOLDCHANGE=0.3), paste0("MBU_spb54_005_", "_enrichKEGG_", "Markers_C5024T", "_l2fc_","0.3" ,".csv"))
write.csv(CalculateKeggEnrichment(Markers_C5024T_vs_A5019G, LOG2FOLDCHANGE=0.3), paste0("MBU_spb54_005_", "_enrichKEGG_", "Markers_C5024T_vs_A5019G", "_l2fc_","0.3" ,".csv"))




PlotKeggPathways <- function(kk_results, RESULTS_TABLE, selected_kegg, NO_OF_PATHWAYS_TO_PLOT= 13, col_UP="red",col_DOWN="blue", plot_title=NA  ){
  PLOT_NAME <- deparse(substitute(kk_results))
  PLOT_NAME <- gsub( "kk_all_", "", PLOT_NAME)
  PLOT_NAME <- gsub( "_", " ", PLOT_NAME)
  enrichKegg_selected <- kk_results[kk_results$ID %in% selected_kegg,]
  enrichKegg_selected$Description <- gsub("endoplasmic reticulum", "ER"  , enrichKegg_selected$Description)
  enrichKegg_selected <- enrichKegg_selected[order(enrichKegg_selected$qvalue),]
  if(nrow(enrichKegg_selected)>NO_OF_PATHWAYS_TO_PLOT){
    enrichKegg_selected <- enrichKegg_selected[c(1:NO_OF_PATHWAYS_TO_PLOT),]
  }
  
  list_up <- list()
  list_down <- list()
  for (i in 1:nrow(enrichKegg_selected)){
    df_tmp <- RESULTS_TABLE[RESULTS_TABLE$gene %in% unlist(strsplit(enrichKegg_selected[i, "geneID"] , split="/")),]
    tmp_up <- length(subset(df_tmp[,2], df_tmp[,2] > 0))
    list_up[[i]] <-tmp_up
    tmp_down <- length(subset(df_tmp[,2], df_tmp[,2] < 0))
    list_down[[i]] <-tmp_down
  }
  
  enrichKegg_selected$genes_UP <- as.numeric(list_up)
  enrichKegg_selected$genes_DOWN <- as.numeric(list_down)
  enrichKegg_selected$genes_DOWN <- -enrichKegg_selected$genes_DOWN
  enrichKegg_molten <- melt(enrichKegg_selected[,c(1,2,6,7,10:11)], id.vars=c("ID","Description","p.adjust", "qvalue") )
  enrichKegg_molten$Description <- gsub( "Chemical","chem.", enrichKegg_molten$Description)
  enrichKegg_molten$Description <- gsub( "reactive oxygen species","ROS", enrichKegg_molten$Description)
  enrichKegg_molten$Description <- gsub( "- multiple diseases","", enrichKegg_molten$Description)
  enrichKegg_molten$Description <- gsub( "Signaling pathways regulating","Sig. path. reg.", enrichKegg_molten$Description)
  
  p_kegg_mlt <- ggplot(enrichKegg_molten, aes(x=reorder(Description, -qvalue), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue))) +
    coord_flip() +  xlab(" ") +
    scale_fill_manual( values = c(col_DOWN, col_UP)) +  ylab("Gene count") + ggtitle(plot_title) +
    ylim(-max(abs(enrichKegg_molten$value)), max(abs(enrichKegg_molten$value))) +
    #ylim(-15, 15) +
    theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_alpha_continuous( range = c(0.5, 1))
  
  return(p_kegg_mlt)
  
}



selected_kegg <- c("mmu03018","mmu05020","mmu03013","mmu04530","mmu05208","mmu05016" ,"mmu04110", "mmu04120",  "mmu00190", "mmu05010",  "mmu01200", "mmu03008", "mmu00020","mmu04141","mmu04144","mmu05022","mmu03010","mmu03040","mmu03050","mmu03015","mmu00480","mmu00010","mmu04550","mmu04810","mmu00970","mmu04520")



plt1 <- PlotKeggPathways(kk_all_Markers_A5019G_vs_WT, Markers_A5019G, selected_kegg, NO_OF_PATHWAYS_TO_PLOT = 13, col_DOWN="dodgerblue4",col_UP="darkolivegreen4", plot_title="WT vs m.5019A>G")
plt2 <- PlotKeggPathways(kk_all_Markers_C5024T_vs_WT, Markers_C5024T, selected_kegg, NO_OF_PATHWAYS_TO_PLOT = 13, col_DOWN="firebrick",col_UP="darkolivegreen4", plot_title="WT vs m.5024C>T")
plt3 <- PlotKeggPathways(kk_all_Markers_C5024T_vs_A5019G, Markers_C5024T_vs_A5019G, selected_kegg, NO_OF_PATHWAYS_TO_PLOT = 13, col_DOWN="firebrick",col_UP="dodgerblue4", plot_title="m.5019A>G vs m.5024C>T")

cowplot::plot_grid(plt1,plt2,plt3, ncol=3)



pdf(paste( Project, "SCTBatchRegr","fig_3A", "___KEGG_barplots_", "markers_l2fc0.3", "v2.pdf", sep="_"), width=22,height=5)
par(bg=NA)
cowplot::plot_grid(plt1,plt2,plt3, ncol=3)
dev.off()
















