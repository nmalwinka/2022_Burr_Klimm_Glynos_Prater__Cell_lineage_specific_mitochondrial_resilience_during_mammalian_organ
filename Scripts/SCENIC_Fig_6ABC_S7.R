#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()

library("SCENIC")
library("ComplexHeatmap")
library("Seurat") 
library("grDevices")
library("dichromat")


Project        <- "MBU_spb54_005__SCENIC"
baseDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005"
setwd(baseDir)
dbDir         <- "/Users/mn367/Documents/MBU-Projects/Databases/Mouse_cisTarget_databases"
list.files(dbDir)

SCENIC_Dir <- paste0(baseDir,"/SCENIC_batch_regressed" )
setwd(SCENIC_Dir)



message("+-------------------------------------------------------------------------------")
message("+                        set up color scheme                                    ")
message("+-------------------------------------------------------------------------------")

celltype_cols <- c( "amnion"="plum1", "mesoderm progenitors"="darkolivegreen1", "neural tube"="plum4","mixed mesoderm"="gold2", "neural crest"="purple4", "mid hindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomitic mesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngeal mesoderm"="royalblue","foregut"="violetred3", "extraembryonic mesoderm"="steelblue3", "cardiac"="firebrick" , "somitic mesoderm"= "darkgreen",  "endothelial"="orange1", "mid hindgut"="violetred2", "blood"="black", "PGCs"="grey10")



message("+-------------------------------------------------------------------------------")
message("+                          load in functions                                    ")
message("+-------------------------------------------------------------------------------")



MakeCustomFeaturePlot <- function(regulon2plot, outdir, colorRampBIAS = 1){
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep="")) 
  regulon_name <- print(regulon2plot)
  matrix.umap$selected_gene <- AUC_data2[,regulon2plot]
  matrix.umap$selected_gene <- as.numeric(as.character(matrix.umap$selected_gene))
  print(head(matrix.umap, 2))
  
  #myPalette <- colorRampPalette((brewer.pal(9, "YlGnBu")), bias = colorRampBIAS)
  ####myPalette <- colorRampPalette(colorschemes$GreentoMagenta.16, bias = colorRampBIAS)   
  #sc <- scale_colour_gradientn(colours = myPalette(100) , limits=c(min(matrix.umap$selected_gene), max(matrix.umap$selected_gene))) 
  ####sc <- scale_color_viridis(option = "D", direction = -1) # c= plasma, d= viridis
  
  myPalette <- colorRampPalette(colorschemes$GreentoMagenta.16, bias = colorRampBIAS)  
  myPalette <-  myPalette(100)
  myPalette <- myPalette[!myPalette %in% c("#FFFDFF" ,"#FFFBFF","#FCFFFC","#FFF9FF")]
  sc <- scale_colour_gradientn(colours = myPalette , limits=c(0, max(matrix.umap$selected_gene))) 
  
  umap.regulon    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point( aes(colour=selected_gene), alpha=0.5, size=0.5) +
    sc +  coord_fixed() +    theme_cowplot(12)  +
    theme(title=element_text(size=16), 
          text=element_text(size=18), 
          axis.text=element_text(size=18), 
          axis.title=element_text(size=18)) +
    xlab("UMAP_1") + ylab("UMAP_2") + ggtitle(paste0(regulon2plot, " regulon activity")) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(title=element_text(size=6), text=element_text(size=6), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))+
    theme_cowplot(12)  + theme(legend.title = element_blank()) + scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL)
  
  
  pdf(paste(outdir, Project, "_FeaturePlot_regulon_", regulon2plot, "_v2.pdf", sep=""), width=5, height=4.5)
  #png(paste(outdir, Project, "_FeaturePlot_regulon_", regulon2plot, ".png", sep=""), width=1500, height=500, type = "cairo")
  par(bg=NA)
  print({ umap.regulon})
  dev.off()
  
}



#max(matrix.umap$UMAP_1) - min(matrix.umap$UMAP_1)
#max(matrix.umap$UMAP_2) - min(matrix.umap$UMAP_2)

MakeCustomFeaturePlotWithZoom <- function(regulon2plot, colorRampBIAS = 1, zoom_celltype = NULL, zoom_level = 2, color_scheme = "A"){
  regulon_name <- print(regulon2plot)
  matrix.umap$selected_gene <- AUC_data2[,regulon2plot]
  matrix.umap$selected_gene <- as.numeric(as.character(matrix.umap$selected_gene))
  print(head(matrix.umap, 2))
  matrix.umap <- matrix.umap[order(matrix.umap$selected_gene),]
  
  zoom_coords_1 <- median(matrix.umap[matrix.umap$celltype == zoom_celltype,]$UMAP_1)
  zoom_coords_2 <- median(matrix.umap[matrix.umap$celltype == zoom_celltype,]$UMAP_2) # + 0.5 # (for Sox9 & Sox10)
  umap1_ylimits <- c(zoom_coords_1 - zoom_level, zoom_coords_1 + zoom_level)
  umap2_ylimits <- c(zoom_coords_2 - zoom_level, zoom_coords_2 + zoom_level)
  
  if(color_scheme == "A") {
    myPalette <- c("#e0ecf4","#bfd3e6", "#9ebcda" ,"#8c96c6" ,"#8c6bb1" ,"#88419d" ,"#810f7c","#4d004b" )
    sc <- scale_colour_gradientn(colours = myPalette , limits=c(0, max(matrix.umap$selected_gene))) 
  } else if (color_scheme == "B"){
    myPalette <- colorRampPalette((brewer.pal(9, "YlGnBu")), bias = colorRampBIAS)
    sc <- scale_colour_gradientn(colours = myPalette(100) , limits=c(0, max(matrix.umap$selected_gene))) 
  } else if (color_scheme == "C"){
    myPalette <- colorRampPalette(colorschemes$GreentoMagenta.16, bias = colorRampBIAS)  
    myPalette <-  myPalette(100)
    myPalette <- myPalette[!myPalette %in% c("#FFFDFF" ,"#FFFBFF","#FCFFFC","#FFF9FF")]
    sc <- scale_colour_gradientn(colours = myPalette , limits=c(0, max(matrix.umap$selected_gene))) 
  } else if (color_scheme == "D"){
    sc <- scale_color_viridis(option = "D", direction = 1) 
  } else {"please provide valid color scheme"}
  
  umap.regulon_WT    <- ggplot(subset(matrix.umap, matrix.umap$mouse == "WT"), aes(x=UMAP_1, y=UMAP_2)) +
    geom_point( aes(colour=selected_gene), alpha=0.6, size=1) +
    sc + xlab("") + ylab("") +  
    coord_fixed(ratio = 1, xlim= c(-13, 15), ylim = c(-16, 12) ) + 
    #scale_y_continuous(breaks=c(-10,0,10)) + scale_x_continuous(breaks=c(-10,0,10))+
    theme_cowplot(12)  + xlab("UMAP_1") + ylab("UMAP_2") +
    geom_rect( mapping=aes(xmin=umap1_ylimits[1], xmax=umap1_ylimits[2], ymin=umap2_ylimits[1], ymax=umap2_ylimits[2]), color="grey35", alpha=0, size = 0.2) + scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL)+ 
    theme(legend.position = c(0.03, 0.18)) + labs(color= paste0(regulon_name, " activity")) 
  
  umap.gene_all    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point( aes(colour=selected_gene), alpha=0.6, size=1) +
    sc + xlab("") + ylab("") +  
    #coord_fixed(ratio = 1, xlim= c(-13, 15), ylim = c(-16, 12) ) + 
    scale_y_continuous(breaks=c(-10,0,10)) + scale_x_continuous(breaks=c(-10,0,10))+
    theme_cowplot(12)  +
    theme(title=element_text(size=16), 
          text=element_text(size=18), 
          axis.text=element_text(size=18), 
          axis.title=element_text(size=18)) +
    xlab("UMAP_1") + ylab("UMAP_2") +
    geom_rect( mapping=aes(xmin=umap1_ylimits[1], xmax=umap1_ylimits[2], ymin=umap2_ylimits[1], ymax=umap2_ylimits[2]), color="grey35", alpha=0, size = 0.2) + scale_y_continuous(breaks=NULL) + scale_x_continuous(breaks=NULL)+ 
    theme(legend.position = c(0.03, 0.13)) + labs(color= paste0(regulon_name, " activity")) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  
  matrix.umap_A5019G <- subset(matrix.umap, matrix.umap$mouse == "A5019G")
  matrix.umap_A5019G <- matrix.umap_A5019G[matrix.umap_A5019G$UMAP_2 >= umap2_ylimits[1] & matrix.umap_A5019G$UMAP_2 <= umap2_ylimits[2], ]
  matrix.umap_A5019G <- matrix.umap_A5019G[matrix.umap_A5019G$UMAP_1 >= umap1_ylimits[1] & matrix.umap_A5019G$UMAP_1 <= umap1_ylimits[2], ]
  
  umap.gene_A5019G_zoom    <- ggplot(matrix.umap_A5019G, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point( aes(colour=selected_gene), alpha=0.8, size=1.5) +
    sc +   xlab("") + ylab("") +  coord_fixed(ratio = 1, expand = FALSE) +
    ylim(umap2_ylimits[1], umap2_ylimits[2]) + xlim(umap1_ylimits[1], umap1_ylimits[2])+
    theme( title=element_text(size=6), text=element_text(size=6), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))+
    theme_cowplot(12) + 
    geom_rect( mapping=aes(xmin=umap1_ylimits[1], xmax=umap1_ylimits[2], ymin=umap2_ylimits[1], ymax=umap2_ylimits[2]), color="grey35", alpha=0, size = 0.5) + scale_y_continuous(breaks=NULL) +
    theme(axis.text.x = element_text(face="bold", color="white", size=1),
          axis.text.y = element_text(face="bold", color="white", size=1), 
          axis.ticks.length = unit(0, "cm")) +
    theme( axis.line = element_line(colour = "grey35",  size = 0.5, linetype = "solid"))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  matrix.umap_C5024T <- subset(matrix.umap, matrix.umap$mouse == "C5024T")
  matrix.umap_C5024T <- matrix.umap_C5024T[matrix.umap_C5024T$UMAP_2 >= umap2_ylimits[1] & matrix.umap_C5024T$UMAP_2 <= umap2_ylimits[2], ]
  matrix.umap_C5024T <- matrix.umap_C5024T[matrix.umap_C5024T$UMAP_1 >= umap1_ylimits[1] & matrix.umap_C5024T$UMAP_1 <= umap1_ylimits[2], ]
  
  umap.gene_C5024T_zoom    <- ggplot(matrix.umap_C5024T, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point( aes(colour=selected_gene), alpha=0.8, size=1.5) +
    sc +  coord_fixed(ratio = 1, expand = FALSE) + xlab("") + ylab("") +  
    ylim(umap2_ylimits[1], umap2_ylimits[2]) + xlim(umap1_ylimits[1], umap1_ylimits[2])+
    theme(title=element_text(size=6), text=element_text(size=6), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))+
    theme_cowplot(12) + 
    geom_rect( mapping=aes(xmin=umap1_ylimits[1], xmax=umap1_ylimits[2], ymin=umap2_ylimits[1], ymax=umap2_ylimits[2]), color="grey35", alpha=0, size = 0.5) + scale_y_continuous(breaks=NULL) +
    theme(axis.text.x = element_text(face="bold", color="white", size=1),
          axis.text.y = element_text(face="bold", color="white", size=1), 
          axis.ticks.length = unit(0, "cm"))+
    theme( axis.line = element_line(colour = "grey35",  size = 0.5, linetype = "solid")) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  matrix.umap_WT <- subset(matrix.umap, matrix.umap$mouse == "WT")
  matrix.umap_WT <- matrix.umap_WT[matrix.umap_WT$UMAP_2 >= umap2_ylimits[1] & matrix.umap_WT$UMAP_2 <= umap2_ylimits[2], ]
  matrix.umap_WT <- matrix.umap_WT[matrix.umap_WT$UMAP_1 >= umap1_ylimits[1] & matrix.umap_WT$UMAP_1 <= umap1_ylimits[2], ]
  
  umap.gene_WT_zoom    <- ggplot(matrix.umap_WT, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point( aes(colour=selected_gene), alpha=0.8, size=1.5) +
    sc +  coord_fixed(ratio = 1, expand = FALSE) + xlab("") + ylab("") +  
    ylim(umap2_ylimits[1], umap2_ylimits[2]) + xlim(umap1_ylimits[1], umap1_ylimits[2])+
    theme(title=element_text(size=6), text=element_text(size=6), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))+
    theme_cowplot(12) + 
    geom_rect( mapping=aes(xmin=umap1_ylimits[1], xmax=umap1_ylimits[2], ymin=umap2_ylimits[1], ymax=umap2_ylimits[2]), color="grey35", alpha=0, size = 0.5) + scale_y_continuous(breaks=NULL) +
    theme(axis.text.x = element_text(face="bold", color="white", size=1),
          axis.text.y = element_text(face="bold", color="white", size=1), 
          axis.ticks.length = unit(0, "cm"))+
    theme( axis.line = element_line(colour = "grey35",  size = 0.5, linetype = "solid"))+ 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  plots_grid_zoom3 <- plot_grid(umap.gene_WT_zoom + theme(legend.position="none"), umap.gene_C5024T_zoom + theme(legend.position="none"), umap.gene_A5019G_zoom + theme(legend.position="none"), align = "hv", ncol = 1, nrow = 3, scale = 1.08, 
                                labels = c('    WT', 'm.5024C>T', 'm.5019A>G'), label_size = 16,
                                label_x = -0.43, label_y = 0.5,
                                hjust = -0.5, vjust = -0.5) # 
  
  plots_grid_final_all <- plot_grid(umap.gene_all , plots_grid_zoom3, ncol = 2, align = "hv", rel_widths = c(2,1), rel_heights = 1, scale = c(0.98,0.96), labels = c(paste0(regulon_name, " regulon activity") , ""), label_size = 20 )
  
  
  pdf(paste( Project, "_FeaturePlot_with_Zoom_","_wSqr_","_regulon_", regulon2plot, "__col", color_scheme, ".pdf", sep=""), width=10, height=6.5)
  par(bg=NA)
  print({ plots_grid_final_all})
  dev.off()
}


MakeCustomFeaturePlotWithZoom("Sox9", zoom_celltype = "neuralCrest", color_scheme = "C")
MakeCustomFeaturePlotWithZoom("Mecom", zoom_celltype = "endothelial", color_scheme = "C")
MakeCustomFeaturePlotWithZoom("Sox10", zoom_celltype = "neuralCrest", color_scheme = "C")










message("+-------------------------------------------------------------------------------")
message("+                             load SCENIC                                       ")
message("+-------------------------------------------------------------------------------")


cellInfo <- read.csv("int/MBU_spb54_005_Flo_SCENIC_cell_info.csv")
head(cellInfo)

rownames(cellInfo) <- cellInfo$X
cellInfo$CellType_Genotype <- paste(cellInfo$celltype, cellInfo$mouse, sep = "_")

cellInfo$CellType_Genotype2 <- gsub( "_", " ", cellInfo$CellType_Genotype)
cellInfo$CellType_Genotype2 <- gsub( "Hind", " hind", cellInfo$CellType_Genotype2)
cellInfo$CellType_Genotype2 <- gsub( "Crest", " crest", cellInfo$CellType_Genotype2)
cellInfo$CellType_Genotype2 <- gsub( "Progenitors", " progenitors", cellInfo$CellType_Genotype2)
cellInfo$CellType_Genotype2 <- gsub( "Tube", " tube", cellInfo$CellType_Genotype2)
unique(cellInfo$CellType_Genotype2)
unique(cellInfo$CellType_Genotype)

scenicOptions <- readRDS("int/scenicOptions_SCT_batch_regrMBU_spb54_005_Flo_SCENIC.Rds")
exprMat  <- readRDS( paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC__SCT_counts_batch_regressed.Rds"))
Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))






message("--------------------------------------------------------------------------------")
message("+                            plot AUC heatmap                                   ")
message("+-------------------------------------------------------------------------------")

selected_TFs <- c("Hand1", "Hoxc6", "Pax6", "Hoxd1", "Pax7", "Pou3f2", "T","Id4","Hoxa5","Six3","Twist1","Six2","Gata3","Ovol","Pitx1","Nkx2−5","Gata4","Nkx3−1","Irf2","Foxa3","Tfap2a","Nkx1−2","Lef1","Gbx2","Foxc2","Foxc1", "Sox9", "Sox10", "Nkx1-2", "Lhx1", "Mecom", "Gata5", "Gata6", "Taf1", "E2f3", "Maz", "Rad21","Cnot3","Bclaf1","E2f5","Klf12" )

selected_TFs <- c("Taf1", "E2f3", "Maz", "Rad21","Cnot3","Bclaf1","E2f5","Klf12")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUCx <- regulonAUC
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

GROUPS <- "CellType_Genotype2"
#GROUPS <- "celltype"
#GROUPS <- "mouseBatch"

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo[[GROUPS]]), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale=TRUE))

head(regulonActivity_byCellType_Scaled,2)
rownames(regulonActivity_byCellType_Scaled)


all_TFs <- as.data.frame(rownames(regulonActivity_byCellType_Scaled))
colnames(all_TFs)[1] <- "TF"
all_TFs$TF_short <- all_TFs$TF
all_TFs$TF_short <- gsub( " \\(.*", "", all_TFs$TF_short)
all_TFs$TF_short <- gsub( "_extended", "", all_TFs$TF_short)



genes_for_heatmap <-  all_TFs$TF
genes_for_heatmap <-  all_TFs[all_TFs$TF_short %in% selected_TFs,]$TF


regulonActivity_byCellType_Scaled2 <- regulonActivity_byCellType_Scaled[rownames(regulonActivity_byCellType_Scaled) %in% genes_for_heatmap,]

colnames(regulonActivity_byCellType_Scaled2) <- gsub( "_", " ", colnames(regulonActivity_byCellType_Scaled2))
colnames(regulonActivity_byCellType_Scaled2) <- gsub( "Hind", " hind", colnames(regulonActivity_byCellType_Scaled2))
colnames(regulonActivity_byCellType_Scaled2) <- gsub( "Mesoderm", " mesoderm", colnames(regulonActivity_byCellType_Scaled2))
colnames(regulonActivity_byCellType_Scaled2) <- gsub( "Crest", " crest", colnames(regulonActivity_byCellType_Scaled2))
colnames(regulonActivity_byCellType_Scaled2) <- gsub( "Progenitors", " progenitors", colnames(regulonActivity_byCellType_Scaled2))
colnames(regulonActivity_byCellType_Scaled2) <- gsub( "Tube", " tube", colnames(regulonActivity_byCellType_Scaled2))
colnames(regulonActivity_byCellType_Scaled2) <- gsub( "C5024T", "m.5024C>T", colnames(regulonActivity_byCellType_Scaled2))
colnames(regulonActivity_byCellType_Scaled2) <- gsub( "A5019G", "m.5019A>G", colnames(regulonActivity_byCellType_Scaled2))

#rownames(regulonActivity_byCellType_Scaled2) <- gsub( " \\(.*", "", rownames(regulonActivity_byCellType_Scaled2))
#rownames(regulonActivity_byCellType_Scaled2) <- gsub( "_extended", "", rownames(regulonActivity_byCellType_Scaled2))

col_split <- gsub( " WT", "", colnames(regulonActivity_byCellType_Scaled2))
col_split <- gsub( " m.5024C>T", "", col_split)
col_split <- gsub( " m.5019A>G", "", col_split)
col_split <- gsub( " &", "", col_split)
ha = columnAnnotation(CellType = col_split, col = list(CellType = celltype_cols))

df <- as.data.frame(regulonActivity_byCellType_Scaled2)
colnames(df) <- gsub( ".* ", "", colnames(df))

f2 = circlize::colorRamp2( c(-2, 0, 2.5), c("blue", "white","red"), space = "RGB") #  colorRampPalette(c("blue","white","red"))(100)

ht2 <- ComplexHeatmap::Heatmap(df, name="Scaled \nregulon \nactivity", col = f2, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean", cluster_columns = F, show_row_names = TRUE, width  = unit(ncol(df)/2, "cm"),height  = unit(nrow(df)/2, "cm"), top_annotation = ha)
ht2

pdf(paste( Project, "Fig7a_SCTBatchRegr", "___Pheatmap_regulonAUC",  GROUPS, "_scaled", "_topregulators", "v5.pdf", sep="_"), width=18,height=10)
par(bg=NA)
ht2
dev.off()



idx <- grep("m.5019A>G", colnames(regulonActivity_byCellType_Scaled2))
df_x <- regulonActivity_byCellType_Scaled2[,-idx]

col_split_df <- data.frame(group= colnames(df_x) )
col_split_df$celltype <- gsub( " WT", "", col_split_df$group)
col_split_df$celltype <- gsub( " m.5024C>T", "", col_split_df$celltype)
col_split_df$celltype <- gsub( " m.5019A>G", "", col_split_df$celltype)
col_split_df$celltype <- gsub( " &", "", col_split_df$celltype)
col_split_df$genotype <- gsub( ".* ", "", col_split_df$group)
mouse_cols <- c( "m.5024C>T"= "firebrick", "WT"="darkolivegreen4", "m.5019A>G"="dodgerblue4")

ha = columnAnnotation(CellType = col_split_df$celltype, Genotype = col_split_df$genotype, col = list(CellType = celltype_cols, Genotype=mouse_cols))

df_xx <- as.data.frame(df_x)
colnames(df_xx) <- col_split_df$celltype

f2 = circlize::colorRamp2( c(-2, 0, 2.5), c("blue", "white","red"), space = "RGB") #  colorRampPalette(c("blue","white","red"))(100)

ht2 <- ComplexHeatmap::Heatmap(df_xx, name="Scaled \nregulon \nactivity", col = f2, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean", cluster_columns = F, show_column_names = T, show_row_names = TRUE, width  = unit(ncol(df)/3.3, "cm"),height  = unit(nrow(df)/2, "cm"), top_annotation = ha)
ht2

setwd(baseDir)
pdf(paste( Project, "Fig7a_SCTBatchRegr", "___Pheatmap_regulonAUC",  GROUPS, "_scaled", "_topGenotypeRegulators", "v2.pdf", sep="_"), width=15,height=5)
par(bg=NA)
ht2
dev.off()




message("--------------------------------------------------------------------------------")
message("+                  plot binarised AUC heatmap                                   ")
message("+-------------------------------------------------------------------------------")

setwd(SCENIC_Dir)

minPerc <- .5

GROUPS <- "CellType_Genotype2"
#GROUPS <- "celltype"
#GROUPS <- "mouse"

binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")

setwd(baseDir)


cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells[[GROUPS]]), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))

binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
colnames(binaryActPerc_subset) <- gsub(" A5019G", " m.5019A>G", colnames(binaryActPerc_subset))
colnames(binaryActPerc_subset) <- gsub(" C5024T", " m.5024C>T", colnames(binaryActPerc_subset))
colnames(binaryActPerc_subset) <- gsub(" WT", " WT", colnames(binaryActPerc_subset))
binaryActPerc_subset <- as.data.frame(binaryActPerc_subset)

colnames(binaryActPerc_subset) <- gsub( "_", " ", colnames(binaryActPerc_subset))
colnames(binaryActPerc_subset) <- gsub( "Hind", " hind", colnames(binaryActPerc_subset))
colnames(binaryActPerc_subset) <- gsub( "Mesoderm", " mesoderm", colnames(binaryActPerc_subset))
colnames(binaryActPerc_subset) <- gsub( "Crest", " crest", colnames(binaryActPerc_subset))
colnames(binaryActPerc_subset) <- gsub( "Progenitors", " progenitors", colnames(binaryActPerc_subset))
colnames(binaryActPerc_subset) <- gsub( "Tube", " tube", colnames(binaryActPerc_subset))

col_split <- gsub( " WT", "", colnames(binaryActPerc_subset))
col_split <- gsub( " m.5024C>T", "", col_split)
col_split <- gsub( " m.5019A>G", "", col_split)
col_split <- gsub(" &", "", col_split)

ha = HeatmapAnnotation(CellType = col_split, col = list(CellType = celltype_cols))

binaryActPerc_subset2 <- binaryActPerc_subset*100

ComplexHeatmap::Heatmap(binaryActPerc_subset2, name="Regulon activity (%)", col = c("white","pink","red"), cluster_columns = TRUE,  bottom_annotation = ha)


pdf(paste(Project, "ComplexHeatmap",  "binaryActPerc_subset", minPerc , GROUPS ,".pdf", sep="_"), width=16,height=25)
par(bg=NA)
ComplexHeatmap::Heatmap(binaryActPerc_subset2, name="Binarised \nRegulon \nactivity (%)", col = c("white","pink","red"), cluster_columns = TRUE , width = unit(23, "cm"), height = unit(50,"cm"),  bottom_annotation = ha, row_names_rot = 0, column_names_rot = 270)
dev.off()


binaryActPerc_t <- as.data.frame(t(binaryActPerc_subset2))
ha = rowAnnotation(CellType = col_split, col = list(CellType = celltype_cols))

pdf(paste(Project, "ComplexHeatmap",  "binaryActPerc_subset", minPerc , GROUPS ,"_t.pdf", sep="_"), width=25,height=16)
par(bg=NA)
ComplexHeatmap::Heatmap(binaryActPerc_t, name="Binarised \nRegulon \nactivity (%)", col = c("white","pink","red"), cluster_columns = TRUE , width = unit(50, "cm"), height = unit(20,"cm"),  right_annotation = ha, row_labels = gsub(".* ", "", rownames(binaryActPerc_t)  )) 
dev.off()



message("--------------------------------------------------------------------------------")
message("+                               RSS plot                                        ")
message("+-------------------------------------------------------------------------------")

unique(cellInfo$celltype)
cellInfo$celltype <- factor(cellInfo$celltype, levels = c("endothelial", "amnion", "cardiac","midHindgut","foregut", "notochord" ,"blood","mesodermProgenitors","presomiticMesoderm", "extraembryonicMesoderm" ,  "mixedMesoderm", "pharyngealMesoderm" , "somiticMesoderm", "placodes","neuralCrest", "neuralTube", "forebrain", "midHindbrain"))

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "celltype"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

for(i in levels(cellInfo$celltype)){
  #pdf(paste( Project, "_plotRSS_oneSet",  "celltype", i, ".pdf", sep="_"), width=5, height=5)
  #par(bg=NA)
  print({plotRSS_oneSet(rss, setName = i)})
  #dev.off()
}

plotRSS_oneSet(rss, setName = "endothelial")



rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "mouse"])
rssPlot <- plotRSS(rss)
rssPlot

cellInfo$mouse2 <- gsub( "C5024T", "Het", cellInfo$mouse)
cellInfo$mouse2 <- gsub( "A5019G", "Het", cellInfo$mouse2)

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "mouse2"])
rssPlot <- plotRSS(rss)
rssPlot





message("--------------------------------------------------------------------------------")
message("+                   Plot specific Regulons on UMAP                              ")
message("+-------------------------------------------------------------------------------")

AUC_data <- as.data.frame(t(regulonAUC@assays@data$AUC))

matrix.umap <- as.data.frame(Embeddings(object=Seurat_obj, reduction="umap"))
matrix.meta <- Seurat_obj@meta.data
matrix.umap$mouse <- matrix.meta$mouse
matrix.umap$celltype <- matrix.meta$celltype
head(matrix.umap)

AUC_data2 <- AUC_data
colnames(AUC_data2) <- gsub(" \\(.*", "", colnames(AUC_data2))
colnames(AUC_data2) <- gsub("_extended", "", colnames(AUC_data2))
length(unique(colnames(AUC_data2))) == length((colnames(AUC_data2)))




MakeCustomFeaturePlot("Gata4")
MakeCustomFeaturePlot("Gata5")
MakeCustomFeaturePlot("Gata6")
MakeCustomFeaturePlot("Foxa3")




















