#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()

library("ggplot2")
library("ggrepel")
library("cowplot")
library("Seurat")  
library("dplyr")
library("ComplexHeatmap")
theme_set(theme_cowplot())



Project        <- "MBU_spb54_005__OXPHOS_GENE_CORR"
baseDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Stephen_Burr/MBU_spb54_005"
setwd(baseDir)

Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))





message("+-------------------------------------------------------------------------------")
message("+                         color schemes                                         ")
message("+-------------------------------------------------------------------------------")

celltype_cols <- c( "amnion"="plum1", "mesoderm progenitors"="darkolivegreen1", "neural tube"="plum4","mixed mesoderm"="gold2", "neural crest"="purple4", "mid hindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomitic mesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngeal mesoderm"="royalblue","foregut"="violetred3", "extraembryonic mesoderm"="steelblue3", "cardiac"="firebrick" , "somitic mesoderm"= "darkgreen",  "endothelial"="orange1", "mid hindgut"="violetred2", "blood"="black")

mouse_cols <- c("WT"="darkolivegreen4", "m.5024C>T"= "firebrick",  "m.5019A>G"="dodgerblue4")





message("+-------------------------------------------------------------------------------")
message("+                         Load in mitocarta                                     ")
message("+-------------------------------------------------------------------------------")

MitoCarta3 <- read.csv("/Users/mn367/Documents/MBU-Projects/Databases/Mouse.MitoCarta3.0_summarised.csv")
MitoCarta3_pathways <- read.csv("/Users/mn367/Documents/MBU-Projects/Databases/Mouse.MitoCarta3.0_pathways.csv")
MitoCarta3_genes_mouse <- unique(MitoCarta3$Symbol)

MitoCarta_oxPhos <- MitoCarta3[grep("OXPHOS", MitoCarta3$MitoCarta3.0_MitoPathways), c("Symbol","EnsemblGeneID","MitoCarta3.0_MitoPathways","Synonyms","Description","MitoCarta3.0_SubMitoLocalization","Tissues")]

MitoCarta_oxPhos$ComplexOxphos <- gsub( "OXPHOS > " , "", MitoCarta_oxPhos$MitoCarta3.0_MitoPathways)
MitoCarta_oxPhos$ComplexOxphos <- gsub( " \\| OXPHOS subunits.*" , "", MitoCarta_oxPhos$ComplexOxphos)
MitoCarta_oxPhos$ComplexOxphos <- gsub( " \\| Metabolism.*" , "", MitoCarta_oxPhos$ComplexOxphos)
MitoCarta_oxPhos$ComplexOxphos <- gsub( " > .*" , "", MitoCarta_oxPhos$ComplexOxphos)


Mito_C1 <- data.frame(gene = unique(unlist(strsplit(unlist(MitoCarta3_pathways[MitoCarta3_pathways$MitoPathway %in% c("Complex I","CI subunits", "CI assembly factors"), ]$Genes), ", "))), category = "Complex_I")
Mito_C2 <- data.frame(gene = unique(unlist(strsplit(unlist(MitoCarta3_pathways[MitoCarta3_pathways$MitoPathway %in% c("Complex II","CII subunits", "CII assembly factors"), ]$Genes), ", "))), category = "Complex_II")
Mito_C3 <- data.frame(gene = unique(unlist(strsplit(unlist(MitoCarta3_pathways[MitoCarta3_pathways$MitoPathway %in% c("Complex III","CIII subunits", "CIII assembly factors"), ]$Genes), ", "))), category = "Complex_III")
Mito_C4 <- data.frame(gene = unique(unlist(strsplit(unlist(MitoCarta3_pathways[MitoCarta3_pathways$MitoPathway %in% c("Complex IV","CIV subunits", "CIV assembly factors"), ]$Genes), ", "))), category = "Complex_IV")
Mito_C5 <- data.frame(gene = unique(unlist(strsplit(unlist(MitoCarta3_pathways[MitoCarta3_pathways$MitoPathway %in% c("Complex V","CV subunits", "CV assembly factors"), ]$Genes), ", "))), category = "Complex_V")
Mito_CytC <- data.frame(gene = unique(unlist(strsplit(unlist(MitoCarta3_pathways[MitoCarta3_pathways$MitoPathway == "Cytochrome C", ]$Genes), ", "))), category = "Cytochrome_C")
Mito_Resp <- data.frame(gene = unique(unlist(strsplit(unlist(MitoCarta3_pathways[MitoCarta3_pathways$MitoPathway == "Respirasome assembly", ]$Genes), ", "))), category = "Respirasome_assembly")

Mito_oxphos_mlt <- rbind(Mito_C1, Mito_C2)
Mito_oxphos_mlt <- rbind(Mito_oxphos_mlt, Mito_C3)
Mito_oxphos_mlt <- rbind(Mito_oxphos_mlt, Mito_C4)
Mito_oxphos_mlt <- rbind(Mito_oxphos_mlt, Mito_C5)
Mito_oxphos_mlt <- rbind(Mito_oxphos_mlt, Mito_CytC)
idx2 <- which(Mito_oxphos_mlt$gene == 'Tmem70' & Mito_oxphos_mlt$category == "Complex_I")
Mito_oxphos_mlt <- Mito_oxphos_mlt[-idx2,]
rownames(Mito_oxphos_mlt) <- Mito_oxphos_mlt$gene





message("+-------------------------------------------------------------------------------")
message("+            Correlation of Mito genes per cell type genotype                   ")
message("+-------------------------------------------------------------------------------")


expr_mat <- GetAssayData(Seurat_obj, assay = "SCT", slot = "data")
expr_mat[1:5,1:5]



selected_genes <- unique(Mito_oxphos_mlt$gene)





calc_GENE_corr_heatmap_for_celltype <- function(Seurat_obj_to_use=Seurat_obj, celltype=NULL, mean_expr_val = 0.2, selected_genes_name = "Oxphos"){
  
  dataset_name <- deparse(substitute(Seurat_obj_to_use))
  expr_mat2 <- expr_mat[, colnames(expr_mat) %in% rownames(Seurat_obj_to_use@meta.data[Seurat_obj_to_use@meta.data$celltype2 == celltype,]) ]
  cor_list <- list()
  count <- 1

  for ( i in c("WT", "m.5024C>T" , "m.5019A>G")){
    mat_tmp <- as.matrix(expr_mat2[rownames(expr_mat2) %in% selected_genes, colnames(expr_mat2) %in% rownames(Seurat_obj_to_use@meta.data[Seurat_obj_to_use@meta.data$mouse == i,]) ])
    genes_to_keep <- data.frame(gene  = rownames(mat_tmp), mean_expression = rowMeans(mat_tmp), sum_expression = rowSums(mat_tmp))
    genes_to_keep <- genes_to_keep[genes_to_keep$mean_expression > mean_expr_val,]

    matrix_mod <- t(scale(t(mat_tmp[rownames(mat_tmp) %in% genes_to_keep$gene,])))
    corr_df <- data.frame(matrix(0, nrow = nrow(matrix_mod), ncol = nrow(matrix_mod)))
    colnames(corr_df) <- rownames(matrix_mod)
    rownames(corr_df) <- rownames(matrix_mod)
    
    for (gene_name in colnames(corr_df)){
      gene <- as.numeric(matrix_mod[gene_name,])
      correlations <- (apply(matrix_mod,1,function(x){cor(gene,x)}))
      corr_df[,colnames(corr_df) == gene_name] <- correlations
    }
    cor_list[[count]] <- corr_df
    count <- count + 1
  }
  
  names(cor_list) <- c("WT","m.5024C>T" , "m.5019A>G")

  idx1 <- which(is.na(cor_list[[1]][,1]))
  idx2 <- which(is.na(cor_list[[2]][,1]))
  idx3 <- which(is.na(cor_list[[3]][,1]))
  idx <- c(idx1, idx2, idx3)
  
  print(paste( "remove genes not shared between genotypes:  ", length(idx)) )
  
  if( length(idx) > 0){
    cor_list[[1]] <- cor_list[[1]][-idx,-idx]
    cor_list[[2]] <- cor_list[[2]][-idx,-idx]
    cor_list[[3]] <- cor_list[[3]][-idx,-idx]
  }
  
  cor_list[[1]] <- cor_list[[1]][rownames(cor_list[[1]]) %in% rownames(cor_list[[2]]),] 
  cor_list[[1]] <- cor_list[[1]][rownames(cor_list[[1]]) %in% rownames(cor_list[[3]]),] 
  cor_list[[2]] <- cor_list[[2]][rownames(cor_list[[2]]) %in% rownames(cor_list[[1]]),] 
  cor_list[[2]] <- cor_list[[2]][rownames(cor_list[[2]]) %in% rownames(cor_list[[3]]),] 
  cor_list[[3]] <- cor_list[[3]][rownames(cor_list[[3]]) %in% rownames(cor_list[[1]]),] 
  cor_list[[3]] <- cor_list[[3]][rownames(cor_list[[3]]) %in% rownames(cor_list[[2]]),] 
  cor_list[[1]] <- cor_list[[1]][,colnames(cor_list[[1]]) %in% colnames(cor_list[[2]])] 
  cor_list[[1]] <- cor_list[[1]][,colnames(cor_list[[1]]) %in% colnames(cor_list[[3]])] 
  cor_list[[2]] <- cor_list[[2]][,colnames(cor_list[[2]]) %in% colnames(cor_list[[1]])] 
  cor_list[[2]] <- cor_list[[2]][,colnames(cor_list[[2]]) %in% colnames(cor_list[[3]])] 
  cor_list[[3]] <- cor_list[[3]][,colnames(cor_list[[3]]) %in% colnames(cor_list[[1]])] 
  cor_list[[3]] <- cor_list[[3]][,colnames(cor_list[[3]]) %in% colnames(cor_list[[2]])] 

  col_anno <- data.frame(gene=rownames(cor_list[[1]]))
  col_anno$Complex <- Mito_oxphos_mlt[match(col_anno$gene , Mito_oxphos_mlt$gene),]$category
  col_anno$Complex[is.na(col_anno$Complex)] <- "nuclear_genes"
  col_anno$genome <- ifelse(col_anno$gene %in% col_anno$gene[grep("mt-", col_anno$gene)] , "MT", "nuclear" )

  Complex_cols <- c("Complex_I"="olivedrab1","Complex_II"="skyblue" ,  "Complex_III"="slategrey",  "Complex_IV"="plum1","Complex_V"="mediumpurple" , "Cytochrome_C"="cornsilk")
  Genome_cols <- c("MT"="blue", "nuclear"="green")

  ha = HeatmapAnnotation(Oxphos = col_anno$Complex, Genome= col_anno$genome, col = list(Oxphos = Complex_cols, Genome=Genome_cols),show_legend=FALSE, annotation_label = "")
  ha_side = rowAnnotation(Oxphos = col_anno$Complex, Genome= col_anno$genome , col = list(Oxphos = Complex_cols, Genome=Genome_cols), show_legend=FALSE, annotation_label = "" )

  f1 = circlize::colorRamp2( c( -1,0, 1), c("blue","white","red"), space = "sRGB") 
  lgd1 = Legend(col_fun = f1, title = "Correlation")
  ht_opt$legend_gap <- unit(1, "cm")
  
  count <- 1
  heatmaps_list <- list()
  heatmaps_list2 <- list()

  for (i in names(cor_list)){
    
    ht_name <- i
    heatmaps_list[[count]] <- ComplexHeatmap::Heatmap(cor_list[[i]], col = f1, column_title = paste(celltype, ht_name),  name=ht_name, cluster_rows = TRUE,  cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE,  cluster_column_slices = FALSE, row_title_side = "left", row_names_side = "left" , row_names_gp = gpar(fontface = "italic"), column_names_gp = gpar(fontface = "italic"), clustering_distance_columns = "spearman", clustering_distance_rows = "spearman", top_annotation = ha, left_annotation = ha_side,  row_title_rot = 0, use_raster=FALSE, heatmap_legend_param = list( title = "Correlation" )) # column_split = col_anno$Complex, row_split = col_anno$Complex , cluster_row_slices = FALSE,
    count <- count + 1
  }

  names(heatmaps_list) <- names(cor_list)

  mean_correlations_list <- list()
  count <- 1
  
  for(i in names(cor_list) ){

    corr_means <- cor_list[[count]]
    
    corr_nuclear <- corr_means[rownames(corr_means) %in% col_anno[col_anno$Complex == "nuclear_genes" ,]$gene , colnames(corr_means) %in% col_anno[col_anno$Complex == "nuclear_genes" ,]$gene ]
    
    corr_mt <- corr_means[rownames(corr_means) %in% col_anno[col_anno$Complex != "nuclear_genes" ,]$gene , colnames(corr_means) %in% col_anno[col_anno$Complex != "nuclear_genes" ,]$gene ]
    
    corr_nuclear_mt <- corr_means[rownames(corr_means) %in% col_anno[col_anno$Complex == "nuclear_genes" ,]$gene , colnames(corr_means) %in% col_anno[col_anno$Complex != "nuclear_genes" ,]$gene ]
    mean_correlations <- data.frame(Gene_group = c("nuclear", "mt", "nuclear_vs_mt"), Mean_correlation=c(mean(rowMeans(corr_nuclear)), mean(rowMeans(corr_mt)), mean(rowMeans(corr_nuclear_mt))))
    print(mean_correlations) 
    mean_correlations_list[[count]] <- mean_correlations
    count <- count + 1
  }

  ht1 <- ComplexHeatmap::Heatmap(cor_list[[1]], col = f1, column_title = "",  name="WT", cluster_rows = TRUE,  cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE,  cluster_column_slices = FALSE, row_title_side = "left", row_names_side = "left" , row_names_gp = gpar(fontface = "italic"), column_names_gp = gpar(fontface = "italic"), clustering_distance_columns = "spearman", clustering_distance_rows = "spearman", top_annotation = ha, left_annotation = ha_side,  row_title_rot = 0, use_raster = FALSE, height = unit(10, "cm") , width = unit(10, "cm"), show_heatmap_legend = FALSE) # column_title = paste(celltype, "WT"),
  
  WT_col_order <- column_order(ht1)
  WT_row_order <- row_order(ht1)

  ht2 <- ComplexHeatmap::Heatmap(cor_list[[2]], col = f1, column_title = "",  name="m.5024C>T", cluster_rows = F,  column_order = WT_col_order, show_row_names = FALSE, show_column_names = FALSE,  cluster_column_slices = FALSE, row_title_side = "left", row_names_side = "left" , row_names_gp = gpar(fontface = "italic"), column_names_gp = gpar(fontface = "italic"), clustering_distance_columns = "spearman", clustering_distance_rows = "spearman", top_annotation = ha, left_annotation = ha_side,  row_title_rot = 0, use_raster = FALSE, show_heatmap_legend = FALSE, height = unit(10, "cm"), width = unit(10, "cm")) # column_title = paste(celltype, "m.5024C>T"), 
  
  ht3 <- ComplexHeatmap::Heatmap(cor_list[[3]], col = f1, column_title = "",  name="m.5019A>G", cluster_rows = F,  column_order = WT_col_order, show_row_names = FALSE, show_column_names = FALSE,  cluster_column_slices = FALSE, row_title_side = "left", row_names_side = "left" , row_names_gp = gpar(fontface = "italic"), column_names_gp = gpar(fontface = "italic"), clustering_distance_columns = "spearman", clustering_distance_rows = "spearman", top_annotation = ha, left_annotation = ha_side,  row_title_rot = 0, use_raster = FALSE, show_heatmap_legend = FALSE,  height = unit(10, "cm"), width = unit(10, "cm")) #  column_title = paste(celltype, "m.5019A>G"),
  
  ht_list = ht1 + ht2 + ht3

  pdf(paste( Project, "ComplexHeatmap_Oxphos", selected_genes_name, "cutoff", mean_expr_val, "", dataset_name , celltype,  "3_heatmaps",  "", ".pdf", sep="_"), width=15,height=5.5)
  par(bg=NA)
  draw(ht1 + ht2 + ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1.5, "cm")) # , annotation_legend_list = lgd1
  dev.off()

  return(list(mean_correlations_list, cor_list, ht_list ))
}



ht_cardiac          <- calc_GENE_corr_heatmap_for_celltype(Seurat_obj,"cardiac", mean_expr_val = 0.2, selected_genes_name="oxphos")
ht_mixed_mesoderm   <- calc_GENE_corr_heatmap_for_celltype(Seurat_obj,"mixed mesoderm", mean_expr_val = 0.2, selected_genes_name="oxphos")
ht_neuralcrest      <- calc_GENE_corr_heatmap_for_celltype(Seurat_obj,"neural crest", mean_expr_val = 0.2, selected_genes_name="oxphos")
ht_forebrain          <- calc_GENE_corr_heatmap_for_celltype(Seurat_obj,"forebrain", mean_expr_val = 0.2, selected_genes_name="oxphos")
ht_foregut          <- calc_GENE_corr_heatmap_for_celltype(Seurat_obj,"foregut", mean_expr_val = 0.2, selected_genes_name="oxphos")










plot_binarised_barplot <- function(corr_list=NULL, celltype=NULL, mean_expr_val=NULL, mean_expr_cutoff=NULL, genes_used=NULL, min_corr = 0.5){
  
  xWT_corr_vals   <- reshape2::melt(cbind(gene = rownames(corr_list[[2]]$WT), corr_list[[2]]$WT))
  x5024_corr_vals <- reshape2::melt(cbind(gene =rownames(corr_list[[2]]$`m.5024C>T`), corr_list[[2]]$`m.5024C>T`))
  x5019_corr_vals <- reshape2::melt(cbind(gene =rownames(corr_list[[2]]$`m.5019A>G`), corr_list[[2]]$`m.5019A>G`))
  gene_to_gene <- paste(xWT_corr_vals$gene, xWT_corr_vals$variable)

  cor_wt <- corr_list[[2]]$WT
  cor_5024 <- corr_list[[2]]$`m.5024C>T`
  cor_5019 <- corr_list[[2]]$`m.5019A>G`
  cor_wt[lower.tri(cor_wt)] <- NA 
  cor_5024[lower.tri(cor_5024)] <- NA 
  cor_5019[lower.tri(cor_5019)] <- NA 
  
  cor_wt[cor_wt == 1] <- NA
  cor_wt[cor_wt == "1"] <- NA
  cor_5024[cor_5024 == 1] <- NA
  cor_5024[cor_5024 == "1"] <- NA
  cor_5019[cor_5019 == 1] <- NA
  cor_5019[cor_5019 == "1"] <- NA
  
  cor_wt[cor_wt > min_corr] <- 1
  cor_wt[cor_wt < -min_corr] <- -1
  cor_wt[cor_wt < min_corr & cor_wt  > -min_corr] <- 0
  cor_5024[cor_5024 > min_corr] <- 1
  cor_5024[cor_5024 < -min_corr] <- -1
  cor_5024[cor_5024 < min_corr & cor_5024  > -min_corr] <- 0
  cor_5019[cor_5019 > min_corr] <- 1
  cor_5019[cor_5019 < -min_corr] <- -1
  cor_5019[cor_5019 < min_corr & cor_5019  > -min_corr] <- 0
  

  df <- data.frame(mouse=c("WT", "m.5024C>T", "m.5019A>G"), sum_binarised_correlations = c( sum(cor_wt, na.rm = TRUE), sum(cor_5024, na.rm = TRUE),sum(cor_5019, na.rm = TRUE) ), positive_correlations = c( sum(cor_wt[cor_wt>0], na.rm = TRUE), sum(cor_5024[cor_5024>0], na.rm = TRUE),sum(cor_5019[cor_5019>0], na.rm = TRUE) ), negative_correlations = c( sum(cor_wt[cor_wt<0], na.rm = TRUE), sum(cor_5024[cor_5024<0], na.rm = TRUE),sum(cor_5019[cor_5019<0], na.rm = TRUE) ))
  
  df$mouse <- factor(df$mouse, levels = names(mouse_cols))
  
  df2 <- reshape2::melt(df[,-2])
  df2$value <- abs(df2$value)
  df2$variable <- gsub( "_", " ", df2$variable)
  df2 <- df2[df2$variable != "negative correlations",]

  plt <- ggplot(data=df2, aes(x=mouse, y=value, fill = variable)) + geom_bar(stat="identity", position=position_dodge()) + xlab("") + ylab("Sum binarised correlations")  + scale_fill_manual(values = c("seagreen", "violetred")) + theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # ggtitle(paste(celltype, ":: \nbinarised corr. cut-off ", min_corr , sep = " "))  + ggtitle(paste(celltype))
  
  pdf(paste( Project, "Oxphos_barplot", celltype, selected_genes_name, "SumBinarisedGeneCorrelations", "cutoff", mean_expr_val, "mean_expr_cutoff", min_corr, "",   "", ".pdf", sep="_"), width=3,height=4)
  par(bg=NA)
  print({plt})
  dev.off()
  
  plt <- ggplot(data=df2, aes(x=mouse, y=value, fill = mouse)) + geom_bar(stat="identity", position=position_dodge()) + xlab("") + ylab("Sum binarised correlations")  + scale_fill_manual(values = mouse_cols) + theme(legend.position="none")
  
  pdf(paste( Project, "Oxphos_barplot_mouseCols", celltype, selected_genes_name, "SumBinarisedGeneCorrelations", "Binarised_CorrCutOff", min_corr, "",   "", ".pdf", sep="_"), width=4,height=4)
  par(bg=NA)
  print({plt})
  dev.off()

  return(plt)  
}

selected_genes_name="oxphos"
min_corr_val <- 0.6

barplot_cardiac <- plot_binarised_barplot(corr_list=ht_cardiac, celltype="cardiac",  min_corr=min_corr_val)
barplot_mixedmeso <- plot_binarised_barplot(corr_list=ht_mixed_mesoderm, celltype="mixed_mesoderm",   min_corr=min_corr_val)
barplot_neuralcrest <- plot_binarised_barplot(corr_list=ht_neuralcrest, celltype="neural_crest",   min_corr=min_corr_val)
barplot_forebrain <- plot_binarised_barplot(corr_list=ht_forebrain, celltype="forebrain",   min_corr=min_corr_val)
barplot_foregut <- plot_binarised_barplot(corr_list=ht_foregut, celltype="foregut",   min_corr=min_corr_val)


















