#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()


##############################################################################################################
#                                                                                                            #
#   Project: MBU_spb54_005_MouseEmbryo                                                                       #
#   Malwina Prater (mn367@cam.ac.uk), 2022                                                                   #
#   MRC MBU, University of Cambridge                                                                         #
#   Script: bulk RNA-seq in mouse MEFs clones - DESeq2 analysis                                              # 
#                                                                                                            #
##############################################################################################################

message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
  library("ggplot2")
  library("ggrepel")
  library("cowplot")
  library("DESeq2")
  library("data.table")
  library("ggdendro")
  library("biomaRt")
  library("dplyr")
  library("RColorBrewer")
  library("clusterProfiler")
  library("org.Mm.eg.db")
  library("enrichR")
})



Project        <- "MBU_spb54_006"
baseDir        <- "/Users/xxx/Documents/xxx/xxx/xxx" # replace with your path
setwd(baseDir)


message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)          


message("+-------------------------------------------------------------------------------")
message("+                           prepare sample table                                ")
message("+-------------------------------------------------------------------------------")

meta_data <- read.csv("Input/RNAseq_Sample_List.csv")
files = list.files(paste0(baseDir, "/Input/STAR_Counts"), "*ReadsPerGene.out.tab$", full.names = T)
seq_contents <- read.csv("Input/SLX-22336.HGNFHDRX2.s_1.contents.csv")

sampleTable <- merge(seq_contents, meta_data, by.x="Sample.name", by.y="Name", all= TRUE)

sampleTable$condition <- gsub("Clone [0-9]{3} ", "", sampleTable$Details)
sampleTable$condition <- gsub("Clone [0-9]{2} ", "", sampleTable$condition)
sampleTable$clone <- gsub(" Galactose", "", sampleTable$Details)
sampleTable$clone <- gsub(" High.*", "", sampleTable$clone)
sampleTable$clone <- gsub(" Low.*", "", sampleTable$clone)
sampleTable$clone <- gsub(" ", "_", sampleTable$clone)
sampleTable$Heteroplasmy2 <- as.numeric(as.character(gsub( "%", "", sampleTable$Heteroplasmy)))
sampleTable$Heteroplasmy_group <- ifelse(sampleTable$Heteroplasmy2 > 50, "highHeteroplasmy", "lowHeteroplasmy")
sampleTable$group <- paste(sampleTable$Heteroplasmy_group, sampleTable$condition, sep = "_")
sampleTable$group <- gsub(" ", "", sampleTable$group)
sampleTable$sample_id <- paste0(sampleTable$Barcode,"_s_1")
sampleTable_second_run <- sampleTable
sampleTable_second_run$sample_id <- paste0(sampleTable_second_run$Barcode,"_s_2")

sampleTable <- rbind(sampleTable, sampleTable_second_run)
sampleTable <- sampleTable[,c("sample_id", "group", "condition", "Heteroplasmy2", "Heteroplasmy_group", "clone", "Barcode","Sample.name")] 

rm(sampleTable_final)
rm(sampleTable_second_run)

sampleTable$Heteroplasmy_group <- factor(sampleTable$Heteroplasmy_group, levels = c("lowHeteroplasmy",  "highHeteroplasmy"))
sampleTable$run <- gsub( ".*_s_", "", sampleTable$sample_id)


message("+-------------------------------------------------------------------------------")
message("+               create count data from STAR counts outputs                      ")
message("+-------------------------------------------------------------------------------")

countData = data.frame(fread(files[1]))[c(1,4)]

# Loop and read the 4th column remaining files
for(i in 2:length(files)) {
  countData = cbind(countData, data.frame(fread(files[i]))[4])
}

# Skip first 4 lines, count data starts on the 5th line
countData = countData[c(5:nrow(countData)),]
colnames(countData) = c("GeneID", gsub(paste0(baseDir,"/Input/STAR_Counts/SLX-22336\\."), "", files))
colnames(countData) = gsub(".r_1.fq.gz_AlignReadsPerGene.out.tab", "", colnames(countData))
colnames(countData) = gsub(".HGNFHDRX2.", "_", colnames(countData))
rownames(countData) = countData$GeneID

countData = countData[,c(2:ncol(countData))]
countData[1:10,1:5]
countData_filt <- countData[rowSums(countData) > 10,]

sampleTable <- sampleTable[order(sampleTable$sample_id),]
sampleTable$sample_id == colnames(countData_filt)





message("+-------------------------------------------------------------------------------")
message("+      Build a DESeqDataSet with Collapsed replicates - just HIGH GLUC          ")
message("+-------------------------------------------------------------------------------")

sampleTable <- sampleTable[sampleTable$condition == "High Glucose DMEM",]
sampleTable_coll <- sampleTable[!duplicated(sampleTable$Barcode) , -c(1,9)]
countData_filt <- countData_filt[, colnames(countData_filt) %in% sampleTable$sample_id]

#write.csv(countData_filt,"raw_counts.csv", quote = FALSE )


dds = DESeqDataSetFromMatrix(countData = countData_filt, colData = sampleTable, design = ~ group )
dds <- collapseReplicates(dds, groupby = sampleTable$Barcode)
dds = DESeq(dds, minReplicatesForReplace=5)

cbind(resultsNames(dds))

norm_counts <- counts(dds, normalized=TRUE)
raw_counts <- counts(dds, normalized=FALSE)
#write.csv(norm_counts,"norm_counts.csv", quote = FALSE )
#write.csv(raw_counts,"raw_counts.csv", quote = FALSE )

#saveRDS(dds, "MBU_spb54_006_outliersRM_dds_coll.Rds")

rld = rlog(dds)
vsd = vst(dds)



message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

elementTextSize <- 14
pca = prcomp(t(assay(rld)))
rv = rowVars(assay(rld))

pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

scores <- data.frame(pca$x[,c(1:4)], sampleTable_coll)
scores$Heteroplasmy <- paste0(scores$Heteroplasmy2, "%")

plt_pca1 <- ggplot(scores, aes(x = PC1, y = PC2, col = condition )) + 
  geom_point(size = 4 , alpha = 0.6) + 
  xlab(pc1lab) + ylab(pc2lab) + theme_classic()+
  scale_colour_manual(name="Condition", values = c("grey", "blue", "green")) +
  theme(text = element_text(size=elementTextSize)) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))+
  coord_fixed(ratio = 1, xlim = c(-40,40), ylim = c(-40,40), expand = TRUE, clip = "on")
plt_pca1


plt_pca2 <- ggplot(scores, aes(x = PC1, y = PC2, col = Heteroplasmy_group, label=Heteroplasmy )) + 
  geom_point(size = 4 , alpha = 0.6) + 
  xlab(pc1lab) + ylab(pc2lab) + theme_classic()+
  scale_colour_manual(name="Heteroplasmy", values = c("darkolivegreen4","firebrick")) +
  theme(text = element_text(size=elementTextSize)) + geom_text_repel() +
  coord_fixed(ratio = 1, xlim = c(-45,45), ylim = c(-45,45), expand = TRUE, clip = "on")
plt_pca2


plt_pca3 <- ggplot(scores, aes(x = PC1, y = PC2, col = group )) + 
  geom_point(size = 4, alpha = 0.7 ) + 
  xlab(pc1lab) + ylab(pc2lab) + theme_classic()+
  scale_colour_manual(name="Group", values = group_cols) +
  theme(text = element_text(size=elementTextSize)) + 
  coord_fixed(ratio = 1, xlim = c(-40,40), ylim = c(-40,40), expand = TRUE, clip = "on")
plt_pca3


plt_pca4 <- ggplot(scores, aes(x = PC1, y = PC2, col = clone )) + 
  geom_point(size = 4, alpha = 0.5 ) + 
  xlab(pc1lab) + ylab(pc2lab) + theme_classic()+
  theme(text = element_text(size=elementTextSize)) + 
  geom_text_repel(aes(label = clone), nudge_x = -1, nudge_y = 0.2, size = 2, max.overlaps = 60) +
  coord_fixed(ratio = 1, xlim = c(-40,40), ylim = c(-40,40), expand = TRUE, clip = "on")
plt_pca4


pdf(paste(Project,"QC", "CollapsedReplicates" , "PCA.pdf", sep="_"), onefile=FALSE, width=15, height=10) 
par(bg=NA)
plot_grid(plt_pca1,plt_pca2,plt_pca3,plt_pca4, ncol = 2, nrow = 2, align = "hv")
dev.off()



pdf(paste(Project,"QC", "CollapsedReplicates" , "PCA_CellRev_v3.pdf", sep="_"), onefile=FALSE, width=5, height=4) 
par(bg=NA)
plt_pca2
dev.off()



message("+-------------------------------------------------------------------------------")
message("+                                PCA explained                                  ")
message("+-------------------------------------------------------------------------------")

ensanno <- ensEMBL2id[,c(1:2)]
ensanno <- ensanno[!duplicated(ensanno),]
loadings                         <- as.data.frame(pca$rotation)
loadings$ensembl_gene_id         <- rownames(loadings)
loadings                         <- merge(loadings, ensanno, by="ensembl_gene_id")

pca.1         <-  loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <-  pca.1[c(1:25),]
pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.2         <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <-  pca.2[c(1:25),]
pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

plt_loadings <- plot_grid(pca.1.25.plot, pca.2.25.plot, labels=c(" ", " "), ncol = 1, nrow = 2)

pdf(paste(Project,"QC", "CollapsedReplicates" , "PCA_loadings.pdf", sep="_"), onefile=FALSE, width=8, height=6) 
par(bg=NA)
plt_loadings
dev.off()




message("+-------------------------------------------------------------------------------")
message("+       Hierarchical clustering        ")
message("+-------------------------------------------------------------------------------")

sampleTable_coll$sample_id <- paste(sampleTable_coll$clone, gsub("DMEM", "", sampleTable_coll$group) )
sampleTable_coll$sample_id <- gsub("Heteroplasmy" , "Het", sampleTable_coll$sample_id)

hclust_matrix <- t(assay(rld))
rownames(hclust_matrix) <- sampleTable_coll[match( rownames(hclust_matrix), sampleTable_coll$Barcode),]$sample_id
sample_distances <- dist(hclust_matrix)
dd.row <- as.dendrogram(hclust(sample_distances), rotate = FALSE, segments = TRUE)
ddata_x <- dendro_data(dd.row)

p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

labs <- label(ddata_x)
labs$group <- gsub(".* ", "", labs$label)
labs$label <- gsub( "Clone", ".                                   Clone", labs$label)
labs


plt_hclust <- p2 + geom_text(data=label(ddata_x),  aes(label=labs$label, x=x, y=-20, colour=labs$group), angle = 90) +  theme(legend.position = "none",  panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),  axis.ticks = element_blank()) + scale_colour_manual(name="Group", values = group_cols2)


pdf(paste(Project,"QC", "CollapsedReplicates" , "hclust.pdf", sep="_"), onefile=FALSE, width=10, height=10) 
par(bg=NA)
plt_hclust
dev.off()



message("+-------------------------------------------------------------------------------")
message("+              Heatmap of the sample-to-sample distances                        ")
message("+-------------------------------------------------------------------------------")

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, vsd$clone, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plt_sample_dist <- pheatmap::pheatmap(sampleDistMatrix,
                                      clustering_distance_rows=sampleDists,
                                      clustering_distance_cols=sampleDists,
                                      col=colors, cellwidth = 10, cellheight = 10, annotation_colors = group_cols)

pdf(paste(Project,"QC", "CollapsedReplicates" , "sample_distances.pdf", sep="_"), onefile=FALSE, width=10, height=10) 
par(bg=NA)
plt_sample_dist
dev.off()



message("+-------------------------------------------------------------------------------")
message("+           more QC  : Effects of transformations on the variance               ")
message("+-------------------------------------------------------------------------------")


rld = rlog(dds)
vsd = vst(dds)

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

plotDispEsts(dds)




message("+-------------------------------------------------------------------------------")
message("+                                 DESeq results                                 ")
message("+-------------------------------------------------------------------------------")

cbind(resultsNames(dds))

res_highGlucose <-  lfcShrink(dds, contrast=c("group","highHeteroplasmy_HighGlucoseDMEM","lowHeteroplasmy_HighGlucoseDMEM"), type = "ashr")

res_highGlucose.ann <- as.data.frame(res_highGlucose[order(res_highGlucose$padj),])

res_highGlucose.ann$external_gene_name <- ensEMBL2id[match( rownames(res_highGlucose.ann), ensEMBL2id$ensembl_gene_id),]$external_gene_name

res_highGlucose.ann <- na.omit(res_highGlucose.ann)
resSig_highGlucose <- res_highGlucose.ann[res_highGlucose.ann$padj < 0.05,]





message("+-------------------------------------------------------------------------------")
message("+ Now creating gene plots.")
message("+-------------------------------------------------------------------------------")

elementTextSize <- 8


makeGeneCountPlot_v2 <- function(gene2plot,outdir) {
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep=""))
  gene_ens <- ensEMBL2id[ensEMBL2id$external_gene_name == gene2plot, ]$ensembl_gene_id 
  t2            <- plotCounts(dds, gene=gene_ens, intgroup=c("group"), normalized=TRUE, returnData=TRUE)
  t2$heteroplasmy <-sampleTable_coll[match(rownames(t2) , sampleTable_coll$Barcode),]$Heteroplasmy_group
  t2$condition    <-sampleTable_coll[match(rownames(t2) , sampleTable_coll$Barcode),]$condition
  t2$heteroplasmy_short <- gsub( "highHeteroplasmy", "HH", t2$heteroplasmy)
  t2$heteroplasmy_short <- gsub( "lowHeteroplasmy", "LH", t2$heteroplasmy_short)
  t2$condition_short <- gsub( "Galactose" , "Gal", t2$condition)
  t2$condition_short <- gsub( "Low Glucose DMEM" , "LG", t2$condition_short)
  t2$condition_short <- gsub( "High Glucose DMEM" , "HG", t2$condition_short)
  
  t2$samples <- paste( t2$heteroplasmy_short , rownames(t2))
  t2$heteroplasmy <- gsub( "Het", " het", t2$heteroplasmy)
  t2$heteroplasmy <- factor(t2$heteroplasmy, levels = c("low heteroplasmy" , "high heteroplasmy" ))
  het_cols <- c( "high heteroplasmy"= "firebrick", "low heteroplasmy"="darkolivegreen4")
  
  plt_box <- ggplot(t2, aes(x=heteroplasmy, y=count, fill=heteroplasmy)) + geom_boxplot() + 
    scale_fill_manual(values = het_cols) +stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +  theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(paste( gene2plot, sep="")) +
    labs(x="", y = "Normalised counts") + ylim(0, max(t2$count)*1.1)
  
  
  plt_bar <- ggplot(t2, aes(x=samples, y=count, fill=heteroplasmy)) + 
    geom_bar(stat="identity", alpha=.75, position=position_dodge()) + 
    scale_fill_manual(values = het_cols) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(paste(Project, " ::: ", gene2plot, sep=""))
  
  pdf(paste(outdir, Project, "-DGE_", gene2plot, "_boxplot.pdf", sep=""),width=3,height=4, onefile=FALSE)
  par(bg=NA)
  print({ plt_box})
  dev.off()
  print(paste("Created plot for", gene2plot), sep=" ")
}

makeGeneCountPlot_v2("Klf12")
makeGeneCountPlot_v2("Bclaf1")
makeGeneCountPlot_v2("Bclaf3")
makeGeneCountPlot_v2("Taf1")
makeGeneCountPlot_v2("E2f3")
makeGeneCountPlot_v2("Maz")



message("+-------------------------------------------------------------------------------")
message("+                                volcano plots                                  ")
message("+-------------------------------------------------------------------------------")

elementTextSize <- 8
#results.df <- res_highGlucose.ann

makeVolcano <- function(results.df){
  res_name <- deparse(substitute(results.df)) 
  
  data <- data.frame(gene = rownames(results.df),
                     symbol = results.df$external_gene_name,
                     padj = results.df$padj,
                     pvalue = -log10(results.df$padj), 
                     lfc = results.df$log2FoldChange)
  data <- na.omit(data)
  
  data <- data %>%  dplyr::mutate(color = ifelse(data$lfc > 1 & data$pvalue > 1.3, yes = "High heteroplasmy", no = ifelse(data$lfc < -1 & data$pvalue > 1.3, yes = "Low heteroplasmy", no = "none")))
  # pvalue = 2 for padj = 0.01
  # pvalue = 1.3 for padj = 0.05
  
  colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
    scale_x_continuous( limits = c(-15, 15)) +
    theme(legend.position = "none") + 
    xlab(expression( title = "log[2] NormCounts (Low heteroplasmy/High heteroplasmy)")) + 
    ylab(expression(-log[10]("adjusted p-value"))) +   ggtitle(label = "" , subtitle = "") +  
    geom_vline(xintercept = -1, colour = "black",linetype="dotted") +  
    geom_vline(xintercept = 1, colour = "black",linetype="dotted") + 
    geom_hline(yintercept = 1.3, colour = "black",linetype="dotted") + 
    annotate(geom = "text", label = "Low heteroplasmy", x = -8, y = 40, size = 5, colour = "black") + 
    annotate(geom = "text", label = "High heteroplasmy", x = 8, y = 40, size = 5, colour = "black") + 
    scale_color_manual(values = c("High heteroplasmy" = "indianred", "Low heteroplasmy" = "slateblue", "none" = "#636363")) +  
    theme(text = element_text(size=elementTextSize)) + theme_classic()
  
  Fig_1_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 1 & padj < 0.005), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"), max.overlaps = 30)
  
  Fig_1_volcano
  
  pdf(paste( Project, "volcano",  res_name, ".pdf", sep="_"), onefile=FALSE, width=7, height=7) 
  par(bg=NA)
  print({Fig_1_volcano})
  dev.off()
}

makeVolcano(res_highGlucose.ann)

plotMA(res_highGlucose, ylim=c(-2,2))




message("+-------------------------------------------------------------------------------")
message("+                     heatmap of DEGs                                           ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 1

rld_df <- as.data.frame(assay(rld))
rld_df_meanCentered <- rld_df - rowMeans(rld_df)   

rld_df_meanCentered$ensembl_gene_id <- rownames(rld_df_meanCentered)
rld_df_meanCentered$external_gene_name <- ensEMBL2id[match(rld_df_meanCentered$ensembl_gene_id, ensEMBL2id$ensembl_gene_id),]$external_gene_name

mat <- rld_df_meanCentered
mat <- rld_df_meanCentered[rld_df_meanCentered$external_gene_name %in% res_highGlucose.ann[res_highGlucose.ann$padj < 0.05 & abs(res_highGlucose.ann$log2FoldChange) > l2fc_cutoff,]$external_gene_name , ]


mat <- mat[!duplicated(mat$external_gene_name),]
mat <- mat[!is.na(mat$external_gene_name),]
rownames(mat) <- mat$external_gene_name
mat <- mat[,-c(ncol(mat), ncol(mat)-1)]

sampleTable_coll <- sampleTable_coll[order(sampleTable_coll$Heteroplasmy2),]
sampleTable_coll <- sampleTable_coll[order(sampleTable_coll$condition),]
mat2 <- mat[, match(sampleTable_coll$Barcode, colnames(mat))]

ht_anno <- data.frame(sample= colnames(mat2))
ht_anno$condition <- sampleTable_coll[match(ht_anno$sample, sampleTable_coll$Barcode),]$condition
ht_anno$group <- sampleTable_coll[match(ht_anno$sample, sampleTable_coll$Barcode),]$group
ht_anno$Heteroplasmy_group <- sampleTable_coll[match(ht_anno$sample, sampleTable_coll$Barcode),]$Heteroplasmy_group
ht_anno$Heteroplasmy <- sampleTable_coll[match(ht_anno$sample, sampleTable_coll$Barcode),]$Heteroplasmy2
ht_anno <- ht_anno[order(ht_anno$Heteroplasmy),]
ht_anno <- ht_anno[order(ht_anno$condition),]

ht_anno$group <- as.character(ht_anno$group)
ht_anno$group <- gsub( "Heteroplasmy_HighGlucoseDMEM", " Heteroplasmy", ht_anno$group)
ht_anno$group <- factor(ht_anno$group, levels = c("low Heteroplasmy","high Heteroplasmy"))
group_cols3 <-  c("low Heteroplasmy"="darkolivegreen4", "high Heteroplasmy"="firebrick" )


ha_top = HeatmapAnnotation(Heteroplasmy=anno_points(ht_anno$Heteroplasmy, pch = 16, size = unit(3, "mm")),  show_annotation_name = TRUE)
ha_bottom = HeatmapAnnotation(Group = ht_anno$group, col = list(Group = group_cols3),  show_annotation_name = TRUE)

mat2 <- mat2[, match( ht_anno$sample, colnames(mat2))]
dim(mat2)

ht1 = Heatmap(as.matrix(mat2),  name = "Bulk RNA-seq",  row_title = "", column_title = "MBU_spb54_006 RNA-seq", show_row_names = F, heatmap_legend_param = list(title = "Expression", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE ,  row_names_side ="right", column_split = ht_anno$condition, bottom_annotation = ha_bottom, top_annotation = ha_top, use_raster = TRUE)# col = f1, width = unit(3, "cm"),
ht1


pdf(paste(Project, "ComplexHeatmap",  "DEGs_just_HighGlucose", "v4.pdf", sep="_"), onefile=FALSE, width=6, height=10)
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                 enrichR                                          ")
message("+-------------------------------------------------------------------------------")

enrichR_DB <- as.data.frame(listEnrichrDbs())
db <- c("WikiPathways_2019_Mouse", "GO_Biological_Process_2018", "KEGG_2019_Mouse","GO_Cellular_Component_2018", "GO_Molecular_Function_2018","Reactome_2016")
websiteLive <- TRUE

l2fc_cutoff <- 1

resSig_df <- resSig_highGlucose

DEGs <- resSig_df[(resSig_df$log2FoldChange) > l2fc_cutoff,]$external_gene_name

enrich_res <- enrichr(DEGs, databases = db )
enrich_WP <- enrich_res[[1]][enrich_res[[1]]$Adjusted.P.value < 0.05,]
enrich_GOBP <- enrich_res[[2]][enrich_res[[2]]$Adjusted.P.value < 0.05,]
enrich_Kegg <- enrich_res[[3]][enrich_res[[3]]$Adjusted.P.value < 0.05,]
enrich_GOCC <- enrich_res[[4]][enrich_res[[4]]$Adjusted.P.value < 0.05,]
enrich_GOMF <- enrich_res[[5]][enrich_res[[5]]$Adjusted.P.value < 0.05,]
enrich_React <- enrich_res[[6]][enrich_res[[6]]$Adjusted.P.value < 0.05,]

enrich_WP$gene_count <- as.numeric(as.character(gsub( "\\/.*", "", enrich_WP$Overlap)))
enrich_WP$Description <- gsub( " WP.*", "", enrich_WP$Term)
enrich_WP$Description <- gsub( " WP.*", "", enrich_WP$Description)
enrich_WP$Description <- gsub( "Signaling", "sig.", enrich_WP$Description)
enrich_WP$Description <- gsub( "signaling", "sig.", enrich_WP$Description)
enrich_WP$Description <- gsub( "Development", "dev.", enrich_WP$Description)
enrich_WP$Description <- gsub( "Pathways", "path.", enrich_WP$Description)
enrich_WP$Description <- gsub( "Pathway", "path.", enrich_WP$Description)
enrich_WP$WP_id <- gsub( ".* WP", "WP", enrich_WP$Term)

selected_terms <- c("WP2074","WP339","WP403","WP1496","WP431","WP431","WP523","WP509","WP113","WP447","WP164","WP3857")
remove_terms <- c("WP403","WP3632","WP2432")
col_DOWN="firebrick"
col_UP="darkolivegreen4"

#enrich_WP <- enrich_WP[enrich_WP$WP_id %in% selected_terms,]
enrich_WP <- enrich_WP[!enrich_WP$WP_id %in% remove_terms,]

p_WP_mlt <- ggplot(enrich_WP, aes(x=reorder(Description, -Adjusted.P.value), y=gene_count, fill=-Adjusted.P.value)) + geom_bar(stat="identity", aes(alpha = -Adjusted.P.value), fill = col_UP) +
  coord_flip() +  xlab(" ") +
  ylab("Gene count") + 
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))


pdf(paste( Project,  "___WP_barplots_HighGlucose", "markers_l2fc1", ".pdf", sep="_"), width=8,height=6) # "_celltype_regulators",
par(bg=NA)
p_WP_mlt
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                           Fig 3A   KEGG                                       ")
message("+-------------------------------------------------------------------------------")



ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl) 


CalculateKeggEnrichment <- function(RESULTS_TABLE, LOG2FOLDCHANGE=1 ){
  Kegg_genes <- na.omit(RESULTS_TABLE)
  Kegg_genes$Gene_ID <- ensEMBL2id[match( Kegg_genes$external_gene_name, ensEMBL2id$external_gene_name),]$entrezgene_id
  Kegg_genes <- na.omit(Kegg_genes)
  Kegg_genes <- Kegg_genes[abs(Kegg_genes$log2FoldChange ) > LOG2FOLDCHANGE,]
  
  foldchanges = Kegg_genes$log2FoldChange
  names(foldchanges) = Kegg_genes$Gene_ID
  foldchanges <- sort(foldchanges, decreasing = T)
  head(foldchanges)
  
  kk_all <- enrichKEGG(names(foldchanges), organism="mmu", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg")
  head(summary(kk_all))
  kk_all <- setReadable(kk_all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  
  kk_down <- enrichKEGG(names(foldchanges[foldchanges<0]), organism="mmu", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg") 
  head(summary(kk_down))
  kk2_down <- setReadable(kk_down, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kk_res_down <- as.data.frame(kk2_down)
  kk_res_down$direction <- "down"
  
  kk_up <- enrichKEGG(names(foldchanges[foldchanges>0]), organism="mmu", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg", ) 
  head(summary(kk_up))
  kk2_up <- setReadable(kk_up, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kk_res_up <- as.data.frame(kk2_up)
  kk_res_up$direction <- "up"
  
  kk_results <- rbind(kk_res_up, kk_res_down)
  return(as.data.frame(kk_all))
}

kk_all_HG <- CalculateKeggEnrichment(resSig_highGlucose, LOG2FOLDCHANGE=0.6)


selected_kegg <- c("mmu04550","mmu03320","mmu04350","mmu04151","mmu04020","mmu04974","mmu04514","mmu04550" ,"mmu04015", "mmu04010",  "mmu00350", "mmu04144",  "mmu04310", "mmu04390", "mmu00260","mmu04060","mmu04979","mmu01521","mmu00480","mmu04923") 



PlotKeggPathways <- function(kk_results, RESULTS_TABLE, selected_kegg, NO_OF_PATHWAYS_TO_PLOT= 13, col_UP="red",col_DOWN="blue"  ){
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
    df_tmp <- RESULTS_TABLE[RESULTS_TABLE$external_gene_name %in% unlist(strsplit(enrichKegg_selected[i, "geneID"] , split="/")),]
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
  enrichKegg_molten$Description <- gsub( "signaling pathway","sig. path.", enrichKegg_molten$Description)
  enrichKegg_molten$Description <- gsub( "Signaling pathways","Sig. path.", enrichKegg_molten$Description)
  
  p_kegg_mlt <- ggplot(enrichKegg_molten, aes(x=reorder(Description, -qvalue), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue))) +
    coord_flip() +  xlab(" ") +
    scale_fill_manual( values = c(col_DOWN, col_UP)) +  ylab("Gene count") +
    ylim(-max(abs(enrichKegg_molten$value)), max(abs(enrichKegg_molten$value))) +
    theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_alpha_continuous( range = c(0.5, 1))
  
  return(p_kegg_mlt)
  
}


pdf(paste( Project,  "___KEGG_barplots_HighGlucose", "markers_l2fc1", ".pdf", sep="_"), width=9,height=4)
par(bg=NA)
PlotKeggPathways(kk_all_HG, resSig_highGlucose, selected_kegg, NO_OF_PATHWAYS_TO_PLOT = 10, col_DOWN="firebrick",col_UP="darkolivegreen4")
dev.off()






















