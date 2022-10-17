#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()


##############################################################################################################
#                                                                                                            #
#   Project: MBU_spb54_005_MouseEmbryo                                                                       #
#   Malwina Prater (mn367@cam.ac.uk), 2022                                                                   #
#   MRC MBU, University of Cambridge                                                                         #
#   Script: scRNA-seq mouse dataset - Gene networks                                                          # 
#                                                                                                            #
##############################################################################################################



message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
  library("STRINGdb")
  library("igraph")
  library("biomaRt")
  library("rWikiPathways")
})



Project        <- "MBU_spb54_005__Enrichment_Network_Analysis"
baseDir        <- "/Users/xxx/Documents/xxx/xxx/xxx" # replace with your path
SCENIC_Dir     <- paste0(baseDir,"/SCENIC" ) # replace with your path
setwd(baseDir)


message("+-------------------------------------------------------------------------------")
message("+                             load in data                                      ")
message("+-------------------------------------------------------------------------------")

Seurat_obj <- readRDS(paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC-_harmony_matrix.su_SCT_batch_regressed_.Rds"))
scenicOptions <- readRDS(paste0(baseDir,"/SCENIC_batch_regressed/int/scenicOptions_SCT_batch_regrMBU_spb54_005_Flo_SCENIC.Rds"))
exprMat  <- readRDS( paste0(baseDir, "/Input/MBU_spb54_005_Flo_SCENIC__SCT_counts_batch_regressed.Rds"))

setwd(SCENIC_Dir)

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
cellInfo <- read.csv("int/MBU_spb54_005_Flo_SCENIC_cell_info.csv")
rownames(cellInfo) <- cellInfo$X
cellInfo$CellType_Genotype <- paste(cellInfo$celltype, cellInfo$mouse, sep = "_")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

setwd(baseDir)



message("+-------------------------------------------------------------------------------")
message("+                       add extra labeling                                      ")
message("+-------------------------------------------------------------------------------")

celltype_cols <- c( "amnion"="plum1", "mesodermProgenitors"="darkolivegreen1", "neuralTube"="plum4","mixedMesoderm"="gold2", "neuralCrest"="purple4", "midHindbrain"="orchid3","notochord"="grey28", "placodes"="red", "presomiticMesoderm"="darkolivegreen3", "forebrain"="mediumpurple2",  "pharyngealMesoderm"="royalblue","foregut"="violetred3", "extraembryonicMesoderm"="steelblue3", "cardiac"="firebrick" , "somiticMesoderm"= "darkgreen",  "endothelial"="orange1", "midHindgut"="violetred2", "blood"="black")


Seurat_obj@meta.data$celltype_genotype <- paste(Seurat_obj@meta.data$celltype, Seurat_obj@meta.data$mouse, sep = "_")
Seurat_obj@meta.data$celltype_genotype <- factor(Seurat_obj@meta.data$celltype_genotype, levels = c(
  "blood_A5019G",
  "amnion_WT", "amnion_C5024T" ,"amnion_A5019G" ,
  "mesodermProgenitors_WT" ,"mesodermProgenitors_C5024T","mesodermProgenitors_A5019G",
  "neuralTube_WT","neuralTube_C5024T" , "neuralTube_A5019G",
  "mixedMesoderm_WT",  "mixedMesoderm_C5024T" ,"mixedMesoderm_A5019G" ,
  "neuralCrest_WT", "neuralCrest_C5024T"  ,"neuralCrest_A5019G", 
  "midHindbrain_WT", "midHindbrain_C5024T", "midHindbrain_A5019G", 
  "notochord_WT","notochord_C5024T", "notochord_A5019G",          
  "placodes_WT", "placodes_C5024T"   ,   "placodes_A5019G",     
  "presomiticMesoderm_WT",  "presomiticMesoderm_C5024T", "presomiticMesoderm_A5019G",  
  "forebrain_WT","forebrain_C5024T" ,"forebrain_A5019G" ,         
  "pharyngealMesoderm_WT", "pharyngealMesoderm_C5024T", "pharyngealMesoderm_A5019G", 
  "foregut_WT", "foregut_C5024T", "foregut_A5019G" ,  
  "extraembryonicMesoderm_WT", "extraembryonicMesoderm_C5024T", "extraembryonicMesoderm_A5019G",
  "cardiac_WT", "cardiac_C5024T", "cardiac_A5019G" ,       
  "somiticMesoderm_WT" , "somiticMesoderm_C5024T" ,"somiticMesoderm_A5019G",    
  "endothelial_WT", "endothelial_C5024T", "endothelial_A5019G",    
  "midHindgut_WT","midHindgut_C5024T", "midHindgut_A5019G" ))
Idents(Seurat_obj) <- Seurat_obj@meta.data$celltype_genotype




message("+-------------------------------------------------------------------------------")
message("+              ensembl annotation and pathways data                             ")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl) # , "hgnc_symbol", "hgnc_id"
gene.data <- getBM(attributes=c('external_gene_name', 'go_id',"name_1006"), mart = ensembl) #, filters = 'go_id', values = 'GO:0007507')

wp.mm.gmt <- rWikiPathways::downloadPathwayArchive(organism="Mus musculus", format = "gmt")
wp2gene <- readPathwayGMT(wp.mm.gmt)
wp2gene$external_gene_name <- ensEMBL2id[match(wp2gene$gene, ensEMBL2id$entrezgene_id),]$external_gene_name



message("+-------------------------------------------------------------------------------")
message("+                     load in data for the plot                                 ")
message("+-------------------------------------------------------------------------------")

Markers_neuralCrest_A5019G <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/neuralCrest_A5019G.csv"), row.names = 1)
Markers_neuralCrest_C5024T <- read.csv(paste0(baseDir,"/Input/celltypeDEGs_Seurat/neuralCrest_C5024T.csv"), row.names = 1)


message("--------------------------------------------------------------------------------")
message("+                              StringDB PPI                                     ")
message("+-------------------------------------------------------------------------------")

# https://www.bioconductor.org/packages/devel/bioc/vignettes/netSmooth/inst/doc/buildingPPIsFromStringDB.html

# 1. getSTRINGdb for human
string_db <- STRINGdb$new(species=10090) # human::9606)
mouse_graph <- string_db$get_graph()

# 2. get edges with high confidence score
edge.scores <- E(mouse_graph)$combined_score
ninetyth.percentile <- quantile(edge.scores, 0.9)
thresh <- data.frame(name='90th percentile',
                     val=ninetyth.percentile)
mouse_graph <- subgraph.edges(mouse_graph,  E(mouse_graph)[combined_score > ninetyth.percentile])

# 3. create adjacency matrix
adj_matrix <- as_adjacency_matrix(mouse_graph)
adj_matrix[1:10,1:10]


# 4. map gene ids to protein ids

### extract protein ids from the human network
protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'),
                      function(x) x[2])

### get protein to gene id mappings
mart_results <- getBM(attributes = c("external_gene_name", "ensembl_peptide_id"), filters = "ensembl_peptide_id", values = protein_ids, mart = ensembl)

### replace protein ids with gene ids
ix <- match(protein_ids, mart_results$ensembl_peptide_id)
ix <- ix[!is.na(ix)]

newnames <- protein_ids
newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
  mart_results[ix, 'external_gene_name']
rownames(adj_matrix) <- newnames
colnames(adj_matrix) <- newnames

ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
nullrows <- Matrix::rowSums(ppi)==0
ppi <- ppi[!nullrows,!nullrows] ## ppi is the network with gene ids

ppi[1:10,1:10]

ppi_2 <- as.data.frame(ppi)
nullrows <- Matrix::rowSums(ppi_2)==0
ppi_2 <- ppi_2[!nullrows,!nullrows]
ppi_2$gene <- rownames(ppi_2)
ppi_molten <- reshape2::melt(ppi_2, id.vars= "gene" )
ppi_molten <- ppi_molten[ppi_molten$value != 0,]
ppi_molten <- ppi_molten[,c(1,2)]
colnames(ppi_molten) <-c("from","to")
ppi_molten <- ppi_molten[!duplicated(ppi_molten$from) | !duplicated(ppi_molten$to),]
ppi_molten$source <-"STRING"





message("--------------------------------------------------------------------------------")
message("+                     iGRAPH PLOT FOR SOX9 SOX10                                ")
message("+-------------------------------------------------------------------------------")

PATHWAY_GENES=wp2gene[wp2gene$wpid == "WP2074",]$external_gene_name
PATHWAY_NAME= "WP2074"

PLOT_NAME = "Markers_neuralCrest_C5024T"
REGULONS = c("Sox10", "Sox9")

MARKER_GENES_C5024T = rownames(Markers_neuralCrest_C5024T)
MARKER_GENES_A5019G = rownames(Markers_neuralCrest_A5019G)


REGULON_NAMES <- stringr::str_c(REGULONS, collapse = "_")
regulonTargetsInfo_subset <- regulonTargetsInfo[regulonTargetsInfo$TF %in% REGULONS,]

edges <- regulonTargetsInfo_subset[,c(1,2)]
colnames(edges) <- c("from","to")

motifs <- unique(as.character(edges$from))
genes <- unique(as.character(edges$to))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes),  
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("pink", length(motifs)), rep("honeydew1", length(genes))))
nodes <- nodes[!duplicated(nodes$id),]
nodes$color <- ifelse(nodes$id %in% MARKER_GENES_C5024T , "dodgerblue" , nodes$color) 
nodes$color <- ifelse(nodes$id %in% MARKER_GENES_A5019G , "green3", nodes$color) 
nodes$color <- ifelse(nodes$id %in% MARKER_GENES_A5019G & nodes$id %in% MARKER_GENES_C5024T, "turquoise2"  , nodes$color)
nodes$cex_color <- ifelse(nodes$id %in% PATHWAY_GENES, "red", "black")

ppi_subset <- as.data.frame(ppi[colnames(ppi) %in% c(nodes$id), rownames(ppi) %in% nodes$id])
nullrows <- Matrix::rowSums(ppi_subset)==0
ppi_subset <- ppi_subset[!nullrows,!nullrows] 
ppi_subset$gene <- rownames(ppi_subset)
ppi_subset_molten <- reshape2::melt(ppi_subset, id.vars= "gene" )
ppi_subset_molten <- ppi_subset_molten[ppi_subset_molten$value != 0,]
ppi_subset_molten <- ppi_subset_molten[,c(1,2)]
colnames(ppi_subset_molten) <-c("from","to")
ppi_subset_molten <- ppi_subset_molten[!duplicated(ppi_subset_molten$from) | !duplicated(ppi_subset_molten$to),]

edges2 <- unique(rbind(edges, ppi_subset_molten ))

net <- graph_from_data_frame(d=edges, vertices=nodes$id, directed=T) 
class(net)
E(net)$arrow.mode <- 0
nodes$shape2 <- ifelse(nodes$shape == "elypse", "circle", "square")
V(net)$shape <- nodes$shape2
nodes$size <- ifelse(nodes$id %in% unique(c(MARKER_GENES_C5024T, MARKER_GENES_A5019G)), 1.3, 0.8)
vertex_size <- ifelse(nodes$id %in% c("Sox9", "Sox10"), 15, 10)
#plot(net, edge.arrow.size=.4,  vertex.size=7, vertex.color=nodes$color, vertex.label.cex = nodes$size, edge.curved=.1, main=REGULON_NAMES)


# remove all non-significant genes
edges3 <- edges2
edges3$blue <- ifelse(edges3$to %in% nodes[nodes$color == "honeydew1",]$id | edges3$from %in% nodes[nodes$color == "honeydew1",]$id, "blue", "other")
edges3 <- edges3[edges3$blue != "blue", c(1,2)]
nodes2 <- nodes[nodes$id %in% edges3$from | nodes$id %in% edges3$to,]



net <- graph_from_data_frame(d=edges3, vertices=nodes2$id, directed=T) 
class(net)
net <- simplify(net, remove.multiple = F, remove.loops = T) 
E(net)$arrow.mode <- 0
nodes2$shape2 <- ifelse(nodes2$shape == "elypse", "circle", "square")
V(net)$shape <- nodes2$shape2
nodes2$size <- ifelse(nodes2$id %in% PATHWAY_GENES, 1.2, 1)
print({plot(net, edge.arrow.size=.4,  vertex.size=7, vertex.color=nodes2$color, vertex.label.color= nodes2$cex_color, vertex.label.cex = nodes2$size, edge.curved=.1, main=REGULON_NAMES) }) 
vertex_size <- ifelse(nodes2$id %in% c("Sox9", "Sox10"), 17, 12)

pdf(paste0( Project, "__", "regulons_", REGULON_NAMES, "__" , PLOT_NAME, "_pathway_in_purple_", PATHWAY_NAME, "_trimmed2_SeuratMarkers.pdf"), onefile=FALSE, width=12, height=12) 
par(bg=NA)
print({plot(net, edge.arrow.size=0.2,  vertex.size=vertex_size, vertex.color=nodes2$color, vertex.label.cex = nodes2$size, vertex.label.color= nodes2$cex_color, edge.curved=0, main="Regulons Sox10 and Sox9") })
dev.off()


































