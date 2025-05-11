# Load necessary packages:
library(Seurat)
library(ggplot2)
library(tidyverse)
library(harmony)
library(celldex)
library(SingleR)
library(pheatmap)
library(gridExtra)

# Read the data:
dirs <- list.dirs( path = "GSE154778_met/metastatic_tissue/data/" , recursive = F , full.names = F)
dirs

# Create seurat objects:
for (x in dirs){
  name <- x
  
  cts <- ReadMtx(mtx = paste0('GSE154778_met/metastatic_tissue/data/', x, '/matrix.mtx.gz'),
                 features = paste0('GSE154778_met/metastatic_tissue/data/' , x, '/features.tsv.gz'),
                 cells = paste0('GSE154778_met/metastatic_tissue/data/', x, '/barcodes.tsv.gz'))
  
  assign(name, CreateSeuratObject(counts = cts))
  
}

# Merge all the seurat objects:
merged_seurat <- merge(MET01, y = c( MET02, MET04, MET05, MET06),
                       add.cell.ids = ls()[3:7],
                       project = 'Cancer')

merged_seurat
view(merged_seurat@meta.data)

# Create separate columns according to data:
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample' , into = c('patient', 'barcode'), sep = '_')
view(merged_seurat@meta.data)

# sanity check:
unique(merged_seurat@meta.data$patient)


# DATA ANALYSIS:

# 1. Quality Control:
# calculate mitochondrial content:
merged_seurat$mitopercent <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')

# Filtering:
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 500 & 
                                   nFeature_RNA > 1000 & nFeature_RNA < 20000 &
                                   mitopercent < 5)

merged_seurat_filtered


# 2. Normalization:
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)


# 3. Feature Selection:
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered, nfeatures = 1000)


# 4. Scaling:
all.genes <- rownames(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered,features = all.genes)


# 5. Dimensionality Reduction:
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered , features = VariableFeatures(object = merged_seurat_filtered))

# Determine dimensionality of dataset:
e <- ElbowPlot(object = merged_seurat_filtered)
# To save elbow plot:
e <- e + theme(plot.background = element_rect(fill = "white", color = NA),
               plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
ggsave("/home/ecgi4/pravallika/PDAC/metastatic_tissue_images/before_int/elbow_plot.png", e, width = 8, height = 6, units = "in")


# Jackstraw analysis:
merged_seurat_filtered <- JackStraw(merged_seurat_filtered, num.replicate = 100)
merged_seurat_filtered <- ScoreJackStraw(merged_seurat_filtered)
JackStrawPlot(merged_seurat_filtered)


# 6. Clustering:
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:10)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.5)


# 7. Visualization:

# UMAP:
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:10)
DimPlot(merged_seurat_filtered, reduction = "umap" , label = TRUE)
u1 <- UMAPPlot(merged_seurat_filtered, group.by = "RNA_snn_res.0.5")

# tsne:
merged_seurat_filtered <- RunTSNE(object = merged_seurat_filtered, dims = 1:10)
DimPlot(merged_seurat_filtered, reduction = "tsne" , label = TRUE)
s1 <- TSNEPlot(merged_seurat_filtered, group.by = "RNA_snn_res.0.5") 

# plot grouped by patient:
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'patient')
p2 <- DimPlot(merged_seurat_filtered ,reduction = 'tsne', group.by = 'patient')


# 8. PERFORM INTEGRATION (Batch effect correction):(optional)

obj.list <- SplitObject(merged_seurat_filtered, split.by = 'patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], nfeatures = 1000)
}

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors)

# Save the integrated seurat object:
saveRDS(seurat.integrated, file = "seurat_object_met_1.rds")

# Read integrated seurat object:
seurat.integrated <- readRDS("seurat_object_met_1.rds")

# 8.1 Scaling:
all.genes <- rownames(seurat.integrated)
seurat.integrated <- ScaleData(object = seurat.integrated, features = all.genes) 

# 8.2 Dimensionality Reduction:
seurat.integrated <- RunPCA(object = seurat.integrated, features = VariableFeatures(object = seurat.integrated)) 

# Determine dimensionality of dataset:
ElbowPlot(object = seurat.integrated)

# 8.3 Clustering:
seurat.integrated <- FindNeighbors(object = seurat.integrated, dims = 1:10)
seurat.integrated <- FindClusters(object = seurat.integrated, resolution = 0.5) 

# 8.4 Visualization:
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:10)
u2 <- DimPlot(seurat.integrated, reduction = "umap" , label = TRUE)
seurat.integrated <- RunTSNE(object = seurat.integrated, dims = 1:10)
s2 <- DimPlot(seurat.integrated, reduction = "tsne" , label = TRUE)

# plot grouped by patient:
p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'patient')
p4 <- DimPlot(seurat.integrated ,reduction = 'tsne', group.by = 'patient')

# To save plots:
ggsave("/home/ecgi4/pravallika/PDAC/metastatic_tissue_images/patient_ts.png", p4, width = 8, height = 6, units = "in")



# 9. DOWNSTREAM ANALYSIS:

# 9.1 FINDING DIFFERENTIALLY EXPRESSED FEATURES:
seurat.integrated.markers <- FindAllMarkers(seurat.integrated )
View(seurat.integrated.markers)


# Sort all the markers with p-value < 0.05 and log2fc < 0.25;
markers <- FindAllMarkers(seurat.integrated, only.pos = TRUE, avg_log2FC = 0.25)
filtered_markers <- markers[markers$p_val_adj < 0.05, ]


# Perform differential expression analysis for a specific cluster: UP-REGULATED GENES
cluster_de_analysis_upregulated <- function(seurat_obj, cluster_id, logfc_threshold = 0.25, p_val_threshold = 0.05) {
  markers <- FindMarkers(seurat_obj, ident.1 = cluster_id, only.pos = TRUE, logfc.threshold = logfc_threshold)
  markers <- markers[markers$p_val_adj < p_val_threshold, ]
  filtered_markers <- markers[order(markers$p_val_adj), ]
  return(filtered_markers)
}
cluster_0_markers_up <- cluster_de_analysis_upregulated(seurat.integrated, cluster_id = 0)
cluster_1_markers_up <- cluster_de_analysis_upregulated(seurat.integrated, cluster_id = 1)
cluster_2_markers_up <- cluster_de_analysis_upregulated(seurat.integrated, cluster_id = 2)
cluster_3_markers_up <- cluster_de_analysis_upregulated(seurat.integrated, cluster_id = 3)
cluster_4_markers_up <- cluster_de_analysis_upregulated(seurat.integrated, cluster_id = 4)
cluster_5_markers_up <- cluster_de_analysis_upregulated(seurat.integrated, cluster_id = 5)
cluster_6_markers_up <- cluster_de_analysis_upregulated(seurat.integrated, cluster_id = 6)


# Perform differential expression analysis for a specific cluster: DOWN-REGULATED GENES
cluster_de_analysis_downregulated <- function(seurat_obj, cluster_id, p_val_threshold = 0.05) {
  markers <- FindMarkers(seurat_obj, ident.1 = cluster_id)
  filtered_markers <- markers[markers$p_val_adj < p_val_threshold & markers$avg_log2FC < 0,]
  filtered_markers <- markers[order(markers$p_val_adj), ]
  filtered_markers <- markers[order( markers$avg_log2FC) ,]
  return(filtered_markers)
}
cluster_0_markers_dr <- cluster_de_analysis_downregulated(seurat.integrated, cluster_id = 0)
cluster_1_markers_dr <- cluster_de_analysis_downregulated(seurat.integrated, cluster_id = 1)
cluster_2_markers_dr <- cluster_de_analysis_downregulated(seurat.integrated, cluster_id = 2)
cluster_3_markers_dr <- cluster_de_analysis_downregulated(seurat.integrated, cluster_id = 3)
cluster_4_markers_dr <- cluster_de_analysis_downregulated(seurat.integrated, cluster_id = 4)
cluster_5_markers_dr <- cluster_de_analysis_downregulated(seurat.integrated, cluster_id = 5)
cluster_6_markers_dr <- cluster_de_analysis_downregulated(seurat.integrated, cluster_id = 6)


# Extract specific marker genes in all clusters:
extract_specific_markers <- function(filtered_markers, genes_of_interest) {
  specific_markers <- filtered_markers[filtered_markers$gene %in% genes_of_interest, c("cluster", "gene", "avg_log2FC" , "p_val_adj")]
  return(specific_markers)
}
# ductal_cells <- c("MMP7", "TSPAN8", "SOX9", "LCN2", "KRT19", "TFF2")
# usage:- ductal_cells_info <- extract_specific_markers(filtered_markers, ductal_cells)


# 9.2 VISUALIZING MARKER GENES:

# Store the marker genes to be visualized in feature set:
ductal_epithelial_cells <- c('EPCAM', 'KRT19', 'MMP7', 'TSPAN8', 'SOX9', 'LCN2')
type_1_ductal <- c('CAPN12' , 'CA12')
type_2_ductal <- c('HMMR' , 'TPX2')
type_3_ductal  <- c('SCGB2A1' , 'IGFBP3')
type_4_ductal  <- c('CTGF' , 'HSPA6')
type_5_ductal  <- c('REG4' , 'CEACAM5')
malignant_ductal_cells <- c('CEACAM1' , 'KRT19', 'CEACAM5', 'CEACAM6')
macrophages <- c('HLA-DRB5', 'NUPR1', 'CD68' , 'G-CSF','BASP1' ,'TNF')
lymphocytes <- c('CD3D', 'CD3G', 'CD3E', 'CD27')
acinar_cells <- c('PRSS1', 'CELA3A')
cell_prolif_markers <- c('MK167','TOP2A')
cell_cycle_markers <- c('CCNB1','CCNB2')
drug_target_genes <- c('CDK1', 'PLK1' , 'AURKA', 'TTK', 'WEE1')
abnormal_cells <- c('MUC1', 'FXYD3','MUC5AC')
normal_ductal_cells <- c('AMBP', 'CFTR', 'MMP7')
tumor_progression_genes <- c('HMGA1', 'FOS', 'KLF5')
features <- c('EPCAM', 'TSPAN8', 'HLA-DRB5','CD68','COL1A1','DCN','CD3D','CD3E', 'KDR', 'VWF')

# Feature_plot:
FeaturePlot(seurat.integrated, features = ductal_epithelial_cells, reduction = "tsne" , label = TRUE, min.cutoff = c(1.5,1.5,1.5,1.5,1.5,1.5))
FeaturePlot(seurat.integrated, features = malignant_ductal_cells, reduction = "tsne" , label = TRUE, min.cutoff = c(0,1,0,0) )
FeaturePlot(seurat.integrated, features = type_1_ductal, reduction = "tsne" , label = TRUE )
FeaturePlot(seurat.integrated, features = type_2_ductal, reduction = "tsne" , label = TRUE, min.cutoff = c(0.5,0.5) )
FeaturePlot(seurat.integrated, features = type_3_ductal, reduction = "tsne" , label = TRUE )
FeaturePlot(seurat.integrated, features = type_4_ductal, reduction = "tsne" , label = TRUE, min.cutoff = c(0,0) )
FeaturePlot(seurat.integrated, features = type_5_ductal, reduction = "tsne" , label = TRUE )
FeaturePlot(seurat.integrated, features = macrophages, reduction = "tsne" , label = TRUE )
FeaturePlot(seurat.integrated, features = lymphocytes, reduction = "tsne" , label = TRUE )
FeaturePlot(seurat.integrated, features = abnormal_cells, reduction = "tsne" , label = TRUE ,min.cutoff = c(0,1,0))
FeaturePlot(seurat.integrated, features = c('WEE1'), reduction = "tsne",label = TRUE )

# Violin_plot:
VlnPlot(seurat.integrated, features= c('TOP2A','CCNB1'), pt.size = 0) 

# Dotplot
e<-DotPlot(object = seurat.integrated, features = c('TOP2A','CCNB1','CCNB2'), cols = c("lightgrey", "red"))
e <- e + theme(plot.background = element_rect(fill = "white", color = NA),
               plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))

# Heatmap
DoHeatmap(object = seurat.integrated, features = drug_target_genes)

# To save the plots:
ggsave("/home/ecgi4/pravallika/PDAC/metastatic_tissue_images/cell_cycle_violin.png", v, width = 8, height = 6, units = "in")


# 9.3 CLUSTER ANNOTATION:

# AUTOMATIC CELL TYPE ANNOTATION:
# Load the reference dataset from celldex:
ref <- celldex::HumanPrimaryCellAtlasData()
view(colData(ref))

# Get the raw count matrix:
met_counts <- GetAssayData(seurat.integrated, layer = 'data')

# Perform SingleR:
pred <- SingleR(test = met_counts,
                ref = ref,
                labels = ref$label.main )

# Add singleR labels to seurat object:
seurat.integrated$singleR.labels <- pred$labels[match(rownames(seurat.integrated@meta.data), rownames(pred))]

# Diagnose the results of singleR:
DimPlot(seurat.integrated, reduction = 'tsne', group.by = 'singleR.labels')
plotScoreHeatmap(pred)
plotDeltaDistribution(pred)
tab <- table(Assigned=pred$labels, Clusters = seurat.integrated$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white', 'blue'))(10))

# MANUAL ANNOTATION:
new.cluster.ids <- c("Type 1 ductal cells",
                     "Type 2 ductal cells",
                     "Type 3 ductal cells",
                     "Type 2 ductal cells",
                     "Type 2 ductal cells",
                     "Macrophages",
                     "T-cells")
names(new.cluster.ids) <- levels(seurat.integrated)
seurat.integrated <- RenameIdents(seurat.integrated, new.cluster.ids)

d <- DimPlot(seurat.integrated, reduction = "tsne", pt.size = 0.2) 

DotPlot(object = seurat.integrated, features = features, cols = c("lightgrey", "red"))

DoHeatmap(object = seurat.integrated, features = features)

labels_of_interest <- c('Epithelial_cells', 'T_cells', 'Tissue_stem_cells', 'Fibroblasts', 'Endothelial_cells', 'Macrophage')
filtered_pred <- subset(pred, labels %in% labels_of_interest)
pred_cell_names <- rownames(filtered_pred)
subset_merged_seurat <- seurat.integrated[, pred_cell_names]
tab <- table(Assigned = filtered_pred$labels, Clusters = subset_merged_seurat$seurat_clusters)
p <- pheatmap(log10(tab + 10), color = colorRampPalette(c('white', 'blue'))(10))

ggsave("/home/ecgi4/pravallika/PDAC/metastatic_tissue_images/ductal_sub_clusters.png", d, width = 8, height = 6, units = "in")


# 9.4 FUNCTIONAL ENRICHMENT ANALYSIS:

# clusterprofiler
# David
# Metascape
# GSVA

# Get top 100 marker genes for each cluster sorted according to adjusted p-value and store:
write_top_genes <- function(markers, output_file, top_genes = 100) {
  top_markers <- head(markers, top_genes)
  write.csv(top_markers, file = output_file, row.names = TRUE)
}

# UP-REGULATED genes:
write_top_genes(cluster_0_markers_up, output_file = "GSE154778_met/met_genes/upregulated/cluster_0_markers.csv")
write_top_genes(cluster_1_markers_up, output_file = "GSE154778_met/met_genes/upregulated/cluster_1_markers.csv")
write_top_genes(cluster_2_markers_up, output_file = "GSE154778_met/met_genes/upregulated/cluster_2_markers.csv")
write_top_genes(cluster_3_markers_up, output_file = "GSE154778_met/met_genes/upregulated/cluster_3_markers.csv")
write_top_genes(cluster_4_markers_up, output_file = "GSE154778_met/met_genes/upregulated/cluster_4_markers.csv")
write_top_genes(cluster_5_markers_up, output_file = "GSE154778_met/met_genes/upregulated/cluster_5_markers.csv")
write_top_genes(cluster_6_markers_up, output_file = "GSE154778_met/met_genes/upregulated/cluster_6_markers.csv")

# DOWN_REGULATED genes:
write_top_genes(cluster_0_markers_dr, output_file = "GSE154778_met/met_genes/downregulated/cluster_0_markers.csv")
write_top_genes(cluster_1_markers_dr, output_file = "GSE154778_met/met_genes/downregulated/cluster_1_markers.csv")
write_top_genes(cluster_2_markers_dr, output_file = "GSE154778_met/met_genes/downregulated/cluster_2_markers.csv")
write_top_genes(cluster_3_markers_dr, output_file = "GSE154778_met/met_genes/downregulated/cluster_3_markers.csv")
write_top_genes(cluster_4_markers_dr, output_file = "GSE154778_met/met_genes/downregulated/cluster_4_markers.csv")
write_top_genes(cluster_5_markers_dr, output_file = "GSE154778_met/met_genes/downregulated/cluster_5_markers.csv")
write_top_genes(cluster_6_markers_dr, output_file = "GSE154778_met/met_genes/downregulated/cluster_6_markers.csv")


# The extracted top marker genes of each cluster are used as input for David analysis;




