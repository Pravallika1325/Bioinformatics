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
dirs <- list.dirs( path = "GSE154778_met/primary_tissue/" , recursive = F , full.names = F)
dirs

# Create seurat objects:
for (x in dirs){
  name <- x
  
  cts <- ReadMtx(mtx = paste0('GSE154778_met/primary_tissue/', x, '/matrix.mtx.gz'),
                 features = paste0('GSE154778_met/primary_tissue/' , x, '/features.tsv.gz'),
                 cells = paste0('GSE154778_met/primary_tissue/', x, '/barcodes.tsv.gz'))
  
  assign(name, CreateSeuratObject(counts = cts))
}

# Merge all the seurat objects:
merged_seurat <- merge(PO1, y = c( PO10, PO2, PO3, PO4, PO5, PO6, PO7, PO8, PO9),
                       add.cell.ids = ls()[43:52],
                       project = 'Cancer')

merged_seurat
view(merged_seurat@meta.data)

# create separate columns according to data:
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

# 5. Dimensionality Reduction
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered , features = VariableFeatures(object = merged_seurat_filtered))

# Determine dimensionality of dataset:
ElbowPlot(object = merged_seurat_filtered)

# Jackstraw analysis:
merged_seurat_filtered <- JackStraw(merged_seurat_filtered, num.replicate = 100)
merged_seurat_filtered <- ScoreJackStraw(merged_seurat_filtered)
JackStrawPlot(merged_seurat_filtered)

# 6. Clustering:
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:10)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.5)

# 7. Visualization:

# UMAP
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:10)
DimPlot(merged_seurat_filtered, reduction = "umap" , label = TRUE)
u1 <- UMAPPlot(merged_seurat_filtered, group.by = "RNA_snn_res.0.5")

# tsne
merged_seurat_filtered <- RunTSNE(object = merged_seurat_filtered, dims = 1:10)
DimPlot(merged_seurat_filtered, reduction = "tsne" , label = TRUE)
s1 <- TSNEPlot(merged_seurat_filtered, group.by = "RNA_snn_res.0.5")

# plot grouped by patient:
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'patient')
p2 <- DimPlot(merged_seurat_filtered ,reduction = 'tsne', group.by = 'patient')

ggsave("/home/ecgi4/pravallika/PDAC/primary_tissue_images/patient_ts.png", p2, width = 8, height = 6, units = "in")


# 8. PERFORM INTEGRATION (Batch effect correction): - seurat object should have PCA data;
# HARMONY:
merged_seurat_filtered_harmony <- RunHarmony(merged_seurat_filtered, "patient")
merged_seurat_filtered_harmony@reductions
h1 <- DimPlot(object = merged_seurat_filtered_harmony, reduction = "harmony", pt.size = .1, group.by = "patient")

# Save the integrated seurat object:
saveRDS(merged_seurat_filtered_harmony, file = "seurat_object_prim_1_all.rds")


# Read integrated seurat object:
merged_seurat_filtered_harmony <- readRDS("seurat_object_prim_1_all.rds")

# 8.1 Clustering:
merged_seurat_filtered_harmony <- FindNeighbors(object = merged_seurat_filtered_harmony, dims = 1:10, reduction = "harmony")
merged_seurat_filtered_harmony <- FindClusters(object = merged_seurat_filtered_harmony, resolution = 0.5)

# 8.2. Visualization
merged_seurat_filtered_harmony <- RunUMAP(object = merged_seurat_filtered_harmony, dims = 1:10, reduction = "harmony")
u2 <- DimPlot(merged_seurat_filtered_harmony, reduction = "umap", label = TRUE, pt.size = .3)
merged_seurat_filtered_harmony <- RunTSNE(object = merged_seurat_filtered_harmony, dims = 1:10, reduction = "harmony")
s2 <- DimPlot(merged_seurat_filtered_harmony, reduction = "tsne", label = TRUE, pt.size = .3)

grid.arrange(u2, s2, ncol(2), nrow(1))


# 9. DOWNSTREAM ANALYSIS:

# 9.1 FINDING DIFFERENTIALLY EXPRESSED FEATURES:
# Join all the layers:
merged_seurat_filtered_harmony <- JoinLayers(merged_seurat_filtered_harmony, 
                                             features = VariableFeatures(merged_seurat_filtered_harmony))
merged_seurat_filtered_harmony.markers <- FindAllMarkers(merged_seurat_filtered_harmony )


# Sort all the markers with p-value < 0.05 and log2fc < 0.25;
markers <- FindAllMarkers(merged_seurat_filtered_harmony, only.pos = TRUE, avg_log2FC = 0.25)
filtered_markers <- markers[markers$p_val_adj < 0.05, ]


# Perform differential expression analysis for a specific cluster:
perform_cluster_de_analysis <- function(seurat_obj, cluster_id, logfc_threshold = 0.25, p_val_threshold = 0.05) {
  markers <- FindMarkers(seurat_obj, ident.1 = cluster_id, only.pos = TRUE, logfc.threshold = logfc_threshold)
  filtered_markers <- markers[markers$p_val_adj < p_val_threshold, ]
  return(filtered_markers)
}
cluster_0_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 0)
cluster_1_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 1)
cluster_2_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 2)
cluster_3_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 3)
cluster_4_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 4)
cluster_5_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 5)
cluster_6_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 6)
cluster_7_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 7)
cluster_8_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 8)
cluster_9_markers <- perform_cluster_de_analysis(merged_seurat_filtered_harmony, cluster_id = 9)


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
fibroblasts <- c('COL1A1', 'ACTA2','SPARC', 'DCN')
macrophages <- c('HLA-DRB5', 'NUPR1', 'CD68' , 'G-CSF')
lymphocytes <- c('CD3D', 'CD3G', 'CD3E')
endothelial_cells <- c('KDR' , 'VWF')
tissue_stem_cells <- c('CD34', 'CD44')
EMT_cells <- c('CDH2', 'SNAI2', 'ZEB1')

cell_cycle_markers <- c('MK167','TOP2A')
cell_proli_markers <- c('CCNB1','CCNB2')
drug_target_genes <- c('CDK1', 'PLK1' , 'AURKA', 'TTK', 'WEE1')
abnormal_cells <- c('MUC1', 'FXYD3','MUC5AC')
normal_ductal_cells <- c('AMBP', 'CFTR', 'MMP7')
tumor_progression_genes <- c('HMGA1', 'FOS', 'KLF5')

# Feature_plot:
FeaturePlot(merged_seurat_filtered_harmony, features = ductal_epithelial_cells, reduction = "tsne" , label = TRUE)
FeaturePlot(merged_seurat_filtered_harmony, features = fibroblasts, reduction = "tsne" , label = TRUE, min.cutoff = c(1,2,3))
FeaturePlot(merged_seurat_filtered_harmony, features = macrophages, reduction = "tsne" , label = TRUE )
FeaturePlot(merged_seurat_filtered_harmony, features = lymphocytes, reduction = "tsne" , label = TRUE )
FeaturePlot(merged_seurat_filtered_harmony, features = endothelial_cells, reduction = "tsne" , label = TRUE )
FeaturePlot(merged_seurat_filtered_harmony, features = abnormal_cells, reduction = "tsne" , label = TRUE, min.cutoff = c(0,2))

FeaturePlot(merged_seurat_filtered_harmony, features = c('WEE1'), reduction = "tsne" , label = TRUE)

ggsave("/home/ecgi4/pravallika/PDAC/primary_tissue_images/cell_cycle_violin.png", v, width = 8, height = 6, units = "in")



# Violin_plot:
features <- c('EPCAM', 'TSPAN8', 'HLA-DRB5','CD68','COL1A1','DCN','CD3D','CD3E', 'KDR', 'VWF')

# Violinplot
VlnPlot(merged_seurat_filtered_harmony, features= c('TOP2A','CCNB1'), pt.size = 0) 

# Dotplot
e<-DotPlot(object = merged_seurat_filtered_harmony, features = c('TOP2A','CCNB1','CCNB2'), cols = c("lightgrey", "red"))
e <- e + theme(plot.background = element_rect(fill = "white", color = NA),
               plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))

# Heatmap
DoHeatmap(object = merged_seurat_filtered_harmony, features = drug_target_genes)


# 9.3 CLUSTER ANNOTATION:

# AUTOMATIC CELL TYPE ANNOTATION:
# Join all the layers:
merged_seurat_filtered_harmony <- JoinLayers(merged_seurat_filtered_harmony, 
                                             features = VariableFeatures(merged_seurat_filtered_harmony))
# Load the reference dataset from celldex:
ref <- celldex::HumanPrimaryCellAtlasData()
view(colData(ref))

# Get the raw count matrix:
met_counts <- GetAssayData(merged_seurat_filtered_harmony, slot = 'data')
# Perform SingleR:
pred <- SingleR(test = met_counts,
                ref = ref,
                labels = ref$label.main )

# Add singleR labels to seurat object:
merged_seurat_filtered_harmony$singleR.labels <- pred$labels[match(rownames(merged_seurat_filtered_harmony@meta.data), rownames(pred))]

# Diagnose the results of singleR:
DimPlot(merged_seurat_filtered_harmony, reduction = 'tsne', group.by = 'singleR.labels')
plotScoreHeatmap(pred)
plotDeltaDistribution(pred)
tab <- table(Assigned=pred$labels, Clusters = merged_seurat_filtered_harmony$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white', 'blue'))(10))

# MANUAL ANNOTATION:
new.cluster.ids <- c("Type 1 ductal cells",
                     "Macrophages",
                     "Fibroblast",
                     "Fibroblast",
                     "Tissue stem cells",
                     "Type 2 ductal cells",
                     "T-cells",
                     "Tissue stem cells",
                     "Tissue stem cells",
                     "Endothelial cells")

names(new.cluster.ids) <- levels(merged_seurat_filtered_harmony)
merged_seurat_filtered_harmony <- RenameIdents(merged_seurat_filtered_harmony, new.cluster.ids)
d <- DimPlot(merged_seurat_filtered_harmony, reduction = "tsne", pt.size = 0.2)
DotPlot(object = merged_seurat_filtered_harmony, features = features, cols = c("lightgrey", "red"))
DoHeatmap(object = merged_seurat_filtered_harmony, features = features)

labels_of_interest <- c('Epithelial_cells', 'T_cells', 'Tissue_stem_cells', 'Fibroblasts', 'Endothelial_cells', 'Macrophage')
filtered_pred <- subset(pred, labels %in% labels_of_interest)
pred_cell_names <- rownames(filtered_pred)
subset_merged_seurat <- merged_seurat_filtered_harmony[, pred_cell_names]
tab <- table(Assigned = filtered_pred$labels, Clusters = subset_merged_seurat$seurat_clusters)
pheatmap(log10(tab + 10), color = colorRampPalette(c('white', 'blue'))(10))


# 9.4 FUNCTIONAL ENRICHMENT ANALYSIS:
# clusterprofiler
# David
# Metascape
# GSVA

# Get top 100 marker genes for each cluster sorted according to adjusted p-value and store:
write_top_genes <- function(markers, output_file, top_genes = 100) {
  markers_sorted <- markers[order(markers$p_val_adj), ]
  top_markers <- head(markers_sorted, top_genes)
  write.csv(top_markers, file = output_file, row.names = TRUE)
}

write_top_genes(cluster_0_markers, output_file = "GSE154778_met/prim_genes/cluster_0_markers.csv")
write_top_genes(cluster_1_markers, output_file = "GSE154778_met/prim_genes/cluster_1_markers.csv")
write_top_genes(cluster_2_markers, output_file = "GSE154778_met/prim_genes/cluster_2_markers.csv")
write_top_genes(cluster_3_markers, output_file = "GSE154778_met/prim_genes/cluster_3_markers.csv")
write_top_genes(cluster_4_markers, output_file = "GSE154778_met/prim_genes/cluster_4_markers.csv")
write_top_genes(cluster_5_markers, output_file = "GSE154778_met/prim_genes/cluster_5_markers.csv")
write_top_genes(cluster_6_markers, output_file = "GSE154778_met/prim_genes/cluster_6_markers.csv")
write_top_genes(cluster_7_markers, output_file = "GSE154778_met/prim_genes/cluster_7_markers.csv")
write_top_genes(cluster_8_markers, output_file = "GSE154778_met/prim_genes/cluster_8_markers.csv")
write_top_genes(cluster_9_markers, output_file = "GSE154778_met/prim_genes/cluster_9_markers.csv")

# The extracted top marker genes of each cluster are used as input for David analysis;



#-------------------
# Save the plots
# feature_plots:
ept <- FeaturePlot(merged_seurat_filtered_harmony, features = c('EPCAM', 'TSPAN8') , reduction = "tsne" , label = TRUE )
mac <- FeaturePlot(merged_seurat_filtered_harmony, features = c('HLA-DRB5','CD68') , reduction = "tsne" , label = TRUE )
fib <- FeaturePlot(merged_seurat_filtered_harmony, features = c('COL1A1','DCN') , reduction = "tsne" , label = TRUE )
lymp <- FeaturePlot(merged_seurat_filtered_harmony, features = c('CD3D','CD3E') , reduction = "tsne" , label = TRUE )
endo <- FeaturePlot(merged_seurat_filtered_harmony, features = c('KDR', 'VWF') , reduction = "tsne" , label = TRUE )
cancer <- FeaturePlot(merged_seurat_filtered_harmony, features = c('MK167','TOP2A','CCNB1','CCNB2') , reduction = "tsne" , label = TRUE )
abnormal <- FeaturePlot(merged_seurat_filtered_harmony, features = c('MUC5AC', 'FXYD3') , reduction = "tsne" , label = TRUE )

# dotplot
dot <- DotPlot(object = merged_seurat_filtered_harmony, features = features, cols = c("lightgrey", "red"))
dot <- dot + theme(plot.background = element_rect(fill = "white", color = NA),
                   plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))

ggsave("/home/ecgi4/pravallika/PDAC/primary_tissue_images/again_annotated_cluster.png", d, width = 8, height = 6, units = "in")


















