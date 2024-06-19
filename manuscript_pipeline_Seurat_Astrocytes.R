## conda activate Seurat
library(Seurat)
library(Matrix)
library(ggplot2)
library(harmony)
set.seed(2022)

project_dir = "/projects/compbio/users/wma36/collaborations/JieXu"
result_dir = file.path(project_dir, 'Manuscript_analysis_figures')
dir.create(result_dir, showWarnings = FALSE)
final_fig_dir = file.path(result_dir, 'Manuscript_final_figures')
dir.create(final_fig_dir, showWarnings = FALSE)
## --------------
# 2. Seurat pipeline - Astrocytes
# subtype resolution, control marker list and GO 
## --------------
## load annotated combined object
annotated_combined_obj = readRDS(file.path(result_dir, 'annotated_Seurat_obj.RDS'))
astrocytes_obj = subset(annotated_combined_obj, subset=annotated_celltype == 'Astrocytes')
astrocytes_obj = astrocytes_obj %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20)
for (resolution in c(0.2)) {
    astrocytes_obj = astrocytes_obj %>% 
        FindClusters(resolution = resolution) %>% 
        identity()
    DimPlot(astrocytes_obj, reduction = "umap", pt.size = .1, label=T)
    ggsave(file.path(result_dir, paste0('Astrocytes_Harmony_integrated_UMAP_louvain', resolution, '.png')),
        width=5, height=5)
    DimPlot(astrocytes_obj, reduction = "umap", pt.size = .1, split.by = 'condition', label=T)
    ggsave(file.path(result_dir, paste0('Astrocytes_Harmony_integrated_UMAP_louvain', resolution, '_condition_split.png')),
        width=10, height=5)
}
saveRDS(astrocytes_obj, file.path(result_dir, 'Astrocytes_Seurat_obj.RDS'))

# proportion difference
annotated_celltype_table = table(astrocytes_obj@meta.data$RNA_snn_res.0.2, astrocytes_obj@meta.data$condition)
annotated_celltype_table = t(t(annotated_celltype_table) / colSums(annotated_celltype_table))
annotated_celltype_df = as.data.frame(annotated_celltype_table)
colnames(annotated_celltype_df) = c('subtype', 'condition', 'prop')
library(ggalluvial)
g = ggplot(annotated_celltype_df, aes(x=condition, y=prop, fill=subtype)) + 
          geom_flow(aes(alluvium = subtype), alpha= .5, color = "white",
          curve_type = "linear", 
          width = .5) +
          geom_col(width = .5, color = "white") +
          scale_y_continuous(NULL, expand = c(0,0)) +
          cowplot::theme_minimal_hgrid() +
          theme(panel.grid.major = element_blank(), 
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank())
ggsave(file.path(result_dir, 'Astrocytes_subtype_proportion_per_condition_lines.pdf'))


library(clusterProfiler)  ## according to the paper, they used the compareCluster with cluster marker genes
library(org.Hs.eg.db)
marker_dir = file.path(result_dir, 'Astrocytes_control_markers')
dir.create(marker_dir, showWarnings = FALSE)

## find cluster specific marker genes for control condition
control_astrocytes_obj = subset(astrocytes_obj, subset=condition == 'control')
# write marker genes for different resolutions
Idents(control_astrocytes_obj) = control_astrocytes_obj$RNA_snn_res.0.2
markers = FindAllMarkers(object = control_astrocytes_obj)
write.csv(markers, file.path(marker_dir, 'Astrocytes_control_cluster_markergenes_resolution_0.2.csv'), quote=F)
for (cluster in unique(astrocytes_obj$RNA_snn_res.0.2)) {
    cell_selection = subset(astrocytes_obj, cells=colnames(astrocytes_obj)[astrocytes_obj@meta.data[, 'RNA_snn_res.0.2'] == cluster])
    cell_selection = SetIdent(cell_selection, value = "condition")
    DEG_cell_selection = FindAllMarkers(cell_selection, test.use = "wilcox", logfc.threshold=0.25, min.pct=0.1, only.pos = F)
    write.csv(DEG_cell_selection, file.path(marker_dir, paste0(gsub('/', '', cluster), '_cluster_DEGs_resolution0.2.csv')), quote=F)
}

## load markers
geneid_list = list()
genename_list = list()
for (cluster in unique(markers$cluster)) {
    cluster_markers = markers[markers$cluster == cluster & markers$p_val_adj < 0.05, ]
    map_res = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, column = "ENTREZID", keytype = "SYMBOL")
    cat('Number of unmapped genes:', sum(is.na(map_res)), 'out of', nrow(cluster_markers), '\n')
    map_res = map_res[!is.na(map_res)]
    geneid_list[[as.character(cluster)]] = unname(map_res)
    genename_list[[as.character(cluster)]] = cluster_markers$gene
}
ck = compareCluster(geneCluster=geneid_list, fun="enrichGO", OrgDb=org.Hs.eg.db)
ck = setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(ck, file.path(marker_dir, 'Astrocytes_control_cluster_markers_enrichGO_resolution0.2.csv'), quote=F)
g = dotplot(ck) + coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(file.path(marker_dir, 'Astrocytes_control_cluster_markers_enrichGO_resolution0.2_dotplot.png'), width=20, height=10)


library(clusterProfiler)  ## according to the paper, they used the compareCluster with cluster marker genes
library(org.Hs.eg.db)
marker_dir = file.path(result_dir, 'Astrocytes_disease_markers')
dir.create(marker_dir, showWarnings = FALSE)

## find cluster specific marker genes for disease condition
disease_astrocytes_obj = subset(astrocytes_obj, subset=condition == 'disease')
# write marker genes for different resolutions
Idents(disease_astrocytes_obj) = disease_astrocytes_obj$RNA_snn_res.0.2
markers = FindAllMarkers(object = disease_astrocytes_obj)
write.csv(markers, file.path(marker_dir, 'Astrocytes_disease_cluster_markergenes_resolution_0.2.csv'), quote=F)
for (cluster in unique(astrocytes_obj$RNA_snn_res.0.2)) {
    cell_selection = subset(astrocytes_obj, cells=colnames(astrocytes_obj)[astrocytes_obj@meta.data[, 'RNA_snn_res.0.2'] == cluster])
    cell_selection = SetIdent(cell_selection, value = "condition")
    DEG_cell_selection = FindAllMarkers(cell_selection, test.use = "wilcox", logfc.threshold=0.25, min.pct=0.1, only.pos = F)
    write.csv(DEG_cell_selection, file.path(marker_dir, paste0(gsub('/', '', cluster), '_cluster_DEGs_resolution0.2.csv')), quote=F)
}

## load markers
geneid_list = list()
genename_list = list()
for (cluster in unique(markers$cluster)) {
    cluster_markers = markers[markers$cluster == cluster & markers$p_val_adj < 0.05, ]
    map_res = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, column = "ENTREZID", keytype = "SYMBOL")
    cat('Number of unmapped genes:', sum(is.na(map_res)), 'out of', nrow(cluster_markers), '\n')
    map_res = map_res[!is.na(map_res)]
    geneid_list[[as.character(cluster)]] = unname(map_res)
    genename_list[[as.character(cluster)]] = cluster_markers$gene
}
ck = compareCluster(geneCluster=geneid_list, fun="enrichGO", OrgDb=org.Hs.eg.db)
ck = setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(ck, file.path(marker_dir, 'Astrocytes_disease_cluster_markers_enrichGO_resolution0.2.csv'), quote=F)
g = dotplot(ck) + coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(file.path(marker_dir, 'Astrocytes_disease_cluster_markers_enrichGO_resolution0.2_dotplot.png'), width=20, height=10)

library(clusterProfiler)  ## according to the paper, they used the compareCluster with cluster marker genes
library(org.Hs.eg.db)
marker_dir = file.path(result_dir, 'Astrocytes_control_disease_markers')
dir.create(marker_dir, showWarnings = FALSE)

# write marker genes for different resolutions
Idents(astrocytes_obj) = astrocytes_obj$RNA_snn_res.0.2
markers = FindAllMarkers(object = astrocytes_obj)
write.csv(markers, file.path(marker_dir, 'Astrocytes_cluster_markergenes_resolution_0.2.csv'), quote=F)
for (cluster in unique(astrocytes_obj$RNA_snn_res.0.2)) {
    cell_selection = subset(astrocytes_obj, cells=colnames(astrocytes_obj)[astrocytes_obj@meta.data[, 'RNA_snn_res.0.2'] == cluster])
    cell_selection = SetIdent(cell_selection, value = "condition")
    DEG_cell_selection = FindAllMarkers(cell_selection, test.use = "wilcox", logfc.threshold=0.25, min.pct=0.1, only.pos = F)
    write.csv(DEG_cell_selection, file.path(marker_dir, paste0(gsub('/', '', cluster), '_cluster_DEGs_resolution0.2.csv')), quote=F)
}

## load markers
geneid_list = list()
genename_list = list()
for (cluster in unique(markers$cluster)) {
    cluster_markers = markers[markers$cluster == cluster & markers$p_val_adj < 0.05, ]
    map_res = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, column = "ENTREZID", keytype = "SYMBOL")
    cat('Number of unmapped genes:', sum(is.na(map_res)), 'out of', nrow(cluster_markers), '\n')
    map_res = map_res[!is.na(map_res)]
    geneid_list[[as.character(cluster)]] = unname(map_res)
    genename_list[[as.character(cluster)]] = cluster_markers$gene
}
ck = compareCluster(geneCluster=geneid_list, fun="enrichGO", OrgDb=org.Hs.eg.db)
ck = setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(ck, file.path(marker_dir, 'Astrocytes_cluster_markers_enrichGO_resolution0.2.csv'), quote=F)
g = dotplot(ck) + coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(file.path(marker_dir, 'Astrocytes_cluster_markers_enrichGO_resolution0.2_dotplot.png'), width=20, height=10)


