## conda activate Seurat
library(monocle3)
library(Seurat)
library(Matrix)
library(ggplot2)
set.seed(2022)

project_dir = "/projects/compbio/users/wma36/collaborations/JieXu"
result_dir = file.path(project_dir, 'Manuscript_analysis_figures')
dir.create(result_dir, showWarnings = FALSE)
final_fig_dir = file.path(result_dir, 'Manuscript_final_figures')
dir.create(final_fig_dir, showWarnings = FALSE)
## --------------
# 1. Seurat pipeline - overall
# Harmony integration, clustering, cell type annotation, cell type proportion, cell type DEGs
## --------------
## get inputFiles
samples = c("p22262-s009_C12N-1", "p22262-s002_C1-2-2", "p22262-s005_IIC1-1", 
            "p22262-s012_FXS-2", "p22262-s003_FXSB2-1", "p22262-s007_soloc3-1")
sample_obj_list = list()
for (sample in samples) {
    input_dir = file.path(project_dir, sample, 'outs', 'filtered_feature_bc_matrix')
    input_mat = readMM(file.path(input_dir, 'matrix.mtx.gz'))
    input_barcodes = scan(file.path(input_dir, 'barcodes.tsv.gz'), what=character())
    input_features = read.csv(file.path(input_dir, 'features.tsv.gz'), header=F, sep='\t')
    dedup_idx = !duplicated(input_features$V2)
    input_mat = input_mat[dedup_idx, ]
    rownames(input_mat) = input_features[dedup_idx, 'V2']
    colnames(input_mat) = input_barcodes

    ## create Seurat object
    sample_obj = CreateSeuratObject(counts=input_mat, project=sample, 
                                    min.cells=0, min.features=0)
    sample_obj_list[[sample]] = sample_obj
}
# run harmony
library(harmony)
total_obj = merge(sample_obj_list[[1]], unlist(sample_obj_list)[2:6],
                    add.cell.ids=names(sample_obj_list)[1:6])
total_obj[["percent.mt"]] = PercentageFeatureSet(total_obj, pattern = "^MT-")  ## already filtered
VlnPlot(total_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(result_dir, 'QuanlityViolinplot.png'))
VlnPlot(total_obj, features='nFeature_RNA') & geom_hline(yintercept=200) & geom_hline(yintercept=8000)
ggsave(file.path(result_dir, 'QuanlityViolinplot_nFeatureRNA.png'))
VlnPlot(total_obj, features='percent.mt') & geom_hline(yintercept=20)
ggsave(file.path(result_dir, 'QuanlityViolinplot_mt-percent.png'))
combined_obj = subset(total_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
combined_obj = Seurat::NormalizeData(combined_obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = combined_obj@var.genes, npcs = 30, verbose = FALSE)
combined_obj@meta.data$condition = ifelse(Idents(combined_obj) %in% samples[1:3], 'control', 'disease')
combined_obj = combined_obj %>% 
        RunHarmony("condition", plot_convergence = TRUE)
DimPlot(object = combined_obj, reduction = "harmony", pt.size = .1, group.by = "condition")
ggsave(file.path(result_dir, 'Harmony_integrated_UMAP_PCA_condition.png'),
    width=5, height=5)
combined_obj = combined_obj %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20)
for (resolution in c(0.8, 1.5)) {
    combined_obj = combined_obj %>% 
        FindClusters(resolution = resolution) %>% 
        identity()
    DimPlot(combined_obj, reduction = "umap", pt.size = .1, label=T)
    ggsave(file.path(result_dir, paste0('Harmony_integrated_UMAP_louvain', resolution, '.png')),
        width=5, height=5)
    DimPlot(combined_obj, reduction = "umap", pt.size = .1, split.by = 'condition', label=T)
    ggsave(file.path(result_dir, paste0('Harmony_integrated_UMAP_louvain', resolution, '_condition_split.png')),
        width=10, height=5)
}
saveRDS(combined_obj, file.path(result_dir, 'Harmony_integrated_Seurat_obj.RDS'))

## Feature plot
marker_dir = file.path(result_dir, 'markers_featureplot')
dir.create(marker_dir, showWarnings=F)
library(readxl)
marker_info = read_excel(file.path(project_dir, 'scRNAseq annatation_marker list 052923.xlsx'),  ## updated marker gene list
                         col_names=c('marker', 'func'))
marker_info = as.data.frame(marker_info)
for (gene_idx in 1:nrow(marker_info)) {
    marker = toupper(marker_info[gene_idx, 'marker'])
    func = marker_info[gene_idx, 'func']
    if (!marker %in% rownames(combined_obj)) {
        cat(marker, 'not found!\n')
        next()
    }
    FeaturePlot(combined_obj, features=marker, label=T) + 
        ggtitle(paste(marker, '(', func, ')'))
    ggsave(file.path(marker_dir, 
                     paste0('Harmony_integrated_', marker, '.png')),
           width=5, height=5)
    FeaturePlot(combined_obj, features=marker, split.by='condition', label=T) & theme(legend.position = c(0.1,0.2)) & ylab(paste(marker, '(', func, ')'))
    ggsave(file.path(marker_dir, 
                     paste0('Harmony_integrated_', marker, '_split.png')),
           width=10, height=5)
}

## Dotplot
order_idx = order(marker_info$func)
marker_info = marker_info[order_idx, ]
features = toupper((marker_info$marker))
features = features[!duplicated(features)]
for (resolution in c(0.8, 1.5)) {
    Idents(combined_obj) = combined_obj@meta.data[, paste0('RNA_snn_res.', resolution)]
    DotPlot(combined_obj, features=features) & coord_flip() & theme_bw()
    ggsave(file.path(marker_dir, paste0('Harmony_integrated_markers_louvain', resolution, '_dotplot.png')), height=20, width=10) 
}


## annotation based on 0.8
combined_obj@meta.data$annotated_celltype = NA
cell_idx = combined_obj@meta.data$RNA_snn_res.0.8 %in% c(2, 4)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'OPCs'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.8 %in% c(0, 1, 5, 13, 21)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Astrocytes'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.8 %in% c(15)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Unknown'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.8 %in% c(16, 19)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Immature neurons'
cell_idx = combined_obj@meta.data$RNA_snn_res.0.8 %in% c(22)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Removed'
# annotate based on 1.5
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(12, 14, 19)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Excitatory neurons'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(25)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Inhibitory neurons'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(28, 29)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Hippocampus NPCs'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(11, 16, 23, 26)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Proliferating radial glia'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(22)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Removed'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(30)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Choroid plexus cells'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(3, 13, 24)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Astrocytes'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(8, 9)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Radial glia'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(10, 15, 32)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'OPCs'
cell_idx = combined_obj@meta.data$RNA_snn_res.1.5 %in% c(34)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Removed'
cell_idx = is.na(combined_obj@meta.data$annotated_celltype)
combined_obj@meta.data[cell_idx, 'annotated_celltype'] = 'Removed'
## do not filter out unknown cell types
removed_cells = combined_obj@meta.data$annotated_celltype == 'Removed'
annotated_combined_obj = combined_obj[, !removed_cells] # 34023 cells out of 34713
umap_coordinates = as.data.frame(annotated_combined_obj[['umap']]@cell.embeddings)
umap_coordinates$UMAP_1 = -umap_coordinates$UMAP_1
umap_coordinates$UMAP_2 = -umap_coordinates$UMAP_2
annotated_combined_obj[['umap_flip']] = CreateDimReducObject(embeddings = as.matrix(umap_coordinates), 
                                        key = "UMAP_flipped_", assay = DefaultAssay(annotated_combined_obj))
DimPlot(annotated_combined_obj, reduction = "umap_flip", group.by='annotated_celltype', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_UMAP.png'),
    width=6, height=5) 
DimPlot(annotated_combined_obj, reduction = "umap_flip", group.by='annotated_celltype', pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_UMAP_condition_split_unlabeled.png'),
    width=10, height=5)
g = DimPlot(annotated_combined_obj, reduction = "umap_flip", group.by='annotated_celltype', pt.size = .1, split.by = 'condition', label=F)
saveRDS(g, file.path(final_fig_dir, 'Fig3A_Seurat_annotation.RDS'))
saveRDS(annotated_combined_obj, file.path(result_dir, 'annotated_Seurat_obj.RDS'))

# use annotated ones to re-run UMAP
annotated_combined_obj = annotated_combined_obj %>% 
    RunUMAP(reduction = "harmony", dims = 1:20)
umap_coordinates = as.data.frame(annotated_combined_obj[['umap']]@cell.embeddings)
umap_coordinates$UMAP_1 = -umap_coordinates$UMAP_1
umap_coordinates$UMAP_2 = -umap_coordinates$UMAP_2
annotated_combined_obj[['umap_flip']] = CreateDimReducObject(embeddings = as.matrix(umap_coordinates), 
                                        key = "UMAP_flipped_", assay = DefaultAssay(annotated_combined_obj))
DimPlot(annotated_combined_obj, reduction = "umap_flip", group.by='annotated_celltype', pt.size = .1, label=F)
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_UMAP_rerun.png'),
    width=6, height=5)
DimPlot(annotated_combined_obj, reduction = "umap_flip", group.by='annotated_celltype', pt.size = .1, split.by = 'condition', label=F)
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_UMAP_condition_split_unlabeled_rerun.png'),
    width=10, height=5)
 
## cell type proportion
annotated_celltype_table = table(annotated_combined_obj@meta.data$annotated_celltype, annotated_combined_obj@meta.data$condition)
annotated_celltype_table = t(t(annotated_celltype_table) / colSums(annotated_celltype_table))
annotated_celltype_df = as.data.frame(annotated_celltype_table)
colnames(annotated_celltype_df) = c('celltype', 'condition', 'prop')
library(ggalluvial)
g = ggplot(annotated_celltype_df, aes(x=condition, y=prop, fill=celltype)) + 
          geom_flow(aes(alluvium = celltype), alpha= .5, color = "white",
          curve_type = "linear", 
          width = .5) +
          geom_col(width = .5, color = "white") +
          scale_y_continuous(NULL, expand = c(0,0)) +
          cowplot::theme_minimal_hgrid() +
          theme(panel.grid.major = element_blank(), 
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank())
saveRDS(g, file.path(final_fig_dir, 'Fig3B_celltype_proportion.RDS'))
ggsave(file.path(result_dir, 'annotated_Harmony_integrated_3v3_celltype_proportion_per_condition_lines.png'))

## === DEG analysis
## per condition
DEG_dir = file.path(result_dir, 'DEGs')
dir.create(DEG_dir, showWarnings=F)
Idents(annotated_combined_obj) = annotated_combined_obj@meta.data$condition
condition_DEGs = FindMarkers(annotated_combined_obj, logfc.threshold=0, ident.1="disease", ident.2="control")  ## default threshold, logfc.threshold = 0.25 & min.pct = 0.1 -> with Wilcox test
condition_DEGs = condition_DEGs[condition_DEGs$p_val_adj < 0.05, ]
write.csv(condition_DEGs, file.path(DEG_dir, 'condition_DEGs.csv'), quote=F)
up_condition_DEGs = condition_DEGs[condition_DEGs$avg_log2FC > 0, ]
write.csv(up_condition_DEGs, file.path(DEG_dir, 'condition_DEGs_UP.csv'), quote=F)
down_condition_DEGs = condition_DEGs[condition_DEGs$avg_log2FC < 0, ]
write.csv(down_condition_DEGs, file.path(DEG_dir, 'condition_DEGs_DOWN.csv'), quote=F)

## per cell type and condition
Idents(annotated_combined_obj) = annotated_combined_obj@meta.data$condition
for (celltype in unique(annotated_combined_obj@meta.data$annotated_celltype)) {
    celltype_combined_obj = subset(annotated_combined_obj, subset = annotated_celltype == as.character(celltype))
    celltype_DEGs = FindMarkers(celltype_combined_obj, logfc.threshold=0, ident.1="disease", ident.2="control")
    celltype_DEGs = celltype_DEGs[celltype_DEGs$p_val_adj < 0.05, ]
    write.csv(celltype_DEGs, file.path(DEG_dir, paste0(gsub('/', '-', celltype), '_condition_DEGs.csv')), quote=F)
    up_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC > 0, ]
    write.csv(up_celltype_DEGs, file.path(DEG_dir, paste0(gsub('/', '-', celltype), '_condition_DEGs_UP.csv')), quote=F)
    down_celltype_DEGs = celltype_DEGs[celltype_DEGs$avg_log2FC < 0, ]
    write.csv(down_celltype_DEGs, file.path(DEG_dir, paste0(gsub('/', '-', celltype), '_condition_DEGs_DOWN.csv')), quote=F)
}

## === GO analysis for cell type specific markers
library(clusterProfiler)  ## according to the paper, they used the compareCluster with cluster marker genes
library(org.Hs.eg.db)
## find cluster specific marker genes for control condition
control_annotated_combined_obj = subset(annotated_combined_obj, subset=condition == 'control')
# write marker genes for different resolutions
Idents(control_annotated_combined_obj) = control_annotated_combined_obj$annotated_celltype
markers = FindAllMarkers(object = control_annotated_combined_obj)
write.csv(markers, file.path(DEG_dir, 'control_celltype_markergenes.csv'), quote=F)
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
#ck = compareCluster(geneCluster=genename_list, fun="enrichGO", OrgDb=org.Hs.eg.db)  ## need gene Entrez ID as input..
ck = compareCluster(geneCluster=geneid_list, fun="enrichGO", OrgDb=org.Hs.eg.db)
ck = setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(ck, file.path(DEG_dir, 'control_celltype_markers_enrichGO.csv'), quote=F)
g = dotplot(ck) + coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(file.path(DEG_dir, 'control_cluster_markers_dotplot_enrichGO.png'), width=20, height=10)

## === GO analysis for cell-type-specific DEGs
library(clusterProfiler)  ## according to the paper, they used the compareCluster with cluster marker genes
library(org.Hs.eg.db)
## load cell type up DEGs
geneid_list = list()
genename_list = list()
for (celltype in unique(annotated_combined_obj@meta.data$annotated_celltype)) {
    celltype_up_markers = read.csv(file.path(DEG_dir, paste0(gsub('/', '-', celltype), '_condition_DEGs_UP.csv')), header=T, row.names=1)
    map_res = mapIds(org.Hs.eg.db, keys = rownames(celltype_up_markers), column = "ENTREZID", keytype = "SYMBOL")
    cat('Number of unmapped genes:', sum(is.na(map_res)), 'out of', nrow(celltype_up_markers), '\n')
    map_res = map_res[!is.na(map_res)]
    geneid_list[[as.character(celltype)]] = unname(map_res)
    genename_list[[as.character(celltype)]] = rownames(celltype_up_markers)
}
ck = compareCluster(geneCluster=geneid_list, fun="enrichGO", OrgDb=org.Hs.eg.db)
ck = setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(ck, file.path(DEG_dir, 'celltype_up_DEGs_enrichGO.csv'), quote=F)
g = dotplot(ck) + coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(file.path(DEG_dir, 'celltype_up_DEGs_dotplot_enrichGO.png'), width=20, height=10)
# load cell type down DEGs
geneid_list = list()
genename_list = list()
for (celltype in unique(annotated_combined_obj@meta.data$annotated_celltype)) {
    celltype_down_markers = read.csv(file.path(DEG_dir, paste0(gsub('/', '-', celltype), '_condition_DEGs_DOWN.csv')), header=T, row.names=1)
    map_res = mapIds(org.Hs.eg.db, keys = rownames(celltype_down_markers), column = "ENTREZID", keytype = "SYMBOL")
    cat('Number of unmapped genes:', sum(is.na(map_res)), 'out of', nrow(celltype_down_markers), '\n')
    map_res = map_res[!is.na(map_res)]
    geneid_list[[as.character(celltype)]] = unname(map_res)
    genename_list[[as.character(celltype)]] = rownames(celltype_down_markers)
}
ck = compareCluster(geneCluster=geneid_list, fun="enrichGO", OrgDb=org.Hs.eg.db)
ck = setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(ck, file.path(DEG_dir, 'celltype_down_DEGs_enrichGO.csv'), quote=F)
g = dotplot(ck) + coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(file.path(DEG_dir, 'celltype_down_DEGs_dotplot_enrichGO.png'), width=20, height=10)


## --- Monocle3 analysis
library(SeuratWrappers)
library(monocle3)
cds = as.cell_data_set(annotated_combined_obj)
cds = cluster_cells(cds, reduction_method='UMAP', cluster_method='louvain')
reducedDims(cds)$UMAP = reducedDims(cds)$UMAP_FLIP
cds = learn_graph(cds)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition") + facet_wrap(~condition, ncol=2)
ggsave(file.path(result_dir, 'monocle3_partitions_by_condition.png'), width=20, height=10)
plot_cells(cds,
           color_cells_by = "annotated_celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE) +
           facet_wrap(~condition, ncol=2)
ggsave(file.path(result_dir, 'monocle3_trajectory_by_condition.png'), width=20, height=10)
## use proliferating radial glia as starting point
get_earliest_principal_node = function(cds, annotated_celltype="Proliferating radial glia"){
  cell_ids <- which(colData(cds)[, "annotated_celltype"] == annotated_celltype)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave(file.path(result_dir, 'monocle3_trajectory_by_pseudotime.png'), width=10, height=10)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + 
           facet_wrap(~condition, ncol=2)
ggsave(file.path(result_dir, 'monocle3_trajectory_by_pseudotime_and_condition.png'), width=20, height=10)

