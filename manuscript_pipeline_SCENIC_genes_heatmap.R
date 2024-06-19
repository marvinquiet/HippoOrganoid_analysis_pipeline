## conda activate Seurat
library(Seurat)
library(Matrix)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
set.seed(2022)
library(harmony)

project_dir = "/projects/compbio/users/wma36/collaborations/JieXu"
result_dir = file.path(project_dir, 'Manuscript_analysis_figures')
dir.create(result_dir, showWarnings = FALSE)
final_fig_dir = file.path(result_dir, 'Manuscript_final_figures')
dir.create(final_fig_dir, showWarnings = FALSE)

SCENIC_output_dir = file.path(result_dir, 'SCENIC_output')
dir.create(SCENIC_output_dir, showWarnings = FALSE)
Astrocytes_obj = readRDS(file.path(result_dir, 'Astrocytes_Seurat_obj.RDS'))
Astrocytes_assignment = read.csv(file.path(SCENIC_output_dir, 'Astrocytes_AUCell_cluster_assignment.csv'), header=T, row.names=1)
common_cells = intersect(colnames(Astrocytes_obj), rownames(Astrocytes_assignment))
Astrocytes_obj = Astrocytes_obj[, common_cells]
Astrocytes_assignment = Astrocytes_assignment[common_cells, ]
Astrocytes_obj$AUC_louvain = Astrocytes_assignment$louvain
table(Astrocytes_obj$AUC_louvain, Astrocytes_obj$RNA_snn_res.0.2)

# load regulons
regulons = c('NKX6-2', 'POU2F1', 'FOS', 'JUNB', 'EGR1', 'SOX9', 'HES5', 'NR2F1', 'FOSB', 'CEBPD', 'JUN')
regulon_dir = file.path(SCENIC_output_dir, 'regulons')
TF_genes = c()
TF_gene_nums = c()
for (regulon in regulons) {
    regulated_genes = scan(file.path(regulon_dir, regulon, paste0(regulon, '_regulon.txt')), what=character())
    TF_genes = c(TF_genes, regulated_genes)
    TF_gene_nums = c(TF_gene_nums, length(regulated_genes))
}

# overlap with Astrocytes DEGs
Astrocytes_DEGs = read.csv(file.path(result_dir, 'DEGs', 'Astrocytes_condition_DEGs.csv'), header=T, row.names=1)
idx = which(Astrocytes_DEGs$p_val_adj < 0.05 & abs(Astrocytes_DEGs$avg_log2FC) > 0.15, ) # 460
Astrocytes_DEGs = Astrocytes_DEGs[idx,]
overlapped_genes = intersect(TF_genes, rownames(Astrocytes_DEGs)) # 106 out of 460 DEGs
# fisher
test_mat = matrix(c(106, 460-106,
                    396-106, 36591-396-460+106), 
                  nrow=2)
print(test_mat)
print(sum(test_mat))
print(fisher.test(test_mat))


Astrocytes_DEGs = Astrocytes_DEGs[overlapped_genes, ]
Astrocytes_DEGs = Astrocytes_DEGs[order(-Astrocytes_DEGs$avg_log2FC), ]
# scale Astrocytes
Astrocytes_obj = ScaleData(Astrocytes_obj, features=rownames(Astrocytes_obj))
# order by TF
regulons = c('CEBPD', 'JUN', 'NR2F1', 'SOX9')
regulon_dir = file.path(SCENIC_output_dir, 'regulons')
genes = c()
gene_nums = c()
TF_nums = c()
for (regulon in regulons) {
    regulated_genes = scan(file.path(regulon_dir, regulon, paste0(regulon, '_regulon.txt')), what=character())
    overlapped_genes = intersect(regulated_genes, rownames(Astrocytes_DEGs))
    overlapped_DEG_FC_df = Astrocytes_DEGs[overlapped_genes, ]
    overlapped_DEG_FC_df = overlapped_DEG_FC_df[order(-overlapped_DEG_FC_df$avg_log2FC), ]
    genes = c(genes, rownames(overlapped_DEG_FC_df))
    gene_nums = c(gene_nums, nrow(overlapped_DEG_FC_df))
    TF_nums = c(TF_nums, rep(regulon, nrow(overlapped_DEG_FC_df)))
}
mat = Astrocytes_obj@assays$RNA@scale.data[genes, ] # 916x2602
ha = Astrocytes_obj$condition
ha_dend = cluster_within_group(mat, ha)
ca = TF_nums
ca_dend = cluster_within_group(t(mat), ca)
library(RColorBrewer)
colors = brewer.pal(n=11, 'RdBu')
col_fun_prop = colorRamp2(c(-2, -0.5, 0, 0.5, 2), 
                            c(colors[11], colors[7], colors[6], colors[5], colors[1]))
pdf(file.path(SCENIC_output_dir, 'Heatmap_ComplexHeatmap_overlapped_orderedbyFC_orderedbyTF_Astrocytes_fourregulons.pdf'), width=6, height=12)
h = Heatmap(mat, name="Gene Exprs", cluster_columns=ha_dend, column_split=2,
            cluster_rows=ca_dend, row_split=length(regulons),
            # cluster_rows=F, row_split=TF_nums,
            col = col_fun_prop, 
            show_column_names=F, row_names_gp = gpar(fontsize = 6),
            heatmap_legend_param = list(title = "Gene Exprs", at = c(-2, -0.5, 0, 0.5, 2)))
#draw(h, merge_legend=T)
draw(h)
dev.off()

regulons = c('NKX6-2', 'POU2F1', 'FOS', 'JUNB', 'EGR1', 'HES5', 'FOSB')
regulon_dir = file.path(SCENIC_output_dir, 'regulons')
genes = c()
gene_nums = c()
TF_nums = c()
for (regulon in regulons) {
    regulated_genes = scan(file.path(regulon_dir, regulon, paste0(regulon, '_regulon.txt')), what=character())
    overlapped_genes = intersect(regulated_genes, rownames(Astrocytes_DEGs))
    overlapped_DEG_FC_df = Astrocytes_DEGs[overlapped_genes, ]
    overlapped_DEG_FC_df = overlapped_DEG_FC_df[order(-overlapped_DEG_FC_df$avg_log2FC), ]
    genes = c(genes, rownames(overlapped_DEG_FC_df))
    gene_nums = c(gene_nums, nrow(overlapped_DEG_FC_df))
    TF_nums = c(TF_nums, rep(regulon, nrow(overlapped_DEG_FC_df)))
}
mat = Astrocytes_obj@assays$RNA@scale.data[genes, ] # 916x2602
mat = as.matrix(mat)
# cluster within each row group
ha = Astrocytes_obj$condition
ha_dend = cluster_within_group(mat, ha)
h = Heatmap(mat, name="Gene Exprs", cluster_columns=ha_dend, column_split=2, cluster_rows=F, row_split=TF_nums,
            show_column_names=F, row_names_gp = gpar(fontsize = 6))
pdf(file.path(SCENIC_output_dir, 'Heatmap_ComplexHeatmap_overlapped_orderedbyFC_orderedbyTF_Astrocytes_sevenregulons.pdf'), width=6, height=18)
draw(h, merge_legend=T)
dev.off()

