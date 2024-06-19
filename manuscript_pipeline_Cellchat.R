## conda activate Seurat
## Follow tutorials from: https://github.com/sqjin/CellChat

library(CellChat)
library(svglite)
library(Seurat)
library(patchwork)
options(stringsAsFactors = FALSE)

set.seed(2022)

CellChatDB = CellChatDB.human
#CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use = CellChatDB
color.use = c('#F8766D', '#D89000', '#39B600', '#A3A500', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC')
patterns_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

project_dir = "/projects/compbio/users/wma36/collaborations/JieXu"
result_dir = file.path(project_dir, 'Manuscript_analysis_figures')
dir.create(result_dir, showWarnings = FALSE)
final_fig_dir = file.path(result_dir, 'Manuscript_final_figures')
dir.create(final_fig_dir, showWarnings = FALSE)
annotated_obj = readRDS(file.path(result_dir, 'annotated_Seurat_obj.RDS'))
# combine neurons
annotated_obj@meta.data$combined_celltype = annotated_obj@meta.data$annotated_celltype
idx = annotated_obj$annotated_celltype %in% c('Immature neurons', 'Inhibitory neurons', 'Excitatory neurons')
annotated_obj@meta.data[idx, 'combined_celltype'] = 'Neurons'
# create cellchat output
cellchat_output_dir = file.path(result_dir, 'Cellchat_output')
dir.create(cellchat_output_dir, showWarnings = FALSE)
run_cellchat = function(Seurat_obj, condition='control') {
    if (condition == "all") {
        obj = Seurat_obj
    } else {
        idx = Seurat_obj$condition == condition
        obj = Seurat_obj[, idx]
        # obj = subset(Seurat_obj, subset= condition == condition) ## this does not work..
    }
    norm_data = GetAssayData(obj, slot='data')
    metadata = obj@meta.data
    cellchat_obj = createCellChat(obj=norm_data, meta=metadata, group.by='combined_celltype')
    cellchat_obj@DB = CellChatDB.use
    # cellchat process
    cellchat_obj = subsetData(cellchat_obj) # This step is necessary even if using the whole database
    cellchat_obj = identifyOverExpressedGenes(cellchat_obj)
    cellchat_obj = identifyOverExpressedInteractions(cellchat_obj)
    cellchat_obj = computeCommunProb(cellchat_obj)
    df.net = subsetCommunication(cellchat_obj)
    write.csv(df.net, file.path(cellchat_output_dir, paste0(condition, '_CellChat_res.csv')), quote=F)
    cellchat_obj = computeCommunProbPathway(cellchat_obj)
    cellchat_obj = aggregateNet(cellchat_obj)
    return(cellchat_obj)
}
#control_cellchat_obj = run_cellchat(annotated_obj, condition="control")
#disease_cellchat_obj = run_cellchat(annotated_obj, condition="disease")
#saveRDS(control_cellchat_obj, file.path(cellchat_output_dir, 'control_cellchat_obj.RDS'))
#saveRDS(disease_cellchat_obj, file.path(cellchat_output_dir, 'disease_cellchat_obj.RDS'))

control_cellchat_obj = readRDS(file.path(cellchat_output_dir, 'control_cellchat_obj.RDS'))
disease_cellchat_obj = readRDS(file.path(cellchat_output_dir, 'disease_cellchat_obj.RDS'))
#' identify glocal pattern
global_dir = file.path(cellchat_output_dir, 'global_cell_comm_patterns')
dir.create(global_dir, showWarnings=F)
library(NMF)
library(ggalluvial)
# interactively select K
# === control, outgoing
color.use = c('#F8766D', '#D89000', '#39B600', '#A3A500', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC')
patterns_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
selectK(control_cellchat_obj, pattern = "outgoing")
ggsave(file.path(global_dir, 'control_outgoing_selectK.png'))
nPatterns = 4
control_cellchat_obj = identifyCommunicationPatterns(control_cellchat_obj, pattern = "outgoing", k = nPatterns)
g = netAnalysis_river(control_cellchat_obj, pattern = "outgoing", 
                      color.use=color.use, color.use.pattern=patterns_color[1:nPatterns], 
                      font.size = 7)
ggsave(file.path(global_dir, 'control_outgoing_comm_patterns_riverplot.pdf'), width=15, height=10)
# === control, incoming
selectK(control_cellchat_obj, pattern = "incoming")
ggsave(file.path(global_dir, 'control_incoming_selectK.png'))
nPatterns = 4
control_cellchat_obj = identifyCommunicationPatterns(control_cellchat_obj, pattern = "incoming", k = nPatterns)
g = netAnalysis_river(control_cellchat_obj, pattern = "incoming",
                      color.use=color.use, color.use.pattern=patterns_color[1:nPatterns],
                      font.size = 7)
ggsave(file.path(global_dir, 'control_incoming_comm_patterns_riverplot.pdf'), width=15, height=10)
# === disease, outgoing
selectK(disease_cellchat_obj, pattern = "outgoing")
ggsave(file.path(global_dir, 'disease_outgoing_selectK.png'))
nPatterns = 4
disease_cellchat_obj = identifyCommunicationPatterns(disease_cellchat_obj, pattern = "outgoing", k = nPatterns)
netAnalysis_river(disease_cellchat_obj, pattern = "outgoing",
                  color.use=color.use, color.use.pattern=patterns_color[1:nPatterns],
                  font.size = 7)
ggsave(file.path(global_dir, 'disease_outgoing_comm_patterns_riverplot.pdf'), width=15, height=10)
# === disease, incoming
selectK(disease_cellchat_obj, pattern = "incoming")
ggsave(file.path(global_dir, 'disease_incoming_selectK.png'))
nPatterns = 4
disease_cellchat_obj = identifyCommunicationPatterns(disease_cellchat_obj, pattern = "incoming", k = nPatterns)
netAnalysis_river(disease_cellchat_obj, pattern = "incoming",
                  color.use=color.use, color.use.pattern=patterns_color[1:nPatterns],
                  font.size = 7)
ggsave(file.path(global_dir, 'disease_incoming_comm_patterns_riverplot.pdf'), width=15, height=10)

#' draw incoming/outgoing dotplot
union_pathways = union(disease_cellchat_obj@netP$pathways, control_cellchat_obj@netP$pathways)
direction = 'incoming'
g = netAnalysis_dot(disease_cellchat_obj, pattern = direction, pathway.show = union_pathways, shape=0, font.size=16, color.use=color.use) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    scale_size_continuous(breaks=c(0.25, 0.50, 0.75), range=c(0.5, 2.5))
ggsave(file.path(global_dir, paste0('disease_', direction, '_comm_patterns_dotplot.pdf')), width=8, height=5)
netAnalysis_dot(control_cellchat_obj, pattern = direction, pathway.show = union_pathways, font.size=16, color.use=color.use) + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_size_continuous(breaks=c(0.25, 0.50, 0.75), range=c(0.5, 2.5))
ggsave(file.path(global_dir, paste0('control_', direction, '_comm_patterns_dotplot.pdf')), width=8, height=5)
direction = 'outgoing'
netAnalysis_dot(disease_cellchat_obj, pattern = direction, pathway.show = union_pathways, shape=0, font.size=16, color.use=color.use) + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_size_continuous(breaks=c(0.25, 0.50, 0.75), range=c(0.5, 2.5))
ggsave(file.path(global_dir, paste0('disease_', direction, '_comm_patterns_dotplot.pdf')), width=8, height=5)
netAnalysis_dot(control_cellchat_obj, pattern = direction, pathway.show = union_pathways, font.size=16, color.use=color.use) + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_size_continuous(breaks=c(0.25, 0.50, 0.75), range=c(0.5, 2.5))
ggsave(file.path(global_dir, paste0('control_', direction, '_comm_patterns_dotplot.pdf')), width=8, height=5)

#' Comparison analysis: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html#load-cellchat-object-of-each-dataset-and-then-merge-together
joint_plot_dir = file.path(cellchat_output_dir, 'disease_vs_control')
dir.create(joint_plot_dir, showWarnings=F)
control_cellchat_obj = netAnalysis_computeCentrality(control_cellchat_obj)
disease_cellchat_obj = netAnalysis_computeCentrality(disease_cellchat_obj)
obj_list = list(control=control_cellchat_obj, disease=disease_cellchat_obj)
merged_obj = mergeCellChat(obj_list, add.names = names(obj_list))
# circle plot normalized by max weight
weight.max = getMaxWeight(obj_list, attribute = c("idents","count"))
pdf(file.path(joint_plot_dir, 'Num_of_interactions_normalized_control.pdf'), width=6, height=6, pointsize=12)
netVisual_circle(obj_list[[1]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                 edge.width.max = 12, title.name = paste0("Number of interactions - ", names(obj_list)[1]),
                 color.use=color.use, margin=0.3, vertex.label.cex=1.3)
dev.off()
pdf(file.path(joint_plot_dir, 'Num_of_interactions_normalized_disease.pdf'), width=6, height=6, pointsize=12)
netVisual_circle(obj_list[[2]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                 edge.width.max = 12, title.name = paste0("Number of interactions - ", names(obj_list)[2]),
                 color.use=color.use, margin=0.3, vertex.label.cex=1.3)
dev.off()

# subset object to neurons, astrocytes, proliferating radial glia
subset_control_cellchat_obj = subsetCellChat(control_cellchat_obj, 
                            idents.use=c('Neurons', 'Astrocytes', 'Proliferating radial glia'))
subset_disease_cellchat_obj = subsetCellChat(disease_cellchat_obj, 
                            idents.use=c('Neurons', 'Astrocytes', 'Proliferating radial glia'))
subset_obj_list = list(control=subset_control_cellchat_obj, disease=subset_disease_cellchat_obj)
weight.max = getMaxWeight(subset_obj_list, attribute = c("idents","count"))
pdf(file.path(joint_plot_dir, 'Num_of_interactions_normalized_control_threetypes.pdf'), width=6, height=6, pointsize=12)
netVisual_circle(subset_obj_list[[1]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                 edge.width.max = 12, title.name = paste0("Number of interactions - ", names(obj_list)[1]),
                 color.use=color.use, margin=0.3, vertex.label.cex=1.3)
dev.off()
pdf(file.path(joint_plot_dir, 'Num_of_interactions_normalized_disease_threetypes.pdf'), width=6, height=6, pointsize=12)
netVisual_circle(subset_obj_list[[2]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                 edge.width.max = 12, title.name = paste0("Number of interactions - ", names(obj_list)[2]),
                 color.use=color.use, margin=0.3, vertex.label.cex=1.3)
dev.off()

#' identify conserved and context-specific signaling pathways
merged_obj = computeNetSimilarityPairwise(merged_obj, type = "functional")
merged_obj = netEmbedding(merged_obj, type = "functional")
merged_obj = netClustering(merged_obj, type = "functional") ## need package future 1.20.1 version
pdf(file.path(joint_plot_dir, 'functional_similarity.pdf'),
    width=8, height=8, pointsize=12)
netVisual_embeddingPairwise(merged_obj, type = "functional", label.size = 5)
dev.off()
rankSimilarity(merged_obj, type = "functional", font.size=12)
ggsave(file.path(joint_plot_dir, 'functional_similarity_rank.pdf'), height=6, width=3)
# Figure H: overall information flow
g = rankNet(merged_obj, mode = "comparison", stacked = T, do.stat = TRUE, color.use=c('#FF6666', '#99CCFF'), font.size=12)
ggsave(file.path(joint_plot_dir, "signaling_pathways_overall_relative_information_flow.pdf"), width=4, height=6)

## dotplot on signaling pathway
pairLR_use = c('LRRC4B_PTPRF', 'NRXN1_NLGN1', 'NRXN1_NLGN2', 'NRXN2_NLGN1', 'NRXN2_NLGN2',
               'NRXN1_NLGN2', 'NRXN2_NLGN2', 'PTN_SDC3', 'SEMA5B_PLXNA1', 'SEMA6D_PLXNA1',
               'EFNB3_EPHB1', 'NCAM1_FGFR1', 'NCAM1_NCAM1', 'NCAM1_NCAM2', 'NCAM1_FGFR1',
               'NCAM1_NCAM1', 'NCAM1_NCAM2', 'DLL3_NOTCH1', 'DLL3_NOTCH2', 'DLL3_NOTCH3',
               'DLL3_NOTCH1')
pairLR_use = as.data.frame(unique(pairLR_use))
colnames(pairLR_use) = c('interaction_name')
bubble_data = netVisual_bubble(merged_obj, sources.use=c('Neurons', 'Astrocytes'), targets.use=c('Astrocytes', 'Neurons', 'Proliferating radial glia'), 
                 comparison=c(1,2), angle.x=45, pairLR.use = pairLR_use, return.data=T)
bubble_comm_data = bubble_data$communication
subset_source_to_target = c("Neurons -> Astrocytes (control)", "Neurons -> Astrocytes (disease)",
                            "Astrocytes -> Neurons (control)", "Astrocytes -> Neurons (disease)",
                            "Neurons -> Proliferating radial glia (control)",
                            "Neurons -> Proliferating radial glia (disease)")
subset_bubble_comm_data = bubble_comm_data[bubble_comm_data$source.target %in% subset_source_to_target, ]
df = subset_bubble_comm_data
g = ggplot(df, aes(x=source.target, y=interaction_name_2, color = prob, size = pval)) + 
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1),
                        axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
g = g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
color.use = RColorBrewer::brewer.pal(n = 10, name = 'Spectral')
color.grid = "grey90"
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
  g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                  breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
} else {
  g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
}
g <- g + theme(text = element_text(size = 10),plot.title = element_text(size=10)) +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
if (length(unique(df$source.target)) > 1) {
  g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
}
if (length(unique(df$interaction_name_2)) > 1) {
  g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
}

dataset.name = names(merged_obj@net)
comparison = c(1,2)
group.names <- paste(rep(levels(df$source), each = length(levels(df$target))), levels(df$target), sep = " -> ")
group.names0 <- group.names
xintercept = seq(0.5+length(dataset.name[comparison]), length(group.names0)*length(dataset.name[comparison]), by = length(dataset.name[comparison]))
g <- g + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = 0.2)

group <- 1:2
names(group) <- c('control', 'disease')
color <- ggPalette(length(unique(group)))
names(color) <- names(group[!duplicated(group)])
color <- color[group]
#names(color) <- dataset.name[comparison]
dataset.name.order <- levels(df$source.target)
dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
xtick.color <- color[dataset.name.order]
g <- g + theme(axis.text.x = element_text(colour = xtick.color))
ggsave(file.path(joint_plot_dir, 'signaling_communication_probs_dotplot.pdf'), width=10, height=6)

