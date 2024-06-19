## conda activate Seurat
## Follow tutorials from: https://github.com/sqjin/CellChat

library(CellChat)
library(svglite)
library(Seurat)
library(patchwork)
options(stringsAsFactors = FALSE)

set.seed(2022)

CellChatDB = CellChatDB.human
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
neurons_annotated_obj = subset(annotated_obj, subset=combined_celltype=='Neurons')
# load astrocytes object
astrocytes_annotated_obj = readRDS(file.path(result_dir, 'Astrocytes_Seurat_obj.RDS'))
merged_obj = merge(neurons_annotated_obj, astrocytes_annotated_obj)
merged_obj$final_celltype = NA
idx = which(merged_obj$combined_celltype == 'Neurons')
merged_obj@meta.data[idx, 'final_celltype'] = 'Neurons'
idx = which(merged_obj$annotated_celltype == 'Astrocytes')
merged_obj@meta.data[idx, 'final_celltype'] = paste0('Astrocytes_', merged_obj@meta.data[idx, 'RNA_snn_res.0.2'])

# create cellchat output
cellchat_output_dir = file.path(result_dir, 'Cellchat_Neurons_Astrocytes_subtypes_output')
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
    cellchat_obj = createCellChat(obj=norm_data, meta=metadata, group.by='final_celltype')
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
#control_cellchat_obj = run_cellchat(merged_obj, condition="control")
#disease_cellchat_obj = run_cellchat(merged_obj, condition="disease")
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
color.use = c('#F8766D', '#D89000', '#39B600', '#00B0F6', '#9590FF', '#E76BF3', '#A3A500')
patterns_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
selectK(control_cellchat_obj, pattern = "outgoing")
ggsave(file.path(global_dir, 'control_outgoing_selectK.png'))
nPatterns = 5
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
nPatterns = 3
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
# --- Astrocytes to Neurons
pairLR_use = c('SEMA6D_PLXNA1', 'SEMA5B_PLXNA1', 'NRXN2_NLGN2', 'NRXN1_NLGN2', 'NRXN2_NLGN1', 'NRXN1_NLGN1', 
               'NRXN2_CLSTN3', 'NRXN1_CLSTN3', 'NRXN2_CLSTN1', 'NRXN1_CLSTN1', 'CADM3_CADM1', 'CADM1_CADM1',
               'NRXN2_ADGRL1', 'NRXN1_ADGRL1')
pairLR_use = as.data.frame(unique(pairLR_use))
colnames(pairLR_use) = c('interaction_name')
# sort pairLR_use
pairLR_use$receptor = sapply(strsplit(pairLR_use$interaction_name, split='_'), '[', 2)
pairLR_use = pairLR_use[order(pairLR_use$receptor), ]
pairLR_use$receptor = NULL

bubble_data = netVisual_bubble(merged_obj, sources.use=c('Astrocytes_0', 'Astrocytes_1', 'Astrocytes_2', 'Astrocytes_3', 'Astrocytes_4', 'Astrocytes_5'), 
                               targets.use=c('Neurons'), 
                               comparison=c(1,2), angle.x=45, pairLR.use = pairLR_use, return.data=T)
bubble_comm_data = bubble_data$communication
plot_pairLR_bubble = function(df){
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
    color <- c('#99CCFF', '#FF6666')
    names(color) <- names(group[!duplicated(group)])
    color <- color[group]
    #names(color) <- dataset.name[comparison]
    dataset.name.order <- levels(df$source.target)
    dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
    dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
    xtick.color <- color[dataset.name.order]
    g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    return(g)
}
g = plot_pairLR_bubble(bubble_comm_data)
ggsave(file.path(joint_plot_dir, 'Astrocytes_to_Neurons_signaling_communication_probs_dotplot.pdf'), width=6, height=6)
export::graph2ppt(x=g, file=file.path(joint_plot_dir, 'Astrocytes_to_Neurons_signaling_communication_probs_dotplot.ppt'), width=6, height=6)

# --- Neurons to Astrocytes
pairLR_use = c('NRXN2_NLGN2', 'NRXN1_NLGN2', 'NRXN2_NLGN1', 'NRXN1_NLGN1', 'NRXN2_CLSTN1', 'NRXN1_CLSTN1',
               'EFNB3_EPHB1', 'CADM3_CADM1', 'CADM1_CADM1', 'TENM2_ADGRL3', 'TENM1_ADGRL3')
pairLR_use = as.data.frame(unique(pairLR_use))
colnames(pairLR_use) = c('interaction_name')
# sort pairLR_use
pairLR_use$receptor = sapply(strsplit(pairLR_use$interaction_name, split='_'), '[', 2)
pairLR_use = pairLR_use[order(pairLR_use$receptor), ]
pairLR_use$receptor = NULL

bubble_data = netVisual_bubble(merged_obj, targets.use=c('Astrocytes_0', 'Astrocytes_1', 'Astrocytes_2', 'Astrocytes_3', 'Astrocytes_4', 'Astrocytes_5'), 
                               sources.use=c('Neurons'), 
                               comparison=c(1,2), angle.x=45, pairLR.use = pairLR_use, return.data=T)
bubble_comm_data = bubble_data$communication
g = plot_pairLR_bubble(bubble_comm_data)
ggsave(file.path(joint_plot_dir, 'Neurons_to_Astrocytes_signaling_communication_probs_dotplot.pdf'), width=6, height=4.5)
export::graph2ppt(x=g, file=file.path(joint_plot_dir, 'Neurons_to_Astrocytes_signaling_communication_probs_dotplot.ppt'), width=6, height=4.5)

#' hierarchical plot for each pathways
pathways = c('EPHB', 'CADM', 'PCDH', 'ADGRL', 'NRXN', 'PTPR', 'NCAM', 'NOTCH',
             'CNTN', 'SEMA5', 'SEMA6')
hierachical_dir = file.path(global_dir, 'hierarchical_plots')
dir.create(hierachical_dir, showWarnings=F)
control_hierachical_dir = file.path(hierachical_dir, 'control')
dir.create(control_hierachical_dir, showWarnings=F)
setwd(control_hierachical_dir)
for (pathway in pathways) {
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    netVisual(control_cellchat_obj, signaling = pathway, layout = "hierarchy", 
              vertex.receiver = seq(1,4), out.format='svg', color.use=color.use, vertex.label.cex=1)
}
disease_hierachical_dir = file.path(hierachical_dir, 'disease')
dir.create(disease_hierachical_dir, showWarnings=F)
setwd(disease_hierachical_dir)
for (pathway in disease_cellchat_obj@netP$pathways) {
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    netVisual(disease_cellchat_obj, signaling = pathway, layout = "hierarchy", 
              vertex.receiver = seq(1,4), out.format='svg', color.use=color.use, vertex.label.cex=1)
}



