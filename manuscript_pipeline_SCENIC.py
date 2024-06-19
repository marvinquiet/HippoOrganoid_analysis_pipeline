# conda activate Seurat
# export HDF5_USE_FILE_LOCKING=FALSE


import scanpy as sc
import os
import glob
import anndata
import numpy as np
import pandas as pd
import loompy as lp

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

project_dir = "/panfs/compbio/users/wma36/collaborations/JieXu"
result_dir = os.path.join(project_dir, 'Manuscript_analysis_figures')
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
final_fig_dir = os.path.join(result_dir, 'Manuscript_final_figures')
if not os.path.exists(final_fig_dir):
    os.makedirs(final_fig_dir)

adata_filepath = os.path.join(result_dir, 'annotated_Seurat_obj.h5ad')
aux_data_dir = os.path.join(project_dir, 'SCENIC_output', 'auxiliary_data')
f_tfs = aux_data_dir+os.sep+'allTFs_hg38.txt'
result_dir = result_dir+os.sep+'SCENIC_output'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# path to unfiltered loom file (this will be created in the optional steps below)
# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = result_dir+os.sep+"filtered_scenic.loom"
# path to pyscenic output
f_pyscenic_output = result_dir+os.sep+"pyscenic_output.loom"

adata = anndata.read_h5ad(adata_filepath)
# combine neurons
replace_neurons = {
        'Immature neurons': 'Neurons',
        'Inhibitory neurons': 'Neurons',
        'Excitatory neurons': 'Neurons'
        }
adata.obs['combined_celltype'] = adata.obs['annotated_celltype'].replace(replace_neurons)

## get highly variable genes
new_adata = anndata.AnnData(X=adata.raw.X, obs=adata.obs, var=adata.raw.var)
sc.pp.highly_variable_genes(new_adata, min_mean=0.0125, max_mean=5, min_disp=0.5) # their cutoff
new_adata = new_adata[:, new_adata.var['highly_variable']]
sc.pp.scale(new_adata, max_value=10)

# combine condition
row_attrs = { 
    "Gene": np.array(new_adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(new_adata.obs.index) ,
    "nGene": np.array( np.sum(abs(new_adata.X.transpose()) > 0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(new_adata.X.transpose(), axis=0)).flatten() ,
}
lp.create(f_loom_path_scenic, new_adata.X.transpose(), row_attrs, col_attrs )
#!pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 20
#pyscenic grn /projects/compbio/users/wma36/collaborations/JieXu/SCENIC_output/filtered_scenic.loom /projects/compbio/users/wma36/collaborations/JieXu/SCENIC_output/auxiliary_data/allTFs_hg38.txt -o /projects/compbio/users/wma36/collaborations/JieXu/SCENIC_output/adj.csv --num_workers 1

adjacencies = pd.read_csv(os.path.join(result_dir, "adj.csv"), index_col=False, sep=',')
adjacencies.head()

f_db_glob = os.path.join(aux_data_dir, '*feather')
f_db_names = ' '.join( glob.glob(f_db_glob) )
f_motif_path = os.path.join(aux_data_dir, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

#!pyscenic ctx adj.tsv \
#   {f_db_names} \
#   --annotations_fname {f_motif_path} \
#   --expression_mtx_fname {f_loom_path_scenic} \
#   --output reg.csv \
#   --mask_dropouts \
#   --num_workers 20

## --- Cellular enrichment
nGenesDetectedPerCell = np.sum(new_adata.X>0, axis=1)
nGenesDetectedPerCell = pd.Series(nGenesDetectedPerCell)
percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])
print(percentiles)
fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
sns.distplot(nGenesDetectedPerCell, norm_hist=False, kde=False, bins='fd')
for i,x in enumerate(percentiles):
    fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
    ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i]*100}%)', color='red', rotation=30, size='x-small',rotation_mode='anchor' )
ax.set_xlabel('# of genes')
ax.set_ylabel('# of cells')
fig.tight_layout()
plt.savefig(os.path.join(result_dir, 'gene_percentiles.png'))
#pyscenic aucell \
#    {f_loom_path_scenic} \
#    reg.csv \
#    --output {f_pyscenic_output} \
#    --num_workers 20
import json
import zlib
import base64
# collect SCENIC AUCell output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
#import umap
## UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( os.path.join(result_dir, "scenic_umap.txt"), sep='\t')
# tSNE
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( os.pah.join(project_dir, 'SCENIC_output', "scenic_tsne.txt"), sep='\t')
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
sig = load_signatures(os.path.join(result_dir, 'reg.csv'))
new_adata = add_scenic_metadata(new_adata, auc_mtx, sig)
from pyscenic.utils import load_motifs
import operator as op
from IPython.display import HTML, display
BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"
def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', 200)
    data = HTML(df.to_html(escape=False))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)
    return data
df_motifs = load_motifs(os.path.join(result_dir, 'reg.csv'))
selected_motifs = df_motifs.index.get_level_values('TF').tolist()
df_motifs_sel = df_motifs.iloc[ [ True if x in selected_motifs else False for x in df_motifs.index.get_level_values('TF') ] ,:]
html_obj = display_logos( df_motifs_sel.sort_values([('Enrichment','NES')], ascending=False))
with open(os.path.join(result_dir, 'motifs.html'), 'w') as f:
    f.write(html_obj.data)


## --- binarize AUC matrix
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize

# scenic output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
## binarize auc mtx
binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers=4 )
binary_mtx.head()
# use pre-set thresholds
AUC_threshold_dir = os.path.join(result_dir, 'AUC_thresholds')
if not os.path.exists(AUC_threshold_dir):
    os.makedirs(AUC_threshold_dir)
auc_thresholds = pd.read_csv(AUC_threshold_dir+os.sep+'AUC_thresholds.csv', index_col=0, header=0)
manual_thresholds = auc_thresholds.to_dict()['0']
binary_mtx, auc_thresholds = binarize( auc_mtx, threshold_overides=manual_thresholds, num_workers=4, seed=2024 )
binary_mtx.head()

# show AUC distributions for all regulons
AUC_threshold_dir = os.path.join(result_dir, 'AUC_thresholds')
if not os.path.exists(AUC_threshold_dir):
    os.makedirs(AUC_threshold_dir)
auc_thresholds.to_csv(AUC_threshold_dir+os.sep+'AUC_thresholds.csv')
r = auc_mtx.columns.tolist()
for r in auc_mtx.columns:
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150, sharey=False)
    sns.distplot(auc_mtx[ r ], ax=ax, norm_hist=True, bins=100)
    ax.plot( [ auc_thresholds[ r ] ]*2, ax.get_ylim(), 'r:')
    ax.title.set_text( r )
    ax.set_xlabel('')
    fig.tight_layout()
    fig.savefig(os.path.join(AUC_threshold_dir, r+'.pdf'), dpi=600, bbox_inches='tight')

# plot binarized heatmap
cats = sorted(list(set(adata.obs['combined_celltype'])))
top_reg = binary_mtx.columns.tolist()
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f
colors = ['#F8766D', '#D89000', '#39B600', '#A3A500', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC']
colorsd = dict( zip( cats, colors ))
colormap = [ colorsd[x] for x in adata.obs['combined_celltype'] ]
sns.set()
sns.set(font_scale=0.4)
fig = palplot( colors, cats, size=1.0)
plt.savefig(os.path.join(result_dir, "heatmap-legend-allRegs.pdf"), dpi=600, bbox_inches = "tight")
sns.set(font_scale=1.2)
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('Custom', ((1,1,1), (0,0,0)), 2)
import random
random.seed(2024)
g = sns.clustermap(binary_mtx, annot=False, square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=0, vmax=1, row_colors=colormap,
    cmap=cmap, figsize=(21, 16), cbar_pos = None) # remove colorbar
#g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig(os.path.join(result_dir, "heatmap-allRegs.pdf"), dpi=600, bbox_inches = "tight")

# --- obtain regulons regulated genes
from pyscenic.utils import modules_from_adjacencies
regulon_dir = result_dir+os.sep+'regulons'
if not os.path.exists(regulon_dir):
    os.makedirs(regulon_dir)
modules = list(modules_from_adjacencies(adjacencies, exprMat.T))
# create a dictionary of regulons:
regulons = {}
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
    regulons[i] =  list(r[r==1].index.values)
lf.close()
for tf in binary_mtx.columns:
    tf = tf.replace('(+)', '')
    tf_dir = regulon_dir+os.sep+tf
    os.makedirs(tf_dir, exist_ok=True)
    tf_mods = [ x for x in modules if x.transcription_factor==tf ]
    for i,mod in enumerate( tf_mods ):
        print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )
    print( f'{tf} regulon: {len(regulons[tf+"(+)"])} genes' )
    for i,mod in enumerate( tf_mods ):
        with open( tf_dir+os.sep+tf+'_module_'+str(i)+'.txt', 'w') as f:
            for item in mod.genes:
                f.write("%s\n" % item)
    with open( tf_dir+os.sep+tf+'_regulon.txt', 'w') as f:
        for item in regulons[tf+'(+)']:
            f.write("%s\n" % item)

# --- extract astrocytes progenitor cells
from pyscenic.export import export2loom, add_scenic_metadata
import operator as op
from cytoolz import compose
from pyscenic.transform import df2regulons
from pyscenic.utils import load_motifs
def derive_regulons(motifs, db_names=('hg38_500bp_up_100bp_down', 
                                 'hg38_10kbp_up_10kbp_down')):
    motifs.columns = motifs.columns.droplevel(0)
    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f
    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        # np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]
    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0) 
                                                                      & ((motifs['Annotation'] == 'gene is directly annotated')
                                                                        | (motifs['Annotation'].str.startswith('gene is orthologous to')
                                                                           & motifs['Annotation'].str.endswith('which is directly annotated for motif')))
                                                                     ])))
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))
df_motifs = load_motifs(os.path.join(result_dir, 'reg.csv'))
regulons = derive_regulons(df_motifs)
add_scenic_metadata(adata, auc_mtx, regulons)
sc.tl.tsne(adata, use_rep='X_aucell')
adata.write_h5ad(result_dir+os.sep+'AUCell_adata.h5ad')

adata = anndata.read_h5ad(result_dir+os.sep+'AUCell_adata.h5ad')
sc.pp.neighbors(adata, use_rep='X_aucell')
sc.tl.louvain(adata, resolution = 3)
sc.pl.tsne(adata, color='louvain')
plt.savefig(result_dir+os.sep+'AUC_tsne-resolution1.png')
embedding_aucell_tsne = pd.DataFrame(adata.obsm['X_tsne'], columns=[['_X', '_Y']], index=adata.obs_names)
sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
# split the figure to individual
plot_images = ['combined_celltype', 'louvain', 'condition']
for per_image in plot_images:
    sc.pl.tsne(adata, color=per_image, ncols=1, use_raw=False, cmap='magma')
    plt.savefig(result_dir+os.sep+per_image+'_AUC_tsne-regulonmarkers.png')

sc.set_figure_params(frameon=False, dpi=150, fontsize=12, figsize=(3.5, 3))
plot_images = ['Regulon(FOS(+))', 'Regulon(JUNB(+))', 'Regulon(NR2F1(+))',
        'Regulon(CEBPD(+))', 'Regulon(JUN(+))', 'Regulon(NKX6-2(+))', 'Regulon(POU2F1(+))',
        'Regulon(EGR1(+))', 'Regulon(SOX9(+))', 'Regulon(PAX6(+))', 'Regulon(HES5(+))', 'Regulon(FOSB(+))',
        'Regulon(ATF3(+))', 'Regulon(CEBPB(+))', 'Regulon(DDIT3(+))']
def custom_formatter(x, pos):
    return f'{x:.2f}'  # Format the numbers to 2 decimal places
for per_image in plot_images:
    ax = sc.pl.tsne(adata, color=per_image, ncols=1, use_raw=False, cmap='magma', show=False)
    colorbar = ax.collections[-1].colorbar
    colorbar.formatter = FuncFormatter(custom_formatter)
    colorbar.update_ticks()  # Update the colorbar with the new formatter
    plt.tight_layout()
    plt.savefig(result_dir+os.sep+per_image+'_AUC_tsne-regulonmarkers.png')
adata_AUCell_cluster.to_csv(result_dir+os.sep+'AUCell_cluster_assignment.csv')

## zoom into astrocytes
astrocytes_adata = adata[adata.obs['combined_celltype'] == 'Astrocytes'].copy()
sc.pp.neighbors(astrocytes_adata, use_rep='X_aucell')
sc.tl.louvain(astrocytes_adata, resolution = 2)
sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
sc.pl.tsne(astrocytes_adata, color=['combined_celltype', 'louvain', 'condition', 'Regulon(FOS(+))', 'Regulon(JUNB(+))', 'Regulon(NR2F1(+))', 
                         'Regulon(CEBPD(+))', 'Regulon(JUN(+))', 'Regulon(NKX6-2(+))', 'Regulon(POU2F1(+))',
			 'Regulon(EGR1(+))', 'Regulon(SOX9(+))', 'Regulon(PAX6(+))', 'Regulon(HES5(+))',
			 'Regulon(FOSB(+))'], ncols=4, use_raw=False)
plt.savefig(result_dir+os.sep+'Astrocytes_AUC_tsne-regulonmarkers.png')
astrocytes_AUCell_cluster = astrocytes_adata.obs[['condition', 'louvain']]
astrocytes_AUCell_cluster.to_csv(result_dir+os.sep+'Astrocytes_AUCell_cluster_assignment.csv')
# obtain regulon AUC
df_obs = astrocytes_adata.obs[astrocytes_adata.obs['louvain'].isin(['2', '8', '11'])]
#df_obs['louvain_condition'] = df_obs['louvain'].astype(str)+'-'+df_obs['condition'].astype(str)
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
#df_scores = df_obs[signature_column_names + ['louvain_condition']]
df_scores = df_obs[signature_column_names + ['condition']]
#df_results = ((df_scores.groupby(by='louvain_condition').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results = ((df_scores.groupby(by='condition').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_heatmap = pd.pivot_table(data=df_results.sort_values('Z', ascending=False),
                                   index='condition', columns='regulon', values='Z')
#df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
fig, ax1 = plt.subplots(1, 1, figsize=(12, 6))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=.7, cbar=False, square=True, linecolor='gray', 
                    cmap="YlGnBu", annot_kws={"size": 6})
ax1.set_ylabel('')
plt.savefig(result_dir+os.sep+'Astrocytes_AUC-Zscore_heatmap_cluster2+8+11.png')

## obtain astrocytes progenitors from all cells
astrocytes_progenitors_adata = adata[adata.obs['louvain']=='3']
sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
sc.pl.tsne(astrocytes_progenitors_adata, color=['combined_celltype', 'louvain', 'condition', 'Regulon(FOS(+))', 'Regulon(JUNB(+))', 'Regulon(NR2F1(+))', 
                         'Regulon(CEBPD(+))', 'Regulon(JUN(+))', 'Regulon(NKX6-2(+))', 'Regulon(POU2F1(+))',
			 'Regulon(EGR1(+))', 'Regulon(SOX9(+))', 'Regulon(PAX6(+))', 'Regulon(HES5(+))',
			 'Regulon(FOSB(+))'], ncols=4, use_raw=False)
plt.savefig(result_dir+os.sep+'AUC_tsne-regulonmarkers.png')


# obtain regulon AUC
df_obs = astrocytes_progenitors_adata.obs
#df_obs['louvain_condition'] = df_obs['louvain'].astype(str)+'-'+df_obs['condition'].astype(str)
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
#df_scores = df_obs[signature_column_names + ['louvain_condition']]
df_scores = df_obs[signature_column_names + ['condition']]
#df_results = ((df_scores.groupby(by='louvain_condition').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results = ((df_scores.groupby(by='condition').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= 0.5].sort_values('Z', ascending=False),
                                   #index='louvain_condition', columns='regulon', values='Z')
                                   index='condition', columns='regulon', values='Z')
#df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
fig, ax1 = plt.subplots(1, 1, figsize=(10, 6))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=.7, cbar=False, square=True, linecolor='gray', 
                    cmap="YlGnBu", annot_kws={"size": 6})
ax1.set_ylabel('')
plt.savefig(result_dir+os.sep+'Astrocytes_progenitors_AUC-Zscore_thres0.5_heatmap_cluster3.png')

df_heatmap = pd.pivot_table(data=df_results.sort_values('Z', ascending=False),
                                   #index='louvain_condition', columns='regulon', values='Z')
                                   index='condition', columns='regulon', values='Z')
#df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
fig, ax1 = plt.subplots(1, 1, figsize=(16, 6))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=.7, cbar=False, square=True, linecolor='gray', 
                    cmap="YlGnBu", annot_kws={"size": 6})
ax1.set_ylabel('')
plt.savefig(result_dir+os.sep+'Astrocytes_progenitors_AUC-Zscore_heatmap_cluster3.png')

