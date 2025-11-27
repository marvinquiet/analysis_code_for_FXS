# conda activate Seurat

import os
import glob
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

project_dir = "/panfs/compbio/users/wma36/collaborations/Yunhee/FMRpolyG_Cortex_10Xmultiomics"
adata_filepath = os.path.join(project_dir, 'Cortex_multiome_Seurat_pFC_200_4000_integration', 'Seurat_integration_annotated_object.h5ad')
result_dir = project_dir+os.sep+'SCENIC_output'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

aux_data_dir = os.path.join(project_dir, 'SCENIC_output', 'auxiliary_data')
f_tfs = aux_data_dir+os.sep+'allTFs_mm.txt'

# path to unfiltered loom file (this will be created in the optional steps below)
# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = result_dir+os.sep+"filtered_scenic.loom"
# path to pyscenic output
f_pyscenic_output = result_dir+os.sep+"pyscenic_output.loom"

adata = anndata.read_h5ad(adata_filepath)
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

#f_db_glob = os.path.join(aux_data_dir, '*feather')
#f_db_names = ' '.join( glob.glob(f_db_glob) )
#f_motif_path = os.path.join(aux_data_dir, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

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
import umap
# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( os.path.join(result_dir, "scenic_umap.txt"), sep='\t')
## tSNE
#tsne = TSNE( n_jobs=20 )
#dr_tsne = tsne.fit_transform( auc_mtx )
#pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( os.pah.join(project_dir, 'SCENIC_output', "scenic_tsne.txt"), sep='\t')
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
# binarize auc mtx
binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers=1 )
binary_mtx.head()

AUC_threshold_dir = os.path.join(result_dir, 'AUC_thresholds')
if not os.path.exists(AUC_threshold_dir):
    os.makedirs(AUC_threshold_dir)
auc_thresholds.to_csv(AUC_threshold_dir+os.sep+'AUC_thresholds.csv')
r = auc_mtx.columns.tolist()
# show AUC distributions for all regulons
for r in auc_mtx.columns:
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150, sharey=False)
    sns.distplot(auc_mtx[ r ], ax=ax, norm_hist=True, bins=100)
    ax.plot( [ auc_thresholds[ r ] ]*2, ax.get_ylim(), 'r:')
    ax.title.set_text( r )
    ax.set_xlabel('')
    fig.tight_layout()
    fig.savefig(os.path.join(AUC_threshold_dir, r+'.pdf'), dpi=600, bbox_inches='tight')
# plot binarized heatmap
cats = sorted(list(set(adata.obs['annotated_celltype'])))
#top_reg = binary_mtx.columns.tolist()
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
colors = sns.color_palette('bright',n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap = [ colorsd[x] for x in adata.obs['annotated_celltype'] ]
sns.set()
sns.set(font_scale=0.4)
fig = palplot( colors, cats, size=1.0)
plt.savefig(os.path.join(result_dir, "heatmap-legend-allRegs.pdf"), dpi=600, bbox_inches = "tight")
sns.set(font_scale=1.2)
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('Custom', ((1,1,1), (0,0,0)), 2)
g = sns.clustermap(binary_mtx, annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=0, vmax=1, row_colors=colormap,
    cmap=cmap, figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig(os.path.join(result_dir, "heatmap-allRegs.pdf"), dpi=600, bbox_inches = "tight")



## split auc_mtx
control_cells = new_adata.obs_names[new_adata.obs['condition'] == 'control'].tolist()
disease_cells = new_adata.obs_names[new_adata.obs['condition'] == 'disease'].tolist()
control_auc_mtx = auc_mtx.loc[control_cells, :]
disease_auc_mtx = auc_mtx.loc[disease_cells, :]
    
control_rss_cellType = regulon_specificity_scores( control_auc_mtx, new_adata.obs.loc[control_cells, 'annotated_celltype'] )
disease_rss_cellType = regulon_specificity_scores( disease_auc_mtx, new_adata.obs.loc[disease_cells, 'annotated_celltype'] )
cats = sorted(list(set(new_adata.obs['annotated_celltype'])))
fig = plt.figure(figsize=(15, 8))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=control_rss_cellType.T[c]
    ax = fig.add_subplot(2,4,num)
    plot_rss(control_rss_cellType, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig(os.path.join(result_dir, 'control_cellType-RSS-top5.pdf'), dpi=600, bbox_inches = "tight")
control_rss_cellType.to_csv(os.path.join(result_dir, 'control_cellType-RSS.csv'))

fig = plt.figure(figsize=(15, 8))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=disease_rss_cellType.T[c]
    ax = fig.add_subplot(2,4,num)
    plot_rss(disease_rss_cellType, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig(os.path.join(result_dir, 'disease_cellType-RSS-top5.pdf'), dpi=600, bbox_inches = "tight")
disease_rss_cellType.to_csv(os.path.join(result_dir, 'disease_cellType-RSS.csv'))

binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers=4 )
binary_mtx.head()
# focus on SOX9 and get regulated genes
fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150, sharey=False)
sns.distplot(auc_mtx[ 'SOX9(+)' ], ax=ax, norm_hist=True, bins=100)
ax.plot( [ auc_thresholds[ 'SOX9(+)' ] ]*2, ax.get_ylim(), 'r:')
ax.title.set_text( 'SOX9(+)' )
ax.set_xlabel('')
fig.tight_layout()
fig.savefig(os.path.join(result_dir, 'SOX9_AUC_distributions.pdf'), dpi=600, bbox_inches='tight')

# modules from network inference output
from pyscenic.utils import modules_from_adjacencies
modules = list(modules_from_adjacencies(adjacencies, exprMat.T))
# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
    regulons[i] =  list(r[r==1].index.values)

tf = 'SOX9'
tf_mods = [ x for x in modules if x.transcription_factor==tf ]
for i,mod in enumerate( tf_mods ):
    print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )
print( f'{tf} regulon: {len(regulons[tf+"(+)"])} genes' )

for i,mod in enumerate( tf_mods ):
    with open( result_dir+os.sep+tf+'_module_'+str(i)+'.txt', 'w') as f:
        for item in mod.genes:
            f.write("%s\n" % item)
            
with open( result_dir+os.sep+tf+'_regulon.txt', 'w') as f:
    for item in regulons[tf+'(+)']:
        f.write("%s\n" % item)


    topreg = []
    for i,c in enumerate(cats):
        topreg.extend(
            list(rss_cellType.T[c].sort_values(ascending=False)[:5].index)
        )
    topreg = list(set(topreg))
    auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
    for col in list(auc_mtx.columns):
        auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
    #auc_mtx_Z.sort_index(inplace=True)
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
    colors = sns.color_palette('bright',n_colors=len(cats) )
    colorsd = dict( zip( cats, colors ))
    colormap = [ colorsd[x] for x in adata.obs['annotated_celltype'] ]
    sns.set()
    sns.set(font_scale=0.4)
    fig = palplot( colors, cats, size=1.0)
    plt.savefig(os.path.join(result_dir, "_cellType-heatmap-legend-top5.pdf"), dpi=600, bbox_inches = "tight")
    sns.set(font_scale=1.2)
    g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
        yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
        cmap="YlGnBu", figsize=(21,16) )
    g.cax.set_visible(True)
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.set_xlabel('')
    plt.savefig(os.path.join(result_dir, "_cellType-heatmap-top5.pdf"), dpi=600, bbox_inches = "tight")


## === @TODO:

binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers=1 )
binary_mtx.head()



# select regulons:
r = [ 'RXRA_(+)', 'NFE2_(+)', 'ETV2_(+)' ]

fig, axs = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)
for i,ax in enumerate(axs):
    sns.distplot(auc_mtx[ r[i] ], ax=ax, norm_hist=True, bins=100)
    ax.plot( [ auc_thresholds[ r[i] ] ]*2, ax.get_ylim(), 'r:')
    ax.title.set_text( r[i] )
    ax.set_xlabel('')
    
fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='large')
fig.text(0.5, -0.01, 'AUC', ha='center', va='center', rotation='horizontal', size='large')

fig.tight_layout()
fig.savefig(os.path.join(project_dir, 'SCENIC_output', 'PBMC10k_cellType-binaryPlot2.pdf'), dpi=600, bbox_inches='tight')


from pyscenic.utils import modules_from_adjacencies
lf = lp.connect( f_pyscenic_output, mode='r', validate=False )
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
    regulons[i] =  list(r[r==1].index.values)
lf.close()
modules = list(modules_from_adjacencies(adjacencies, exprMat))


tf = 'EBF1'
tf_mods = [ x for x in modules if x.transcription_factor==tf ]

for i,mod in enumerate( tf_mods ):
    print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )
print( f'{tf} regulon: {len(regulons[tf+"_(+)"])} genes' )


for i,mod in enumerate( tf_mods ):
    with open( os.path.join(project_dir, 'SCENIC_output', tf+'_module_'+str(i)+'.txt'), 'w') as f:
        for item in mod.genes:
            f.write("%s\n" % item)
            
with open( project_dir, 'SCENIC_output', tf+'_regulon.txt'), 'w') as f:
    for item in regulons[tf+'_(+)']:
        f.write("%s\n" % item)
