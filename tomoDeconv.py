#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import sys
import importlib
import scipy as sci

sys.path.append(".")


# In[2]:


import random


# In[5]:


#select replicate indices and save their mean
def removeReplicas(adata):
    #find genes with a '.' in their name
    noreplica = adata.var.loc[~adata.var.index.str.contains("\\."),:].index
    replicas = adata.var.loc[adata.var.index.str.contains("\\."),:].index
    uniques = set()
    new_var_names = set()
    for x in replicas: 
        uniques.add(x.split(".")[0])
        new_var_names.add(x.split(".")[0])
    new_var_names.update(noreplica.tolist())
    adata_new = sc.AnnData(np.zeros((len(adata.obs_names),len(new_var_names))))
    adata_new.var_names = list(new_var_names)
    adata_new.obs_names = adata.obs_names
    adata_new[:,noreplica].X = adata[:,noreplica].X
    adata_new.obs = adata.obs
    for x in uniques:
        ind = adata.var.loc[adata.var.index.str.contains(x),:].index
        if(len(ind)!=1):
            adata_new[:,x].X = np.mean(adata[:,ind].X,axis=1)
        else:
            adata_new[:,x].X =  adata[:,ind[0]].X
    return adata_new


# # Read Tomo Data

# In[6]:


tomoData = pd.read_csv('../../Data/Tomo1_DR11_count_table.csv',index_col=0) #Tomo1_count_table.csv


# In[7]:


tomoaData_normalized = sc.AnnData(tomoData.T)
tomoaData_normalized = removeReplicas(tomoaData_normalized)


# In[8]:


tomoaData_dataframe = pd.DataFrame(data=tomoaData_normalized.X.T,index=tomoaData_normalized.var_names,columns=tomoaData_normalized.obs_names)
tomoaData_dataframe.to_csv(r'../write/cleaned_ids.csv')


# In[9]:


tomoaData_dataframe.shape


# In[10]:


sc.pp.normalize_per_cell(tomoaData_normalized, counts_per_cell_after=1e4)
tomoData_norm = pd.DataFrame(data=tomoaData_normalized.X.T,columns = tomoaData_normalized.obs_names, index=tomoaData_normalized.var_names)


# In[11]:


import seaborn as sb
import matplotlib.pyplot as plt
counts = pd.DataFrame(data=tomoData.T.sum(1),index= tomoData.columns)
raw_counts = counts
counts = counts.T
counts = counts[counts>1000]
counts = counts.iloc[0].apply(pd.to_numeric)
counts = counts.dropna()
plt.hist(counts, color='blue')
plt.show()


# In[12]:


counts_proc = counts.drop(['74'], axis=0)
counts_proc.shape


# In[13]:


sc.set_figure_params(dpi=200, dpi_save=300, vector_friendly=True, frameon=False)


# In[14]:


plt.rcParams['figure.figsize'] = (6, 4)


# In[15]:


print(counts_proc.shape)
fig, ax = plt.subplots()
x_pos = [i for i, _ in enumerate(counts_proc.index)]
print(x_pos)
plt.xlabel('section')
plt.ylabel('# of counts')
ax.plot(x_pos,counts_proc)


# # Read single-cell data

# In[16]:


file_1 = '../../Data/all_h5_transfer/H5_Dr11_cr31_ffbm.h5'
file_2 = '../../Data/all_h5_transfer/H6_Dr11_cr31_ffbm.h5'
file_3 = '../../Data/all_h5_transfer/H7_Dr11_cr31_ffbm.h5'
file_4 = '../../Data/all_h5_transfer/H8a_Dr11_cr31_ffbm.h5'
file_5 = '../../Data/all_h5_transfer/H8v_Dr11_cr31_ffbm.h5'

adata_1 = sc.read_10x_h5(file_1)
adata_2 = sc.read_10x_h5(file_2)
adata_3 = sc.read_10x_h5(file_3)
adata_4 = sc.read_10x_h5(file_4)
adata_5 = sc.read_10x_h5(file_5)


# In[17]:


adata_1.var_names_make_unique()
adata_2.var_names_make_unique()
adata_3.var_names_make_unique()
adata_4.var_names_make_unique()
adata_5.var_names_make_unique()


# In[18]:


adata=adata_1.concatenate([adata_2,adata_3,adata_4,adata_5],
                                      batch_key='sample',
                                      batch_categories=['H5','H6','H7','H8a','H8v'])
adata.X = adata.X.todense()


# ## Single-cell data analysis

# In[19]:


#read annotations
annotations = pd.read_csv('../../Data/final_metadata_Tmacromerged_2.csv',index_col=0)


# In[20]:


indices = [x for x in annotations.index if not x.startswith('Hr')]
annotations = annotations.loc[indices,:]


# In[21]:


adata.obs_names = [str(adata.obs.loc[x,'sample'])+'_'+str(x.split('-', 1)[0]) for x in adata.obs_names]


# In[22]:


annotations.shape


# In[23]:


anno_drop = annotations.index.difference(adata.obs_names)
annotations = annotations.drop(anno_drop)
adata_proc = adata[annotations.index]
adata_proc


# In[25]:


for s in annotations.columns.values:
    adata_proc.obs[s] = annotations[s].tolist()


# In[26]:


adata_proc = removeReplicas(adata_proc)


# In[27]:


for s in annotations.columns.values:
    adata_proc.obs[s] = annotations[s].tolist()


# In[28]:


#select shared genes between sc and tomo data
adata_proc = adata_proc[:,adata_proc.var_names.intersection(tomoData_norm.index)]
tomoData_norm = tomoData_norm.loc[adata_proc.var_names.intersection(tomoData_norm.index),:]


# In[29]:


adata_norm = sc.pp.normalize_per_cell(adata_proc, counts_per_cell_after=1e4,copy=True)
adata_proc = sc.pp.log1p(adata_norm, copy=True)
sc.pp.highly_variable_genes(adata_proc, flavor='seurat', n_top_genes=5000)


# In[30]:


sc.tl.pca(adata_proc)
adata_proc.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat


# In[31]:


sc.pp.neighbors(adata_proc, n_neighbors=30)
sc.tl.umap(adata_proc)


# In[33]:


sc.pl.umap(adata_proc, color='Cell_type')


# # Deconvolution

# In[37]:


adata_norm.var['highly_variable'] = adata_proc.var['highly_variable']


# In[38]:


adata_norm_ag = adata_norm.copy()
adata_norm_ag = adata_norm_ag[adata_norm_ag.obs['Cell_type'].isin(['Epicardium (Atrium)','Epicardium (Ventricle)','Endocardium (Atrium)','Cardiomyocytes (Ventricle)','Endocardium (Ventricle)','Cardiomyocytes (Atrium)','Smooth muscle cells','Macrophages','Fibroblasts (const.)','T-cells'])]


# In[39]:


sc.pl.umap(adata_proc[adata_norm_ag.obs_names], color='Cell_type')


# In[40]:


import autogenes as ag
ngenes = 400 #400
ngen=5000
centroids = ag.init(adata_norm_ag,use_highly_variable=True,celltype_key='Cell_type')
ag.optimize(ngen=ngen,seed=0,nfeatures=ngenes,mode='fixed')


# In[41]:


ag.plot(index=0)


# In[42]:


pareto =20
selection = ag.select(index=pareto)

centroids_sc_pareto = pd.DataFrame(centroids[:,selection].X.T,index=centroids[:,selection].var_names,columns=centroids[:,selection].obs_names)


# In[43]:


import seaborn as sns
corr = pd.DataFrame(data = np.corrcoef(centroids_sc_pareto.T), columns = centroids_sc_pareto.columns, index = centroids_sc_pareto.columns)
mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    sns_plot =sns.clustermap(np.abs(corr),cmap=sns.color_palette("GnBu", 1000), robust=True)


# In[44]:


import seaborn as sns
subTypes = pd.DataFrame
subTypes = centroids_sc_pareto.columns
type_pal = sns.husl_palette(centroids_sc_pareto.columns.size, s=0.7)
lut = dict(zip(centroids_sc_pareto.columns.unique(), type_pal))
row_colors = subTypes.map(lut)
sns_plot = sns.clustermap(centroids_sc_pareto.T, cmap="mako", robust=True) #, cmap="mako", robust=True, row_cluster = False)
#sns_plot.figure.savefig(figdir+"heatmap_GA.png",dpi=300)


# In[45]:


coef_nnls = ag.deconvolve(tomoData_norm.T, model='nnls') #nusvr
coef_nnls.shape


# In[46]:


def normalize_proportions(data,copy):
    if copy==True:
        data_copy = data.copy()
    else:
        data_copy = data
    data_copy[data_copy < 0] = 0
    for raw in data_copy.index:
        sum = data_copy.loc[raw].sum()
        if sum!=0:
            data_copy.loc[raw] = np.divide(data_copy.loc[raw],sum)
        else:
            data_copy.loc[raw] = 0
    return data_copy


# In[47]:


proportions_nnls = pd.DataFrame(data=coef_nnls.T, index= centroids_sc_pareto.columns, columns = tomoData_norm.columns)
proportions_nnls_norm = normalize_proportions(proportions_nnls.T,copy=True).T


# In[48]:


proportions_nnls_norm=proportions_nnls_norm.astype(float)


# In[49]:


import seaborn as sns; sns.set()
ax = sns.heatmap(proportions_nnls_norm.astype(float),
                vmin = 0, vmax = 1, xticklabels=False)
ax_fig = ax.get_figure()#.savefig('Tomo1_deconvolution_highandepi_400feat_5000gen.pdf')


# In[68]:


#ax_fig.savefig('Tomo1_deconvolution_highandepi_400feat_5000gen.pdf', bbox_inches='tight')


# In[120]:


x_pos = [i for i, _ in enumerate(proportions_nnls_norm.T.index)]
for cell in proportions_nnls_norm.T.columns:
    plt.plot(x_pos,proportions_nnls_norm.T.loc[:,cell])
    plt.xlabel('section')
    plt.ylabel('proportion')
    plt.title(cell)
    plt.show()


# ## Validation 

# In[50]:


import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams
from matplotlib.colors import ColorConverter
import pandas as pd
from pandas import unique, isnull
from scipy.sparse import issparse
from scvelo.plotting.utils import is_categorical, interpret_colorkey, savefig_or_show


def heatmap(adata, var_names, tkey='pseudotime', xkey='Ms', color_map='magma', col_color=None, n_convolve=30,
            standard_scale=0, sort=True, colorbar=None, col_cluster=False, row_cluster=False, figsize=(10, 5),
            font_scale=None, show=True, save=None, ax=None, dpi=200, **kwargs):
    """    Plot time series for genes as heatmap.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    var_names: `str`,  list of `str`
        Names of variables to use for the plot.
    tkey: `str` (default: `'pseudotime'`)
        Observation key to extract time data from.
    xkey: `str` (default: `'Ms'`)
        Layer key to extract count data from.
    color_map: `str` (default: `'viridis'`)
        String denoting matplotlib color map.
    col_color: `str` or `None` (default: `None`)
        String denoting matplotlib color map to use along the columns.
    n_convolve: `int` or `None` (default: `30`)
        If `int` is given, data is smoothed by convolution along the x-axis with kernel size n_convolve.
    standard_scale : `int` or `None` (default: `0`)
        Either 0 (rows) or 1 (columns). Whether or not to standardize that dimension, meaning for each row or column,
        subtract the minimum and divide each by its maximum.
    sort: `bool` (default: `True`)
        Wether to sort the expression values given by xkey.
    colorbar: `bool` or `None` (default: `None`)
        Whether to show colorbar.
    {row,col}_cluster : bool, optional
        If True, cluster the {rows, columns}.
    figsize: tuple (default: `(7,5)`)
        Figure size.
    show: `bool`, optional (default: `None`)
        Show the plot, do not return axis.
    save: `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the default filename.
        Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
    ax: `matplotlib.Axes`, optional (default: `None`)
        A matplotlib axes object. Only works if plotting a single component.
    Returns
    -------
        If `show==False` a `matplotlib.Axis`
    """

    import seaborn as sns
    var_names = [name for name in var_names if name in adata.var_names]

    time = adata.obs[tkey].values
    time = time[np.isfinite(time)]

    df = pd.DataFrame(adata[:, var_names].layers[xkey][np.argsort(time)], columns=var_names)
    print(var_names)
    if n_convolve is not None:
        weights = np.ones(n_convolve) / n_convolve
        for i, gene in enumerate(var_names):
            try:
                df[gene] = np.convolve(df[gene].values, weights, mode='same')
            except: pass

    if sort:
        max_sort = np.argsort(np.argmax(df.values, axis=0))
        df = pd.DataFrame(df.values[:, max_sort], columns=df.columns[max_sort])
    if col_color is not None: 
        col_color_0 = interpret_colorkey(adata, col_color)[np.argsort(time)]
        #col_color_pdt = interpret_colorkey(adata, col_color[1])[np.argsort(time)]
    if font_scale is not None: sns.set(font_scale=font_scale)
    #print(interpret_colorkey)
    cm = sns.clustermap(df.T, col_colors=col_color_0, col_cluster=False, row_cluster=row_cluster, cmap=color_map,
                        xticklabels=False, standard_scale=standard_scale, figsize=figsize, **kwargs)
    if not colorbar: cm.cax.set_visible(False)
    savefig_or_show('heatmap', save=save, show=show,dpi=dpi)
    if not show: return cm
    pl.show()
    return cm


# In[51]:


adata_log_markers = adata_proc[:,centroids_sc_pareto.index].copy()
adata_log_markers = adata_log_markers[adata_norm_ag.obs_names]
sc.tl.rank_genes_groups(adata_log_markers, groupby='Cell_type', key_added='rank_subset')
mg = []
clusters = ['Cardiomyocytes (Ventricle)','Endocardium (Ventricle)','Epicardium (Ventricle)','Cardiomyocytes (Atrium)','Endocardium (Atrium)','Epicardium (Atrium)','Smooth muscle cells','Macrophages','Fibroblasts (const.)','T-cells']
markers_subset = pd.DataFrame(columns = clusters)
list_mg = []
for clust in clusters:
    mg=adata_log_markers.uns['rank_subset']['names'][clust][:10]
    print(clust, mg)
    list_mg.extend(mg)
    markers_subset.loc[:,clust] = mg


# In[52]:


markers = []
for clust in clusters:
    markers = list(markers + list(adata_log_markers.uns['rank_subset']['names'][clust][:10]))


# In[53]:


markers_epicardium = []
for clust in ['Epicardium (Ventricle)','Epicardium (Atrium)']:
    markers_epicardium = list(markers_epicardium + list(adata_log_markers.uns['rank_subset']['names'][clust][:10]))


# In[54]:


adata_log_markers = adata_proc[:,list(dict.fromkeys(markers))].copy()
adata_log_markers = adata_log_markers[adata_norm_ag.obs_names]
#adata_log_markers = adata_log_markers[adata_log_markers.obs['louvain_subset'].isin(['Apical_non','Apical_polyp','Basal_non','Basal_polyp','Ciliated','Glandular'])]
adata_log_markers = adata_log_markers[adata_log_markers.obs.loc[adata_log_markers.obs["Cell_type"].str.lower().sort_values().index].index]
adata_log_markers.obs['order'] = np.linspace(1,len(adata_log_markers.obs_names),len(adata_log_markers.obs_names)).astype('int')
adata_log_markers.layers['counts'] = adata_log_markers.X


# In[55]:


print(clusters)


# In[56]:


sc.pl.umap(adata_proc, color='Cell_type')


# In[57]:


clusters


# In[182]:


heatmap(adata_log_markers,
        var_names=adata_log_markers.var_names,
        xkey='counts', tkey='order', n_convolve=10,sort=False,
        col_color='Cell_type',yticklabels=1, font_scale=.5, figsize=(3, 8))#,save="_all_sc_pareto.png")


# In[193]:


#Epicardium markers are only expressed in Atrium which validates the results
tomoData_markers = tomoData_norm.loc[markers_epicardium,:]
ax = sns.clustermap(tomoData_markers.astype(float),yticklabels=1,figsize=(3, 4), col_cluster=False, row_cluster=False)

