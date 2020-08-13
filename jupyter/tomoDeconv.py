#!/usr/bin/env python
# coding: utf-8

# # Dependencies

# In[1]:


import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import sys
import importlib
import scipy as sci
import seaborn as sb
import matplotlib.pyplot as plt

sys.path.append(".")

import pickle
from autogenes import AutoGenes


# In[2]:


import random


# In[3]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '1')
get_ipython().run_line_magic('aimport', 'autogenes')


# In[4]:


importlib.reload(autogenes)


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


# In[6]:


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


# In[7]:


mito_genes = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/mito.genes.vs.txt',
                        header = None,
                        names = ['Remove_1', 'Remove_2', 'Gene'])
mito_genes = mito_genes.drop(columns = ['Remove_1', 'Remove_2'])


# In[8]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors_2.csv', 
                               index_col = 0)
cell_type_colors.loc['Erythrocytes'] = '#681b0e'


# In[9]:


#read annotations
annotations_nonery = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv', 
                                 index_col = 0)
annotations_ery = pd.read_csv('../../Data/final_erythrocytes.csv', index_col=0)
annotations_ery.columns = ['orig.ident', 'Cell_type']
annotations = pd.concat([annotations_nonery[['orig.ident', 'Cell_type']], annotations_ery])
indices = [x for x in annotations.index if not x.startswith('Hr')]
annotations = annotations.loc[indices,:]
annotations = annotations[~annotations.index.duplicated()]


# In[10]:


#annotations.Cell_type.value_counts()


# In[11]:


highcount_celltypes = annotations.Cell_type.value_counts()[annotations.Cell_type.value_counts() > 1000]
select_celltypes = highcount_celltypes.index.tolist() #+ ['Epicardium (Atrium)', 'Epicardium (Ventricle)']
annotations = annotations[annotations.Cell_type.isin(select_celltypes)]


# # Read tomo data, remove mito reads and duplicates

# In[12]:


tomoData = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Tomo1_DR11_count_table.csv',index_col=0)


# In[13]:


#tomoData = tomoData.drop(index = list(set(mito_genes.Gene.tolist()) & set(tomoData.index.tolist())))


# In[14]:


tomoaData_noreplicates = sc.AnnData(tomoData.T)
tomoaData_noreplicates = removeReplicas(tomoaData_noreplicates)


# # Read sc data, remove mito reads and duplicates

# In[15]:


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


# In[16]:


adata_1.var_names_make_unique()
adata_2.var_names_make_unique()
adata_3.var_names_make_unique()
adata_4.var_names_make_unique()
adata_5.var_names_make_unique()


# In[17]:


adata=adata_1.concatenate([adata_2,adata_3,adata_4,adata_5],
                                      batch_key='sample',
                                      batch_categories=['H5','H6','H7','H8a','H8v'])
adata.X = adata.X.todense()


# In[18]:


#non_mito_genes = [name for name in adata.var_names if name not in mito_genes.Gene.tolist()]
#adata = adata[:, non_mito_genes]


# In[19]:


adata.obs_names = [str(adata.obs.loc[x,'sample'])+'_'+str(x.split('-', 1)[0]) for x in adata.obs_names]


# In[20]:


anno_drop = annotations.index.difference(adata.obs_names)
annotations = annotations.drop(anno_drop)
adata = adata[annotations.index]
for s in annotations.columns.values:
    adata.obs[s] = annotations[s].tolist()


# In[21]:


adata = removeReplicas(adata)


# # Intersect genes in tomo and single cell datasets

# In[22]:


adata = adata[:,adata.var_names.intersection(tomoaData_noreplicates.var_names)]
tomoaData_noreplicates = tomoaData_noreplicates[:,adata.var_names.intersection(tomoaData_noreplicates.var_names)]


# # Merge and normalize tomo-sections; determine highly-variable genes

# In[23]:


counts = pd.DataFrame(data=tomoData.T.sum(1),index= tomoData.columns)
#raw_counts = counts
#counts = counts.T
#counts = counts[counts>1000]
#counts = counts.iloc[0].apply(pd.to_numeric)
#counts = counts.dropna()
#plt.hist(counts, color='blue')
#plt.show()
#counts[41:60]


# In[24]:


sc.set_figure_params(dpi=200, dpi_save=300, vector_friendly=True, frameon=False)
#print(counts_proc.shape)
plt.rcParams['figure.figsize'] = (6, 4)
fig, ax = plt.subplots()
x_pos = [i for i, _ in enumerate(counts.index)]
#print(x_pos)
plt.xlabel('section')
plt.ylabel('# of counts')
ax.plot(x_pos,counts)


# In[25]:


# Add up low-read columns from tomoData for tomosectionmerge: if reads of column in tomosectionmerge
# sum up to < n (default 20k), add next tomoData column to current tomosectionmerge column; 
# if not, make new tomosectionmerge column from tomoData column.
tomo_norep = pd.DataFrame(data = tomoaData_noreplicates.X.T,
                         columns = tomoaData_noreplicates.obs_names,
                         index = tomoaData_noreplicates.var_names)
m = 0
tomosectionmerge = tomo_norep.drop(tomo_norep.columns, axis = 1)
sectionmergewidth = pd.DataFrame(columns=['width'])
for c in tomo_norep.columns:
    if(c==tomo_norep.columns[0]):
        tomosectionmerge[m] = tomo_norep[c]
        sectionmergewidth.loc[m] = 1
    elif(tomosectionmerge[m].sum(axis=0) < 20000):
        tomosectionmerge[m] = tomosectionmerge[m] + tomo_norep[c]
        sectionmergewidth.loc[m] +=1
    else:
        m = m + 1
        tomosectionmerge[m] = tomo_norep[c]
        sectionmergewidth.loc[m] = 1


# In[26]:


tomoaData_normalized = sc.AnnData(tomosectionmerge.T)
sc.pp.normalize_per_cell(tomoaData_normalized, counts_per_cell_after=1e4)
tomoData_norm = pd.DataFrame(data=tomoaData_normalized.X.T,
                             columns = tomoaData_normalized.obs_names, 
                             index=tomoaData_normalized.var_names)


# # Process single cell data to clusters; determine highly-variable genes

# In[27]:


adata_norm = sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4,copy=True)
adata_proc = sc.pp.log1p(adata_norm, copy=True)
sc.pp.highly_variable_genes(adata_proc, flavor='seurat', n_top_genes=4000)
#sc.pp.scale(adata_proc)


# In[28]:


sc.tl.pca(adata_proc)
adata_proc.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
#sc.pl.pca_scatter(adata_proc, save="_PCA.png" )


# In[29]:


sc.pp.neighbors(adata_proc, n_neighbors=30)
sc.tl.umap(adata_proc)


# In[30]:


sc.pl.umap(adata_proc, color='Cell_type')#,
#           palette = cell_type_colors.loc[adata_proc.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[ ]:


# Plot gene expression
gene = 'hgh1'
expression_plot = tomoData_norm.loc[gene]
x_pos = [i for i, _ in enumerate(expression_plot.index)]
fig, ax = plt.subplots()
plt.xlabel('section')
plt.ylabel('Normalized expression')
ax.plot(x_pos,expression_plot)
plt.title(gene)
#plt.savefig('tomo1_normalized_expression_' + gene + '.png')


# # Calculate centroids

# In[31]:


adata_norm.obs = adata_proc.obs


# In[32]:


clusters = list(set(adata_norm.obs['Cell_type']))
sc_mean = pd.DataFrame(index=adata_norm.var_names,columns=clusters)
for cluster in clusters:
    sc_part = adata_norm[adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]].X.T
    sc_mean[cluster] = pd.DataFrame(sc_part.mean(axis=1),index=adata_norm.var_names)
centroids_sc = sc_mean
print(centroids_sc.shape)


# In[33]:


# Centroids constrained to single cell or tomo highly variable genes
centroids_sc_hv = centroids_sc.loc[adata_proc[:,adata_proc.var['highly_variable']==True].var_names,:]
#centroids_tomo_hv = centroids_sc.loc[tomo_loghv.var[tomo_loghv.var['highly_variable']==True].index.array]


# # AutoGeneS

# In[34]:


ag = AutoGenes(centroids_sc_hv.T)
ag.run(ngen=2000,seed=0,nfeatures=100,mode='fixed')


# In[ ]:


#ag_tomo = AutoGenes(centroids_tomo_hv.T)
#ag_tomo.run(ngen=2000,seed=0,nfeatures=100,mode='fixed')


# In[ ]:


#ag.plot(size='large',weights=(1,-1))


# In[35]:


pareto = ag.pareto
#pareto = ag_tomo.pareto


# In[36]:


#ag._AutoGenes__fitness_matrix()


# In[38]:


hv_genes = adata_proc.var[adata_proc.var['highly_variable']==True].index
#hv_genes
tomoData_norm.loc[hv_genes].sum(axis=1).sum(axis=0)
#adata_proc.var[adata_proc.var['highly_variable']==True].index


# In[39]:


# Plot sum of tomodata counts for each pareto point.
pareto_cordist = ag._AutoGenes__fitness_matrix()
#pareto_cordist[3, 2]
# Initialize dataframe: index pareto points, column for correlation, column for distance, column for tomocounts
pareto_results = pd.DataFrame(index= range(0, len(pareto)), columns = ['Correlation', 'Distance', 'Tomocounts'])
# Loop over pareto points and fill dataframe.
for p in pareto_results.index:
    #print(p)
    pareto_results.Correlation[p] = pareto_cordist[p, 1]
    pareto_results.Distance[p] = pareto_cordist[p, 0]
    pareto_results.Tomocounts[p] = tomoData_norm.loc[centroids_sc_hv[pareto[p]].index].sum(axis=1).sum(axis=0)
    #print(pareto_cordist[p, 0])
#    tomoData_norm.loc[centroids_sc_hv[pareto[330]].index].sum(axis=1).sum(axis=0)
#pareto_results


# In[42]:


plt.plot(pareto_results.Correlation, pareto_results.Tomocounts)
plt.show()


# In[43]:


pareto_results[pareto_results.Tomocounts > 14000]


# In[ ]:


#import inspect
#inspect.getsourcelines(ag.plot)


# In[ ]:


#len(pareto)


# In[44]:


centroids_sc_pareto = centroids_sc_hv[pareto[154]] # Also used 291. #ag.select(weights=(1,-1)) -->422
#centroids_sc_pareto = centroids_tomo_hv[pareto[242]]
#centroids_sc_pareto.shape


# In[50]:


add_genes = ['rbp4', 'klf2b']
for g in add_genes:
    if(g not in centroids_sc_pareto.index):
        centroids_sc_pareto = centroids_sc_pareto.append(centroids_sc_hv.loc[g])


# In[51]:


#centroids_sc_pareto['Fibroblast (col11a1a)']


# In[52]:


corr = pd.DataFrame(data = np.corrcoef(centroids_sc_pareto.T), 
                    columns = centroids_sc_pareto.columns, 
                    index = centroids_sc_pareto.columns)
mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True
#corr
with sb.axes_style("white"):
    sns_plot =sb.clustermap(np.abs(corr),cmap=sb.color_palette("GnBu", 1000), robust=True, yticklabels=1)
   #sns_plot.savefig("tomodiffgenes_over1000counts_paretomax_type_type_heatmap_GA.png",dpi=300)


# In[53]:


#subTypes = pd.DataFrame
#subTypes = centroids_sc_pareto.columns
#type_pal = sb.husl_palette(centroids_sc_pareto.columns.size, s=0.7)
#lut = dict(zip(centroids_sc_pareto.columns.unique(), type_pal))
#row_colors = subTypes.map(lut)
#sns_plot = sb.clustermap(centroids_sc_pareto.T, cmap="mako", robust=True) #, cmap="mako", robust=True, row_cluster = False)
#sns_plot.figure.savefig(figdir+"heatmap_GA.png",dpi=300)


# # Regression

# In[54]:


tomoData_hv = tomoData_norm.loc[centroids_sc_hv.index,:] # Tomo-expression of highly-variable sc genes
#tomoData_hv = tomoData_norm.loc[centroids_tomo_hv.index,:] # Tomo-expression of highly-variable tomo genes
tomoData_pareto = tomoData_norm.loc[centroids_sc_pareto.index,:] # Tomo-expression of pareto genes


# In[55]:


# centroids_sc_hv: single cell type centroid expression in highly variable genes.
# centroids_sc_pareto: sc centroid expression in pareto genes.
# tomoData_hv: expression of sc highly variable genes in normalized tomodata
# tomoData_pareto: expression of pareto genes in normalized tomodata
from sklearn.svm import NuSVR
from sklearn import linear_model
#ag.bulk_deconvolution(tomoData_norm.loc[centroids_sc_hv.index,:].T,model=NuSVR)
proportions_NuSVR = pd.DataFrame(index=tomoData_pareto.columns,columns=centroids_sc_pareto.columns)
# proportions_NuSVR is tomo-section x celltype dataframe
# Fit each tomo-section in terms of its pareto genes to the cell types, using a support vector regression.
for s in tomoData_pareto.columns:
    model = NuSVR(nu=0.3,C=0.3,kernel='linear')
    model.fit(centroids_sc_pareto,tomoData_pareto.loc[:,s])
    proportions_NuSVR.loc[s] = model.coef_[0]
# Remove negative coefficients and normalize so each section adds up to 1
proportions_NuSVR_norm = normalize_proportions(proportions_NuSVR, copy = True)


# In[56]:


x_pos = [i for i, _ in enumerate(proportions_NuSVR_norm.index)]
for cell in ['Cardiomyocytes A', 'Cardiomyocytes V',
             'Endocardium (A)', 'Endocardium (V)',
             'Fibroblast', 'Macrophages',
           'Erythrocytes', 'Smooth muscle cells']: #proportions_NuSVR_norm.columns:             
    plt.plot(x_pos,proportions_NuSVR_norm.loc[:,cell])
    plt.xlabel('section')
    plt.ylabel('proportion')
    plt.title(cell)
    #plt.savefig(cell + 'tomodiffgenes_over1000counts_paretomax_NuSVR_AVonly_final_celltypes_tomo.png')
    plt.show()
# Atrial expression is very spiky; ventricular epicardium is in the atrium. Cardiomyocytes A are very 
# small proportion in atrium, cardiomyocytes V are large proportion in ventricle (but also very abundant in atrium).
# ttn2 atrial cardiomyocytes are mostly ubiquitous, apart from one section. 


# In[58]:


x_pos = [i for i, _ in enumerate(proportions_NuSVR_norm.index)]
barbottom = [0 for i, _ in enumerate(proportions_NuSVR_norm.index)]
#cell = 'Cardiomyocytes V'
for cell in ['Cardiomyocytes A', 'Cardiomyocytes V',
             'Endocardium (A)', 'Endocardium (V)',
             'Fibroblast', 'Macrophages',
             'Erythrocytes', 'Smooth muscle cells']: #proportions_NuSVR_norm.columns:             
    plt.bar(x_pos, proportions_NuSVR_norm.loc[:,cell], color=cell_type_colors.color.loc[cell], 
            #width=sectionmergewidth.loc[x_pos].width.tolist(), 
            #width = 1,
            bottom = barbottom, label=cell)
    barbottom = barbottom + proportions_NuSVR_norm.loc[:,cell]
plt.xlabel('section')
#plt.legend(loc="upper left", fontsize = 'small')
plt.ylabel('proportion')
plt.title('Stacking')
#plt.savefig('High_tomo_over10counts_pareto62_200genes_NuSVR_tomo_bar.png')
plt.savefig('High_tomo_over1000counts_pareto62+addgenes_100genes_incl_mito_NuSVR_tomo_bar.png')
plt.show()


# In[ ]:


# Fit each tomo-section in terms of its pareto genes to the cell types, using non-negative least squares.
proportions_nnls = pd.DataFrame(index=tomoData_hv.columns,columns=centroids_sc_pareto.columns)
for s in tomoData_hv.columns:
    proportions_nnls.loc[s] = sci.optimize.nnls(centroids_sc_pareto,tomoData_pareto.loc[:,s])[0]
    #proportions_nnls.loc[s] = model.coef_[0]
proportions_nnls_norm = normalize_proportions(proportions_nnls, copy = True)


# In[ ]:


#cell_type_colors[cell_type_colors['Cell.type']==cell]
#type(cell_type_colors)
#cell_type_colors.color.loc[cell]
#barbottom
sectionmergewidth.loc[x_pos].width.tolist()


# In[ ]:


#x_pos = [i for i, _ in enumerate(proportions_nnls_norm.drop(['25'],axis=0).index)] # Why drop section 25?
x_pos = [i for i, _ in enumerate(proportions_nnls_norm.index)]
for cell in ['Cardiomyocytes A', 'Cardiomyocytes V',
             'Endocardium (A)', 'Endocardium (V)',
             'Fibroblast', 'Macrophages',
            'Erythrocytes', 'Smooth muscle cells']: #proportions_nnls_norm.columns:
    #plt.plot(x_pos,proportions_nnls_norm.drop(['25'],axis=0).loc[:,cell])
    plt.plot(x_pos,proportions_nnls_norm.loc[:,cell])
    plt.xlabel('section')
    plt.ylabel('proportion')
    plt.title(cell)
    #plt.savefig(cell + 'celltype_counts_over_1000_2000ngen__500feature_pareto200_nnls_AVonly_final_celltypes_tomo.png')
    plt.show()
# Cardiomyocytes A is in atrium but quite sparse, CM V mostly ventricular but also quite present in atrium.
# CM ttn2 A is in atrium (spiky but clearly), CM ttn2 V is mostly in ventricle but also spikes in atrium.
# Endo 1 A is spiky in both but more in atrium, endo 1 (V) is only in atrium
# Endo 2 A is present in one section, endo 2 V is ubiquitous with higher baseline in ventricle.
# Epi A is spiky but atrial, epi V is spiky and atrial (but low proportion)

