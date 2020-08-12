
# coding: utf-8

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


mito_genes = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/mito.genes.vs.txt',
                        header = None,
                        names = ['Remove_1', 'Remove_2', 'Gene'])
mito_genes = mito_genes.drop(columns = ['Remove_1', 'Remove_2'])


# In[7]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors_2.csv', 
                               index_col = 0)
cell_type_colors.loc['Erythrocytes'] = '#681b0e'


# In[34]:


#read annotations
annotations_nonery = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv', 
                                 index_col = 0)
annotations_ery = pd.read_csv('../../Data/final_erythrocytes.csv', index_col=0)
annotations_ery.columns = ['orig.ident', 'Cell_type']
annotations = pd.concat([annotations_nonery[['orig.ident', 'Cell_type']], annotations_ery])
indices = [x for x in annotations.index if not x.startswith('Hr')]
annotations = annotations.loc[indices,:]
annotations = annotations[~annotations.index.duplicated()]


# In[9]:


#highcount_celltypes = annotations.Cell_type.value_counts()[annotations.Cell_type.value_counts() > 1000]
#select_celltypes = highcount_celltypes.index.tolist() #+ ['Epicardium (Atrium)', 'Epicardium (Ventricle)']


# In[10]:


#annotations = annotations[annotations.Cell_type.isin(select_celltypes)]


# # Read tomo data, remove mito reads and duplicates

# In[11]:


tomoData = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Tomo1_DR11_count_table.csv',index_col=0)


# In[12]:


tomoData = tomoData.drop(index = list(set(mito_genes.Gene.tolist()) & set(tomoData.index.tolist())))


# In[31]:


tomoaData_noreplicates = sc.AnnData(tomoData.T)
tomoaData_noreplicates = removeReplicas(tomoaData_noreplicates)


# # Read sc data, remove mito reads and duplicates

# In[14]:


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


# In[15]:


adata_1.var_names_make_unique()
adata_2.var_names_make_unique()
adata_3.var_names_make_unique()
adata_4.var_names_make_unique()
adata_5.var_names_make_unique()


# In[16]:


adata=adata_1.concatenate([adata_2,adata_3,adata_4,adata_5],
                                      batch_key='sample',
                                      batch_categories=['H5','H6','H7','H8a','H8v'])
adata.X = adata.X.todense()


# In[26]:


non_mito_genes = [name for name in adata.var_names if name not in mito_genes.Gene.tolist()]
adata = adata[:, non_mito_genes]


# In[35]:


adata.obs_names = [str(adata.obs.loc[x,'sample'])+'_'+str(x.split('-', 1)[0]) for x in adata.obs_names]


# In[36]:


anno_drop = annotations.index.difference(adata.obs_names)
annotations = annotations.drop(anno_drop)
adata = adata[annotations.index]
for s in annotations.columns.values:
    adata.obs[s] = annotations[s].tolist()


# In[32]:


adata = removeReplicas(adata)


# # Intersect genes in tomo and single cell datasets

# In[42]:


adata = adata[:,adata.var_names.intersection(tomoaData_noreplicates.var_names)]
tomoaData_noreplicates = tomoaData_noreplicates[:,adata.var_names.intersection(tomoaData_noreplicates.var_names)]


# # Merge and normalize tomo-sections; determine highly-variable genes

# In[ ]:


#tomoaData_normalized = sc.AnnData(tomoData.T)
#print(tomoaData_normalized)
#tomoaData_normalized = removeReplicas(tomoaData_normalized)
#print(tomoaData_normalized)


# In[ ]:


#tomoaData_dataframe = pd.DataFrame(data=tomoaData_normalized.X.T,
#                                   index=tomoaData_normalized.var_names,
#                                   columns=tomoaData_normalized.obs_names)
#tomoaData_dataframe.to_csv(r'../write/cleaned_ids.csv')


# In[ ]:


#tomoaData_dataframe.shape


# In[43]:


sc.pp.normalize_per_cell(tomoaData_noreplicates, counts_per_cell_after=1e4)
tomoData_norm = pd.DataFrame(data=tomoaData_noreplicates.X.T,
                             columns = tomoaData_noreplicates.obs_names, 
                             index=tomoaData_noreplicates.var_names)


# In[46]:


counts = pd.DataFrame(data=tomoData.T.sum(1),index= tomoData.columns)
#raw_counts = counts
#counts = counts.T
#counts = counts[counts>1000]
#counts = counts.iloc[0].apply(pd.to_numeric)
#counts = counts.dropna()
#plt.hist(counts, color='blue')
#plt.show()
counts


# In[ ]:


counts_proc = counts.drop(['74'], axis=0)
counts_proc.shape


# In[ ]:


sc.set_figure_params(dpi=200, dpi_save=300, vector_friendly=True, frameon=False)


# In[ ]:


plt.rcParams['figure.figsize'] = (6, 4)


# In[ ]:


#print(counts_proc.shape)
fig, ax = plt.subplots()
x_pos = [i for i, _ in enumerate(counts_proc.index)]
#print(x_pos)
plt.xlabel('section')
plt.ylabel('# of counts')
ax.plot(x_pos,counts_proc)


# In[ ]:


# Create dictionary for new columns: if reads of dictionary element are < n (default 10k), add next
# column to current dictionary list; if not, make new dictionary list and add column to it.
#m = 0
#tomosectionmerge = {0:[]}
#for c in tomoData.columns:
#    if(tomoData[tomosectionmerge[m]].sum(axis=1).sum(axis=0) < 20000):
#        tomosectionmerge[m] = tomosectionmerge[m] + [c]
#    else:
#        m = m + 1
#        tomosectionmerge.update({m: [c]})


# In[ ]:


#tomoData[tomoData.columns[0]]
#tomosectionmerge


# In[ ]:


# Add up low-read columns from tomoData for tomosectionmerge: if reads of column in tomosectionmerge
# sum up to < n (default 20k), add next tomoData column to current tomosectionmerge column; 
# if not, make new tomosectionmerge column from tomoData column.
tomo_norep = pd.DataFrame(data = tomoaData_noreplicates.X.T,
                         columns = tomoaData_noreplicates.obs_names,
                         index = tomoaData_noreplicates.var_names)
#tomo_norep
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


# In[ ]:


#sectionmergewidth


# In[ ]:


#tomoaData_dataframe = pd.DataFrame(data=tomoaData_normalized.X.T,index=tomoaData_normalized.var_names,columns=tomoaData_normalized.obs_names)
#tomoaData_dataframe.to_csv(r'../write/cleaned_ids.csv')


# In[ ]:


tomosectionmerge.shape


# In[ ]:


tomoaData_normalized = sc.AnnData(tomosectionmerge.T)
sc.pp.normalize_per_cell(tomoaData_normalized, counts_per_cell_after=1e4)
tomoData_norm = pd.DataFrame(data=tomoaData_normalized.X.T,columns = tomoaData_normalized.obs_names, index=tomoaData_normalized.var_names)


# In[ ]:


#tomosectionmerge[m] = tomoData[c]
#tomosectionmerge.iloc[m].sum(axis=0)
#c = tomoData.columns[1]
#tomosectionmerge[m] = tomosectionmerge[m] + tomoData[c]
#tomosectionmerge.sum(axis=0)


# In[ ]:


# Sum up sections to make new data object
#HERE
#summed_tomo = tomoData.drop[columns=range(0,20)]
#for r in tomosectionmerge:
#    print(r)
#    print(tomosectionmerge[r])
#    print(tomoData[tomosectionmerge[r]].sum(axis=1).sum(axis=0))
#    summed_tomo[r] = tomoData.iloc[:,rangey[r]].sum(axis=1)
#tomoData.iloc[:,range(4)] #.sum(axis=1)


# In[ ]:


#counts_proc[0:50]


# In[ ]:


#adata_proc.var_names.intersection(tomoData_norm.index)
#HR_ps_1 = HR_ps_1[:, all_genes_but_RFP]
tomo_test = tomoaData_normalized[:, adata_proc.var_names.intersection(tomoData_norm.index)]


# In[ ]:


sc.pp.log1p(tomo_test)
sc.pp.highly_variable_genes(tomo_test, flavor='seurat', n_top_genes=4000)


# # Process single cell data to clusters; determine highly-variable genes

# In[ ]:


adata_norm = sc.pp.normalize_per_cell(adata_proc, counts_per_cell_after=1e4,copy=True)
adata_proc = sc.pp.log1p(adata_norm, copy=True)
sc.pp.highly_variable_genes(adata_proc, flavor='seurat', n_top_genes=4000)
#sc.pp.scale(adata_proc)


# In[ ]:


sc.tl.pca(adata_proc)
adata_proc.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pl.pca_scatter(adata_proc, save="_PCA.png" )


# In[ ]:


sc.pp.neighbors(adata_proc, n_neighbors=30)
sc.tl.umap(adata_proc)


# In[ ]:


sc.pl.umap(adata_proc, color='Cell_type')


# In[ ]:


sc.tl.rank_genes_groups(adata_proc, groupby = 'Cell_type')


# In[ ]:


sc.pl.rank_genes_groups_heatmap(adata_proc, n_genes = 5, save = 'marker_genes_over1000.png')


# In[ ]:


sc.tl.rank_genes_groups(adata_proc, groups = ['Endocardium (A)'], 
                        reference = 'Endocardium (V)', groupby = 'Cell_type')


# In[ ]:


sc.pl.rank_genes_groups_heatmap(adata_proc, n_genes = 20)


# In[ ]:


sc.pl.rank_genes_groups(adata_proc)


# In[ ]:


# Plot gene expression
gene = 'bzw1b'
expression_plot = tomoData_norm.loc[gene]
x_pos = [i for i, _ in enumerate(expression_plot.index)]
#print(x_pos)
fig, ax = plt.subplots()
plt.xlabel('section')
plt.ylabel('Normalized expression')
ax.plot(x_pos,expression_plot)
plt.title(gene)
#plt.savefig('tomo1_normalized_expression_' + gene + '.png')


# # Calculate centroids

# In[ ]:


adata_norm.obs = adata_proc.obs


# In[ ]:


clusters = list(set(adata_norm.obs['Cell_type']))
sc_mean = pd.DataFrame(index=adata_norm.var_names,columns=clusters)
for cluster in clusters:
    sc_part = adata_norm[adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]].X.T
    sc_mean[cluster] = pd.DataFrame(sc_part.mean(axis=1),index=adata_norm.var_names)
centroids_sc = sc_mean
print(centroids_sc.shape)


# In[ ]:


#cluster = "Cardiomyocytes A"
#adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]
#adata_norm[adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]] #.X #.T
#adata_norm.obs.index[adata_norm.obs.index.duplicated()]


# In[ ]:


centroids_sc_hv = centroids_sc.loc[adata_proc[:,adata_proc.var['highly_variable']==True].var_names,:]


# In[ ]:


centroids_tomo_hv = centroids_sc.loc[tomo_test.var[tomo_test.var['highly_variable']==True].index.array]


# In[ ]:


#tomo_test.var


# In[ ]:


#centroids_sc_hv = centroids_sc_hv.drop(['Dead cells', 'Neuronal cells'], axis=1)
#centroids_sc_hv = centroids_sc_hv.drop(['Neuronal cells'], axis=1)

centroids_sc_hv.shape


# In[ ]:


centroids_sc_hv.columns


# ## AutoGeneS

# In[ ]:


ag = AutoGenes(centroids_sc_hv.T)
ag.run(ngen=2000,seed=0,nfeatures=500,mode='fixed')


# In[ ]:


ag_tomo = AutoGenes(centroids_tomo_hv.T)
ag_tomo.run(ngen=2000,seed=0,nfeatures=100,mode='fixed')


# In[ ]:


#import pickle
#with open('ag_5000gen_seed0_nfeatures300.pickle', 'wb') as handle:
#    pickle.dump(ag, handle)

#with open('ag_5000gen_seed0_nfeatures300.pickle', 'rb') as handle:
#    b = pickle.load(handle)
    
#print(ag == b)


# In[ ]:


ag.plot(size='large',weights=(1,-1))


# In[ ]:


ag_tomo.plot(size='large',weights=(1,-1))


# In[ ]:


#pareto = ag.pareto
pareto = ag_tomo.pareto


# In[ ]:


#dir(ag)


# In[ ]:


ag._AutoGenes__fitness_matrix()


# In[ ]:


#import inspect
#inspect.getsourcelines(ag.plot)


# In[ ]:


len(pareto)


# In[ ]:


#centroids_sc_pareto = centroids_sc_hv[pareto[262]] # Also used 291. #ag.select(weights=(1,-1)) -->422
centroids_sc_pareto = centroids_tomo_hv[pareto[242]]
#centroids_sc_pareto.shape


# In[ ]:


import seaborn as sns
corr = pd.DataFrame(data = np.corrcoef(centroids_sc_pareto.T), 
                    columns = centroids_sc_pareto.columns, 
                    index = centroids_sc_pareto.columns)
mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    sns_plot =sns.clustermap(np.abs(corr),cmap=sns.color_palette("GnBu", 1000), robust=True, yticklabels=1)
    #sns_plot.savefig("tomodiffgenes_over1000counts_paretomax_type_type_heatmap_GA.png",dpi=300)


# In[ ]:


import seaborn as sns
subTypes = pd.DataFrame
subTypes = centroids_sc_pareto.columns
type_pal = sns.husl_palette(centroids_sc_pareto.columns.size, s=0.7)
lut = dict(zip(centroids_sc_pareto.columns.unique(), type_pal))
row_colors = subTypes.map(lut)
sns_plot = sns.clustermap(centroids_sc_pareto.T, cmap="mako", robust=True) #, cmap="mako", robust=True, row_cluster = False)
#sns_plot.figure.savefig(figdir+"heatmap_GA.png",dpi=300)


# ## Regression

# In[ ]:


#tomoData_hv = tomoData_norm.loc[centroids_sc_hv.index,:] # Tomo-expression of highly-variable sc genes
tomoData_hv = tomoData_norm.loc[centroids_tomo_hv.index,:] # Tomo-expression of highly-variable tomo genes
tomoData_pareto = tomoData_norm.loc[centroids_sc_pareto.index,:] # Tomo-expression of pareto genes


# In[ ]:


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


# In[ ]:


from sklearn.svm import NuSVR
from sklearn import linear_model
#ag.bulk_deconvolution(tomoData_norm.loc[centroids_sc_hv.index,:].T,model=NuSVR)
proportions_NuSVR = pd.DataFrame(index=tomoData_hv.columns,columns=centroids_sc_pareto.columns)
# proportions_NuSVR is tomo-section x celltype dataframe
# Fit each tomo-section in terms of its pareto genes to the cell types, using a support vector regression.
for s in tomoData_hv.columns:
    model = NuSVR(nu=0.3,C=0.3,kernel='linear')
    model.fit(centroids_sc_pareto,tomoData_pareto.loc[:,s])
    proportions_NuSVR.loc[s] = model.coef_[0]
# Remove negative coefficients and normalize so each section adds up to 1
proportions_NuSVR_norm = normalize_proportions(proportions_NuSVR, copy = True)


# In[ ]:


# Fit each tomo-section in terms of its pareto genes to the cell types, using non-negative least squares.
proportions_nnls = pd.DataFrame(index=tomoData_hv.columns,columns=centroids_sc_pareto.columns)
for s in tomoData_hv.columns:
    proportions_nnls.loc[s] = sci.optimize.nnls(centroids_sc_pareto,tomoData_pareto.loc[:,s])[0]
    #proportions_nnls.loc[s] = model.coef_[0]
proportions_nnls_norm = normalize_proportions(proportions_nnls, copy = True)


# In[ ]:


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


# In[ ]:


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
plt.legend(loc="upper left", fontsize = 'small')
plt.ylabel('proportion')
plt.title('Stacking')
plt.show()


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

