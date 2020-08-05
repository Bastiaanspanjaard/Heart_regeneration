#!/usr/bin/env python
# coding: utf-8

# In[11]:


import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import sys
import importlib
import scipy as sci

sys.path.append(".")

import pickle
from autogenes import AutoGenes


# In[12]:


import random


# In[13]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '1')
get_ipython().run_line_magic('aimport', 'autogenes')


# In[14]:


importlib.reload(autogenes)


# In[15]:


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

# In[9]:


#tomoData = pd.read_csv('../Data/Tomo1_DR11_count_table.csv',index_col=0) #Tomo1_count_table.csv
tomoData = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Tomo1_DR11_count_table.csv',index_col=0)
print(tomoData.shape)
tomoData.head()


# In[10]:


tomoaData_normalized = sc.AnnData(tomoData.T)
print(tomoaData_normalized)
tomoaData_normalized = removeReplicas(tomoaData_normalized)
print(tomoaData_normalized)


# In[16]:


tomoaData_dataframe = pd.DataFrame(data=tomoaData_normalized.X.T,index=tomoaData_normalized.var_names,columns=tomoaData_normalized.obs_names)
#tomoaData_dataframe.to_csv(r'../write/cleaned_ids.csv')


# In[17]:


tomoaData_dataframe.shape


# In[18]:


sc.pp.normalize_per_cell(tomoaData_normalized, counts_per_cell_after=1e4)
tomoData_norm = pd.DataFrame(data=tomoaData_normalized.X.T,columns = tomoaData_normalized.obs_names, index=tomoaData_normalized.var_names)


# In[267]:


tomo_test


# In[19]:


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


# In[20]:


counts_proc = counts.drop(['74'], axis=0)
counts_proc.shape


# In[21]:


sc.set_figure_params(dpi=200, dpi_save=300, vector_friendly=True, frameon=False)


# In[22]:


plt.rcParams['figure.figsize'] = (6, 4)


# In[23]:


#print(counts_proc.shape)
fig, ax = plt.subplots()
x_pos = [i for i, _ in enumerate(counts_proc.index)]
#print(x_pos)
plt.xlabel('section')
plt.ylabel('# of counts')
ax.plot(x_pos,counts_proc)


# In[24]:


import os
print(os.getcwd())


# # Read single-cell data

# In[90]:


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


# In[91]:


adata_1.var_names_make_unique()
adata_2.var_names_make_unique()
adata_3.var_names_make_unique()
adata_4.var_names_make_unique()
adata_5.var_names_make_unique()


# In[191]:


adata=adata_1.concatenate([adata_2,adata_3,adata_4,adata_5],
                                      batch_key='sample',
                                      batch_categories=['H5','H6','H7','H8a','H8v'])
adata.X = adata.X.todense()


# In[261]:


adata


# In[192]:


#test = adata[:,adata.var_names.str.startswith('NC_002333.')].X
#np.mean(test,axis=0)


# In[193]:


#tomoData.index.intersection(adata.var_names)


# ## AutoGeneS

# In[194]:


#read annotations
# Update annotation?
#annotations = pd.read_csv('../../Data/all.hearts.all.cells.all.sub.sept03.csv',index_col=0)
#annotations = pd.read_csv('../../Data/celltypes_zoom_allcells.csv',index_col=0)
#annotations = pd.read_csv('../data/all.hearts.all.cells.all.sub.sept03.csv',index_col=0)
#annotations_nonery = pd.read_csv('../../Data/final_metadata_Tmacromerged.csv', index_col = 0)
annotations_nonery = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv', index_col = 0)
annotations_ery = pd.read_csv('../../Data/final_erythrocytes.csv', index_col=0)


# In[195]:


print(annotations_nonery.shape)
annotations_nonery.head()


# In[196]:


annotations_ery.columns = ['orig.ident', 'Cell_type']
#print(annotations_ery.shape)
#annotations_ery.head()


# In[197]:


annotations = pd.concat([annotations_nonery[['orig.ident', 'Cell_type']], annotations_ery])


# In[198]:


indices = [x for x in annotations.index if not x.startswith('Hr')]
annotations = annotations.loc[indices,:]
annotations = annotations[~annotations.index.duplicated()]
#non_dupe_indices = annotations.index[~annotations.index.duplicated()]
#annotations = annotations.loc[non_dupe_indices, :]


# In[199]:


highcount_celltypes = annotations.Cell_type.value_counts()[annotations.Cell_type.value_counts() > 1000]


# In[200]:


annotations = annotations[annotations.Cell_type.isin(highcount_celltypes.index)]


# In[201]:


adata.obs_names = [str(adata.obs.loc[x,'sample'])+'_'+str(x.split('-', 1)[0]) for x in adata.obs_names]


# In[202]:


anno_drop = annotations.index.difference(adata.obs_names)
annotations = annotations.drop(anno_drop)
adata_proc = adata[annotations.index]


# In[203]:


adata_proc
#annotations
#annotations.index[~annotations.index.duplicated()]


# In[204]:


for s in annotations.columns.values:
    #adata_proc.obs[s] = annotations_filter[s].tolist()
    adata_proc.obs[s] = annotations[s].tolist()


# In[242]:


#adata_proc = removeReplicas(adata_proc)


# In[205]:


#select shared genes between sc and tomo data
adata_proc = adata_proc[:,adata_proc.var_names.intersection(tomoData_norm.index)]
tomoData_norm = tomoData_norm.loc[adata_proc.var_names.intersection(tomoData_norm.index),:]


# In[294]:


#adata_proc.var_names.intersection(tomoData_norm.index)
#HR_ps_1 = HR_ps_1[:, all_genes_but_RFP]
tomo_test = tomoaData_normalized[:, adata_proc.var_names.intersection(tomoData_norm.index)]


# In[295]:


sc.pp.log1p(tomo_test)
sc.pp.highly_variable_genes(tomo_test, flavor='seurat', n_top_genes=4000)


# In[296]:


tomo_test


# In[243]:


print(adata_proc)
tomoData_norm.shape


# In[247]:


tomo_test = tomoData_norm
tomo_test


# In[246]:


sc.pp.log1p(tomo_test)


# In[207]:


adata_norm = sc.pp.normalize_per_cell(adata_proc, counts_per_cell_after=1e4,copy=True)
adata_proc = sc.pp.log1p(adata_norm, copy=True)
sc.pp.highly_variable_genes(adata_proc, flavor='seurat', n_top_genes=4000)
#sc.pp.scale(adata_proc)


# In[208]:


sc.tl.pca(adata_proc)
adata_proc.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pl.pca_scatter(adata_proc, save="_PCA.png" )


# In[209]:


sc.pp.neighbors(adata_proc, n_neighbors=30)
sc.tl.umap(adata_proc)


# In[210]:


adata_proc


# In[211]:


#sc.pl.tsne(adata_proc, color='big.Cell.type')
#sc.pl.umap(adata_proc, color='big.Cell.type')
#sc.pl.tsne(adata_proc, color='orig.ident')
#sc.pl.umap(adata_proc, color='orig.ident')
sc.pl.umap(adata_proc, color='Cell_type')


# # Deconvolution

# In[212]:


adata_norm.obs = adata_proc.obs


# In[213]:


clusters = list(set(adata_norm.obs['Cell_type']))
sc_mean = pd.DataFrame(index=adata_norm.var_names,columns=clusters)
for cluster in clusters:
    sc_part = adata_norm[adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]].X.T
    sc_mean[cluster] = pd.DataFrame(sc_part.mean(axis=1),index=adata_norm.var_names)
centroids_sc = sc_mean
print(centroids_sc.shape)


# In[214]:


#cluster = "Cardiomyocytes A"
#adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]
#adata_norm[adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]] #.X #.T
#adata_norm.obs.index[adata_norm.obs.index.duplicated()]


# In[215]:


centroids_sc_hv = centroids_sc.loc[adata_proc[:,adata_proc.var['highly_variable']==True].var_names,:]


# In[298]:


centroids_tomo_hv = centroids_sc.loc[tomo_test.var[tomo_test.var['highly_variable']==True].index.array]


# In[217]:


#centroids_sc_hv = centroids_sc_hv.drop(['Dead cells', 'Neuronal cells'], axis=1)
#centroids_sc_hv = centroids_sc_hv.drop(['Neuronal cells'], axis=1)

centroids_sc_hv.shape


# In[ ]:


centroids_sc_hv.columns


# ## AutoGeneS

# In[219]:


ag = AutoGenes(centroids_sc_hv.T)
ag.run(ngen=2000,seed=0,nfeatures=500,mode='fixed')


# In[299]:


ag_tomo = AutoGenes(centroids_tomo_hv.T)
ag_tomo.run(ngen=2000,seed=0,nfeatures=500,mode='fixed')


# In[44]:


#import pickle
#with open('ag_5000gen_seed0_nfeatures300.pickle', 'wb') as handle:
#    pickle.dump(ag, handle)

#with open('ag_5000gen_seed0_nfeatures300.pickle', 'rb') as handle:
#    b = pickle.load(handle)
    
#print(ag == b)


# In[220]:


ag.plot(size='large',weights=(1,-1))


# In[300]:


ag_tomo.plot(size='large',weights=(1,-1))


# In[301]:


#pareto = ag.pareto
pareto = ag_tomo.pareto


# In[222]:


#dir(ag)


# In[223]:


ag._AutoGenes__fitness_matrix()


# In[224]:


#import inspect
#inspect.getsourcelines(ag.plot)


# In[302]:


len(pareto)


# In[304]:


#centroids_sc_pareto = centroids_sc_hv[pareto[200]] # Also used 291. #ag.select(weights=(1,-1)) -->422
centroids_sc_pareto = centroids_tomo_hv[pareto[345]]
centroids_sc_pareto.shape


# In[305]:


import seaborn as sns
corr = pd.DataFrame(data = np.corrcoef(centroids_sc_pareto.T), 
                    columns = centroids_sc_pareto.columns, 
                    index = centroids_sc_pareto.columns)
mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    sns_plot =sns.clustermap(np.abs(corr),cmap=sns.color_palette("GnBu", 1000), robust=True, yticklabels=1)
    sns_plot.savefig("tomodiffgenes_over1000counts_paretomax_type_type_heatmap_GA.png",dpi=300)


# In[306]:


import seaborn as sns
subTypes = pd.DataFrame
subTypes = centroids_sc_pareto.columns
type_pal = sns.husl_palette(centroids_sc_pareto.columns.size, s=0.7)
lut = dict(zip(centroids_sc_pareto.columns.unique(), type_pal))
row_colors = subTypes.map(lut)
sns_plot = sns.clustermap(centroids_sc_pareto.T, cmap="mako", robust=True) #, cmap="mako", robust=True, row_cluster = False)
#sns_plot.figure.savefig(figdir+"heatmap_GA.png",dpi=300)


# ## Regression

# In[307]:


#tomoData_hv = tomoData_norm.loc[centroids_sc_hv.index,:] # Tomo-expression of highly-variable sc genes
tomoData_hv = tomoData_norm.loc[centroids_tomo_hv.index,:] # Tomo-expression of highly-variable tomo genes
tomoData_pareto = tomoData_norm.loc[centroids_sc_pareto.index,:] # Tomo-expression of pareto genes


# In[252]:


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


# In[308]:


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


# In[253]:


# Fit each tomo-section in terms of its pareto genes to the cell types, using non-negative least squares.
proportions_nnls = pd.DataFrame(index=tomoData_hv.columns,columns=centroids_sc_pareto.columns)
for s in tomoData_hv.columns:
    proportions_nnls.loc[s] = sci.optimize.nnls(centroids_sc_pareto,tomoData_pareto.loc[:,s])[0]
    #proportions_nnls.loc[s] = model.coef_[0]
proportions_nnls_norm = normalize_proportions(proportions_nnls, copy = True)


# In[309]:


x_pos = [i for i, _ in enumerate(proportions_NuSVR_norm.index)]
for cell in ['Cardiomyocytes A', 'Cardiomyocytes V',
             'Endocardium (A)', 'Endocardium (V)',
             'Fibroblast', 'Macrophages',
            'Erythrocytes', 'Smooth muscle cells']: #proportions_NuSVR_norm.columns:             
    plt.plot(x_pos,proportions_NuSVR_norm.loc[:,cell])
    plt.xlabel('section')
    plt.ylabel('proportion')
    plt.title(cell)
    plt.savefig(cell + 'tomodiffgenes_over1000counts_paretomax_NuSVR_AVonly_final_celltypes_tomo.png')
    plt.show()
# Atrial expression is very spiky; ventricular epicardium is in the atrium. Cardiomyocytes A are very 
# small proportion in atrium, cardiomyocytes V are large proportion in ventricle (but also very abundant in atrium).
# ttn2 atrial cardiomyocytes are mostly ubiquitous, apart from one section. 


# In[254]:


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
    plt.savefig(cell + 'celltype_counts_over_1000_2000ngen__500feature_pareto200_nnls_AVonly_final_celltypes_tomo.png')
    plt.show()
# Cardiomyocytes A is in atrium but quite sparse, CM V mostly ventricular but also quite present in atrium.
# CM ttn2 A is in atrium (spiky but clearly), CM ttn2 V is mostly in ventricle but also spikes in atrium.
# Endo 1 A is spiky in both but more in atrium, endo 1 (V) is only in atrium
# Endo 2 A is present in one section, endo 2 V is ubiquitous with higher baseline in ventricle.
# Epi A is spiky but atrial, epi V is spiky and atrial (but low proportion)


# In[90]:


#proportions_nnls_norm.drop(['25'],axis=0)


# In[79]:


#proportions_nnls_norm['Cardiomyocytes V']


# In[173]:


cell = 'Erythrocytes'
plt.plot(x_pos,proportions_nnls_norm.drop(['25'],axis=0).loc[:,cell])
plt.xlabel('section')
plt.ylabel('proportion')
plt.title(cell)
#plt.show()
plt.savefig('Plot_try.png')


# In[104]:


plt.plot([0, 1, 2, 3, 4], [0, 3, 5, 9, 11])
plt.xlabel('Months')
plt.ylabel('Books Read')
#plt.show()
plt.savefig('books_read.png')


# In[57]:


proportions_NuSVR


# In[52]:


proportions_NuSVR


# In[1]:


x_pos = [i for i, _ in enumerate(proportions_nnls_norm.drop(['25'],axis=0).index)]
for cell in proportions_nnls_norm.columns:
    plt.plot(x_pos,proportions_nnls_norm.drop(['25'],axis=0).loc[:,cell])
    plt.xlabel('section')
    plt.ylabel('proportion')
    plt.title(cell)
    plt.show()


# # Subsampling

# In[53]:


#deconvolving subsamples and get a confidence interval 


# In[56]:


tomoData_Resample = pd.read_csv('../data/Tomo1_Resample.csv',index_col=0) #Tomo1_count_table.csv
tomoaData_Resample_norm = sc.AnnData(tomoData_Resample.T)
sc.pp.normalize_per_cell(tomoaData_Resample_norm, counts_per_cell_after=1e4)
tomoData_Resample_norm = pd.DataFrame(data=tomoaData_Resample_norm.X.T,columns = tomoaData_Resample_norm.obs_names, index=tomoaData_Resample_norm.var_names)


# In[58]:


tomoData_Resample_norm.shape


# In[59]:


tomoData_Resample_norm.columns


# In[60]:


proportions_nnls_resample = pd.DataFrame(index=tomoData_Resample_norm.columns,columns=centroids_sc_pareto.columns)
for s in tomoData_Resample_norm.columns:
    proportions_nnls_resample.loc[s] = sci.optimize.nnls(centroids_sc_pareto,tomoData_Resample_norm.loc[centroids_sc_pareto.index,s])[0]
proportions_nnls_norm = normalize_proportions(proportions_nnls_resample, copy = True)


# In[141]:


proportions_nnls_resample


# In[70]:


proportions_nnls_resample.loc[:,"section"]=proportions_nnls_resample.index.str[6:8]
sns.set(rc={'figure.figsize':(20,8)})
sns.set(style="whitegrid")
for celltype in proportions_nnls_resample.columns:
    if celltype!='section':
        #sns.barplot(x="section", y=celltype, data=proportions_nnls_resample, capsize=.2)
        sns.catplot(x="section", y=celltype, data=proportions_nnls_resample,
                kind="bar",color='red')
        #sns.boxplot(data=proportions_nnls_resample,x="section",y=celltype,color='red')


# In[175]:


def plot_pdf(data):
    mu, std = sci.stats.norm.fit(data)

    # Plot the histogram.
    print(np.count_nonzero(data))
    if np.count_nonzero(data)!=0 && mu!=np.nan && std!=nan:
        print(data)
        plt.hist(np.array(data), bins=10, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = sci.stats.norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)

    plt.show()
    return


# In[63]:


mus = pd.DataFrame(columns=proportions_nnls_resample.columns,index=tomoData.columns)
stds = pd.DataFrame(columns=proportions_nnls_resample.columns,index=tomoData.columns)
for s in tomoData.columns:
    section_samples = proportions_nnls_resample.loc[[x for x in proportions_nnls_resample.index if 'X'+str(s) in x],:]
    for celltype in proportions_nnls_resample.columns:
        #print(sci.stats.t.pdf(section_samples.loc[:,celltype],df=1))
        data = section_samples.loc[:,celltype]
        mu, std = sci.stats.norm.fit(list(data))
        mus.loc[s,celltype] = mu
        stds.loc[s,celltype] = std
        #plot_pdf(data)


# In[64]:


for cell in mus.columns:
    plt.plot(x_pos,mus.loc[:,cell])
    print(cell)
    plt.show()


# In[66]:


mus.fillna(value=np.nan, inplace=True)
stds.fillna(value=np.nan, inplace=True)


# In[67]:


with sns.axes_style("white"):
    sns.heatmap(mus,cmap=sns.color_palette("GnBu", 1000))


# In[68]:


with sns.axes_style("white"):
    sns.heatmap(stds,cmap=sns.color_palette("GnBu", 1000))


# In[ ]:


fig, ax = plt_.subplots()
ax.set_title('standard error (%)')
es = [x*100 for x in es]
es_raw = [x*100 for x in es_raw]
es_mg = [x*100 for x in es_mg]
bplot = ax.boxplot([es,es_raw,es_mg], patch_artist=True, showfliers=False)
colors = ['red', 'blue', 'green']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MOO-based genes', 'HV genes','marker genes'])
plt.show()


# In[ ]:




