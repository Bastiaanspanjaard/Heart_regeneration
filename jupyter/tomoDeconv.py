
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


# # Read Tomo Data

# In[6]:


#tomoData = pd.read_csv('../Data/Tomo1_DR11_count_table.csv',index_col=0) #Tomo1_count_table.csv
tomoData = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Tomo2_DR11_count_table.csv',index_col=0)
print(tomoData.shape)
tomoData.head()


# In[7]:


tomoaData_normalized = sc.AnnData(tomoData.T)
print(tomoaData_normalized)
tomoaData_normalized = removeReplicas(tomoaData_normalized)
print(tomoaData_normalized)


# In[8]:


tomoaData_dataframe = pd.DataFrame(data=tomoaData_normalized.X.T,index=tomoaData_normalized.var_names,columns=tomoaData_normalized.obs_names)
#tomoaData_dataframe.to_csv(r'../write/cleaned_ids.csv')


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


# In[16]:


#print(counts_proc.shape)
fig, ax = plt.subplots()
x_pos = [i for i, _ in enumerate(counts_proc.index)]
#print(x_pos)
plt.xlabel('section')
plt.ylabel('# of counts')
ax.plot(x_pos,counts_proc)


# In[17]:


import os
print(os.getcwd())


# # Read single-cell data

# In[18]:


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


# In[19]:


adata_1.var_names_make_unique()
adata_2.var_names_make_unique()
adata_3.var_names_make_unique()
adata_4.var_names_make_unique()
adata_5.var_names_make_unique()


# In[94]:


adata=adata_1.concatenate([adata_2,adata_3,adata_4,adata_5],
                                      batch_key='sample',
                                      batch_categories=['H5','H6','H7','H8a','H8v'])
adata.X = adata.X.todense()


# In[20]:


#test = adata[:,adata.var_names.str.startswith('NC_002333.')].X
#np.mean(test,axis=0)


# In[21]:


#tomoData.index.intersection(adata.var_names)


# ## AutoGeneS

# In[88]:


#read annotations
# Update annotation?
#annotations = pd.read_csv('../../Data/all.hearts.all.cells.all.sub.sept03.csv',index_col=0)
#annotations = pd.read_csv('../../Data/celltypes_zoom_allcells.csv',index_col=0)
#annotations = pd.read_csv('../data/all.hearts.all.cells.all.sub.sept03.csv',index_col=0)
annotations_nonery = pd.read_csv('../../Data/final_metadata_Tmacromerged.csv', index_col = 0)
annotations_ery = pd.read_csv('../../Data/final_erythrocytes.csv', index_col=0)


# In[89]:


print(annotations_nonery.shape)
annotations_nonery.head()


# In[90]:


annotations_ery.columns = ['orig.ident', 'Cell_type']
#print(annotations_ery.shape)
#annotations_ery.head()


# In[91]:


annotations = pd.concat([annotations_nonery[['orig.ident', 'Cell_type']], annotations_ery])


# In[92]:


indices = [x for x in annotations.index if not x.startswith('Hr')]
annotations = annotations.loc[indices,:]


# In[95]:


adata.obs_names = [str(adata.obs.loc[x,'sample'])+'_'+str(x.split('-', 1)[0]) for x in adata.obs_names]


# In[98]:


anno_drop = annotations.index.difference(adata.obs_names)
annotations = annotations.drop(anno_drop)
adata_proc = adata[annotations.index]


# In[99]:


for s in annotations.columns.values:
    #adata_proc.obs[s] = annotations_filter[s].tolist()
    adata_proc.obs[s] = annotations[s].tolist()


# In[100]:


adata_proc = removeReplicas(adata_proc)


# In[101]:


#select shared genes between sc and tomo data
adata_proc = adata_proc[:,adata_proc.var_names.intersection(tomoData_norm.index)]
tomoData_norm = tomoData_norm.loc[adata_proc.var_names.intersection(tomoData_norm.index),:]


# In[102]:


print(adata_proc)
tomoData_norm.shape


# In[103]:


adata_norm = sc.pp.normalize_per_cell(adata_proc, counts_per_cell_after=1e4,copy=True)
adata_proc = sc.pp.log1p(adata_norm, copy=True)
sc.pp.highly_variable_genes(adata_proc, flavor='seurat', n_top_genes=4000)
#sc.pp.scale(adata_proc)


# In[104]:


sc.tl.pca(adata_proc)
adata_proc.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pl.pca_scatter(adata_proc, save="_PCA.png" )


# In[105]:


sc.pp.neighbors(adata_proc, n_neighbors=30)
sc.tl.umap(adata_proc)


# In[106]:


adata_proc


# In[107]:


#sc.pl.tsne(adata_proc, color='big.Cell.type')
#sc.pl.umap(adata_proc, color='big.Cell.type')
#sc.pl.tsne(adata_proc, color='orig.ident')
#sc.pl.umap(adata_proc, color='orig.ident')
sc.pl.umap(adata_proc, color='Cell_type')


# # Deconvolution

# In[108]:


adata_norm.obs = adata_proc.obs


# In[109]:


clusters = list(set(adata_norm.obs['Cell_type']))
sc_mean = pd.DataFrame(index=adata_norm.var_names,columns=clusters)
for cluster in clusters:
    sc_part = adata_norm[adata_norm.obs_names[adata_norm.obs['Cell_type']==cluster]].X.T
    sc_mean[cluster] = pd.DataFrame(sc_part.mean(axis=1),index=adata_norm.var_names)
centroids_sc = sc_mean
print(centroids_sc.shape)


# In[113]:


centroids_sc_hv = centroids_sc.loc[adata_proc[:,adata_proc.var['highly_variable']==True].var_names,:]


# In[115]:


centroids_sc_hv = centroids_sc_hv.drop(['Dead cells', 'Neuronal cells'], axis=1)
centroids_sc_hv.shape


# In[116]:


centroids_sc_hv.columns


# ## AutoGeneS

# In[ ]:


ag = AutoGenes(centroids_sc_hv.T)
ag.run(ngen=2000,seed=0,nfeatures=300,mode='fixed')


# In[99]:


#import pickle
#with open('ag_5000gen_seed0_nfeatures300.pickle', 'wb') as handle:
#    pickle.dump(ag, handle)

#with open('ag_5000gen_seed0_nfeatures300.pickle', 'rb') as handle:
#    b = pickle.load(handle)
    
#print(ag == b)


# In[87]:


ag.plot(size='large',weights=(1,-1))


# In[54]:


pareto = ag.pareto


# In[55]:


len(pareto)


# In[57]:


centroids_sc_pareto = centroids_sc_hv[pareto[240]]#ag.select(weights=(1,-1)) -->422
centroids_sc_pareto.shape


# In[59]:


import seaborn as sns
corr = pd.DataFrame(data = np.corrcoef(centroids_sc_pareto.T), columns = centroids_sc_pareto.columns, index = centroids_sc_pareto.columns)
mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    sns_plot =sns.clustermap(np.abs(corr),cmap=sns.color_palette("GnBu", 1000), robust=True, yticklabels=1)


# In[61]:


import seaborn as sns
subTypes = pd.DataFrame
subTypes = centroids_sc_pareto.columns
type_pal = sns.husl_palette(centroids_sc_pareto.columns.size, s=0.7)
lut = dict(zip(centroids_sc_pareto.columns.unique(), type_pal))
row_colors = subTypes.map(lut)
sns_plot = sns.clustermap(centroids_sc_pareto.T, cmap="mako", robust=True) #, cmap="mako", robust=True, row_cluster = False)
#sns_plot.figure.savefig(figdir+"heatmap_GA.png",dpi=300)


# ## Regression

# In[62]:


tomoData_hv = tomoData_norm.loc[centroids_sc_hv.index,:] # Tomo-expression of highly-variable genes
tomoData_pareto = tomoData_norm.loc[centroids_sc_pareto.index,:] # Tomo-expression of pareto genes


# In[63]:


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


# In[80]:


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


# In[81]:


# Fit each tomo-section in terms of its pareto genes to the cell types, using non-negative least squares.
proportions_nnls = pd.DataFrame(index=tomoData_hv.columns,columns=centroids_sc_pareto.columns)
for s in tomoData_hv.columns:
    proportions_nnls.loc[s] = sci.optimize.nnls(centroids_sc_pareto,tomoData_pareto.loc[:,s])[0]
    #proportions_nnls.loc[s] = model.coef_[0]
proportions_nnls_norm = normalize_proportions(proportions_nnls, copy = True)


# In[84]:


proportions_NuSVR_norm.columns


# In[92]:


x_pos = [i for i, _ in enumerate(proportions_NuSVR_norm.index)]
for cell in ['Cardiomyocytes A', 'Cardiomyocytes V', 
             'Cardiomyocytes (ttn.2) A', 'Cardiomyocytes (ttn.2) V', 
             'Endocardium 1 (A)', 'Endocardium 1 (V)',
             'Endocardium 2 (A)', 'Endocardium 2 (V)',
             'Epicardium (Atrium)','Epicardium (Ventricle)']: #proportions_NuSVR_norm.columns:             
    plt.plot(x_pos,proportions_NuSVR_norm.loc[:,cell])
    plt.xlabel('section')
    plt.ylabel('proportion')
    plt.title(cell)
    plt.savefig(cell + '_2000gen_pareto240_NuSVR_AVonly_final_celltypes_tomo.png')
    plt.show()
# Atrial expression is very spiky; ventricular epicardium is in the atrium. Cardiomyocytes A are very 
# small proportion in atrium, cardiomyocytes V are large proportion in ventricle (but also very abundant in atrium).
# ttn2 atrial cardiomyocytes are mostly ubiquitous, apart from one section. 


# In[91]:


#x_pos = [i for i, _ in enumerate(proportions_nnls_norm.drop(['25'],axis=0).index)] # Why drop section 25?
x_pos = [i for i, _ in enumerate(proportions_nnls_norm.index)]
for cell in ['Cardiomyocytes A', 'Cardiomyocytes V', 
             'Cardiomyocytes (ttn.2) A', 'Cardiomyocytes (ttn.2) V', 
             'Endocardium 1 (A)', 'Endocardium 1 (V)',
             'Endocardium 2 (A)', 'Endocardium 2 (V)',
             'Epicardium (Atrium)','Epicardium (Ventricle)']: #proportions_nnls_norm.columns:
    #plt.plot(x_pos,proportions_nnls_norm.drop(['25'],axis=0).loc[:,cell])
    plt.plot(x_pos,proportions_nnls_norm.loc[:,cell])
    plt.xlabel('section')
    plt.ylabel('proportion')
    plt.title(cell)
    plt.savefig(cell + '_2000gen_pareto240_nnls_AVonly_final_celltypes_tomo.png')
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


# In[105]:


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

