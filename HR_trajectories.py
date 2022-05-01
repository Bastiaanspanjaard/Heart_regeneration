#!/usr/bin/env python
# coding: utf-8

# # Dependencies and parameters

# In[9]:


from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))
get_ipython().run_line_magic('matplotlib', 'inline')
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
from sklearn import preprocessing
import seaborn as sns
import h5py
results_file = './write/HR_trajectories.h5ad'
sc.set_figure_params(dpi_save = 300)
import scvelo as scv
import cellrank as cr
scv.settings.set_figure_params('scvelo', dpi_save = 300)


# In[10]:


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


# In[11]:


def GetDiffGenes(adata):
    # Return differentially expressed genes (p < 0.01, logfoldchange > 1) per cell type in a AnnData -
    # requires calculation of differentially expressed genes first.
    result = pd.DataFrame({})
    for ctype in adata.obs['Cell_type'].cat.categories:
        ctype_d = {'Gene' : adata.uns['rank_genes_groups']['names'][ctype],
                   'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][ctype],
                   'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][ctype],
                   'Cell_type':  np.repeat(ctype, len(adata.var.index))}
        ctype_degdf = pd.DataFrame(data=ctype_d)
        result = result.append(ctype_degdf[(ctype_degdf.pvals_adj < 0.01) & (ctype_degdf.logfoldchanges > 1)], ignore_index=True)
    return result


# In[12]:


def CountSecretome(adata, secretome):
    # Return secretome expression per cell type - mean and standard error of the mean.
    secretome_ind = [adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index]

    result = pd.DataFrame({'Cell_type' : adata.obs['Cell_type'].cat.categories,
                          'Secretome' : np.repeat(-1, len(adata.obs['Cell_type'].cat.categories)),
                          'SEM' : np.repeat(-1, len(adata.obs['Cell_type'].cat.categories))})
    for ctype in adata.obs['Cell_type'].cat.categories:
        ctype_index = [x for x in range(0, len(adata.obs) - 1) if adata.obs.Cell_type[x] == ctype]
        if(len(ctype_index) == 0):
            continue
        secretome_slice = adata.X[:,secretome_ind]
        result.Secretome[result['Cell_type'] == ctype] = (sum(secretome_slice[ctype_index, :].sum(axis = 1))/len(ctype_index))[0,0]
        result.SEM[result['Cell_type'] == ctype] = np.std(secretome_slice[ctype_index, :].sum(axis = 1))/np.sqrt(len(ctype_index))
    result = result[result.Secretome != -1]
    return(result)


# In[13]:


def ExpressedSecretome(adata, secretome):
    # Return mean expression and z-score expression (scaled over genes within a cell type) for secretome genes
    # per cell type in an AnnData object adata.
    secretome_ind = [adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index]

    secretome_slice = adata.X[:,secretome_ind]
    secretome_dense_slice = adata.X[:,secretome_ind].todense()
    secretome_dense_slice_z = preprocessing.scale(secretome_dense_slice)

    ctype_averages = pd.DataFrame(index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index])
    ctype_z_averages = pd.DataFrame(index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index])
    for ctype in adata.obs['Cell_type'].cat.categories:
        ctype_index = [x for x in range(0, len(adata.obs) - 1) if adata.obs.Cell_type[x] == ctype]
        if(len(ctype_index) == 0):
            continue
        ctype_average = pd.DataFrame({ctype : np.squeeze(np.asarray(secretome_slice[ctype_index, :].mean(axis = 0)))},
                                  index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index]) # axis = 0 gets us gene-level averages. #.sum(axis = 1))/len(ctype_index))[0,0]
        ctype_averages = ctype_averages.join(ctype_average)
        
        ctype_z_average = pd.DataFrame({ctype : secretome_dense_slice_z[ctype_index, :].mean(axis = 0)},
                                    index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index])
        ctype_z_averages = ctype_z_averages.join(ctype_z_average)
    return(ctype_averages, ctype_z_averages)


# In[14]:


annotations = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv', index_col = 0)


# In[15]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors_2.csv', index_col = 0)


# In[16]:


epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (cfd)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (cxcl12a)', 'Fibroblasts (mpeg1.1)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)',
                    'Valve fibroblasts', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[17]:


connected_epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)','Fibroblasts (cxcl12a)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[18]:


local_connected_epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells', 'Fibroblasts (cxcl12a)']


# In[19]:


connected_endofibro_types = ['Endocardium (frzb)', 'Endocardium (Ventricle)', 
                            'Fibroblasts (nppc)', 'Fibroblasts (spock3)']


# In[20]:


mito_genes = [line.rstrip('\n').rsplit(',')[2] for line in open('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/mito.genes.vs.txt')]


# In[21]:


HR_setnames = pd.DataFrame({'batch': ['0', '1', '2', '3', '4', 
                                      '5', '6', '7', '8', '9', 
                                      '10', '11', '12','13', '14', 
                                      '15', '16', '17', '18', '19', 
                                      '20', '21', '22','23', '24', 
                                      '25', '26', '27', '28', '29', 
                                      '30', '31', '32','33', '34', 
                                      '35', '36', '37', '38', '39', 
                                      '40', '41', '42'],
                            'heart': ['H5', 'H6', 'H7', 'H8a', 'H8v', 
                                      'Hr1', 'Hr2a', 'Hr2b', 'Hr3', 'Hr4', 
                                      'Hr5', 'Hr6a', 'Hr6v', 'Hr7a', 'Hr7v', 
                                      'Hr8', 'Hr9', 'Hr10', 'Hr11', 'Hr12', 
                                      'Hr13', 'Hr14', 'Hr15', 'Hr16', 'Hr17', 
                                      'Hr18', 'Hr19', 'Hr20', 'Hr21', 'Hr22', 
                                      'Hr23', 'Hr24', 'Hr25', 'Hr26', 'Hr27', 
                                      'Hr28', 'Hr29', 'Hr30', 'Hr31', 'Hr32', 
                                      'Hr33', 'Hr34', 'Hr35'],
                             'dpi':  ['0', '0', '0', '0', '0', 
                                      '7', '7', '7', '30', '30', 
                                      '60', '7', '7', '7', '7', 
                                      '7', '7', '3', '3', '3', 
                                      '7', '7', '7', '15', '15', 
                                      '15', '30', '30', '30', '3', 
                                      '3', '3', '3', '3', '3', 
                                      '3', '3', '7', '7', '7', 
                                      '7', '3', '3'],
                           'inhib':  ['no', 'no', 'no', 'no', 'no', 
                                      'no', 'no', 'no', 'no', 'no', 
                                      'no', 'no', 'no', 'no', 'no', 
                                      'no', 'no', 'no', 'no', 'no', 
                                      'no', 'no', 'no', 'no', 'no', 
                                      'no', 'no', 'no', 'no', 'no', 
                                      'no', 'no', 'no', 'no', 'no', 
                                      'DMSO', 'IWR1', 'DMSO', 'IWR1', 'IWR1', 
                                      'IWR1', 'IWR1', 'IWR1']})
HR_setnames.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/HR_setnames.csv')


# # Load and merge RNA velocity datasets

# In[22]:


# H5_Rv_data = scv.read('../Data/RNAvelo/H5_v3Dr11.loom', cache=True)
# H5_Rv_data.var_names_make_unique()
# H6_Rv_data = scv.read('../Data/RNAvelo/H6_v3Dr11.loom', cache=True)
# H6_Rv_data.var_names_make_unique()
# H7_Rv_data = scv.read('../Data/RNAvelo/H7_v3Dr11.loom', cache=True)
# H7_Rv_data.var_names_make_unique()
# H8a_Rv_data = scv.read('../Data/RNAvelo/H8a_v3Dr11.loom', cache=True)
# H8a_Rv_data.var_names_make_unique()
# H8v_Rv_data = scv.read('../Data/RNAvelo/H8v_v3Dr11.loom', cache=True)
# H8v_Rv_data.var_names_make_unique()
# Hr1_Rv_data = scv.read('../Data/RNAvelo/Hr1_v3Dr11.loom', cache=True)
# Hr1_Rv_data.var_names_make_unique()
# Hr2a_Rv_data = scv.read('../Data/RNAvelo/Hr2a_v3Dr11.loom', cache=True)
# Hr2a_Rv_data.var_names_make_unique()
# Hr2b_Rv_data = scv.read('../Data/RNAvelo/Hr2b_v3Dr11.loom', cache=True)
# Hr2b_Rv_data.var_names_make_unique()
# Hr3_Rv_data = scv.read('../Data/RNAvelo/Hr3_v3Dr11.loom', cache=True)
# Hr3_Rv_data.var_names_make_unique()
# Hr4_Rv_data = scv.read('../Data/RNAvelo/Hr4_v3Dr11.loom', cache=True)
# Hr4_Rv_data.var_names_make_unique()
# Hr5_Rv_data = scv.read('../Data/RNAvelo/Hr5_v3Dr11.loom', cache=True)
# Hr5_Rv_data.var_names_make_unique()
# Hr6a_Rv_data = scv.read('../Data/RNAvelo/Hr6a_v3Dr11.loom', cache=True)
# Hr6a_Rv_data.var_names_make_unique()
# Hr6v_Rv_data = scv.read('../Data/RNAvelo/Hr6v_v3Dr11.loom', cache=True)
# Hr6v_Rv_data.var_names_make_unique()
# Hr7a_Rv_data = scv.read('../Data/RNAvelo/Hr7a_v3Dr11.loom', cache=True)
# Hr7a_Rv_data.var_names_make_unique()
# Hr7v_Rv_data = scv.read('../Data/RNAvelo/Hr7v_v3Dr11.loom', cache=True)
# Hr7v_Rv_data.var_names_make_unique()
# Hr8_Rv_data = scv.read('../Data/RNAvelo/Hr8_v3Dr11.loom', cache=True)
# Hr8_Rv_data.var_names_make_unique()
# Hr9_Rv_data = scv.read('../Data/RNAvelo/Hr9_v3Dr11.loom', cache=True)
# Hr9_Rv_data.var_names_make_unique()
# Hr10_Rv_data = scv.read('../Data/RNAvelo/Hr10_v3Dr11.loom', cache=True)
# Hr10_Rv_data.var_names_make_unique()
# Hr11_Rv_data = scv.read('../Data/RNAvelo/Hr11_v3Dr11.loom', cache=True)
# Hr11_Rv_data.var_names_make_unique()
# Hr12_Rv_data = scv.read('../Data/RNAvelo/Hr12_v3Dr11.loom', cache=True)
# Hr12_Rv_data.var_names_make_unique()
# Hr13_Rv_data = scv.read('../Data/RNAvelo/Hr13_v3Dr11.loom', cache=True)
# Hr13_Rv_data.var_names_make_unique()
# Hr14_Rv_data = scv.read('../Data/RNAvelo/Hr14_v3Dr11.loom', cache=True)
# Hr14_Rv_data.var_names_make_unique()
# Hr15_Rv_data = scv.read('../Data/RNAvelo/Hr15_v3Dr11.loom', cache=True)
# Hr15_Rv_data.var_names_make_unique()
# Hr16_Rv_data = scv.read('../Data/RNAvelo/Hr16_v3Dr11.loom', cache=True)
# Hr16_Rv_data.var_names_make_unique()
# Hr17_Rv_data = scv.read('../Data/RNAvelo/Hr17_v3Dr11.loom', cache=True)
# Hr17_Rv_data.var_names_make_unique()
# Hr18_Rv_data = scv.read('../Data/RNAvelo/Hr18_v3Dr11.loom', cache=True)
# Hr18_Rv_data.var_names_make_unique()
# Hr19_Rv_data = scv.read('../Data/RNAvelo/Hr19_v3Dr11.loom', cache=True)
# Hr19_Rv_data.var_names_make_unique()
# Hr20_Rv_data = scv.read('../Data/RNAvelo/Hr20_v3Dr11.loom', cache=True)
# Hr20_Rv_data.var_names_make_unique()
# Hr21_Rv_data = scv.read('../Data/RNAvelo/Hr21_v3Dr11.loom', cache=True)
# Hr21_Rv_data.var_names_make_unique()
# Hr22_Rv_data = scv.read('../Data/RNAvelo/Hr22_v3Dr11.loom', cache=True)
# Hr22_Rv_data.var_names_make_unique()
# Hr23_Rv_data = scv.read('../Data/RNAvelo/Hr23_v3Dr11.loom', cache=True)
# Hr23_Rv_data.var_names_make_unique()
# Hr24_Rv_data = scv.read('../Data/RNAvelo/Hr24_v3Dr11.loom', cache=True)
# Hr24_Rv_data.var_names_make_unique()
# Hr25_Rv_data = scv.read('../Data/RNAvelo/Hr25_v3Dr11.loom', cache=True)
# Hr25_Rv_data.var_names_make_unique()
# Hr26_Rv_data = scv.read('../Data/RNAvelo/Hr26_v3Dr11.loom', cache=True)
# Hr26_Rv_data.var_names_make_unique()
# Hr27_Rv_data = scv.read('../Data/RNAvelo/Hr27_v3Dr11.loom', cache=True)
# Hr27_Rv_data.var_names_make_unique()
# Hr28_Rv_data = scv.read('../Data/RNAvelo/Hr28_v3Dr11.loom', cache=True)
# Hr28_Rv_data.var_names_make_unique()
# Hr29_Rv_data = scv.read('../Data/RNAvelo/Hr29_v3Dr11.loom', cache=True)
# Hr29_Rv_data.var_names_make_unique()
# Hr30_Rv_data = scv.read('../Data/RNAvelo/Hr30_v3Dr11.loom', cache=True)
# Hr30_Rv_data.var_names_make_unique()
# Hr31_Rv_data = scv.read('../Data/RNAvelo/Hr31_v3Dr11.loom', cache=True)
# Hr31_Rv_data.var_names_make_unique()
# Hr32_Rv_data = scv.read('../Data/RNAvelo/Hr32_v3Dr11.loom', cache=True)
# Hr32_Rv_data.var_names_make_unique()
# Hr33_Rv_data = scv.read('../Data/RNAvelo/Hr33_v3Dr11.loom', cache=True)
# Hr33_Rv_data.var_names_make_unique()
# Hr34_Rv_data = scv.read('../Data/RNAvelo/Hr34_v3Dr11.loom', cache=True)
# Hr34_Rv_data.var_names_make_unique()
# Hr35_Rv_data = scv.read('../Data/RNAvelo/Hr35_v3Dr11.loom', cache=True)
# Hr35_Rv_data.var_names_make_unique()


# In[23]:


# HR_Rv =\
#     H5_Rv_data.concatenate(H6_Rv_data, H7_Rv_data, H8a_Rv_data, H8v_Rv_data,
#                         Hr1_Rv_data, Hr2a_Rv_data, Hr2b_Rv_data, Hr3_Rv_data, Hr4_Rv_data, Hr5_Rv_data, 
#                         Hr6a_Rv_data, Hr6v_Rv_data, Hr7a_Rv_data, Hr7v_Rv_data, Hr8_Rv_data, Hr9_Rv_data, Hr10_Rv_data,
#                         Hr11_Rv_data, Hr12_Rv_data, Hr13_Rv_data, Hr14_Rv_data, Hr15_Rv_data,
#                         Hr16_Rv_data, Hr17_Rv_data, Hr18_Rv_data, Hr19_Rv_data, Hr20_Rv_data,
#                         Hr21_Rv_data, Hr22_Rv_data, Hr23_Rv_data, Hr24_Rv_data, Hr25_Rv_data,
#                         Hr26_Rv_data, Hr27_Rv_data, Hr28_Rv_data, Hr29_Rv_data, Hr30_Rv_data,
#                         Hr31_Rv_data, Hr32_Rv_data, Hr33_Rv_data, Hr34_Rv_data, Hr35_Rv_data)
# HR_Rv.shape


# In[24]:


# HR_Rv = sc.read('./write/HR_Rv.h5ad')


# In[25]:


# HR_Rv.obs = HR_Rv.obs.reset_index().merge(HR_setnames, how="inner").set_index('index')


# In[26]:


# # Rename cells to match cell names in annotation file
# HR_Rv.obs_names = [str(HR_Rv.obs.loc[x,'heart'])+'_'+str(x.split(':', 1)[1])[0:16] for x in HR_Rv.obs_names]
# # Drop annotations that are not in the single-cell object
# anno_drop_Rv = annotations.index.difference(HR_Rv.obs_names)
# annotations_Rv = annotations.drop(anno_drop_Rv)


# In[27]:


# HR_Rv_filter = HR_Rv[annotations_Rv.index]


# In[28]:


# HR_Rv_filter.obs['Cell_type'] = annotations_Rv['Cell_type'].tolist()


# In[29]:


#HR_Rv_filter.write('./write/HR_Rv_filter.h5ad')
HR_Rv_filter = sc.read('./write/HR_Rv_filter.h5ad')


# In[30]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]


# In[31]:


HR_Rv_3dpi_winh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi_winh.var_names if not name == 'RFP']
HR_Rv_3dpi_winh = HR_Rv_3dpi_winh[:, all_genes_but_RFP]


# In[32]:


HR_Rv_7dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]


# In[33]:


HR_Rv_7dpi_winh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi_winh.var_names if not name == 'RFP']
HR_Rv_7dpi_winh = HR_Rv_7dpi_winh[:, all_genes_but_RFP]


# # Naive trajectory analysis at 3dpi

# In[16]:


HR_Rv_3dpi_epifibro = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(epifibro_types)]
sc.pp.filter_genes(HR_Rv_3dpi_epifibro, min_cells=3)


# In[17]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_epifibro.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_epifibro.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_epifibro[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_epifibro.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_epifibro, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_3dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[18]:


sc.pp.normalize_total(HR_Rv_3dpi_epifibro, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi_epifibro)


# In[19]:


sc.pp.regress_out(HR_Rv_3dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[20]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_epifibro)
sc.tl.pca(HR_Rv_3dpi_epifibro)
sc.external.pp.bbknn(HR_Rv_3dpi_epifibro, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_epifibro)


# In[21]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color='batch',
          title = '3dpi epicardium and fibroblasts')


# In[22]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epifibro.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '', legend_loc = 'none', frameon=False) #save = 'subset1_3dpi_umap.png', 


# In[23]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '3dpi epicardium and fibroblasts')


# In[24]:


sc.tl.paga(HR_Rv_3dpi_epifibro, groups='Cell_type')


# In[26]:


sc.pl.paga(HR_Rv_3dpi_epifibro, threshold=0.3, show=True, frameon=False) #,
        #   labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2`)#, save = 'subset1_3dpi_paga.png')


# # Naive trajectory analysis at 7dpi

# In[14]:


epiendofibro_types = ['Fibroblasts (const.)', 'Fibroblasts (cfd)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                      'Fibroblasts (cxcl12a)', 'Fibroblasts (mpeg1.1)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)',
                      'Valve fibroblasts', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
                      'Fibroblasts (proliferating)', 'Perivascular cells', 'Endocardium (frzb)', 'Endocardium (Ventricle)', 
                      'Endocardium (Atrium)']


# In[15]:


HR_Rv_7dpi_epifibro = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(epiendofibro_types)]
sc.pp.filter_genes(HR_Rv_7dpi_epifibro, min_cells=3)


# In[16]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_epifibro.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_epifibro.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_epifibro[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_epifibro.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_epifibro, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[17]:


sc.pp.normalize_total(HR_Rv_7dpi_epifibro, target_sum=1e4)
sc.pp.log1p(HR_Rv_7dpi_epifibro)


# In[18]:


sc.pp.regress_out(HR_Rv_7dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[19]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_epifibro)
sc.tl.pca(HR_Rv_7dpi_epifibro)
sc.external.pp.bbknn(HR_Rv_7dpi_epifibro, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_epifibro)


# In[20]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color='batch',
          title = '7dpi epicardium and fibroblasts')


# In[21]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_epifibro.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '', legend_loc = 'none', frameon=False) #, save = 'subset1_7dpi_umap.png'


# In[210]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '7dpi epicardium and fibroblasts')


# In[22]:


sc.tl.paga(HR_Rv_7dpi_epifibro, groups='Cell_type')


# In[36]:


sc.pl.paga(HR_Rv_7dpi_epifibro, threshold=0.2, show=True, frameon=False,
           min_edge_width = 3, max_edge_width = 3,
           node_size_scale = 2, node_size_power = 0,
          labels = ['', '', '', '', 
                    '', '', '', '', 
                    '', '', '', '', 
                    '', '', '', ''], #), , save = 'subset1_7dpi_paga.png')
           save = 'epiendo_7dpi.pdf')


# # Trajectories in 3dpi epicardial connected niche

# In[16]:


HR_Rv_3dpi_epiconn = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(local_connected_epifibro_types)] # connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_3dpi_epiconn, min_cells=3)


# In[17]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_epiconn.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_epiconn.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_epiconn[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_epiconn.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_epiconn, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_3dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[18]:


scv.pp.normalize_per_cell(HR_Rv_3dpi_epiconn, counts_per_cell_after=1e4) #sc.pp.normalize_total(HR_Rv_3dpi_epiconn, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi_epiconn)


# In[19]:


sc.pp.regress_out(HR_Rv_3dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[20]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_epiconn)
sc.tl.pca(HR_Rv_3dpi_epiconn)
sc.external.pp.bbknn(HR_Rv_3dpi_epiconn, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_epiconn)


# In[21]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='batch',
          title = '3dpi epicardial niche')


# In[22]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[23]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color=['n_genes_by_counts'], 
          title = '3dpi epicardial niche')


# In[24]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color=['npr3'])


# In[25]:


sc.tl.leiden(HR_Rv_3dpi_epiconn, resolution=2)


# In[26]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[27]:


sc.tl.paga(HR_Rv_3dpi_epiconn, groups='leiden')


# In[28]:


sc.pl.paga(HR_Rv_3dpi_epiconn, show=True, node_size_scale = 2)


# In[29]:


sc.tl.umap(HR_Rv_3dpi_epiconn, init_pos = 'paga')


# In[30]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='batch',
          title = '3dpi epicardial niche')


# In[31]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[32]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='tcf21')


# In[33]:


sc.pl.violin(HR_Rv_3dpi_epiconn, keys = 'tcf21', groupby = 'Cell_type')


# In[34]:


scv.pp.moments(HR_Rv_3dpi_epiconn)


# In[36]:


scv.tl.velocity(HR_Rv_3dpi_epiconn, mode='stochastic')


# In[37]:


scv.tl.velocity_graph(HR_Rv_3dpi_epiconn)


# In[38]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
                                 #save = '_scvelo_local_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[39]:


# Stochastic
scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'Cell_type')#, save = '_scvelo_epicardial_niche_3dpi_20n_cell_types_umap_pub.png')


# In[40]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'leiden')


# In[38]:


scv.tl.velocity_confidence(HR_Rv_3dpi_epiconn)


# In[39]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_3dpi_epiconn, c=keys, cmap='coolwarm', perc=[5, 95])#,
#                                save = 'strength_coherence_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[41]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6', 'nrg1']
scv.pl.velocity(HR_Rv_3dpi_epiconn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
#                save = 'gene_velocity_local_epicardial_niche_3dpi_umap.png')


# # Trajectories in 3dpi epicardial connected regular and WNT-inhibited niche

# ## Both datasets on same UMAP & clusters

# In[46]:


#HR_Rv_3dpi_epiconn_winh = HR_Rv_3dpi_winh[HR_Rv_3dpi_winh.obs['Cell_type'].isin(local_connected_epifibro_types)] #connected_epifibro_types)]
HR_Rv_3dpi_epiconn_all = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3')]
HR_Rv_3dpi_epiconn_all = HR_Rv_3dpi_epiconn_all[HR_Rv_3dpi_epiconn_all.obs['Cell_type'].isin(local_connected_epifibro_types)]
all_genes_but_RFP = [name for name in HR_Rv_3dpi_epiconn_all.var_names if not name == 'RFP']
HR_Rv_3dpi_epiconn_all = HR_Rv_3dpi_epiconn_all[:, all_genes_but_RFP]
sc.pp.filter_genes(HR_Rv_3dpi_epiconn_all, min_cells=3)


# In[47]:


HR_Rv_3dpi_epiconn_all


# In[48]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_epiconn_all.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_epiconn_all.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_epiconn_all[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_epiconn_all.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_epiconn_all, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_3dpi_epiconn_all, ['total_counts', 'n_genes_by_counts', 'percent_mito'], 
             jitter=0.4, multi_panel=True, size = 0.1)


# In[49]:


sc.pp.normalize_per_cell(HR_Rv_3dpi_epiconn_all, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_3dpi_epiconn_all)


# In[50]:


sc.pp.regress_out(HR_Rv_3dpi_epiconn_all, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[51]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_epiconn_all)
sc.tl.pca(HR_Rv_3dpi_epiconn_all)
sc.external.pp.bbknn(HR_Rv_3dpi_epiconn_all, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_epiconn_all)


# In[52]:


sc.pl.umap(HR_Rv_3dpi_epiconn_all, color='batch',
          title = '3dpi epicardial niche')


# In[53]:


sc.pl.umap(HR_Rv_3dpi_epiconn_all, color='inhib',
          title = '3dpi epicardial niche')


# In[54]:


sc.pl.umap(HR_Rv_3dpi_epiconn_all, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn_all.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epicardial niche')


# In[55]:


sc.tl.leiden(HR_Rv_3dpi_epiconn_all, resolution=2)


# In[56]:


sc.pl.umap(HR_Rv_3dpi_epiconn_all, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[57]:


sc.tl.paga(HR_Rv_3dpi_epiconn_all, groups='leiden')


# In[58]:


sc.pl.paga(HR_Rv_3dpi_epiconn_all, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[59]:


sc.tl.umap(HR_Rv_3dpi_epiconn_all, init_pos = 'paga')


# ## RNA velocity on both and separate afterwards

# In[88]:


scv.pp.moments(HR_Rv_3dpi_epiconn_all)
scv.tl.velocity(HR_Rv_3dpi_epiconn_all, mode='stochastic')
scv.tl.velocity_graph(HR_Rv_3dpi_epiconn_all)


# In[115]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn_all, title = '', 
                                 color = 'Cell_type', legend_loc = 'none', size = 100,
                                save = '_scvelo_local_epicardial_niche_3dpi_all_20n_cell_types_umap.png')


# In[91]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn_all, title = 'RNA velocity epicardium 3 dpi', 
                                 color = 'inhib', legend_loc = 'right margin')


# In[114]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn_all[HR_Rv_3dpi_epiconn_all.obs['inhib'] != 'IWR1'], title = '', 
                                 color = 'Cell_type', legend_loc = 'none', size = 100,
                                save = '_scvelo_local_epicardial_niche_3dpi_all_ctrl_20n_cell_types_umap.png')


# In[116]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn_all[HR_Rv_3dpi_epiconn_all.obs['inhib'] == 'IWR1'], title = '', 
                                 color = 'Cell_type', legend_loc = 'none', size = 100,
                                 save = '_scvelo_local_epicardial_niche_3dpi_all_winh_20n_cell_types_umap.png')


# ## Separate control and Wnt-inhibited and do RNA velocity for both separately

# In[60]:


HR_Rv_3dpi_epiconn_all_subctrl = HR_Rv_3dpi_epiconn_all[HR_Rv_3dpi_epiconn_all.obs['inhib'] != 'IWR1']
HR_Rv_3dpi_epiconn_all_subwinh = HR_Rv_3dpi_epiconn_all[HR_Rv_3dpi_epiconn_all.obs['inhib'] == 'IWR1']


# In[68]:


sc.pl.umap(HR_Rv_3dpi_epiconn_all, color='inhib')


# In[67]:


sc.pl.umap(HR_Rv_3dpi_epiconn_all_subctrl, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn_all_subctrl.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          size = 19, title = '3dpi epicardial niche, only control')


# In[66]:


sc.pl.umap(HR_Rv_3dpi_epiconn_all_subwinh, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn_all_subwinh.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          size = 19, title = '3dpi epicardial niche, only Wnt-inhibited')


# In[69]:


scv.pp.moments(HR_Rv_3dpi_epiconn_all_subctrl)
scv.tl.velocity(HR_Rv_3dpi_epiconn_all_subctrl, mode='stochastic')
scv.tl.velocity_graph(HR_Rv_3dpi_epiconn_all_subctrl)


# In[70]:


scv.pp.moments(HR_Rv_3dpi_epiconn_all_subwinh)
scv.tl.velocity(HR_Rv_3dpi_epiconn_all_subwinh, mode='stochastic')
scv.tl.velocity_graph(HR_Rv_3dpi_epiconn_all_subwinh)


# In[71]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn_all_subctrl, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')


# In[72]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn_all_subwinh, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')


# In[29]:


scv.pp.moments(HR_Rv_3dpi_epiconn_winh)


# In[30]:


scv.tl.velocity(HR_Rv_3dpi_epiconn_winh, mode='stochastic')


# In[31]:


scv.tl.velocity_graph(HR_Rv_3dpi_epiconn_winh)


# In[34]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn_winh, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
                                 #save = '_scvelo_local_epicardial_niche_3dpi_winhib_20n_cell_types_umap.png')


# # Trajectories in 7dpi epicardial connected niche

# In[53]:


HR_Rv_7dpi_epiconn = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(local_connected_epifibro_types)] #connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_7dpi_epiconn, min_cells=3)


# In[54]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_epiconn.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_epiconn.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_epiconn[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_epiconn.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_epiconn, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[55]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_epiconn, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_epiconn)


# In[56]:


sc.pp.regress_out(HR_Rv_7dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[57]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_epiconn)
sc.tl.pca(HR_Rv_7dpi_epiconn)
sc.external.pp.bbknn(HR_Rv_7dpi_epiconn, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_epiconn)


# In[58]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='batch',
          title = '7dpi epicardial niche')


# In[59]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color=['total_counts', 'n_genes_by_counts', 'percent_mito'],
          title = '7dpi epicardial niche')


# In[60]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epicardial niche')


# In[61]:


sc.tl.leiden(HR_Rv_7dpi_epiconn, resolution=2)


# In[62]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[63]:


sc.tl.paga(HR_Rv_7dpi_epiconn, groups='leiden')


# In[64]:


sc.pl.paga(HR_Rv_7dpi_epiconn, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[65]:


sc.tl.umap(HR_Rv_7dpi_epiconn, init_pos = 'paga')


# In[66]:


scv.pp.moments(HR_Rv_7dpi_epiconn)


# In[67]:


scv.tl.velocity(HR_Rv_7dpi_epiconn, mode='stochastic')


# In[68]:


scv.tl.velocity_graph(HR_Rv_7dpi_epiconn)


# In[69]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_epiconn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
                                 #save = '_scvelo_local_epicardial_niche_7dpi_noinhib_20n_cell_types_umap.png')


# In[71]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_epiconn, title = '', 
                                 color = 'leiden')


# In[70]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6','nrg1']
scv.pl.velocity(HR_Rv_7dpi_epiconn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
                #save = 'gene_velocity_local_epicardial_niche_7dpi_noinhib_umap.png')


# # Trajectories in 7dpi connected WNT-inhibited epicardium

# In[111]:


HR_Rv_7dpi_epiconn_winh = HR_Rv_7dpi_winh[HR_Rv_7dpi_winh.obs['Cell_type'].isin(local_connected_epifibro_types)] #connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_7dpi_epiconn_winh, min_cells=3)


# In[112]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_epiconn_winh.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_epiconn_winh.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_epiconn_winh[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_epiconn_winh.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_epiconn_winh, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_epiconn_winh, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[113]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_epiconn_winh, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_epiconn_winh)


# In[114]:


sc.pp.regress_out(HR_Rv_7dpi_epiconn_winh, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[115]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_epiconn_winh)
sc.tl.pca(HR_Rv_7dpi_epiconn_winh)
sc.external.pp.bbknn(HR_Rv_7dpi_epiconn_winh, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_epiconn_winh)


# In[116]:


sc.pl.umap(HR_Rv_7dpi_epiconn_winh, color='batch',
          title = '7dpi epicardial niche')


# In[117]:


sc.pl.umap(HR_Rv_7dpi_epiconn_winh, color=['total_counts', 'n_genes_by_counts', 'percent_mito'],
          title = '7dpi epicardial niche')


# In[118]:


sc.pl.umap(HR_Rv_7dpi_epiconn_winh, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_epiconn_winh.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epicardial niche')


# In[119]:


sc.tl.leiden(HR_Rv_7dpi_epiconn_winh, resolution=2)


# In[120]:


sc.pl.umap(HR_Rv_7dpi_epiconn_winh, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[121]:


sc.tl.paga(HR_Rv_7dpi_epiconn_winh, groups='leiden')


# In[122]:


sc.pl.paga(HR_Rv_7dpi_epiconn_winh, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[123]:


sc.tl.umap(HR_Rv_7dpi_epiconn_winh, init_pos = 'paga')


# In[124]:


scv.pp.moments(HR_Rv_7dpi_epiconn_winh)


# In[125]:


scv.tl.velocity(HR_Rv_7dpi_epiconn_winh, mode='stochastic')


# In[126]:


scv.tl.velocity_graph(HR_Rv_7dpi_epiconn_winh)


# In[130]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_epiconn_winh, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
                                 #save = '_scvelo_local_epicardial_niche_7dpi_winhib_20n_cell_types_umap.png')


# In[128]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6','nrg1']
scv.pl.velocity(HR_Rv_7dpi_epiconn_winh, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
                #save = 'gene_velocity_local_epicardial_niche_7dpi_noinhib_umap.png')


# # Endocardial niche for deep injuries at 7dpi

# In[16]:


all_endocardial_lineage = ['Endocardium (frzb)', 'Endocardium (Ventricle)', 
                            'Fibroblasts (nppc)']


# In[17]:


deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9'] # All samples with > 50 nppc fibroblasts.
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(all_endocardial_lineage)]#connected_endofibro_types)]
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['heart'].isin(deep_injury_libraries)]
scv.pp.remove_duplicate_cells(HR_Rv_7dpi_deep_endo)
sc.pp.filter_genes(HR_Rv_7dpi_deep_endo, min_cells=3)


# In[18]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_deep_endo.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_deep_endo.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_deep_endo[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_deep_endo.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_deep_endo, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[19]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_deep_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_deep_endo)


# In[20]:


sc.pp.regress_out(HR_Rv_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[21]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_deep_endo)
sc.tl.pca(HR_Rv_7dpi_deep_endo)
sc.external.pp.bbknn(HR_Rv_7dpi_deep_endo, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_deep_endo)


# In[22]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi deep injury endocardial niche')


# In[23]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '7dpi deep injury endocardial niche')


# In[24]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='batch',
          title = '7dpi deep injury endocardial niche')


# In[26]:


sc.tl.leiden(HR_Rv_7dpi_deep_endo, resolution=2)


# In[ ]:


sc.tl.paga(HR_Rv_7dpi_deep_endo, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_7dpi_deep_endo, show=True, color = 'Cell_type', node_size_scale = 2)#, 
    #       labels = ['', '', '', '', '', 
   #                  '', '', '', '', '', 
   #                  '', '', '', '', '', '', ''], 
  #         node_size_scale = 2)#, save = 'scvelo_deep_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[ ]:


sc.tl.umap(HR_Rv_7dpi_deep_endo, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'leiden')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'batch')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'Cell_type')


# In[ ]:


endmt_genes = ["acta2", "cdh5", "pecam1", 
                 "vwf", "cdh2", "vim", "icn", "fn1a",
                 "fn1b", "fap", "cnn1b", "tagln",
                 "vcana", "vcanb", "snai1a", "snai1b", "snai2", "snai3"]


# In[73]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = ['pecam1', 'cdh5', 'nppc'])


# In[74]:


scv.pp.moments(HR_Rv_7dpi_deep_endo)


# In[75]:


scv.tl.velocity(HR_Rv_7dpi_deep_endo, mode='stochastic')


# In[76]:


scv.tl.velocity_graph(HR_Rv_7dpi_deep_endo)


# In[78]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, basis='umap', title = '', 
                                 color = 'Cell_type',
                                save = 'gene_velocity_local_endocardial_deep_niche_7dpi_umap.png')


# In[42]:


scv.pl.velocity_embedding(HR_Rv_7dpi_deep_endo, basis='umap', title = '', scale = 0.3,
                          color = 'Cell_type', legend_loc = 'none')#,
#                         save = 'gene_velocity_endoventricular_deep_niche_7dpi_umap.png')


# In[ ]:


scv.tl.velocity_confidence(HR_Rv_7dpi_deep_endo)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_deep_endo, c=keys, cmap='coolwarm', perc=[5, 95])#,
#                                save = 'strength_coherence_scvelo_deep_endo_niche_7dpi_20n_cell_types_umap.png')


# # Endocardial niche inhibited and uninhibited at 7 dpi

# In[98]:


HR_Rv_7dpi_endoconn_all = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7')]
HR_Rv_7dpi_endoconn_all = HR_Rv_7dpi_endoconn_all[HR_Rv_7dpi_endoconn_all.obs['Cell_type'].isin(all_endocardial_lineage)]
all_genes_but_RFP = [name for name in HR_Rv_7dpi_endoconn_all.var_names if not name == 'RFP']
HR_Rv_7dpi_endoconn_all = HR_Rv_7dpi_endoconn_all[:, all_genes_but_RFP]
sc.pp.filter_genes(HR_Rv_7dpi_endoconn_all, min_cells=3)


# In[99]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_endoconn_all.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_endoconn_all.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_endoconn_all[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_endoconn_all.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_endoconn_all, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_endoconn_all, ['total_counts', 'n_genes_by_counts', 'percent_mito'], 
             jitter=0.4, multi_panel=True, size = 0.1)


# In[100]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_endoconn_all, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_endoconn_all)
sc.pp.regress_out(HR_Rv_7dpi_endoconn_all, ['total_counts', 'n_genes_by_counts', 'percent_mito'])
sc.pp.highly_variable_genes(HR_Rv_7dpi_endoconn_all)
sc.tl.pca(HR_Rv_7dpi_endoconn_all)
sc.external.pp.bbknn(HR_Rv_7dpi_endoconn_all, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_endoconn_all)
sc.tl.leiden(HR_Rv_7dpi_endoconn_all, resolution=2)
sc.tl.paga(HR_Rv_7dpi_endoconn_all, groups='leiden')


# In[101]:


sc.pl.paga(HR_Rv_7dpi_endoconn_all, show=True, color = 'Cell_type', node_size_scale = 2)
sc.tl.umap(HR_Rv_7dpi_endoconn_all, init_pos = 'paga')


# In[102]:


sc.pl.umap(HR_Rv_7dpi_endoconn_all, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_endoconn_all.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi endocardial niche')


# ## RNA velocity on both and separate afterwards

# In[103]:


scv.pp.moments(HR_Rv_7dpi_endoconn_all)
scv.tl.velocity(HR_Rv_7dpi_endoconn_all, mode='stochastic')
scv.tl.velocity_graph(HR_Rv_7dpi_endoconn_all)


# In[119]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endoconn_all, title = '', 
                                 color = 'Cell_type', legend_loc = 'none', size = 100,
                                save = '_scvelo_local_endocardial_niche_7dpi_all_20n_cell_types_umap.png')


# In[120]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endoconn_all[HR_Rv_7dpi_endoconn_all.obs['inhib'] != 'IWR1'], title = '', 
                                 color = 'Cell_type', legend_loc = 'none', size = 100,
                                save = '_scvelo_local_endocardial_niche_7dpi_all_ctrl_20n_cell_types_umap.png')


# In[121]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endoconn_all[HR_Rv_7dpi_endoconn_all.obs['inhib'] == 'IWR1'], title = '', 
                                 color = 'Cell_type', legend_loc = 'none', size = 100,
                                save = '_scvelo_local_endocardial_niche_7dpi_all_winh_20n_cell_types_umap.png')


# # Secretome analysis

# How much of the transcriptome is part of the secretome for each cell type?

# ## Load secretome genes

# In[34]:


secretome = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Alliance_secretome_gene_names_noDRduplicates.scsv', sep = ';')
secretome = secretome.rename(columns={'DR_name': 'external_gene_name'})
secretome = secretome[secretome['external_gene_name'].notna()]
secretome = (secretome[['external_gene_name']]).drop_duplicates()


# ## Calculate differentially expressed genes per timepoint

# In[35]:


HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '0') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)
sc.pp.log1p(HR_Rv_ctrl)
sc.tl.rank_genes_groups(HR_Rv_ctrl, groupby = 'Cell_type')
dg_ctrl = GetDiffGenes(HR_Rv_ctrl)


# In[36]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi)
sc.tl.rank_genes_groups(HR_Rv_3dpi, groupby = 'Cell_type')
dg_3dpi = GetDiffGenes(HR_Rv_3dpi)


# In[40]:


HR_Rv_7dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi, target_sum=1e4)
sc.pp.log1p(HR_Rv_7dpi)
sc.tl.rank_genes_groups(HR_Rv_7dpi, groupby = 'Cell_type')
dg_7dpi = GetDiffGenes(HR_Rv_7dpi)


# ## Annotate genes with secretome and count numbers differentially expressed

# In[41]:


dg_ctrl_counts = pd.DataFrame({'DG_count' : dg_ctrl.Cell_type.value_counts()})
diff_secretome_ctrl = dg_ctrl[dg_ctrl.Gene.isin(secretome.external_gene_name)]
ds_ctrl_counts = pd.DataFrame({'DSec_count' : diff_secretome_ctrl.Cell_type.value_counts()})
diff_ctrl_counts = pd.concat([dg_ctrl_counts, ds_ctrl_counts], axis = 1, join='inner')
diff_ctrl_counts['DSec_ratio'] = diff_ctrl_counts['DSec_count']/diff_ctrl_counts['DG_count']


# In[42]:


dg_3dpi_counts = pd.DataFrame({'DG_count' : dg_3dpi.Cell_type.value_counts()})
diff_secretome_3dpi = dg_3dpi[dg_3dpi.Gene.isin(secretome.external_gene_name)]
ds_3dpi_counts = pd.DataFrame({'DSec_count' : diff_secretome_3dpi.Cell_type.value_counts()})
diff_3dpi_counts = pd.concat([dg_3dpi_counts, ds_3dpi_counts], axis = 1, join='inner')
diff_3dpi_counts['DSec_ratio'] = diff_3dpi_counts['DSec_count']/diff_3dpi_counts['DG_count']


# In[43]:


dg_7dpi_counts = pd.DataFrame({'DG_count' : dg_7dpi.Cell_type.value_counts()})
diff_secretome_7dpi = dg_7dpi[dg_7dpi.Gene.isin(secretome.external_gene_name)]
ds_7dpi_counts = pd.DataFrame({'DSec_count' : diff_secretome_7dpi.Cell_type.value_counts()})
diff_7dpi_counts = pd.concat([dg_7dpi_counts, ds_7dpi_counts], axis = 1, join='inner')
diff_7dpi_counts['DSec_ratio'] = diff_7dpi_counts['DSec_count']/diff_7dpi_counts['DG_count']


# ## Secretome expression per cell type

# In[44]:


HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '0') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)
secretome_ctrl = CountSecretome(HR_Rv_ctrl, secretome)
secretome_averages_ctrl, secretome_z_averages_ctrl = ExpressedSecretome(HR_Rv_ctrl, secretome)


# In[30]:


#secretome_ctrl


# In[21]:


#secretome_old_ctrl = CountSecretome(HR_Rv_ctrl, secretome_old)
#secretome_old_averages_ctrl, secretome_old_z_averages_ctrl = ExpressedSecretome(HR_Rv_ctrl, secretome_old)


# In[108]:


#secretome_z_averages_ctrl


# In[109]:


#secretome_old_z_averages_ctrl


# In[45]:


secretome_averages_ctrl_save = secretome_averages_ctrl[secretome_averages_ctrl.max(axis = 1) > 10]
#secretome_averages_ctrl_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_ctrl.csv')
secretome_z_averages_ctrl_save = secretome_z_averages_ctrl[secretome_z_averages_ctrl.max(axis = 1) > 2]
#secretome_z_averages_ctrl_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_ctrl.csv')


# In[27]:


#secretome_old_averages_ctrl_save = secretome_old_averages_ctrl[secretome_old_averages_ctrl.max(axis = 1) > 10]
#secretome_old_z_averages_ctrl_save = secretome_old_z_averages_ctrl[secretome_old_z_averages_ctrl.max(axis = 1) > 2]


# In[110]:


#secretome_z_averages_ctrl_save


# In[111]:


#secretome_old_z_averages_ctrl_save


# In[112]:


#import matplotlib.pyplot as plt

# secretome_ctrl_plot = secretome_ctrl[secretome_ctrl.Cell_type.isin(epifibro_types)] 
# secretome_ctrl_plot = secretome_ctrl_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
# secretome_ctrl_plot = secretome_ctrl_plot.sort_values(by = 'Secretome', ascending = False)

# pl.bar(x = np.arange(len(secretome_ctrl_plot)),
#        height = secretome_ctrl_plot.Secretome/100,
#        yerr = 3 * secretome_ctrl_plot.SEM/100, capsize = 2,
#        color = secretome_ctrl_plot.color)
# pl.xlabel('Cell type')
# pl.xticks(np.arange(len(secretome_ctrl_plot)), secretome_ctrl_plot.Cell_type, rotation=270)
# pl.ylabel('Secretome (%)')
# pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_ctrl_fibroniche.pdf', bbox_inches = 'tight', transparent = True)
# pl.show()


# In[ ]:


# START PLOT IMPROVEMENT


# In[31]:


# HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '0') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
# all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
# HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
# sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)
#secretome_ctrl = CountSecretome(HR_Rv_ctrl, secretome)


# In[47]:


secretome_ind = [HR_Rv_ctrl.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in HR_Rv_ctrl.var.index]
secretome_df = HR_Rv_ctrl.obs.loc[:, ['Cell_type']]
secretome_df['Secretome'] = np.nan
secretome_slice = HR_Rv_ctrl.X[:,secretome_ind]
for ctype in HR_Rv_ctrl.obs['Cell_type'].cat.categories:
    ctype_index = HR_Rv_ctrl.obs[HR_Rv_ctrl.obs['Cell_type'] == ctype].index
    ctype_index_numeric = [x for x in range(0, len(HR_Rv_ctrl.obs)) if HR_Rv_ctrl.obs.Cell_type[x] == ctype]
    if(len(ctype_index_numeric) == 0):
        continue
    #print(len(ctype_index_numeric))
    secretome_df.Secretome.loc[ctype_index] = np.array((secretome_slice[ctype_index_numeric, :].sum(axis = 1)))[:,0]/100


# In[126]:


# for ctype in HR_Rv_ctrl.obs['Cell_type'].cat.categories:
#     ctype_index = HR_Rv_ctrl.obs[HR_Rv_ctrl.obs['Cell_type'] == ctype].index
#     ctype_index_numeric = [x for x in range(0, len(HR_Rv_ctrl.obs)) if HR_Rv_ctrl.obs.Cell_type[x] == ctype]
#     if(len(ctype_index_numeric) == 0):
#         continue
#     #print(len(ctype_index_numeric))
#     secretome_df.Secretome.loc[ctype_index] = np.array((secretome_slice[ctype_index_numeric, :].sum(axis = 1)))[:,0]/100


# In[164]:


# celltype_order = (secretome_ctrl.sort_values(by = 'Secretome', ascending = False)).Cell_type
# fig, ax = pl.subplots(figsize=(12,4))
# g = sns.violinplot(ax=ax, x="Cell_type", y="Secretome", data=secretome_df, scale = 'width', inner=None,
#               order = celltype_order, palette = cell_type_colors.loc[celltype_order].color)
# g.set_xticklabels(g.get_xticklabels(), rotation=270)
# pl.show()


# In[163]:


# celltype_order = (secretome_ctrl.sort_values(by = 'Secretome', ascending = False)).Cell_type
# fig, ax = pl.subplots(figsize=(12,4))
# g = sns.pointplot(ax=ax, x="Cell_type", y="Secretome", data=secretome_df,
#                   capsize=.4, ci=99, join=False,
#               order = celltype_order, color = 'black')#, palette = cell_type_colors.loc[celltype_order].color)
# #ax = sns.pointplot(x="day", y="tip", data=tips, ci=68)
# g.set_xticklabels(g.get_xticklabels(), rotation=270)
# pl.show()


# In[48]:


celltype_order = (secretome_ctrl.sort_values(by = 'Secretome', ascending = False)).Cell_type
fig, ax = pl.subplots(figsize=(12,4))
g = sns.violinplot(ax=ax, x="Cell_type", y="Secretome", data=secretome_df, scale = 'width', inner=None,
              order = celltype_order, palette = cell_type_colors.loc[celltype_order].color)
pl.setp(ax.collections, alpha=.5)
g = sns.pointplot(ax=g, x="Cell_type", y="Secretome", data=secretome_df,
                  capsize=.4, ci=99, join=False,
              order = celltype_order, color = 'black')
g.set_xticklabels(g.get_xticklabels(), rotation=270)
g.spines['top'].set_visible(False)
g.spines['right'].set_visible(False)
pl.show()


# In[ ]:


# END PLOT IMPROVEMENTS


# In[19]:


secretome_averages_ctrl_plot = np.log1p(secretome_averages_ctrl[secretome_averages_ctrl.max(axis = 1) > 10])
secretome_averages_ctrl_plot = secretome_averages_ctrl_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)         


# In[21]:


sns.clustermap(secretome_averages_ctrl_plot, method = 'ward', xticklabels = True, yticklabels = True)
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_ctrl_high_genes.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[22]:


secretome_z_averages_ctrl_plot = secretome_z_averages_ctrl[secretome_z_averages_ctrl.max(axis = 1) > 2]
secretome_z_averages_ctrl_plot[secretome_z_averages_ctrl_plot > 5] = 5
secretome_z_averages_ctrl_plot = secretome_z_averages_ctrl_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)

sns.clustermap(secretome_z_averages_ctrl_plot[secretome_z_averages_ctrl_plot.max(axis = 1) > 2], method = 'ward', xticklabels = True, yticklabels = True, figsize = (10,25))
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_ctrl_high_z.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[18]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi, target_sum=1e4)
secretome_3dpi = CountSecretome(HR_Rv_3dpi, secretome)
secretome_averages_3dpi, secretome_z_averages_3dpi = ExpressedSecretome(HR_Rv_3dpi, secretome)


# In[103]:


#secretome_old_3dpi = CountSecretome(HR_Rv_3dpi, secretome_old)
#secretome_old_averages_3dpi, secretome_old_z_averages_3dpi = ExpressedSecretome(HR_Rv_3dpi, secretome_old)


# In[94]:


#adata = HR_Rv_3dpi
#secretome_ind = [adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index]
#secretome_ind = np.unique([adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index])

#secretome_slice = adata.X[:,secretome_ind]
#secretome_dense_slice = adata.X[:,secretome_ind].todense()
#secretome_dense_slice_z = preprocessing.scale(secretome_dense_slice)
#secretome_dense_slice_z
#ctype_averages = pd.DataFrame(index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index])
#ctype_z_averages = pd.DataFrame(index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index])
#ctype_averages
#for ctype in adata.obs['Cell_type'].cat.categories:
#    ctype_index = [x for x in range(0, len(adata.obs) - 1) if adata.obs.Cell_type[x] == ctype]
#    if(len(ctype_index) == 0):
#        continue
#    ctype_average = pd.DataFrame({ctype : np.squeeze(np.asarray(secretome_slice[ctype_index, :].mean(axis = 0)))},
#                                  index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index]) # axis = 0 gets us gene-level averages. #.sum(axis = 1))/len(ctype_index))[0,0]
#    ctype_averages = ctype_averages.join(ctype_average)
#        
#    ctype_z_average = pd.DataFrame({ctype : secretome_dense_slice_z[ctype_index, :].mean(axis = 0)},
#                                    index = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index])
#    ctype_z_averages = ctype_z_averages.join(ctype_z_average)


# In[114]:


#len(secretome_ind)


# In[115]:


#secretome_dense_slice.shape


# In[116]:


#sec_genes = [x for x in np.array(secretome.external_gene_name) if x in adata.var.index]


# In[117]:


#sec_genes[1:10]


# In[118]:


#adata.var.loc[sec_genes[1:10]]


# In[119]:


#raw_sec_index = [adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index]


# In[120]:


#raw_sec_index[1:10]


# In[121]:


#adata.var.iloc[raw_sec_index[1:10]]


# In[122]:


#adata.var.iloc[secretome_ind[1:10]]


# In[123]:


#adata.var


# In[124]:


#adata.var.iloc[adata.var.index.get_loc('adprhl1')]


# In[125]:


#secretome_old_ind[1:10]


# In[126]:


#[x for x in np.array(secretome_old.external_gene_name) if x in adata.var.index]


# In[96]:


#adata = HR_Rv_3dpi
#secretome_old_ind = [adata.var.index.get_loc(x) for x in np.array(secretome_old.external_gene_name) if x in adata.var.index]
#secretome_old_ind = np.unique([adata.var.index.get_loc(x) for x in np.array(secretome_old.external_gene_name) if x in adata.var.index])

#secretome_old_slice = adata.X[:,secretome_old_ind]
#secretome_old_dense_slice = adata.X[:,secretome_old_ind].todense()
#secretome_old_dense_slice_z = preprocessing.scale(secretome_old_dense_slice)

#ctype_old_averages = pd.DataFrame(index = [x for x in np.array(secretome_old.external_gene_name) if x in adata.var.index])
#ctype_old_z_averages = pd.DataFrame(index = [x for x in np.array(secretome_old.external_gene_name) if x in adata.var.index])
#ctype_old_averages
#for ctype_old in adata.obs['Cell_type'].cat.categories:
#    ctype_old_index = [x for x in range(0, len(adata.obs) - 1) if adata.obs.Cell_type[x] == ctype_old]
#    if(len(ctype_old_index) == 0):
#        continue
#    ctype_old_average = pd.DataFrame({ctype_old : np.squeeze(np.asarray(secretome_old_slice[ctype_old_index, :].mean(axis = 0)))},
#                                  index = [x for x in np.array(secretome_old.external_gene_name) if x in adata.var.index]) # axis = 0 gets us gene-level averages. #.sum(axis = 1))/len(ctype_old_index))[0,0]
#    ctype_old_averages = ctype_old_averages.join(ctype_old_average)
        
#    ctype_old_z_average = pd.DataFrame({ctype_old : secretome_old_dense_slice_z[ctype_old_index, :].mean(axis = 0)},
#                                    index = [x for x in np.array(secretome_old.external_gene_name) if x in adata.var.index])
#    ctype_old_z_averages = ctype_old_z_averages.join(ctype_old_z_average)


# In[127]:


#ctype_averages.loc['c1qa']


# In[128]:


#ctype_old_averages.loc['c1qa']


# In[129]:


#secretome_z_averages_3dpi


# In[130]:


#secretome_old_z_averages_3dpi


# In[35]:


#secretome_old_averages_3dpi_save = secretome_old_averages_3dpi[secretome_old_averages_3dpi.max(axis = 1) > 10]
#secretome_old_z_averages_3dpi_save = secretome_old_z_averages_3dpi[secretome_old_z_averages_3dpi.max(axis = 1) > 2]


# In[131]:


secretome_averages_3dpi_save = secretome_averages_3dpi[secretome_averages_3dpi.max(axis = 1) > 10]
secretome_averages_3dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_3dpi.csv')
secretome_z_averages_3dpi_save = secretome_z_averages_3dpi[secretome_z_averages_3dpi.max(axis = 1) > 2]
secretome_z_averages_3dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_3dpi.csv')


# In[132]:


#secretome_averages_3dpi.loc['c1qa']


# In[133]:


#secretome_averages_3dpi.loc['c1qa']


# In[134]:


#secretome_averages_3dpi.loc['tuft1b']


# In[135]:


#secretome_old_averages_3dpi.loc['tuft1b']


# In[136]:


#secretome_old_averages_3dpi.loc['tuft1b']


# In[137]:


import matplotlib.pyplot as plt

secretome_3dpi_plot = secretome_3dpi[secretome_3dpi.Cell_type.isin(epifibro_types)]
secretome_3dpi_plot = secretome_3dpi_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_3dpi_plot = secretome_3dpi_plot.sort_values(by = 'Secretome', ascending = False)

plt.bar(x = np.arange(len(secretome_3dpi_plot)),
       height = secretome_3dpi_plot.Secretome/100,
       yerr = 3 * secretome_3dpi_plot.SEM/100, capsize = 2,
       color = secretome_3dpi_plot.color)
plt.xlabel('Cell type')
plt.xticks(np.arange(len(secretome_3dpi_plot)), secretome_3dpi_plot.Cell_type, rotation=270)
plt.ylabel('Secretome (%)')
plt.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_fibroniche.pdf', bbox_inches = 'tight', transparent = True)
plt.show()


# In[196]:


secretome_averages_3dpi_plot = np.log1p(secretome_averages_3dpi[secretome_averages_3dpi.max(axis = 1) > 10])
secretome_averages_3dpi_plot = secretome_averages_3dpi_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)         


# In[215]:


sns.clustermap(secretome_averages_3dpi_plot, method = 'ward', xticklabels = True, yticklabels = True)
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_high_genes.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[216]:


#secretome_z_averages_3dpi_plot = secretome_z_averages_3dpi
secretome_z_averages_3dpi_plot = secretome_z_averages_3dpi[secretome_z_averages_3dpi.max(axis = 1) > 2]
secretome_z_averages_3dpi_plot[secretome_z_averages_3dpi_plot > 5] = 5
secretome_z_averages_3dpi_plot = secretome_z_averages_3dpi_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)

sns.clustermap(secretome_z_averages_3dpi_plot[secretome_z_averages_3dpi_plot.max(axis = 1) > 2], method = 'ward', xticklabels = True, yticklabels = True, figsize = (10,25))
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_high_z.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[24]:


HR_Rv_7dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi, target_sum=1e4)
secretome_7dpi = CountSecretome(HR_Rv_7dpi, secretome)
secretome_averages_7dpi, secretome_z_averages_7dpi = ExpressedSecretome(HR_Rv_7dpi, secretome)


# In[141]:


secretome_averages_7dpi_save = secretome_averages_7dpi[secretome_averages_7dpi.max(axis = 1) > 10]
secretome_averages_7dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_7dpi.csv')
secretome_z_averages_7dpi_save = secretome_z_averages_7dpi[secretome_z_averages_7dpi.max(axis = 1) > 2]
secretome_z_averages_7dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_7dpi.csv')


# In[142]:


import matplotlib.pyplot as plt

secretome_7dpi_plot = secretome_7dpi #[secretome_7dpi.Cell_type.isin(epifibro_types)]
secretome_7dpi_plot = secretome_7dpi_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_7dpi_plot = secretome_7dpi_plot.sort_values(by = 'Secretome', ascending = False)

pl.bar(x = np.arange(len(secretome_7dpi_plot)),
       height = secretome_7dpi_plot.Secretome/100,
       yerr = 3 * secretome_7dpi_plot.SEM/100, capsize = 2,
       color = secretome_7dpi_plot.color)
pl.xlabel('Cell type')
pl.xticks(np.arange(len(secretome_7dpi_plot)), secretome_7dpi_plot.Cell_type, rotation=270)
pl.ylabel('Secretome (%)')
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_7dpi.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[25]:


secretome_averages_7dpi_plot = np.log1p(secretome_averages_7dpi[secretome_averages_7dpi.max(axis = 1) > 10])
secretome_averages_7dpi_plot = secretome_averages_7dpi_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)         


# In[28]:


sns.clustermap(secretome_averages_7dpi_plot, method = 'ward', xticklabels = True, yticklabels = True)
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_7dpi_high_genes.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[30]:


secretome_z_averages_7dpi_plot = secretome_z_averages_7dpi[secretome_z_averages_7dpi.max(axis = 1) > 2]
secretome_z_averages_7dpi_plot[secretome_z_averages_7dpi_plot > 5] = 5
secretome_z_averages_7dpi_plot = secretome_z_averages_7dpi_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)

sns.clustermap(secretome_z_averages_7dpi_plot[secretome_z_averages_7dpi_plot.max(axis = 1) > 2], method = 'ward', xticklabels = True, yticklabels = True, figsize = (10,25))
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_7dpi_high_z.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[26]:


HR_Rv_3dpi_inh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi_inh.var_names if not name == 'RFP']
HR_Rv_3dpi_inh = HR_Rv_3dpi_inh[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi_inh, target_sum=1e4)
secretome_3dpi_inh = CountSecretome(HR_Rv_3dpi_inh, secretome)
secretome_averages_3dpi_inh, secretome_z_averages_3dpi_inh = ExpressedSecretome(HR_Rv_3dpi_inh, secretome)


# In[144]:


secretome_averages_3dpi_inh_save = secretome_averages_3dpi_inh[secretome_averages_3dpi_inh.max(axis = 1) > 10]
secretome_averages_3dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_3dpi_inh.csv')
secretome_z_averages_3dpi_inh_save = secretome_z_averages_3dpi_inh[secretome_z_averages_3dpi_inh.max(axis = 1) > 2]
secretome_z_averages_3dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_3dpi_inh.csv')


# In[145]:


import matplotlib.pyplot as plt

secretome_3dpi_inh_plot = secretome_3dpi_inh #[secretome_3dpi_inh.Cell_type.isin(epifibro_types)]
secretome_3dpi_inh_plot = secretome_3dpi_inh_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_3dpi_inh_plot = secretome_3dpi_inh_plot.sort_values(by = 'Secretome', ascending = False)

plt.bar(x = np.arange(len(secretome_3dpi_inh_plot)),
       height = secretome_3dpi_inh_plot.Secretome/100,
       yerr = 3 * secretome_3dpi_inh_plot.SEM/100, capsize = 2,
       color = secretome_3dpi_inh_plot.color)
plt.xlabel('Cell type')
plt.xticks(np.arange(len(secretome_3dpi_inh_plot)), secretome_3dpi_inh_plot.Cell_type, rotation=270)
plt.ylabel('Secretome (%)')
plt.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_inh.pdf', bbox_inches = 'tight', transparent = True)
plt.show()


# In[31]:


secretome_averages_3dpi_inh_plot = np.log1p(secretome_averages_3dpi_inh[secretome_averages_3dpi_inh.max(axis = 1) > 10])
secretome_averages_3dpi_inh_plot = secretome_averages_3dpi_inh_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)         


# In[32]:


sns.clustermap(secretome_averages_3dpi_inh_plot, method = 'ward', xticklabels = True, yticklabels = True)
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_inh_high_genes.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[33]:


secretome_z_averages_3dpi_inh_plot = secretome_z_averages_3dpi_inh[secretome_z_averages_3dpi_inh.max(axis = 1) > 2]
secretome_z_averages_3dpi_inh_plot[secretome_z_averages_3dpi_inh_plot > 5] = 5
secretome_z_averages_3dpi_inh_plot = secretome_z_averages_3dpi_inh_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)

sns.clustermap(secretome_z_averages_3dpi_inh_plot[secretome_z_averages_3dpi_inh_plot.max(axis = 1) > 2], method = 'ward', xticklabels = True, yticklabels = True, figsize = (10,25))
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_inh_high_z.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[27]:


HR_Rv_7dpi_inh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi_inh.var_names if not name == 'RFP']
HR_Rv_7dpi_inh = HR_Rv_7dpi_inh[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi_inh, target_sum=1e4)
secretome_7dpi_inh = CountSecretome(HR_Rv_7dpi_inh, secretome)
secretome_averages_7dpi_inh, secretome_z_averages_7dpi_inh = ExpressedSecretome(HR_Rv_7dpi_inh, secretome)


# In[147]:


secretome_averages_7dpi_inh_save = secretome_averages_7dpi_inh[secretome_averages_7dpi_inh.max(axis = 1) > 10]
secretome_averages_7dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_7dpi_inh.csv')
secretome_z_averages_7dpi_inh_save = secretome_z_averages_7dpi_inh[secretome_z_averages_7dpi_inh.max(axis = 1) > 2]
secretome_z_averages_7dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_7dpi_inh.csv')


# In[214]:


import matplotlib.pyplot as plt

secretome_7dpi_inh_plot = secretome_7dpi_inh #[secretome_7dpi_inh.Cell_type.isin(epifibro_types)]
secretome_7dpi_inh_plot = secretome_7dpi_inh_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_7dpi_inh_plot = secretome_7dpi_inh_plot.sort_values(by = 'Secretome', ascending = False)

plt.bar(x = np.arange(len(secretome_7dpi_inh_plot)),
       height = secretome_7dpi_inh_plot.Secretome/100,
       yerr = 3 * secretome_7dpi_inh_plot.SEM/100, capsize = 2,
       color = secretome_7dpi_inh_plot.color)
plt.xlabel('Cell type')
plt.xticks(np.arange(len(secretome_7dpi_inh_plot)), secretome_7dpi_inh_plot.Cell_type, rotation=270)
plt.ylabel('Secretome (%)')
plt.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_7dpi_inh_test.pdf', bbox_inches = 'tight', transparent = True)
plt.show()


# In[34]:


secretome_averages_7dpi_inh_plot = np.log1p(secretome_averages_7dpi_inh[secretome_averages_7dpi_inh.max(axis = 1) > 10])
secretome_averages_7dpi_inh_plot = secretome_averages_7dpi_inh_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)         


# In[35]:


sns.clustermap(secretome_averages_7dpi_inh_plot, method = 'ward', xticklabels = True, yticklabels = True)
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_7dpi_inh_high_genes.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[37]:


secretome_z_averages_7dpi_inh_plot = secretome_z_averages_7dpi_inh[secretome_z_averages_7dpi_inh.max(axis = 1) > 2]
secretome_z_averages_7dpi_inh_plot[secretome_z_averages_7dpi_inh_plot > 5] = 5
secretome_z_averages_7dpi_inh_plot = secretome_z_averages_7dpi_inh_plot.drop(['Myelin cells', 'Neuronal cells'], axis = 1)

sns.clustermap(secretome_z_averages_7dpi_inh_plot[secretome_z_averages_7dpi_inh_plot.max(axis = 1) > 2], method = 'ward', xticklabels = True, yticklabels = True, figsize = (10,25))
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_7dpi_inh_high_z.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[127]:


# Prep data to plot together
secretome_ctrl_plot = CountSecretome(HR_Rv_ctrl, secretome)
secretome_ctrl_plot = secretome_ctrl_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_ctrl_plot = secretome_ctrl_plot.sort_values(by = 'Secretome', ascending = False)

secretome_3dpi_plot = CountSecretome(HR_Rv_3dpi, secretome)
secretome_3dpi_plot = secretome_3dpi_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_3dpi_plot = secretome_3dpi_plot.sort_values(by = 'Secretome', ascending = False)

secretome_7dpi_plot = CountSecretome(HR_Rv_7dpi, secretome)
secretome_7dpi_plot = secretome_7dpi_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_7dpi_plot = secretome_7dpi_plot.sort_values(by = 'Secretome', ascending = False)

secretome_3dpi_inh_plot = CountSecretome(HR_Rv_3dpi_inh, secretome)
secretome_3dpi_inh_plot = secretome_3dpi_inh_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_3dpi_inh_plot = secretome_3dpi_inh_plot.sort_values(by = 'Secretome', ascending = False)

secretome_7dpi_inh_plot = CountSecretome(HR_Rv_7dpi_inh, secretome)
secretome_7dpi_inh_plot = secretome_7dpi_inh_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_7dpi_inh_plot = secretome_7dpi_inh_plot.sort_values(by = 'Secretome', ascending = False)

secretome_array = [secretome_ctrl_plot, secretome_3dpi_plot, secretome_7dpi_plot, 
                   secretome_3dpi_inh_plot, secretome_7dpi_inh_plot]
annotation_array = ['ctrl', '3dpi', '7dpi', '3dpi IWR1', '7dpi IWR1']


# In[217]:


# Plot all secretome together - all cell types
cell_types_drop = ['Neuronal cells', 'Myelin cells']
plot_all_secretomes = [x[~x.Cell_type.isin(cell_types_drop)] for x in secretome_array]

fig, axs = plt.subplots(1, len(plot_all_secretomes), sharey='row', figsize=(30,10), gridspec_kw={'hspace': 0.1, 'wspace': 0.1})
for i in np.arange(len(plot_all_secretomes)):
    plot_data = plot_all_secretomes[i]
    condition = annotation_array[i]
    axs[i].bar(x = np.arange(len(plot_data)),
               height = plot_data.Secretome/100,
               yerr = plot_data.SEM/100, capsize = 2,
               color = plot_data.color)
    axs[i].set_xticks(np.arange(len(plot_data))) 
    axs[i].set_xticklabels(plot_data.Cell_type, rotation=270)

for ax, col in zip(axs, annotation_array):
    ax.set_title(col)

for ax, row in zip(axs, 'Secretome (%)'):
    ax.set_ylabel('Secretome (%)')

for ax in axs.flat:
    ax.label_outer()

plt.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_all_conditions.pdf', 
            bbox_inches = 'tight', transparent = True)


# In[218]:


# Plot all secretome together - epiniche
cell_types_drop = ['Neuronal cells', 'Myelin cells']
plot_all_secretomes = [x[x.Cell_type.isin(epifibro_types)] for x in secretome_array]
#secretome_3dpi.Cell_type.isin(epifibro_types)

fig, axs = plt.subplots(1, len(plot_all_secretomes), sharey='row', figsize=(20,10), gridspec_kw={'hspace': 0.1, 'wspace': 0.1})
for i in np.arange(len(plot_all_secretomes)):
    plot_data = plot_all_secretomes[i]
    condition = annotation_array[i]
    axs[i].bar(x = np.arange(len(plot_data)),
               height = plot_data.Secretome/100,
               yerr = plot_data.SEM/100, capsize = 2,
               color = plot_data.color)
    axs[i].set_xticks(np.arange(len(plot_data))) 
    axs[i].set_xticklabels(plot_data.Cell_type, rotation=270)

for ax, col in zip(axs, annotation_array):
    ax.set_title(col)

for ax, row in zip(axs, 'Secretome (%)'):
    ax.set_ylabel('Secretome (%)')

for ax in axs.flat:
    ax.label_outer()

plt.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_all_conditions_epiniche.pdf', 
            bbox_inches = 'tight', transparent = True)


# In[135]:


plot_all_secretomes[0]


# In[ ]:


# Plot all secretome together - selected cell types


# In[27]:


secretome_averages_3dpi.loc[['wif1', 'dkk3a', 'dkk3b']]


# In[26]:


secretome_z_averages_7dpi.loc['wif1']

