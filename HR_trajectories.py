#!/usr/bin/env python
# coding: utf-8

# # Dependencies and parameters

# In[1]:


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
scv.settings.set_figure_params('scvelo', dpi_save = 300)


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


connected_epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)',
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

# In[ ]:


H5_Rv_data = scv.read('../Data/RNAvelo/H5_v3Dr11.loom', cache=True)
H5_Rv_data.var_names_make_unique()
H6_Rv_data = scv.read('../Data/RNAvelo/H6_v3Dr11.loom', cache=True)
H6_Rv_data.var_names_make_unique()
H7_Rv_data = scv.read('../Data/RNAvelo/H7_v3Dr11.loom', cache=True)
H7_Rv_data.var_names_make_unique()
H8a_Rv_data = scv.read('../Data/RNAvelo/H8a_v3Dr11.loom', cache=True)
H8a_Rv_data.var_names_make_unique()
H8v_Rv_data = scv.read('../Data/RNAvelo/H8v_v3Dr11.loom', cache=True)
H8v_Rv_data.var_names_make_unique()
Hr1_Rv_data = scv.read('../Data/RNAvelo/Hr1_v3Dr11.loom', cache=True)
Hr1_Rv_data.var_names_make_unique()
Hr2a_Rv_data = scv.read('../Data/RNAvelo/Hr2a_v3Dr11.loom', cache=True)
Hr2a_Rv_data.var_names_make_unique()
Hr2b_Rv_data = scv.read('../Data/RNAvelo/Hr2b_v3Dr11.loom', cache=True)
Hr2b_Rv_data.var_names_make_unique()
Hr3_Rv_data = scv.read('../Data/RNAvelo/Hr3_v3Dr11.loom', cache=True)
Hr3_Rv_data.var_names_make_unique()
Hr4_Rv_data = scv.read('../Data/RNAvelo/Hr4_v3Dr11.loom', cache=True)
Hr4_Rv_data.var_names_make_unique()
Hr5_Rv_data = scv.read('../Data/RNAvelo/Hr5_v3Dr11.loom', cache=True)
Hr5_Rv_data.var_names_make_unique()
Hr6a_Rv_data = scv.read('../Data/RNAvelo/Hr6a_v3Dr11.loom', cache=True)
Hr6a_Rv_data.var_names_make_unique()
Hr6v_Rv_data = scv.read('../Data/RNAvelo/Hr6v_v3Dr11.loom', cache=True)
Hr6v_Rv_data.var_names_make_unique()
Hr7a_Rv_data = scv.read('../Data/RNAvelo/Hr7a_v3Dr11.loom', cache=True)
Hr7a_Rv_data.var_names_make_unique()
Hr7v_Rv_data = scv.read('../Data/RNAvelo/Hr7v_v3Dr11.loom', cache=True)
Hr7v_Rv_data.var_names_make_unique()
Hr8_Rv_data = scv.read('../Data/RNAvelo/Hr8_v3Dr11.loom', cache=True)
Hr8_Rv_data.var_names_make_unique()
Hr9_Rv_data = scv.read('../Data/RNAvelo/Hr9_v3Dr11.loom', cache=True)
Hr9_Rv_data.var_names_make_unique()
Hr10_Rv_data = scv.read('../Data/RNAvelo/Hr10_v3Dr11.loom', cache=True)
Hr10_Rv_data.var_names_make_unique()
Hr11_Rv_data = scv.read('../Data/RNAvelo/Hr11_v3Dr11.loom', cache=True)
Hr11_Rv_data.var_names_make_unique()
Hr12_Rv_data = scv.read('../Data/RNAvelo/Hr12_v3Dr11.loom', cache=True)
Hr12_Rv_data.var_names_make_unique()
Hr13_Rv_data = scv.read('../Data/RNAvelo/Hr13_v3Dr11.loom', cache=True)
Hr13_Rv_data.var_names_make_unique()
Hr14_Rv_data = scv.read('../Data/RNAvelo/Hr14_v3Dr11.loom', cache=True)
Hr14_Rv_data.var_names_make_unique()
Hr15_Rv_data = scv.read('../Data/RNAvelo/Hr15_v3Dr11.loom', cache=True)
Hr15_Rv_data.var_names_make_unique()
Hr16_Rv_data = scv.read('../Data/RNAvelo/Hr16_v3Dr11.loom', cache=True)
Hr16_Rv_data.var_names_make_unique()
Hr17_Rv_data = scv.read('../Data/RNAvelo/Hr17_v3Dr11.loom', cache=True)
Hr17_Rv_data.var_names_make_unique()
Hr18_Rv_data = scv.read('../Data/RNAvelo/Hr18_v3Dr11.loom', cache=True)
Hr18_Rv_data.var_names_make_unique()
Hr19_Rv_data = scv.read('../Data/RNAvelo/Hr19_v3Dr11.loom', cache=True)
Hr19_Rv_data.var_names_make_unique()
Hr20_Rv_data = scv.read('../Data/RNAvelo/Hr20_v3Dr11.loom', cache=True)
Hr20_Rv_data.var_names_make_unique()
Hr21_Rv_data = scv.read('../Data/RNAvelo/Hr21_v3Dr11.loom', cache=True)
Hr21_Rv_data.var_names_make_unique()
Hr22_Rv_data = scv.read('../Data/RNAvelo/Hr22_v3Dr11.loom', cache=True)
Hr22_Rv_data.var_names_make_unique()
Hr23_Rv_data = scv.read('../Data/RNAvelo/Hr23_v3Dr11.loom', cache=True)
Hr23_Rv_data.var_names_make_unique()
Hr24_Rv_data = scv.read('../Data/RNAvelo/Hr24_v3Dr11.loom', cache=True)
Hr24_Rv_data.var_names_make_unique()
Hr25_Rv_data = scv.read('../Data/RNAvelo/Hr25_v3Dr11.loom', cache=True)
Hr25_Rv_data.var_names_make_unique()
Hr26_Rv_data = scv.read('../Data/RNAvelo/Hr26_v3Dr11.loom', cache=True)
Hr26_Rv_data.var_names_make_unique()
Hr27_Rv_data = scv.read('../Data/RNAvelo/Hr27_v3Dr11.loom', cache=True)
Hr27_Rv_data.var_names_make_unique()
Hr28_Rv_data = scv.read('../Data/RNAvelo/Hr28_v3Dr11.loom', cache=True)
Hr28_Rv_data.var_names_make_unique()
Hr29_Rv_data = scv.read('../Data/RNAvelo/Hr29_v3Dr11.loom', cache=True)
Hr29_Rv_data.var_names_make_unique()
Hr30_Rv_data = scv.read('../Data/RNAvelo/Hr30_v3Dr11.loom', cache=True)
Hr30_Rv_data.var_names_make_unique()
Hr31_Rv_data = scv.read('../Data/RNAvelo/Hr31_v3Dr11.loom', cache=True)
Hr31_Rv_data.var_names_make_unique()
Hr32_Rv_data = scv.read('../Data/RNAvelo/Hr32_v3Dr11.loom', cache=True)
Hr32_Rv_data.var_names_make_unique()
Hr33_Rv_data = scv.read('../Data/RNAvelo/Hr33_v3Dr11.loom', cache=True)
Hr33_Rv_data.var_names_make_unique()
Hr34_Rv_data = scv.read('../Data/RNAvelo/Hr34_v3Dr11.loom', cache=True)
Hr34_Rv_data.var_names_make_unique()
Hr35_Rv_data = scv.read('../Data/RNAvelo/Hr35_v3Dr11.loom', cache=True)
Hr35_Rv_data.var_names_make_unique()


# In[ ]:


HR_Rv =    H5_Rv_data.concatenate(H6_Rv_data, H7_Rv_data, H8a_Rv_data, H8v_Rv_data,
                        Hr1_Rv_data, Hr2a_Rv_data, Hr2b_Rv_data, Hr3_Rv_data, Hr4_Rv_data, Hr5_Rv_data, 
                        Hr6a_Rv_data, Hr6v_Rv_data, Hr7a_Rv_data, Hr7v_Rv_data, Hr8_Rv_data, Hr9_Rv_data, Hr10_Rv_data,
                        Hr11_Rv_data, Hr12_Rv_data, Hr13_Rv_data, Hr14_Rv_data, Hr15_Rv_data,
                        Hr16_Rv_data, Hr17_Rv_data, Hr18_Rv_data, Hr19_Rv_data, Hr20_Rv_data,
                        Hr21_Rv_data, Hr22_Rv_data, Hr23_Rv_data, Hr24_Rv_data, Hr25_Rv_data,
                        Hr26_Rv_data, Hr27_Rv_data, Hr28_Rv_data, Hr29_Rv_data, Hr30_Rv_data,
                        Hr31_Rv_data, Hr32_Rv_data, Hr33_Rv_data, Hr34_Rv_data, Hr35_Rv_data)
HR_Rv.shape


# In[ ]:


HR_Rv = sc.read('./write/HR_Rv.h5ad')


# In[ ]:


HR_Rv.obs = HR_Rv.obs.reset_index().merge(HR_setnames, how="inner").set_index('index')


# In[ ]:


# Rename cells to match cell names in annotation file
HR_Rv.obs_names = [str(HR_Rv.obs.loc[x,'heart'])+'_'+str(x.split(':', 1)[1])[0:16] for x in HR_Rv.obs_names]
# Drop annotations that are not in the single-cell object
anno_drop_Rv = annotations.index.difference(HR_Rv.obs_names)
annotations_Rv = annotations.drop(anno_drop_Rv)


# In[ ]:


HR_Rv_filter = HR_Rv[annotations_Rv.index]


# In[ ]:


HR_Rv_filter.obs['Cell_type'] = annotations_Rv['Cell_type'].tolist()


# In[2]:


#HR_Rv_filter.write('./write/HR_Rv_filter.h5ad')
HR_Rv_filter = sc.read('./write/HR_Rv_filter.h5ad')


# In[11]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]


# In[12]:


HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]


# # Naive trajectory analysis at 3dpi

# In[ ]:


HR_Rv_3dpi_epifibro = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(epifibro_types)]
sc.pp.filter_genes(HR_Rv_3dpi_epifibro, min_cells=3)


# In[ ]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_epifibro.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_epifibro.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_epifibro[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_epifibro.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_epifibro, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_3dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[ ]:


sc.pp.normalize_total(HR_Rv_3dpi_epifibro, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi_epifibro)


# In[ ]:


sc.pp.regress_out(HR_Rv_3dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_epifibro)
sc.tl.pca(HR_Rv_3dpi_epifibro)
sc.external.pp.bbknn(HR_Rv_3dpi_epifibro, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_epifibro)


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color='batch',
          title = '3dpi epicardium and fibroblasts')


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epifibro.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '', save = 'subset1_3dpi_umap.png', legend_loc = 'none', frameon=False)


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '3dpi epicardium and fibroblasts')


# In[ ]:


sc.tl.paga(HR_Rv_3dpi_epifibro, groups='Cell_type')


# In[ ]:


sc.pl.paga(HR_Rv_3dpi_epifibro, threshold=0.3, show=True, frameon=False,
           labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_3dpi_paga.png')


# # Naive trajectory analysis at 7dpi

# In[ ]:


HR_Rv_7dpi_epifibro = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(epifibro_types)]
sc.pp.filter_genes(HR_Rv_7dpi_epifibro, min_cells=3)


# In[ ]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_epifibro.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_epifibro.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_epifibro[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_epifibro.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_epifibro, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[ ]:


sc.pp.normalize_total(HR_Rv_7dpi_epifibro, target_sum=1e4)
sc.pp.log1p(HR_Rv_7dpi_epifibro)


# In[ ]:


sc.pp.regress_out(HR_Rv_7dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_epifibro)
sc.tl.pca(HR_Rv_7dpi_epifibro)
sc.external.pp.bbknn(HR_Rv_7dpi_epifibro, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_epifibro)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color='batch',
          title = '7dpi epicardium and fibroblasts')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_epifibro.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '', save = 'subset1_7dpi_umap.png', legend_loc = 'none', frameon=False)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '7dpi epicardium and fibroblasts')


# In[ ]:


sc.tl.paga(HR_Rv_7dpi_epifibro, groups='Cell_type')


# In[ ]:


sc.pl.paga(HR_Rv_7dpi_epifibro, threshold=0.3, show=True, frameon=False,
           labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_7dpi_paga.png')


# # Trajectories in 3dpi epicardial connected niche

# In[13]:


HR_Rv_3dpi_epiconn = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_3dpi_epiconn, min_cells=3)


# In[14]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_epiconn.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_epiconn.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_epiconn[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_epiconn.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_epiconn, percent_top=None, log1p=True, inplace=True)
#HR_Rv_3dpi_epiconn.obs['n_counts'] = HR_Rv_3dpi_epiconn.X.sum(axis=1)
sc.pl.violin(HR_Rv_3dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[15]:


sc.pp.normalize_total(HR_Rv_3dpi_epiconn, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi_epiconn)


# In[16]:


sc.pp.regress_out(HR_Rv_3dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[17]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_epiconn)
sc.tl.pca(HR_Rv_3dpi_epiconn)
sc.external.pp.bbknn(HR_Rv_3dpi_epiconn, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_epiconn)


# In[18]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='batch',
          title = '3dpi epicardial niche')


# In[19]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[20]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color=['n_genes_by_counts'], 
          title = '3dpi epicardial niche')


# In[21]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color=['npr3'])


# In[21]:


sc.tl.leiden(HR_Rv_3dpi_epiconn, resolution=2)


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='leiden', legend_loc='on data', legend_fontsize='x-large')#,
#          palette = ['#92aad4', '#E9D723',  '#7b7f51', '#92aad4',
#                    '#E9D723', '#E9D723', '#E9D723', '#c6ba83',
#                    '#E9D723', '#c6ba83', '#525566', '#CE3A39',
#                    '#E9D723', '#CE3A39', '#e1e3d9', '#E9D723'])


# In[ ]:


sc.tl.paga(HR_Rv_3dpi_epiconn, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_3dpi_epiconn, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[ ]:


sc.tl.umap(HR_Rv_3dpi_epiconn, init_pos = 'paga')


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='batch',
          title = '3dpi epicardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='tcf21')


# In[ ]:


sc.pl.violin(HR_Rv_3dpi_epiconn, keys = 'tcf21', groupby = 'Cell_type')


# In[ ]:


scv.pp.moments(HR_Rv_3dpi_epiconn)


# In[ ]:


scv.tl.velocity(HR_Rv_3dpi_epiconn, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_3dpi_epiconn)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
                                 #save = '_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'Cell_type')


# In[ ]:


scv.tl.velocity_confidence(HR_Rv_3dpi_epiconn)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_3dpi_epiconn, c=keys, cmap='coolwarm', perc=[5, 95])#,
#                                save = 'strength_coherence_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[ ]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6', 'nrg1']
scv.pl.velocity(HR_Rv_3dpi_epiconn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
#                save = 'gene_velocity_epicardial_niche_3dpi_umap.png')


# ## Which genes drive the transitions?

# In[ ]:





# In[ ]:


sc.pl.paga(HR_Rv_3dpi_epiconn, show=True, color = 'Cell_type', node_size_scale = 2)


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color=['leiden', 'Cell_type'], legend_loc='on data', legend_fontsize='x-large')


# For each cell of type 1, extract the maximum transition probability to type 2. The mean and distribution of these probabilities tell us how likely a source type 1 is for type 2. We compare these maximum transition probabilities to transition probabilities of type 1 to 'anything but type 1+2' to estimate whether the transition probabilities from type 1 to type 2 are significantly higher than a null hypothesis.
# The questions we want to answer with this are:  
# 1) Is there a difference between transitions from Epi V to col11/col12 at 3dpi and at 7dpi?  
# 2) Is the seeming transition from constitutional fibroblasts to epicardium as meaningful as that from constitutional fibroblasts to col11/col12 at 3dpi and 7dpi? At 3dpi, the constitutional fibroblast clusters involved are 2 and 7; col11/col12 are 4,5,11,13; epicardium is 0,3,10

# In[ ]:


print(HR_Rv_3dpi_epiconn.uns['velocity_graph'].max(axis=1).max(axis=0))


# In[50]:


# Epicardium to col11/col12
leiden_1 = ['0', '3', '10']
leiden_2 = ['4', '5', '11', '13']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if not HR_Rv_3dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2, sharex=True)
fig.suptitle('3dpi', y = 1.05)
ax1.set_title('Epicardium to col11/col12')
ax2.set_title('Epicardium to rest')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
pl.savefig('figures/Epi_col1112_transitions_3dpi.png')
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[100]:


# Col11/col12 to epicardium
leiden_1 = ['4', '5', '11', '13']
leiden_2 = ['0', '3', '10']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if not HR_Rv_3dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[101]:


# Const. fibroblasts to col11/col12 fibroblasts
leiden_1 = ['2', '7']
leiden_2 = ['4', '5', '11', '13']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if not HR_Rv_3dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[102]:


# Col11/col12 to cfibro
leiden_1 = ['4', '5', '11', '13']
leiden_2 = ['2', '7']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if not HR_Rv_3dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[51]:


# Const. fibroblasts to epicardium
leiden_1 = ['2', '7']
leiden_2 = ['0', '3', '10']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if not HR_Rv_3dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2, sharex=True)
#fig.suptitle('Horizontally stacked subplots')
ax1.set_title('Fibroblasts to epicardium')
ax2.set_title('Fibroblasts to rest')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
pl.savefig('figures/Fibroc_to_epi_transitions_3dpi.png')
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[104]:


# Epicardium to const. fibroblasts
leiden_1 = ['0', '3', '10']
leiden_2 = ['2', '7']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if not HR_Rv_3dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[44]:


from scipy import stats
print(stats.ttest_ind(max_transitions_1to2, max_transitions_1tonot2, equal_var=False))
test_2 = stats.ttest_ind(max_transitions_1to2[max_transitions_1to2 != 0], max_transitions_1tonot2, equal_var=False)
print(test_2[1])


# In[86]:


# Determine frequency of preferred transitions
from scipy import sparse
vgraph = HR_Rv_3dpi_epiconn.uns['velocity_graph'] > 0.1
# For each cell, determine its favorite transition cell
#sparse.csr_argmax(HR_Rv_3dpi_epiconn.uns['velocity_graph'])
favorite_cell = HR_Rv_3dpi_epiconn.uns['velocity_graph'].argmax(axis = 0) # Why am I seeing cell 0 overrepresented here? Is this maybe the favorite for cells that have very low
# likelihoods? Maybe if transition probability is 0 to everything, this is the cell selected?
# Aggregate favorite transitions over cell types (from and to) to give a matrix of #transitions from cell type A (row) to cell type B (column)
#HR_Rv_3dpi_epiconn.obs.Cell_type.iloc(favorite_cell) # HERE! This does not work
# Normalize by number of from and to cells 
# (normalize how exactly? what is the background model here? Random connections, feels like it should be divided by both number of to and from cells)
# Note we can also downsample.


# In[107]:


vgraph = HR_Rv_3dpi_epiconn.uns['velocity_graph']
print(vgraph[1])
vgraph[1].argmax()


# In[102]:


vgraph[vgraph < 0.1] = 0
vgraph.eliminate_zeros()
favorite_cell = vgraph.argmax(axis = 0)
unique, counts = np.unique(favorite_cell, return_counts=True, axis = 1)
dict(zip(unique[0], counts))
# Lol that doesn't help - now we have 460 instead of 305.


# In[85]:


#favorite_cell
unique, counts = np.unique(favorite_cell, return_counts=True, axis = 1)
#counts.argmin()
dict(zip(unique[0], counts))


# # Trajectories in 7dpi epicardial connected niche

# In[68]:


HR_Rv_7dpi_epiconn = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_7dpi_epiconn, min_cells=3)


# In[69]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_epiconn.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_epiconn.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_epiconn[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_epiconn.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_epiconn, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[70]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_epiconn, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_epiconn)


# In[71]:


sc.pp.regress_out(HR_Rv_7dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[72]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_epiconn)
sc.tl.pca(HR_Rv_7dpi_epiconn)
sc.external.pp.bbknn(HR_Rv_7dpi_epiconn, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_epiconn)


# In[73]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='batch',
          title = '7dpi epicardial niche')


# In[74]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color=['total_counts', 'n_genes_by_counts', 'percent_mito'],
          title = '7dpi epicardial niche')


# In[75]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epicardial niche')


# In[76]:


sc.tl.leiden(HR_Rv_7dpi_epiconn, resolution=2)


# In[77]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[78]:


sc.tl.paga(HR_Rv_7dpi_epiconn, groups='leiden')


# In[79]:


sc.pl.paga(HR_Rv_7dpi_epiconn, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[80]:


sc.tl.umap(HR_Rv_7dpi_epiconn, init_pos = 'paga')


# In[81]:


scv.pp.moments(HR_Rv_7dpi_epiconn)


# In[82]:


scv.tl.velocity(HR_Rv_7dpi_epiconn, mode='stochastic')


# In[83]:


scv.tl.velocity_graph(HR_Rv_7dpi_epiconn)


# In[84]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_epiconn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
#                                 save = '_scvelo_epicardial_niche_7dpi_20n_cell_types_umap.png')


# In[85]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6','nrg1']
scv.pl.velocity(HR_Rv_7dpi_epiconn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
#                save = 'gene_velocity_epicardial_niche_7dpi_umap.png')


# In[86]:


sc.pl.paga(HR_Rv_7dpi_epiconn, show=True, color = 'Cell_type', node_size_scale = 2)


# In[87]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color=['leiden', 'Cell_type'], legend_loc='on data', legend_fontsize='x-large')


# In[88]:


print(HR_Rv_7dpi_epiconn.uns['velocity_graph'].max(axis=1).max(axis=0))


# Epi: 1, 2, 5, 10; col11/col12: 4, 13; cfibro: 7

# In[105]:


# Epicardium to col11/col12
leiden_1 = ['1', '2', '5', '10']
leiden_2 = ['4', '13']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[106]:


# Col11/col12 to epicardium
leiden_1 = ['4', '13']
leiden_2 = ['1', '2', '5', '10']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[107]:


# Const. fibroblasts to col11/col12 fibroblasts
leiden_1 = ['7']
leiden_2 = ['4', '13']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[108]:


# Col11/col12 to cfibro
leiden_1 = ['4', '13']
leiden_2 = ['7']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[109]:


# Const. fibroblasts to epicardium
leiden_1 = ['7']
leiden_2 = ['1', '2', '5', '10']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[110]:


# Epicardium to const. fibroblasts
leiden_1 = ['1', '2', '5', '10']
leiden_2 = ['7']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[95]:


# Epi to col11/col12 fibroblasts
leiden_1 = ['1', '2', '5', '10']
leiden_2 = ['4', '13']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2))
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2))


# In[98]:


from scipy import stats
print(stats.ttest_ind(max_transitions_1to2, max_transitions_1tonot2, equal_var=False))
test_2 = stats.ttest_ind(max_transitions_1to2[max_transitions_1to2 != 0], max_transitions_1tonot2, equal_var=False)
print(test_2[1])


# In[96]:


# Const. fibroblasts to col11/col12 fibroblasts
leiden_1 = ['7']
leiden_2 = ['4', '13']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2))
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2))


# In[97]:


# Const. fibroblasts to epicardium - closer than col11/col12 fibroblasts.
leiden_1 = ['7']
leiden_2 = ['1', '2', '5', '10']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if HR_Rv_7dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_leiden_2_index = [i for i in range(len(HR_Rv_7dpi_epiconn.obs.leiden)) if not HR_Rv_7dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
transitions_1tonot2 = (HR_Rv_7dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2))
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2))


# In[ ]:


HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'].isin(['0', '3'])) & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]


# In[ ]:


HR_Rv_ctrl_epiconn = HR_Rv_ctrl[HR_Rv_ctrl.obs['Cell_type'].isin(connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_ctrl, min_cells=3)


# In[ ]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_ctrl.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_ctrl.obs['percent_mito'] = np.sum(
    HR_Rv_ctrl[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_ctrl.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_ctrl, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_ctrl, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[ ]:


sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)
sc.pp.log1p(HR_Rv_ctrl)


# In[ ]:


sc.pp.regress_out(HR_Rv_ctrl, ['total_counts', 'n_genes_by_counts'])


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_ctrl)
sc.tl.pca(HR_Rv_ctrl)
sc.external.pp.bbknn(HR_Rv_ctrl, batch_key='batch')
sc.tl.umap(HR_Rv_ctrl)


# In[ ]:


sc.pp.filter_genes(HR_Rv_ctrl_epiconn, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_ctrl_epiconn, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_ctrl_epiconn)
sc.pp.highly_variable_genes(HR_Rv_ctrl_epiconn)
sc.tl.pca(HR_Rv_ctrl_epiconn)
sc.external.pp.bbknn(HR_Rv_ctrl_epiconn, batch_key='batch')#, neighbors_within_batch = 3, n_pcs = 20)
sc.tl.umap(HR_Rv_ctrl_epiconn)


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_ctrl_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          )#save = 'Epicardial_niche_control_3dpi_Celltypes.png')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epiconn, color='dpi',
          )#save = 'Epicardial_niche_control_3dpi_timepoints.png')


# In[ ]:


sc.tl.leiden(HR_Rv_ctrl_epiconn, resolution=2)


# In[ ]:


sc.tl.paga(HR_Rv_ctrl_epiconn, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_ctrl_epiconn, show=True, color ='Cell_type')
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', ''], node_size_scale = 2, 
#           save = 'scvelo_epicardial_niche_control_20n_leiden_diffmap.png')


# In[ ]:


sc.pl.paga(HR_Rv_ctrl_epiconn, show=True, color ='dpi')


# In[ ]:


sc.tl.umap(HR_Rv_ctrl_epiconn, init_pos = 'paga')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epiconn, color='Cell_type',
          save = 'Epicardial_niche_control_3dpi_Celltypes.png')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epiconn, color='dpi',
          save = 'Epicardial_niche_control_3dpi_timepoints.png')


# In[ ]:


scv.pp.moments(HR_Rv_ctrl_epiconn)
scv.tl.velocity(HR_Rv_ctrl_epiconn, mode='stochastic')
scv.tl.velocity_graph(HR_Rv_ctrl_epiconn)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_ctrl_epiconn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none',
                                save = 'Epicardial_niche_control_3dpi_Celltypes_velo.png')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_ctrl_epiconn, title = '', scale = 0.2,
                                 color = 'Cell_type', legend_loc = 'none')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_ctrl_epiconn[HR_Rv_ctrl_epiconn.obs['dpi'] == '0'], title = '', scale = 0.2,
                                color = 'dpi', legend_loc = 'none',
                         save = 'Epicardial_niche_control_3dpi_control_velo.png')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_ctrl_epiconn[HR_Rv_ctrl_epiconn.obs['dpi'] == '3'], title = '', scale = 0.2,
                                color = 'dpi', legend_loc = 'none',
                         save = 'Epicardial_niche_control_3dpi_3dpi_velo.png')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epiconn, color='leiden',
          save = 'Epicardial_niche_control_3dpi_leidenclusters.png')


# In[ ]:


HR_Rv_ctrl_epiconn[HR_Rv_ctrl_epiconn.obs['dpi'] == '3'].obs.leiden.value_counts()/sum(HR_Rv_ctrl_epiconn[HR_Rv_ctrl_epiconn.obs['dpi'] == '3'].obs.leiden.value_counts())


# In[ ]:


HR_Rv_ctrl_epiconn[HR_Rv_ctrl_epiconn.obs['dpi'] == '0'].obs.leiden.value_counts()/sum(HR_Rv_ctrl_epiconn[HR_Rv_ctrl_epiconn.obs['dpi'] == '0'].obs.leiden.value_counts())


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epiconn[HR_Rv_ctrl_epiconn.obs['dpi'] == '0'], color = 'Cell_type',
          save = 'Epicardial_niche_control_3dpi_Celltypes_control_only.png')


# # Endocardial niche for deep injuries at 7dpi

# In[22]:


deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9'] # All samples with > 50 nppc fibroblasts.
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_endofibro_types)]
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['heart'].isin(deep_injury_libraries)]
scv.pp.remove_duplicate_cells(HR_Rv_7dpi_deep_endo)
sc.pp.filter_genes(HR_Rv_7dpi_deep_endo, min_cells=3)


# In[23]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_deep_endo.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_deep_endo.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_deep_endo[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_deep_endo.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_deep_endo, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[24]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_deep_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_deep_endo)


# In[25]:


sc.pp.regress_out(HR_Rv_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[26]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_deep_endo)
sc.tl.pca(HR_Rv_7dpi_deep_endo)
sc.external.pp.bbknn(HR_Rv_7dpi_deep_endo, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_deep_endo)


# In[27]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi deep injury endocardial niche')


# In[28]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '7dpi deep injury endocardial niche')


# In[29]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='batch',
          title = '7dpi deep injury endocardial niche')


# In[30]:


sc.tl.leiden(HR_Rv_7dpi_deep_endo, resolution=2)


# In[31]:


sc.tl.paga(HR_Rv_7dpi_deep_endo, groups='leiden')


# In[32]:


sc.pl.paga(HR_Rv_7dpi_deep_endo, show=True, color = 'Cell_type', node_size_scale = 2)#, 
    #       labels = ['', '', '', '', '', 
   #                  '', '', '', '', '', 
   #                  '', '', '', '', '', '', ''], 
  #         node_size_scale = 2)#, save = 'scvelo_deep_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[33]:


sc.tl.umap(HR_Rv_7dpi_deep_endo, init_pos='paga')


# In[34]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'leiden')


# In[35]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'batch')


# In[36]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'Cell_type')


# In[37]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = ['pecam1', 'cdh5', 'nppc'])


# In[38]:


scv.pp.moments(HR_Rv_7dpi_deep_endo)


# In[39]:


scv.tl.velocity(HR_Rv_7dpi_deep_endo, mode='stochastic')


# In[40]:


scv.tl.velocity_graph(HR_Rv_7dpi_deep_endo)


# In[41]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, basis='umap', title = '', 
                                 color = 'Cell_type',
                                save = 'gene_velocity_endoventricular_deep_niche_7dpi_umap.png')


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


# # Test transition probabilities on endocrine pancreas dataset

# In[52]:


adata = scv.datasets.pancreas()


# In[53]:


scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


# In[54]:


scv.tl.velocity(adata)


# In[55]:


scv.tl.velocity_graph(adata)


# In[56]:


scv.pl.velocity_embedding_stream(adata, basis='umap')


# In[57]:


print(adata.uns['velocity_graph'].max(axis=1).max(axis=0))


# In[58]:


adata.obs


# clusters: 'Ductal' to 'Ngn3 low EP', 'Ngn3 high EP' to 'Pre-endocrine'

# In[119]:


# Ductal to Ngn3 low EP
clusters_1 = ['Ductal']
clusters_2 = ['Ngn3 low EP']
# Find (numeric) indices for both types
clusters_1_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_1] 
clusters_2_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_2]
# Subset velocity graph for both types
transitions_1to2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, clusters_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_clusters_2_index = [i for i in range(len(adata.obs.clusters)) if not adata.obs.clusters[i] in (clusters_1 + clusters_2)]
transitions_1tonot2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, not_clusters_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[59]:


# Ngn3 low EP to Ductal
clusters_1 = ['Ngn3 low EP']
clusters_2 = ['Ductal']
# Find (numeric) indices for both types
clusters_1_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_1] 
clusters_2_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_2]
# Subset velocity graph for both types
transitions_1to2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, clusters_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_clusters_2_index = [i for i in range(len(adata.obs.clusters)) if not adata.obs.clusters[i] in (clusters_1 + clusters_2)]
transitions_1tonot2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, not_clusters_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[60]:


# 'Ngn3 high EP' to 'Pre-endocrine'
clusters_1 = ['Ngn3 high EP']
clusters_2 = ['Pre-endocrine']
# Find (numeric) indices for both types
clusters_1_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_1] 
clusters_2_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_2]
# Subset velocity graph for both types
transitions_1to2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, clusters_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_clusters_2_index = [i for i in range(len(adata.obs.clusters)) if not adata.obs.clusters[i] in (clusters_1 + clusters_2)]
transitions_1tonot2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, not_clusters_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2, sharex=True)
ax1.set_title('Ngn3 high EP to Pre-endocrine')
ax2.set_title('Ngn3 high EP to rest')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
pl.savefig('figures/Ngn3highEP_to_Pre-endo_transitions_3dpi.png')
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[124]:


# 'Pre-endocrine' to 'Ngn3 high EP'
clusters_1 = ['Pre-endocrine']
clusters_2 = ['Ngn3 high EP']
# Find (numeric) indices for both types
clusters_1_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_1] 
clusters_2_index = [i for i in range(len(adata.obs.clusters)) if adata.obs.clusters[i] in clusters_2]
# Subset velocity graph for both types
transitions_1to2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, clusters_2_index]
# Take row maxima
max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
not_clusters_2_index = [i for i in range(len(adata.obs.clusters)) if not adata.obs.clusters[i] in (clusters_1 + clusters_2)]
transitions_1tonot2 = (adata.uns['velocity_graph'])[clusters_1_index,:][:, not_clusters_2_index]
# Take row maxima
max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
fig, (ax1, ax2) = pl.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')
ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
ax2.hist(max_transitions_1tonot2)
print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# # Cardiomyocyte states at 7dpi

# In[169]:


HR_Rv_3dpi.obs['Cell_type'].value_counts()/HR_Rv_3dpi.obs['Cell_type'].value_counts().sum()


# In[170]:


CM_types = ['Cardiomyocytes (Ventricle)', 'Cardiomyocytes (ttn.2)', 'Cardiomyocytes (proliferating)']


# In[171]:


HR_Rv_3dpi_CM = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(CM_types)]
sc.pp.filter_genes(HR_Rv_3dpi_CM, min_cells=3)


# In[172]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_CM.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_CM.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_CM[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_CM.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_CM, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_3dpi_CM, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[173]:


sc.pp.normalize_per_cell(HR_Rv_3dpi_CM, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_3dpi_CM)


# In[174]:


sc.pp.regress_out(HR_Rv_3dpi_CM, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[175]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_CM)
sc.tl.pca(HR_Rv_3dpi_CM)
sc.external.pp.bbknn(HR_Rv_3dpi_CM, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_CM)


# In[176]:


sc.pl.umap(HR_Rv_3dpi_CM, color='batch',
          title = '3dpi Cardiomyocytes')


# In[177]:


sc.pl.umap(HR_Rv_3dpi_CM, color=['total_counts', 'n_genes_by_counts', 'percent_mito'],
          title = '3dpi Cardiomyocytes')


# In[178]:


sc.pl.umap(HR_Rv_3dpi_CM, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_CM.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi Cardiomyocytes')


# In[179]:


sc.tl.leiden(HR_Rv_3dpi_CM, resolution=2)


# In[180]:


sc.pl.umap(HR_Rv_3dpi_CM, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[181]:


sc.tl.paga(HR_Rv_3dpi_CM, groups='leiden')


# In[182]:


sc.pl.paga(HR_Rv_3dpi_CM, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[183]:


sc.tl.umap(HR_Rv_3dpi_CM, init_pos = 'paga')


# In[184]:


scv.pp.moments(HR_Rv_3dpi_CM)


# In[185]:


scv.tl.velocity(HR_Rv_3dpi_CM, mode='stochastic')


# In[186]:


scv.tl.velocity_graph(HR_Rv_3dpi_CM)


# In[187]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_CM, title = '', 
                                 color = 'Cell_type')#, legend_loc = 'none')#,
#                                 save = '_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# # Secretome analysis

# How much of the transcriptome is part of the secretome for each cell type?

# ## Load secretome genes

# In[3]:


#secretome = pd.read_csv('~/Documents/Projects/heart_Bo/Data/1604312481_danio_rerio/Secretome_gene_names.csv', sep = '\t')
#secretome = secretome[secretome['external_gene_name'].notna()]
#secretome = (secretome[['external_gene_name']]).drop_duplicates()
#secretome_nocol = secretome[~secretome['external_gene_name'].str.contains('^col[0-9]', na=False)]


# Remove NAs and duplicates

# In[5]:


#secretome = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Human_secretome_translated.tsv', sep = '\;')
#secretome = secretome[secretome['external_gene_name'].notna()]
#secretome = (secretome[['external_gene_name']]).drop_duplicates()


# In[3]:


secretome = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Alliance_secretome_gene_names_noDRduplicates.scsv', sep = ';')
#secretome = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Alliance_secretome_gene_names.scsv', sep = '\t')
secretome = secretome.rename(columns={'DR_name': 'external_gene_name'})
secretome = secretome[secretome['external_gene_name'].notna()]
secretome = (secretome[['external_gene_name']]).drop_duplicates()


# In[4]:


secretome


# ## Calculate differentially expressed genes per timepoint

# In[12]:


def GetDiffGenes(adata):
    result = pd.DataFrame({})
    for ctype in adata.obs['Cell_type'].cat.categories:
        ctype_d = {'Gene' : adata.uns['rank_genes_groups']['names'][ctype],
                   'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][ctype],
                   'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][ctype],
                   'Cell_type':  np.repeat(ctype, len(adata.var.index))}
        ctype_degdf = pd.DataFrame(data=ctype_d)
        result = result.append(ctype_degdf[(ctype_degdf.pvals_adj < 0.01) & (ctype_degdf.logfoldchanges > 1)], ignore_index=True)
    return result


# In[13]:


HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '0') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)
sc.pp.log1p(HR_Rv_ctrl)
sc.tl.rank_genes_groups(HR_Rv_ctrl, groupby = 'Cell_type')
dg_ctrl = GetDiffGenes(HR_Rv_ctrl)


# In[14]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi)
sc.tl.rank_genes_groups(HR_Rv_3dpi, groupby = 'Cell_type')
dg_3dpi = GetDiffGenes(HR_Rv_3dpi)


# In[15]:


HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi, target_sum=1e4)
sc.pp.log1p(HR_Rv_7dpi)
sc.tl.rank_genes_groups(HR_Rv_7dpi, groupby = 'Cell_type')
dg_7dpi = GetDiffGenes(HR_Rv_7dpi)


# ## Annotate genes with secretome and count numbers differentially expressed

# In[29]:


dg_ctrl_counts = pd.DataFrame({'DG_count' : dg_ctrl.Cell_type.value_counts()})
diff_secretome_ctrl = dg_ctrl[dg_ctrl.Gene.isin(secretome.external_gene_name)]
ds_ctrl_counts = pd.DataFrame({'DSec_count' : diff_secretome_ctrl.Cell_type.value_counts()})
diff_ctrl_counts = pd.concat([dg_ctrl_counts, ds_ctrl_counts], axis = 1, join='inner')
diff_ctrl_counts['DSec_ratio'] = diff_ctrl_counts['DSec_count']/diff_ctrl_counts['DG_count']


# In[30]:


dg_3dpi_counts = pd.DataFrame({'DG_count' : dg_3dpi.Cell_type.value_counts()})
diff_secretome_3dpi = dg_3dpi[dg_3dpi.Gene.isin(secretome.external_gene_name)]
ds_3dpi_counts = pd.DataFrame({'DSec_count' : diff_secretome_3dpi.Cell_type.value_counts()})
diff_3dpi_counts = pd.concat([dg_3dpi_counts, ds_3dpi_counts], axis = 1, join='inner')
diff_3dpi_counts['DSec_ratio'] = diff_3dpi_counts['DSec_count']/diff_3dpi_counts['DG_count']


# In[31]:


dg_7dpi_counts = pd.DataFrame({'DG_count' : dg_7dpi.Cell_type.value_counts()})
diff_secretome_7dpi = dg_7dpi[dg_7dpi.Gene.isin(secretome.external_gene_name)]
ds_7dpi_counts = pd.DataFrame({'DSec_count' : diff_secretome_7dpi.Cell_type.value_counts()})
diff_7dpi_counts = pd.concat([dg_7dpi_counts, ds_7dpi_counts], axis = 1, join='inner')
diff_7dpi_counts['DSec_ratio'] = diff_7dpi_counts['DSec_count']/diff_7dpi_counts['DG_count']


# In[32]:


diff_ctrl_counts


# In[33]:


diff_3dpi_counts


# In[34]:


diff_7dpi_counts


# ## Secretome expression per cell type

# In[5]:


def CountSecretome(adata, secretome):
    secretome_ind = np.unique([adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index])
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


# In[6]:


def ExpressedSecretome(adata, secretome):
    secretome_ind = np.unique([adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index])
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


# In[7]:


HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '0') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)
secretome_ctrl = CountSecretome(HR_Rv_ctrl, secretome)
secretome_averages_ctrl, secretome_z_averages_ctrl = ExpressedSecretome(HR_Rv_ctrl, secretome)


# In[37]:


secretome_averages_ctrl_save = secretome_averages_ctrl[secretome_averages_ctrl.max(axis = 1) > 10]
secretome_averages_ctrl_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_ctrl.csv')
secretome_z_averages_ctrl_save = secretome_z_averages_ctrl[secretome_z_averages_ctrl.max(axis = 1) > 2]
secretome_z_averages_ctrl_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_ctrl.csv')


# In[22]:


#import matplotlib.pyplot as plt

secretome_ctrl_plot = secretome_ctrl[secretome_ctrl.Cell_type.isin(epifibro_types)] 
secretome_ctrl_plot = secretome_ctrl_plot.join(cell_type_colors, on = 'Cell_type', how = 'inner')
secretome_ctrl_plot = secretome_ctrl_plot.sort_values(by = 'Secretome', ascending = False)

pl.bar(x = np.arange(len(secretome_ctrl_plot)),
       height = secretome_ctrl_plot.Secretome/100,
       yerr = 3 * secretome_ctrl_plot.SEM/100, capsize = 2,
       color = secretome_ctrl_plot.color)
pl.xlabel('Cell type')
pl.xticks(np.arange(len(secretome_ctrl_plot)), secretome_ctrl_plot.Cell_type, rotation=270)
pl.ylabel('Secretome (%)')
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_ctrl_fibroniche.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[11]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi, target_sum=1e4)
secretome_3dpi = CountSecretome(HR_Rv_3dpi, secretome)
secretome_averages_3dpi, secretome_z_averages_3dpi = ExpressedSecretome(HR_Rv_3dpi, secretome)


# In[40]:


secretome_averages_3dpi_save = secretome_averages_3dpi[secretome_averages_3dpi.max(axis = 1) > 10]
secretome_averages_3dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_3dpi.csv')
secretome_z_averages_3dpi_save = secretome_z_averages_3dpi[secretome_z_averages_3dpi.max(axis = 1) > 2]
secretome_z_averages_3dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_3dpi.csv')


# In[23]:


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


# In[42]:


sns.clustermap(np.log1p(secretome_averages_3dpi[secretome_averages_3dpi.max(axis = 1) > 10]), method = 'ward', xticklabels = True, yticklabels = True)
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_high_genes.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[43]:


secretome_z_averages_3dpi_plot = secretome_z_averages_3dpi
secretome_z_averages_3dpi_plot[secretome_z_averages_3dpi_plot > 5] = 5

sns.clustermap(secretome_z_averages_3dpi_plot[secretome_z_averages_3dpi_plot.max(axis = 1) > 2], method = 'ward', xticklabels = True, yticklabels = True, figsize = (10,25))
pl.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_3dpi_high_z.pdf', bbox_inches = 'tight', transparent = True)
pl.show()


# In[44]:


HR_Rv_7dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi, target_sum=1e4)
secretome_7dpi = CountSecretome(HR_Rv_7dpi, secretome)
secretome_averages_7dpi, secretome_z_averages_7dpi = ExpressedSecretome(HR_Rv_7dpi, secretome)


# In[45]:


secretome_averages_7dpi_save = secretome_averages_7dpi[secretome_averages_7dpi.max(axis = 1) > 10]
secretome_averages_7dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_7dpi.csv')
secretome_z_averages_7dpi_save = secretome_z_averages_7dpi[secretome_z_averages_7dpi.max(axis = 1) > 2]
secretome_z_averages_7dpi_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_7dpi.csv')


# In[46]:


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


# In[47]:


HR_Rv_3dpi_inh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi_inh.var_names if not name == 'RFP']
HR_Rv_3dpi_inh = HR_Rv_3dpi_inh[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi_inh, target_sum=1e4)
secretome_3dpi_inh = CountSecretome(HR_Rv_3dpi_inh, secretome)
secretome_averages_3dpi_inh, secretome_z_averages_3dpi_inh = ExpressedSecretome(HR_Rv_3dpi_inh, secretome)


# In[48]:


secretome_averages_3dpi_inh_save = secretome_averages_3dpi_inh[secretome_averages_3dpi_inh.max(axis = 1) > 10]
secretome_averages_3dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_3dpi_inh.csv')
secretome_z_averages_3dpi_inh_save = secretome_z_averages_3dpi_inh[secretome_z_averages_3dpi_inh.max(axis = 1) > 2]
secretome_z_averages_3dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_3dpi_inh.csv')


# In[49]:


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


# In[50]:


HR_Rv_7dpi_inh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi_inh.var_names if not name == 'RFP']
HR_Rv_7dpi_inh = HR_Rv_7dpi_inh[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi_inh, target_sum=1e4)
secretome_7dpi_inh = CountSecretome(HR_Rv_7dpi_inh, secretome)
secretome_averages_7dpi_inh, secretome_z_averages_7dpi_inh = ExpressedSecretome(HR_Rv_7dpi_inh, secretome)


# In[51]:


secretome_averages_7dpi_inh_save = secretome_averages_7dpi_inh[secretome_averages_7dpi_inh.max(axis = 1) > 10]
secretome_averages_7dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_7dpi_inh.csv')
secretome_z_averages_7dpi_inh_save = secretome_z_averages_7dpi_inh[secretome_z_averages_7dpi_inh.max(axis = 1) > 2]
secretome_z_averages_7dpi_inh_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_7dpi_inh.csv')


# In[52]:


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
plt.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Secretome_Alliance_conversion_nodup_7dpi_inh.pdf', bbox_inches = 'tight', transparent = True)
plt.show()


# In[ ]:




