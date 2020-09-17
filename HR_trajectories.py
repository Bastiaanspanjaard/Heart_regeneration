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
import h5py
results_file = './write/HR_trajectories.h5ad'
sc.set_figure_params(dpi_save = 300)
import scvelo as scv
scv.settings.set_figure_params('scvelo', dpi_save = 300)


# In[2]:


annotations = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv', index_col = 0)


# In[3]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors_2.csv', index_col = 0)


# In[4]:


epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (cfd)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (cxcl12a)', 'Fibroblasts (mpeg1.1)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)',
                    'Valve fibroblasts', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[5]:


connected_epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)','Fibroblasts (cxcl12a)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[6]:


connected_epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells', 'Fibroblasts (cxcl12a)']


# In[7]:


connected_endofibro_types = ['Endocardium (frzb)', 'Endocardium (Ventricle)', 
                            'Fibroblasts (nppc)', 'Fibroblasts (spock3)']


# In[8]:


mito_genes = [line.rstrip('\n').rsplit(',')[2] for line in open('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/mito.genes.vs.txt')]


# In[10]:


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


# In[11]:


#HR_Rv_filter.write('./write/HR_Rv_filter.h5ad')
HR_Rv_filter = sc.read('./write/HR_Rv_filter.h5ad')


# In[ ]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]


# In[ ]:


HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]


# # Naive trajectory analysis at 3dpi

# In[196]:


HR_Rv_3dpi_epifibro = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(epifibro_types)]
sc.pp.filter_genes(HR_Rv_3dpi_epifibro, min_cells=3)


# In[197]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_epifibro.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_epifibro.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_epifibro[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_epifibro.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_epifibro, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_3dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[198]:


sc.pp.normalize_total(HR_Rv_3dpi_epifibro, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi_epifibro)


# In[199]:


sc.pp.regress_out(HR_Rv_3dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[200]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_epifibro)
sc.tl.pca(HR_Rv_3dpi_epifibro)
sc.external.pp.bbknn(HR_Rv_3dpi_epifibro, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_epifibro)


# In[201]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color='batch',
          title = '3dpi epicardium and fibroblasts')


# In[202]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epifibro.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '', save = 'subset1_3dpi_umap.png', legend_loc = 'none', frameon=False)


# In[203]:


sc.pl.umap(HR_Rv_3dpi_epifibro, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '3dpi epicardium and fibroblasts')


# In[204]:


sc.tl.paga(HR_Rv_3dpi_epifibro, groups='Cell_type')


# In[205]:


sc.pl.paga(HR_Rv_3dpi_epifibro, threshold=0.3, show=True, frameon=False,
           labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_3dpi_paga.png')


# # Naive trajectory analysis at 7dpi

# In[186]:


HR_Rv_7dpi_epifibro = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(epifibro_types)]
sc.pp.filter_genes(HR_Rv_7dpi_epifibro, min_cells=3)


# In[187]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_epifibro.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_epifibro.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_epifibro[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_epifibro.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_epifibro, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[188]:


sc.pp.normalize_total(HR_Rv_7dpi_epifibro, target_sum=1e4)
sc.pp.log1p(HR_Rv_7dpi_epifibro)


# In[189]:


sc.pp.regress_out(HR_Rv_7dpi_epifibro, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[190]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_epifibro)
sc.tl.pca(HR_Rv_7dpi_epifibro)
sc.external.pp.bbknn(HR_Rv_7dpi_epifibro, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_epifibro)


# In[191]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color='batch',
          title = '7dpi epicardium and fibroblasts')


# In[192]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_epifibro.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '', save = 'subset1_7dpi_umap.png', legend_loc = 'none', frameon=False)


# In[193]:


sc.pl.umap(HR_Rv_7dpi_epifibro, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '7dpi epicardium and fibroblasts')


# In[194]:


sc.tl.paga(HR_Rv_7dpi_epifibro, groups='Cell_type')


# In[195]:


sc.pl.paga(HR_Rv_7dpi_epifibro, threshold=0.3, show=True, frameon=False,
           labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_7dpi_paga.png')


# # Trajectories in 3dpi epicardial connected niche

# In[12]:


HR_Rv_3dpi_epiconn = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_3dpi_epiconn, min_cells=3)


# In[13]:


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


# In[21]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='batch',
          title = '3dpi epicardial niche')


# In[23]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[24]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color=['n_genes_by_counts'], 
          title = '3dpi epicardial niche')


# In[25]:


sc.tl.leiden(HR_Rv_3dpi_epiconn, resolution=2)


# In[26]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='leiden', legend_loc='on data', legend_fontsize='x-large')#,
#          palette = ['#92aad4', '#E9D723',  '#7b7f51', '#92aad4',
#                    '#E9D723', '#E9D723', '#E9D723', '#c6ba83',
#                    '#E9D723', '#c6ba83', '#525566', '#CE3A39',
#                    '#E9D723', '#CE3A39', '#e1e3d9', '#E9D723'])


# In[27]:


sc.tl.paga(HR_Rv_3dpi_epiconn, groups='leiden')


# In[28]:


sc.pl.paga(HR_Rv_3dpi_epiconn, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[29]:


sc.tl.umap(HR_Rv_3dpi_epiconn, init_pos = 'paga')


# In[30]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='batch',
          title = '3dpi epicardial niche')


# In[31]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[32]:


scv.pp.moments(HR_Rv_3dpi_epiconn)


# In[33]:


scv.tl.velocity(HR_Rv_3dpi_epiconn, mode='stochastic')


# In[38]:


scv.tl.velocity_graph(HR_Rv_3dpi_epiconn)


# In[49]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none',
                                 save = '_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[45]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'Cell_type')


# In[47]:


scv.tl.velocity_confidence(HR_Rv_3dpi_epiconn)


# In[48]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_3dpi_epiconn, c=keys, cmap='coolwarm', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[206]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6', 'nrg1']
scv.pl.velocity(HR_Rv_3dpi_epiconn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300,
                save = 'gene_velocity_epicardial_niche_3dpi_umap.png')


# # Trajectories in 7dpi epicardial connected niche

# In[143]:


HR_Rv_7dpi_epiconn = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_7dpi_epiconn, min_cells=3)


# In[144]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_epiconn.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_epiconn.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_epiconn[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_epiconn.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_epiconn, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[145]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_epiconn, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_epiconn)


# In[146]:


sc.pp.regress_out(HR_Rv_7dpi_epiconn, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[147]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_epiconn)
sc.tl.pca(HR_Rv_7dpi_epiconn)
sc.external.pp.bbknn(HR_Rv_7dpi_epiconn, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_epiconn)


# In[148]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='batch',
          title = '7dpi epicardial niche')


# In[150]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color=['total_counts', 'n_genes_by_counts', 'percent_mito'],
          title = '7dpi epicardial niche')


# In[149]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_epiconn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epicardial niche')


# In[151]:


sc.tl.leiden(HR_Rv_7dpi_epiconn, resolution=2)


# In[152]:


sc.pl.umap(HR_Rv_7dpi_epiconn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[153]:


sc.tl.paga(HR_Rv_7dpi_epiconn, groups='leiden')


# In[154]:


sc.pl.paga(HR_Rv_7dpi_epiconn, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[155]:


sc.tl.umap(HR_Rv_7dpi_epiconn, init_pos = 'paga')


# In[156]:


scv.pp.moments(HR_Rv_7dpi_epiconn)


# In[157]:


scv.tl.velocity(HR_Rv_7dpi_epiconn, mode='stochastic')


# In[158]:


scv.tl.velocity_graph(HR_Rv_7dpi_epiconn)


# In[159]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_epiconn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none',
                                 save = '_scvelo_epicardial_niche_7dpi_20n_cell_types_umap.png')


# In[207]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6','nrg1']
scv.pl.velocity(HR_Rv_7dpi_epiconn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300,
                save = 'gene_velocity_epicardial_niche_7dpi_umap.png')


# # Trajectories in control epicardial connected niche

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

# In[209]:


deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9'] # All samples with > 50 nppc fibroblasts.
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_endofibro_types)]
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['heart'].isin(deep_injury_libraries)]
scv.pp.remove_duplicate_cells(HR_Rv_7dpi_deep_endo)
sc.pp.filter_genes(HR_Rv_7dpi_deep_endo, min_cells=3)


# In[210]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_deep_endo.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_deep_endo.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_deep_endo[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_deep_endo.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_deep_endo, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[211]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_deep_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_deep_endo)


# In[212]:


sc.pp.regress_out(HR_Rv_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[213]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_deep_endo)
sc.tl.pca(HR_Rv_7dpi_deep_endo)
sc.external.pp.bbknn(HR_Rv_7dpi_deep_endo, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_deep_endo)


# In[214]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi deep injury endocardial niche')


# In[215]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '7dpi deep injury endocardial niche')


# In[216]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='batch',
          title = '7dpi deep injury endocardial niche')


# In[217]:


sc.tl.leiden(HR_Rv_7dpi_deep_endo, resolution=2)


# In[218]:


sc.tl.paga(HR_Rv_7dpi_deep_endo, groups='leiden')


# In[219]:


sc.pl.paga(HR_Rv_7dpi_deep_endo, show=True, color = 'Cell_type', node_size_scale = 2)#, 
    #       labels = ['', '', '', '', '', 
   #                  '', '', '', '', '', 
   #                  '', '', '', '', '', '', ''], 
  #         node_size_scale = 2)#, save = 'scvelo_deep_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[220]:


sc.tl.umap(HR_Rv_7dpi_deep_endo, init_pos='paga')


# In[221]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'leiden')


# In[222]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'batch')


# In[223]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = 'Cell_type')


# In[224]:


scv.pp.moments(HR_Rv_7dpi_deep_endo)


# In[225]:


scv.tl.velocity(HR_Rv_7dpi_deep_endo, mode='stochastic')


# In[226]:


scv.tl.velocity_graph(HR_Rv_7dpi_deep_endo)


# In[227]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, basis='umap', title = '', 
                                 color = 'Cell_type',
                                save = 'gene_velocity_endoventricular_deep_niche_7dpi_umap.png')


# In[138]:


scv.pl.velocity_embedding(HR_Rv_7dpi_deep_endo, basis='umap', title = '', scale = 0.3,
                          color = 'Cell_type', legend_loc = 'none')#,
#                         save = 'gene_velocity_endoventricular_deep_niche_7dpi_umap.png')


# In[141]:


scv.tl.velocity_confidence(HR_Rv_7dpi_deep_endo)


# In[142]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_deep_endo, c=keys, cmap='coolwarm', perc=[5, 95])#,
#                                save = 'strength_coherence_scvelo_deep_endo_niche_7dpi_20n_cell_types_umap.png')

