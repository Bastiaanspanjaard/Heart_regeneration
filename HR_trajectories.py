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


trajectory_subset = ['Fibroblasts (const.)', 'Fibroblasts (cfd)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (cxcl12a)', 'Fibroblasts (mpeg1.1)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)',
                    'Valve fibroblasts', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[5]:


connected_3dpi = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)','Fibroblasts (cxcl12a)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[6]:


connected_7dpi = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells', 'Fibroblasts (cxcl12a)']


# In[7]:


endo_7dpi = ['Endocardium (Atrium)', 'Endocardium (frzb)', 'Endocardium (Ventricle)', 'Fibroblasts (const.)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)', 'Valve fibroblasts']


# In[8]:


mito_genes = [line.rstrip('\n').rsplit(',')[2] for line in open('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/mito.genes.vs.txt')]


# In[9]:


#fibro_colors = cell_type_colors[['color', 'setFibro']]
#fibro_colors = fibro_colors[fibro_colors.notna().setFibro]
#fibro_colors = fibro_colors.set_index('setFibro')


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


# In[ ]:


#HR_setnames[HR_setnames.dpi == '30']


# # Load and annotate single-cell data

# In[ ]:


H5_data = sc.read_10x_h5('../Data/all_h5_transfer/H5_Dr11_cr31_ffbm.h5')
H5_data.var_names_make_unique()
H6_data = sc.read_10x_h5('../Data/all_h5_transfer/H6_Dr11_cr31_ffbm.h5')
H6_data.var_names_make_unique()
H7_data = sc.read_10x_h5('../Data/all_h5_transfer/H7_Dr11_cr31_ffbm.h5')
H7_data.var_names_make_unique()
H8a_data = sc.read_10x_h5('../Data/all_h5_transfer/H8a_Dr11_cr31_ffbm.h5')
H8a_data.var_names_make_unique()
H8v_data = sc.read_10x_h5('../Data/all_h5_transfer/H8v_Dr11_cr31_ffbm.h5')
H8v_data.var_names_make_unique()
Hr1_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr1_Dr11_cr31_ffbm.h5')
Hr1_data.var_names_make_unique()
Hr2a_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr2a_Dr11_cr31_ffbm.h5')
Hr2a_data.var_names_make_unique()
Hr2b_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr2b_Dr11_cr31_ffbm.h5')
Hr2b_data.var_names_make_unique()
Hr3_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr3_Dr11_cr31_ffbm.h5')
Hr3_data.var_names_make_unique()
Hr4_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr4_Dr11_cr31_ffbm.h5')
Hr4_data.var_names_make_unique()
Hr5_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr5_Dr11_cr31_ffbm.h5')
Hr5_data.var_names_make_unique()
Hr6a_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr6a_Dr11_cr31_ffbm.h5')
Hr6a_data.var_names_make_unique()
Hr6v_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr6v_Dr11_cr31_ffbm.h5')
Hr6v_data.var_names_make_unique()
Hr7a_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr7a_Dr11_cr31_ffbm.h5')
Hr7a_data.var_names_make_unique()
Hr7v_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr7v_Dr11_cr31_ffbm.h5')
Hr7v_data.var_names_make_unique()
Hr8_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr8_Dr11_cr31_ffbm.h5')
Hr8_data.var_names_make_unique()
Hr9_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr9_Dr11_cr31_ffbm.h5')
Hr9_data.var_names_make_unique()
Hr10_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr10_Dr11_cr31_ffbm.h5')
Hr10_data.var_names_make_unique()
Hr11_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr11_Dr11_cr31_ffbm.h5')
Hr11_data.var_names_make_unique()
Hr12_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr12_Dr11_cr31_ffbm.h5')
Hr12_data.var_names_make_unique()
Hr13_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr13_Dr11_cr31_ffbm.h5')
Hr13_data.var_names_make_unique()
Hr14_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr14_Dr11_cr31_ffbm.h5')
Hr14_data.var_names_make_unique()
Hr15_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr15_Dr11_cr31_ffbm.h5')
Hr15_data.var_names_make_unique()
Hr16_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr16_Dr11_cr31_ffbm.h5')
Hr16_data.var_names_make_unique()
Hr17_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr17_Dr11_cr31_ffbm.h5')
Hr17_data.var_names_make_unique()
Hr18_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr18_Dr11_cr31_ffbm.h5')
Hr18_data.var_names_make_unique()
Hr19_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr19_Dr11_cr31_ffbm.h5')
Hr19_data.var_names_make_unique()
Hr20_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr20_Dr11_cr31_ffbm.h5')
Hr20_data.var_names_make_unique()
Hr21_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr21_Dr11_cr31_ffbm.h5')
Hr21_data.var_names_make_unique()
Hr22_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr22_Dr11_cr31_ffbm.h5')
Hr22_data.var_names_make_unique()
Hr23_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr23_Dr11_cr31_ffbm.h5')
Hr23_data.var_names_make_unique()
Hr24_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr24_Dr11_cr31_ffbm.h5')
Hr24_data.var_names_make_unique()
Hr25_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr25_Dr11_cr31_ffbm.h5')
Hr25_data.var_names_make_unique()
Hr26_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr26_Dr11_cr31_ffbm.h5')
Hr26_data.var_names_make_unique()
Hr27_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr27_Dr11_cr31_ffbm.h5')
Hr27_data.var_names_make_unique()
Hr28_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr28_Dr11_cr31_ffbm.h5')
Hr28_data.var_names_make_unique()
Hr29_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr29_Dr11_cr31_ffbm.h5')
Hr29_data.var_names_make_unique()
Hr30_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr30_Dr11_cr31_ffbm.h5')
Hr30_data.var_names_make_unique()
Hr31_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr31_Dr11_cr31_ffbm.h5')
Hr31_data.var_names_make_unique()
Hr32_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr32_Dr11_cr31_ffbm.h5')
Hr32_data.var_names_make_unique()
Hr33_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr33_Dr11_cr31_ffbm.h5')
Hr33_data.var_names_make_unique()
Hr34_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr34_Dr11_cr31_ffbm.h5')
Hr34_data.var_names_make_unique()
Hr35_data = sc.read_10x_h5('../Data/all_h5_transfer/Hr35_Dr11_cr31_ffbm.h5')
Hr35_data.var_names_make_unique()


# In[ ]:


HR =    H5_data.concatenate(H6_data, H7_data, H8a_data, H8v_data,
                        Hr1_data, Hr2a_data, Hr2b_data, Hr3_data, Hr4_data, Hr5_data, 
                        Hr6a_data, Hr6v_data, Hr7a_data, Hr7v_data, Hr8_data, Hr9_data, Hr10_data,
                        Hr11_data, Hr12_data, Hr13_data, Hr14_data, Hr15_data,
                        Hr16_data, Hr17_data, Hr18_data, Hr19_data, Hr20_data,
                        Hr21_data, Hr22_data, Hr23_data, Hr24_data, Hr25_data,
                        Hr26_data, Hr27_data, Hr28_data, Hr29_data, Hr30_data,
                        Hr31_data, Hr32_data, Hr33_data, Hr34_data, Hr35_data)
HR.shape


# In[ ]:


HR_obs = HR.obs
HR.obs = HR.obs.reset_index().merge(HR_setnames, how="inner").set_index('index')


# In[ ]:


# Rename cells so the cell names correspond to the ones in the annotation file
HR.obs_names = [str(HR.obs.loc[x,'heart'])+'_'+str(x.split('-', 1)[0]) for x in HR.obs_names]


# In[ ]:


# Drop annotations that are not in the single-cell object
anno_drop = annotations.index.difference(HR.obs_names)
annotations = annotations.drop(anno_drop)


# In[ ]:


HR_filter = HR[annotations.index]
HR_filter


# In[ ]:


HR_filter.obs['Cell_type'] = annotations['Cell_type'].tolist()


# In[ ]:


#HR_filter.write('./write/HR_filter.h5ad')


# # Naive trajectory analysis at 3dpi

# In[ ]:


HR_ps_1 = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '3']
all_genes_but_RFP = [name for name in HR_ps_1.var_names if not name == 'RFP']
HR_ps_1 = HR_ps_1[:, all_genes_but_RFP]


# In[ ]:


sc.pp.filter_genes(HR_ps_1, min_cells=3)
HR_ps_1_norm = sc.pp.normalize_per_cell(HR_ps_1, counts_per_cell_after=1e4,copy=True)
HR_ps_1 = sc.pp.log1p(HR_ps_1_norm, copy=True)


# In[ ]:


HR_3d_sub1_cbr = HR_ps_1[HR_ps_1.obs['Cell_type'].isin(trajectory_subset)]
#HR_3d_sub1_cbr = sc.read('./write/HR_3dpi_subset1.h5ad')


# In[ ]:


sc.pp.highly_variable_genes(HR_3d_sub1_cbr)
sc.tl.pca(HR_3d_sub1_cbr)
#sc.pp.neighbors(HR_3d_sub1_cbr, n_neighbors=30)
sc.external.pp.bbknn(HR_3d_sub1_cbr, batch_key='batch')
sc.tl.umap(HR_3d_sub1_cbr)
# The order here is very important, to my surprise. If I run bbknn before neighbors, I do not get batch correction at all.


# In[ ]:


sc.pl.umap(HR_3d_sub1_cbr, color='batch',
          title = 'Subset 1, 3dpi', save = 'subset1_3dpi_umap_batch.png')


# In[ ]:


sc.pl.umap(HR_3d_sub1_cbr, color='Cell_type', palette = cell_type_colors.loc[HR_3d_sub1_cbr.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Subset 1, 3dpi')#, save = 'subset1_3dpi_umap.png')


# In[ ]:


sc.tl.paga(HR_3d_sub1_cbr, groups='Cell_type')


# In[ ]:


sc.pl.paga(HR_3d_sub1_cbr, threshold=0.3, show=True, frameon=False,labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_3dpi_paga.png')


# In[ ]:


sc.tl.umap(HR_3d_sub1_cbr, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_3d_sub1_cbr, color='Cell_type', title = '', legend_loc = 'none', frameon=False, save = 'subset1_3dpi_paga.png')


# In[ ]:


sc.tl.diffmap(HR_3d_sub1_cbr)


# In[ ]:


sc.pp.neighbors(HR_3d_sub1_cbr, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_3d_sub1_cbr)


# In[ ]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color='Cell_type', title = 'Subset 1, 3dpi', save = 'subset1_3dpi_20n_diffmap.png')


# In[ ]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color = 'batch', title = 'Subset 1, 3dpi', save = 'subset1_3dpi_20n_diffmap_batches.png')


# In[ ]:


sc.tl.paga(HR_3d_sub1_cbr, groups='Cell_type')


# In[ ]:


sc.pl.paga(HR_3d_sub1_cbr, threshold=0.3, show=True, labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_3dpi_20n_diffmap.png')


# In[ ]:


sc.tl.draw_graph(HR_3d_sub1_cbr, init_pos='paga')


# In[ ]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color = 'Cell_type', title = 'Subset 1, 3dpi', save = 'subset1_3dpi_20n_diffmap_pagainit.png')


# We can see in the PAGA that everything is connected, even cell types such as fibroblast-like cells and spock3 and nppc-fibroblasts that we know are not connected.  
# Maybe this has to do with the coarse grained clustering?

# In[ ]:


sc.tl.leiden(HR_3d_sub1_cbr, resolution=0.5)


# In[ ]:


sc.pl.umap(HR_3d_sub1_cbr, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.pl.umap(HR_3d_sub1_cbr, color='Cell_type')


# In[ ]:


sc.tl.paga(HR_3d_sub1_cbr, groups='leiden')


# In[ ]:


sc.tl.draw_graph(HR_3d_sub1_cbr, init_pos='paga')


# In[ ]:


sc.pl.paga(HR_3d_sub1_cbr, threshold = 0.2, color=['leiden'])


# In[ ]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color = ['Cell_type'])


# In[ ]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# After some fiddling with parameters, we manage to disconnect the fibroblast-like cells although they are still connected to another cluster of fibroblasts. However, the spock3 and nppc-fibroblasts remain connected, the perivascular cells are also disconnected, and a part of the atrial epicardium is now also disconnected.

# # Naive trajectory analysis at 7dpi

# In[ ]:


HR_7d = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_7d.var_names if not name == 'RFP']
HR_7d = HR_7d[:, all_genes_but_RFP]


# In[ ]:


sc.pp.filter_genes(HR_7d, min_cells=3)
HR_7d_norm = sc.pp.normalize_per_cell(HR_7d, counts_per_cell_after=1e4,copy=True)
HR_7d = sc.pp.log1p(HR_7d_norm, copy=True)


# In[ ]:


HR_7d_sub1_cbr = HR_7d[HR_7d.obs['Cell_type'].isin(trajectory_subset)]
#HR_3d_sub1_cbr = sc.read('./write/HR_3dpi_subset1.h5ad')


# In[ ]:


sc.pp.highly_variable_genes(HR_7d_sub1_cbr)
sc.tl.pca(HR_7d_sub1_cbr)
#sc.pp.neighbors(HR_7d_sub1_cbr, n_neighbors=30)
sc.external.pp.bbknn(HR_7d_sub1_cbr, batch_key='batch')
sc.tl.umap(HR_7d_sub1_cbr)
# The order here is very important, to my surprise. If I run bbknn before neighbors, I do not get batch correction at all.


# In[ ]:


sc.pl.umap(HR_7d_sub1_cbr, color='batch',
          title = 'Subset 1, 7dpi', save = 'subset1_7dpi_umap_batch.png')


# In[ ]:


sc.pl.umap(HR_7d_sub1_cbr, color='Cell_type', palette = cell_type_colors.loc[HR_7d_sub1_cbr.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Subset 1, 7dpi')#, save = 'subset1_7dpi_umap.png')


# In[ ]:


sc.tl.paga(HR_7d_sub1_cbr, groups='Cell_type')


# In[ ]:


sc.pl.paga(HR_7d_sub1_cbr, threshold=0.3, show=True, frameon=False,labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_7dpi_paga.png')


# In[ ]:


sc.tl.umap(HR_7d_sub1_cbr, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_7d_sub1_cbr, color='Cell_type', title = '', legend_loc = 'none', frameon=False, save = 'subset1_7dpi_paga.png')


# # RNA velocity

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


#HR_Rv.write('./write/HR_Rv.h5ad')


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
HR_Rv_filter


# In[ ]:


HR_Rv_filter.obs['Cell_type'] = annotations_Rv['Cell_type'].tolist()


# In[ ]:


#HR_Rv_filter.obs['Cell_type'].value_counts()
#HR_Rv_7dpi_deep_endo.obs['Cell_type'].unique()


# In[11]:


#HR_Rv_filter.write('./write/HR_Rv_filter.h5ad')
HR_Rv_filter = sc.read('./write/HR_Rv_filter.h5ad')


# In[ ]:


#del(HR_Rv_filter)


# # Trajectories in 3dpi epicardial connected niche

# In[12]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
HR_Rv_3dpi_conn = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(connected_3dpi)]
sc.pp.filter_genes(HR_Rv_3dpi_conn, min_cells=3)


# In[13]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_conn.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_conn.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_conn[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_conn.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_conn, percent_top=None, log1p=True, inplace=True)
#HR_Rv_3dpi_conn.obs['n_counts'] = HR_Rv_3dpi_conn.X.sum(axis=1)
sc.pl.violin(HR_Rv_3dpi_conn, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[14]:


#HR_Rv_3dpi_conn.obs


# In[15]:


sc.pp.normalize_total(HR_Rv_3dpi_conn, target_sum=1e4)
sc.pp.log1p(HR_Rv_3dpi_conn)


# In[16]:


sc.pp.regress_out(HR_Rv_3dpi_conn, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[17]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_conn)
sc.tl.pca(HR_Rv_3dpi_conn)
#sc.pp.neighbors(HR_Rv_3dpi_conn, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_3dpi_conn, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_conn)


# In[54]:


sc.pl.pca(HR_Rv_3dpi_conn, color = 'Cell_type', components=['1,2', '3,4'], legend_loc='none',
         palette = cell_type_colors.loc[HR_Rv_3dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[55]:


sc.pl.pca(HR_Rv_3dpi_conn, color = 'Cell_type', components=['5,6', '7,8'], legend_loc='none',
         palette = cell_type_colors.loc[HR_Rv_3dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[47]:


sc.pl.pca(HR_Rv_3dpi_conn, color = 'n_genes_by_counts', components=['1,2', '3,4'], legend_loc='none')


# In[17]:


sc.pl.umap(HR_Rv_3dpi_conn, color='batch',
          title = '3dpi epicardial niche')


# In[18]:


#sc.pl.umap(HR_Rv_3dpi_conn, color='batch',
#          title = '3dpi epicardial niche')


# In[19]:


sc.pl.umap(HR_Rv_3dpi_conn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[20]:


sc.pl.umap(HR_Rv_3dpi_conn, color=['n_genes_by_counts'], palette = cell_type_colors.loc[HR_Rv_3dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')


# In[ ]:


#sc.pl.umap(HR_Rv_3dpi_conn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
#          title = '3dpi epicardial niche')
#cell_type_colors.loc[HR_Rv_7dpi_endo.obs.Cell_type.cat.categories.tolist()].color.tolist()


# In[ ]:


#cell_type_colors


# In[ ]:


#sc.tl.diffmap(HR_Rv_3dpi_conn)
#sc.pp.neighbors(HR_Rv_3dpi_conn, n_neighbors=20, use_rep='X_diffmap')
#sc.tl.draw_graph(HR_Rv_3dpi_conn)


# In[21]:


sc.tl.leiden(HR_Rv_3dpi_conn, resolution=2)


# In[22]:


sc.pl.umap(HR_Rv_3dpi_conn, color='leiden', legend_loc='on data', legend_fontsize='x-large')#,
#          palette = ['#92aad4', '#E9D723',  '#7b7f51', '#92aad4',
#                    '#E9D723', '#E9D723', '#E9D723', '#c6ba83',
#                    '#E9D723', '#c6ba83', '#525566', '#CE3A39',
#                    '#E9D723', '#CE3A39', '#e1e3d9', '#E9D723'])


# In[23]:


sc.tl.paga(HR_Rv_3dpi_conn, groups='leiden')


# In[24]:


sc.pl.paga(HR_Rv_3dpi_conn, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[25]:


sc.tl.umap(HR_Rv_3dpi_conn, init_pos = 'paga')


# In[26]:


sc.pl.umap(HR_Rv_3dpi_conn, color='batch',
          title = '3dpi epicardial niche')


# In[27]:


sc.pl.umap(HR_Rv_3dpi_conn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche')
#cell_type_colors.loc[HR_Rv_7dpi_endo.obs.Cell_type.cat.categories.tolist()].color.tolist()


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_conn, color=['col1a1a', 'col11a1a', 'tcf21', 'tbx18'])


# In[ ]:


#sc.pl.umap(HR_Rv_3dpi_conn, color=['col1a1a', 'col11a1a', 'tcf21', 'tbx18'])


# In[ ]:


#sc.tl.draw_graph(HR_Rv_3dpi_conn, init_pos='paga')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_3dpi_conn, color = 'leiden', title = '3dpi epicardial niche', save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_3dpi_conn, color = ['Cell_type'], title = '3dpi epicardial niche', save = '_scvelo_connected_niche_3dpi_20n_cell_types_diffmap.png')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_3dpi_conn, ncols = 4,
#                 color = ['col1a1a', 'tbx18', 'tcf21', 'fn1a', 'fn1b',
#                          'notch3', 'pdgfrb', 'mpeg1.1', 'cfd', 
#                          'col11a1a', 'col12a1a', 
#                          'postnb', 'pcna', 'aldh1a2', 'stra6'],
#                save = '_scvelo_epicardial_niche_3dpi_20n_marker_genes_diffmap.png')


# Slightly different but very similar. I like the scanpy colors better for the heatmap.

# # RNA velocity on 3dpi connected

# In[28]:


#scv.pp.filter_and_normalize(HR_Rv_3dpi_conn)
scv.pp.moments(HR_Rv_3dpi_conn)


# In[29]:


scv.tl.velocity(HR_Rv_3dpi_conn, mode='stochastic')


# In[30]:


HR_Rv_3dpi_conn.layers['velocity'] # Array (cells) of arrays (gene velocities); 


# In[31]:


HR_Rv_3dpi_conn.layers['velocity'][0].max()


# In[38]:


HR_Rv_3dpi_conn


# In[34]:


scv.tl.velocity_graph(HR_Rv_3dpi_conn)


# In[35]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, basis = 'pca', color = 'Cell_type', legend_loc = 'none', components=('1,2')) # This needs velocity_graph to be computed


# In[41]:


get_ipython().run_line_magic('pinfo2', 'scv.tl.compute_velocity_embedding')


# In[ ]:


#scv.pl.velocity_embedding(HR_Rv_3dpi_conn, basis='draw_graph_fa')


# In[ ]:


#scv.pl.velocity_embedding_grid(HR_Rv_3dpi_conn, basis='draw_graph_fa', color = 'Cell_type')


# In[ ]:


#scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, basis='draw_graph_fa', title = '', 
#                                 color = 'Cell_type', legend_loc = 'none',
#                                save = '_scvelo_epicardial_niche_3dpi_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
#                                 save = '_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, title = '', 
                                 color = 'Cell_type')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_3dpi_conn, title = '', legend_loc = 'none', scale = 0.3,figsize=[16,9],
                                 color = 'Cell_type')


# In[ ]:


scv.tl.velocity_confidence(HR_Rv_3dpi_conn)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_3dpi_conn, c=keys, cmap='coolwarm', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# In[13]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6', 'nrg1', 'tcf21', 'tbx18']
scv.pl.velocity(HR_Rv_3dpi_conn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
#                save = 'gene_velocity_epicardial_niche_3dpi_umap.png')


# In[14]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, title = '', 
                                 color = 'leiden')


# In[19]:


sc.pl.pca(HR_Rv_3dpi_conn, color = 'Cell_type')


# In[93]:


scv.pl.velocity_embedding(HR_Rv_3dpi_conn, basis = 'pca', color = 'Cell_type', scale = 0.3, figsize=[16,9], components = '21,22')


# In[92]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, basis = 'pca', color = 'Cell_type', legend_loc = 'none', components=('1,2', '21,22'))


# In[23]:


sc.pl.pca(HR_Rv_3dpi_conn, color = 'Cell_type', components=['1,2', '3,4'], legend_loc='none')


# In[31]:


sc.pl.pca(HR_Rv_3dpi_conn, color = 'n_genes_by_counts', components=['1,2', '3,4'], legend_loc='none')


# In[ ]:





# In[ ]:


sc.tl.draw_graph(HR_Rv_3dpi_conn, init_pos = 'paga')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, title = '', basis='draw_graph_fa',
                                 color = 'Cell_type', legend_loc = 'none')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_3dpi_conn, color = 'Cell_type', basis='draw_graph_fa', figsize=(16,9), arrow_size=5, legend_loc='right')


# In[ ]:


#HR_Rv_3dpi_conn.write('./write/HR_Rv_3dpi_conn_regressed_filter.h5ad')
#HR_Rv_filter = sc.read('./write/HR_Rv_filter.h5ad')
del(HR_Rv_3dpi_conn)


# In[ ]:


scv.tl.rank_velocity_genes(HR_Rv_3dpi_conn, groupby='leiden', min_corr=.3)
leidenvelocities_HR_3dpi_conn = scv.DataFrame(HR_Rv_3dpi_conn.uns['rank_velocity_genes']['names'])


# In[ ]:


leidenvelocities_HR_3dpi_conn['2'].tolist()


# In[ ]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6', 'nrg1']
for gene in epi_velo_genes:
    scv.pl.velocity(HR_Rv_3dpi_conn, gene, basis='umap', ncols=1, fontsize=16,
                    save = gene+'_velocity_epicardial_niche_3dpi_20n_umap.png')


# In[ ]:


scv.pl.velocity(HR_Rv_3dpi_conn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300,
                save = 'gene_velocity_epicardial_niche_3dpi_umap.png')


# In[ ]:


scv.tl.paga(HR_Rv_3dpi_conn, groups='leiden')


# In[ ]:


#sc.tl.umap(HR_Rv_3dpi_conn, init_pos='paga')


# In[42]:


scv.tl.paga(HR_Rv_3dpi_conn, groups='leiden')
transitions_HR_3dpi_conn = scv.get_df(HR_Rv_3dpi_conn, 'paga/transitions_confidence', precision=2).T
transitions_HR_3dpi_conn.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[43]:


scv.pl.paga(HR_Rv_3dpi_conn, basis='umap', color = 'Cell_type', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, legend_loc = 'on data', title = '3dpi epicardial niche')#, 
  #         save = '_paga_scvelo_epicardial_niche_3dpi_20n_cell_types_diffmap.png')


# In[45]:


import louvain
scv.tl.velocity_clusters(HR_Rv_3dpi_conn)


# In[47]:


scv.pl.umap(HR_Rv_3dpi_conn, color='velocity_clusters', legend_loc='right')


# In[ ]:


sc.tl.paga(HR_Rv_3dpi_conn, groups='leiden', use_rna_velocity=True)


# In[ ]:


sc.pl.paga(HR_Rv_3dpi_conn, show=True, node_size_scale = 2)#,


# In[ ]:


scv.tl.score_genes_cell_cycle(HR_Rv_3dpi_conn)
scv.pl.scatter(HR_Rv_3dpi_conn, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95]) #, ]


# Intron/exon proportions per cluster - this is very strange; why would it be 9% across the board? It's also 9% for each sample. I hope this is just the plotting calculation going wrong and not an outcome of, say, the batch integration.

# In[ ]:


scv.pl.proportions(HR_Rv_3dpi_conn, groupby = 'leiden_anno')


# How are the speed and coherence?

# In[ ]:


scv.tl.velocity_confidence(HR_Rv_3dpi_conn)


# In[12]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_3dpi_conn, c=keys, cmap='coolwarm', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_3dpi_20n_cell_types_umap.png')


# Which genes drive the transitions?

# In[ ]:


scv.pl.velocity_graph(HR_Rv_3dpi_conn, color = 'Cell_type')


# In[ ]:


scv.pl.scatter(HR_Rv_3dpi_conn, c='col11a1a', vkey = 'velocity', 
               cmap='coolwarm', perc=[5, 95])


# In[ ]:


scv.pl.scatter(HR_Rv_3dpi_conn, 'col11a1a', c = ['Cell_type', 'leiden', 'velocity'])


# In[ ]:


scv.pl.velocity(HR_Rv_3dpi_conn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300,
                save = 'gene_velocity_epicardial_niche_3dpi_umap.png')


# In[11]:


#HR_Rv_3dpi_conn.write('./write/HR_Rv_3dpi_conn.h5ad')
HR_Rv_3dpi_conn = sc.read('./write/HR_Rv_3dpi_conn.h5ad')


# In[ ]:


#scv.tl.rank_velocity_genes(HR_Rv_3dpi_conn, groupby='leiden_anno', min_corr=.3)


# In[ ]:


#HR_3dpi_conn_velogenes = scv.DataFrame(HR_Rv_3dpi_conn.uns['rank_velocity_genes']['names'])
#HR_3dpi_conn_velogenes.head()


# In[ ]:


#scv.pl.scatter(HR_Rv_3dpi_conn, HR_3dpi_conn_velogenes['FB 1'][:5], ylabel='FB 1', color = 'leiden_anno', add_outline='FB 1, FB 3, FB 5, EpiC V')


# # Confirm velocity in single sample and different visualizations

# In[ ]:


HR_Rv_3dpi_conn.obs.heart.value_counts()


# ## Hr10

# In[ ]:


Hr10_Rv = scv.read('../Data/RNAvelo/Hr10_v3Dr11.loom', cache=True)
Hr10_Rv.var_names_make_unique()


# In[ ]:


# Rename cells to match cell names in annotation file
Hr10_Rv.obs['heart'] = "Hr10"
Hr10_Rv.obs_names = [str(Hr10_Rv.obs.loc[x,'heart'])+'_'+str(x.split(':', 1)[1])[0:16] for x in Hr10_Rv.obs_names]
# Drop annotations that are not in the single-cell object
anno_drop_Hr10 = annotations.index.difference(Hr10_Rv.obs_names)
annotations_Hr10 = annotations.drop(anno_drop_Hr10)
Hr10_Rv_filter = Hr10_Rv[annotations_Hr10.index]
Hr10_Rv_filter.obs['Cell_type'] = annotations_Hr10['Cell_type'].tolist()


# In[ ]:


Hr10_Rv_filter.obs['Cell_type'].value_counts()


# In[ ]:


all_genes_but_RFP = [name for name in Hr10_Rv_filter.var_names if not name == 'RFP']
Hr10_Rv_filter = Hr10_Rv_filter[:, all_genes_but_RFP]
Hr10_Rv_epi = Hr10_Rv_filter[Hr10_Rv_filter.obs['Cell_type'].isin(connected_3dpi)]
sc.pp.filter_genes(Hr10_Rv_epi, min_cells=3)
sc.pp.normalize_per_cell(Hr10_Rv_epi, counts_per_cell_after=1e4)
sc.pp.log1p(Hr10_Rv_epi)
sc.pp.highly_variable_genes(Hr10_Rv_epi)
sc.tl.pca(Hr10_Rv_epi)
sc.pp.neighbors(Hr10_Rv_epi, n_neighbors=30)
sc.tl.umap(Hr10_Rv_epi)


# In[ ]:


sc.pl.umap(Hr10_Rv_epi, color='Cell_type', palette = cell_type_colors.loc[Hr10_Rv_epi.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche in Hr10')


# In[ ]:


sc.tl.diffmap(Hr10_Rv_epi)
sc.pp.neighbors(Hr10_Rv_epi, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(Hr10_Rv_epi)
sc.tl.leiden(Hr10_Rv_epi, resolution=0.4)


# In[ ]:


Hr10_Rv_epi = Hr10_Rv_epi[Hr10_Rv_epi.obs.leiden != '15']


# In[ ]:


sc.pl.umap(Hr10_Rv_epi, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(Hr10_Rv_epi, groups='leiden')


# In[ ]:


sc.pl.paga(Hr10_Rv_epi, show=True, 
           labels = ['', '', '', '', '', 
                     '', '', '', '', '', 
                     '', '', '', '', '',
                    ''], node_size_scale = 2, 
           save = 'scvelo_epicardial_niche_Hr10_20n_leiden_diffmap.png')


# In[ ]:


sc.tl.draw_graph(Hr10_Rv_epi, init_pos='paga')


# In[ ]:


sc.pl.draw_graph(Hr10_Rv_epi, color = 'leiden', title = 'Hr10 epicardial niche', save = 'scvelo_epicardial_niche_Hr10_20n_leiden_diffmap.png')


# In[ ]:


sc.pl.draw_graph(Hr10_Rv_epi, color = ['Cell_type'], title = 'Hr10 epicardial niche', save = '_scvelo_connected_niche_Hr10_20n_cell_types_diffmap.png')


# In[ ]:


scv.pp.moments(Hr10_Rv_epi)
scv.tl.velocity(Hr10_Rv_epi, mode='stochastic')
scv.tl.velocity_graph(Hr10_Rv_epi)


# In[ ]:


scv.pl.velocity_embedding_stream(Hr10_Rv_epi, basis='draw_graph_fa', title = 'Hr10 epicardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_epicardial_niche_Hr10_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(Hr10_Rv_epi, title = 'Hr10 epicardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                 save = '_scvelo_epicardial_niche_Hr10_20n_cell_types_umap.png')


# Intron/exon proportions per cluster - this is very strange; why would it be 9% across the board? It's also 9% for each sample. I hope this is just the plotting calculation going wrong and not an outcome of, say, the batch integration.

# How are the speed and coherence?

# In[ ]:


scv.tl.velocity_confidence(Hr10_Rv_epi)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(Hr10_Rv_epi, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_Hr10_20n_cell_types_diffmap.png')


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(Hr10_Rv_epi, c=keys, cmap='coolwarm', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_Hr10_20n_cell_types_umap.png')


# ## Hr26

# In[ ]:


Hr26_Rv = scv.read('../Data/RNAvelo/Hr26_v3Dr11.loom', cache=True)
Hr26_Rv.var_names_make_unique()


# In[ ]:


# Rename cells to match cell names in annotation file
Hr26_Rv.obs['heart'] = "Hr26"
Hr26_Rv.obs_names = [str(Hr26_Rv.obs.loc[x,'heart'])+'_'+str(x.split(':', 1)[1])[0:16] for x in Hr26_Rv.obs_names]
# Drop annotations that are not in the single-cell object
anno_drop_Hr26 = annotations.index.difference(Hr26_Rv.obs_names)
annotations_Hr26 = annotations.drop(anno_drop_Hr26)
Hr26_Rv_filter = Hr26_Rv[annotations_Hr26.index]
Hr26_Rv_filter.obs['Cell_type'] = annotations_Hr26['Cell_type'].tolist()


# In[ ]:


Hr26_Rv_filter.obs['Cell_type'].value_counts()


# In[ ]:


all_genes_but_RFP = [name for name in Hr26_Rv_filter.var_names if not name == 'RFP']
Hr26_Rv_filter = Hr26_Rv_filter[:, all_genes_but_RFP]
Hr26_Rv_epi = Hr26_Rv_filter[Hr26_Rv_filter.obs['Cell_type'].isin(connected_3dpi)]
sc.pp.filter_genes(Hr26_Rv_epi, min_cells=3)
sc.pp.normalize_per_cell(Hr26_Rv_epi, counts_per_cell_after=1e4)
sc.pp.log1p(Hr26_Rv_epi)
sc.pp.highly_variable_genes(Hr26_Rv_epi)
sc.tl.pca(Hr26_Rv_epi)
sc.pp.neighbors(Hr26_Rv_epi, n_neighbors=30)
sc.tl.umap(Hr26_Rv_epi)


# In[ ]:


sc.pl.umap(Hr26_Rv_epi, color='Cell_type', palette = cell_type_colors.loc[Hr26_Rv_epi.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi epicardial niche in Hr26')


# In[ ]:


sc.tl.diffmap(Hr26_Rv_epi)
sc.pp.neighbors(Hr26_Rv_epi, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(Hr26_Rv_epi)
sc.tl.leiden(Hr26_Rv_epi, resolution=0.4)


# In[ ]:


sc.pl.umap(Hr26_Rv_epi, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(Hr26_Rv_epi, groups='leiden')


# In[ ]:


sc.pl.paga(Hr26_Rv_epi, show=True, 
           labels = ['', '', '', '', '', 
                     '', '', '', '', '', 
                     '', '', '', '', ''], node_size_scale = 2, 
           save = 'scvelo_epicardial_niche_Hr26_20n_leiden_diffmap.png')


# In[ ]:


sc.tl.draw_graph(Hr26_Rv_epi, init_pos='paga')


# In[ ]:


sc.pl.draw_graph(Hr26_Rv_epi, color = 'leiden', title = 'Hr26 epicardial niche', save = 'scvelo_epicardial_niche_Hr26_20n_leiden_diffmap.png')


# In[ ]:


sc.pl.draw_graph(Hr26_Rv_epi, color = ['Cell_type'], title = 'Hr26 epicardial niche', save = '_scvelo_connected_niche_Hr26_20n_cell_types_diffmap.png')


# In[ ]:


scv.pp.moments(Hr26_Rv_epi)
scv.tl.velocity(Hr26_Rv_epi, mode='stochastic')
scv.tl.velocity_graph(Hr26_Rv_epi)


# In[ ]:


scv.pl.velocity_embedding_stream(Hr26_Rv_epi, basis='draw_graph_fa', title = 'Hr26 epicardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_epicardial_niche_Hr26_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(Hr26_Rv_epi, title = 'Hr26 epicardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                 save = '_scvelo_epicardial_niche_Hr26_20n_cell_types_umap.png')


# How are the speed and coherence?

# In[ ]:


scv.tl.velocity_confidence(Hr26_Rv_epi)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(Hr26_Rv_epi, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_Hr26_20n_cell_types_diffmap.png')


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(Hr26_Rv_epi, c=keys, cmap='coolwarm', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_Hr26_20n_cell_types_umap.png')


# # Epi niche at control

# In[ ]:


#HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
#all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
#HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
#HR_Rv_3dpi_conn = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(connected_3dpi)]
#sc.pp.filter_genes(HR_Rv_3dpi_conn, min_cells=3)


# In[ ]:


# Find mito_genes in dataset.
#mito_in_index = list(set(HR_Rv_3dpi_conn.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
#HR_Rv_3dpi_conn.obs['percent_mito'] = np.sum(
#    HR_Rv_3dpi_conn[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_conn.X, axis=1)
#sc.pp.calculate_qc_metrics(HR_Rv_3dpi_conn, percent_top=None, log1p=True, inplace=True)
#HR_Rv_3dpi_conn.obs['n_counts'] = HR_Rv_3dpi_conn.X.sum(axis=1)
#sc.pl.violin(HR_Rv_3dpi_conn, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
#             jitter=0.4, multi_panel=True, size = 0.1)


# In[ ]:


#HR_Rv_3dpi_conn.obs


# In[ ]:


# HERE
HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'].isin(['0', '3'])) & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
HR_Rv_ctrl_epi = HR_Rv_ctrl[HR_Rv_ctrl.obs['Cell_type'].isin(connected_3dpi)]
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


#HR_Rv_ctrl_epi


# In[ ]:


sc.pp.filter_genes(HR_Rv_ctrl_epi, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_ctrl_epi, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_ctrl_epi)
sc.pp.highly_variable_genes(HR_Rv_ctrl_epi)
sc.tl.pca(HR_Rv_ctrl_epi)
sc.external.pp.bbknn(HR_Rv_ctrl_epi, batch_key='batch')#, neighbors_within_batch = 3, n_pcs = 20)
#sc.pp.neighbors(HR_Rv_ctrl_epi, n_neighbors=30)
sc.tl.umap(HR_Rv_ctrl_epi)


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epi, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_ctrl_epi.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          )#save = 'Epicardial_niche_control_3dpi_Celltypes.png')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epi, color='dpi',
          )#save = 'Epicardial_niche_control_3dpi_timepoints.png')


# In[ ]:


#sc.pl.umap(HR_Rv_ctrl_epi, color='batch',
 #         title = 'Epicardial niche at control')


# In[ ]:


#sc.pl.umap(HR_Rv_ctrl_epi, color=['tcf21', 'col1a1a', 'tbx18'])


# In[ ]:


#sc.tl.diffmap(HR_Rv_ctrl_epi)
#sc.pp.neighbors(HR_Rv_ctrl_epi, n_neighbors=20, use_rep='X_diffmap')
#sc.tl.draw_graph(HR_Rv_ctrl_epi)
sc.tl.leiden(HR_Rv_ctrl_epi, resolution=2)


# In[ ]:


#sc.pl.umap(HR_Rv_ctrl_epi, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_ctrl_epi, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_ctrl_epi, show=True, color ='Cell_type')
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', ''], node_size_scale = 2, 
#           save = 'scvelo_epicardial_niche_control_20n_leiden_diffmap.png')


# In[ ]:


sc.pl.paga(HR_Rv_ctrl_epi, show=True, color ='dpi')


# In[ ]:


sc.tl.umap(HR_Rv_ctrl_epi, init_pos = 'paga')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epi, color='Cell_type',
          save = 'Epicardial_niche_control_3dpi_Celltypes.png')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epi, color='dpi',
          save = 'Epicardial_niche_control_3dpi_timepoints.png')


# In[ ]:


#sc.tl.draw_graph(HR_Rv_ctrl_epi, init_pos='paga')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_ctrl_epi, color = 'leiden', title = 'Control epicardial niche', 
#                 save = 'scvelo_epicardial_niche_control_20n_leiden_diffmap.png')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_ctrl_epi, color = ['Cell_type'], title = 'Control epicardial niche', 
#                 save = '_scvelo_connected_niche_control_20n_cell_types_diffmap.png')


# In[ ]:


scv.pp.moments(HR_Rv_ctrl_epi)
scv.tl.velocity(HR_Rv_ctrl_epi, mode='stochastic')
scv.tl.velocity_graph(HR_Rv_ctrl_epi)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_ctrl_epi, title = '', 
                                 color = 'Cell_type', legend_loc = 'none',
                                save = 'Epicardial_niche_control_3dpi_Celltypes_velo.png')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_ctrl_epi, title = '', scale = 0.2,
                                 color = 'Cell_type', legend_loc = 'none')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_ctrl_epi[HR_Rv_ctrl_epi.obs['dpi'] == '0'], title = '', scale = 0.2,
                                color = 'dpi', legend_loc = 'none',
                         save = 'Epicardial_niche_control_3dpi_control_velo.png')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_ctrl_epi[HR_Rv_ctrl_epi.obs['dpi'] == '3'], title = '', scale = 0.2,
                                color = 'dpi', legend_loc = 'none',
                         save = 'Epicardial_niche_control_3dpi_3dpi_velo.png')


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epi, color='leiden',
          save = 'Epicardial_niche_control_3dpi_leidenclusters.png')


# In[ ]:


HR_Rv_ctrl_epi[HR_Rv_ctrl_epi.obs['dpi'] == '3'].obs.leiden.value_counts()/sum(HR_Rv_ctrl_epi[HR_Rv_ctrl_epi.obs['dpi'] == '3'].obs.leiden.value_counts())


# In[ ]:


HR_Rv_ctrl_epi[HR_Rv_ctrl_epi.obs['dpi'] == '0'].obs.leiden.value_counts()/sum(HR_Rv_ctrl_epi[HR_Rv_ctrl_epi.obs['dpi'] == '0'].obs.leiden.value_counts())


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epi[HR_Rv_ctrl_epi.obs['dpi'] == '0'], color = 'Cell_type',
          save = 'Epicardial_niche_control_3dpi_Celltypes_control_only.png')


# In[ ]:


scv.tl.score_genes_cell_cycle(HR_Rv_ctrl_epi)
scv.pl.scatter(HR_Rv_ctrl_epi, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])


# In[ ]:


sc.pp.calculate_qc_metrics(HR_Rv_ctrl_epi, percent_top=None, log1p=False, inplace=True)
#HR_Rv_ctrl_epi.obs


# In[ ]:


HR_Rv_ctrl_epi.obs


# In[ ]:


sc.pl.umap(HR_Rv_ctrl_epi, color = ['n_genes_by_counts', 'n_counts'])


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_ctrl_epi, basis='draw_graph_fa', title = '', 
                                 color = 'Cell_type', legend_loc = 'none',
                                save = '_scvelo_epicardial_niche_control_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_ctrl_epi, title = '', 
                                 color = 'Cell_type', legend_loc = 'none',
                                 save = '_scvelo_epicardial_niche_control_20n_cell_types_umap.png')


# How are the speed and coherence?

# In[ ]:


scv.tl.velocity_confidence(HR_Rv_ctrl_epi)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_ctrl_epi, c=keys, cmap='coolwarm', basis='umap', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_control3dpi.png')


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_ctrl_epi, c=keys, cmap='coolwarm', perc=[5, 95],
                                save = 'strength_coherence_scvelo_epicardial_niche_Hr26_20n_cell_types_umap.png')


# # Same niche at 7dpi

# In[ ]:


HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
HR_Rv_7dpi_conn = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_7dpi_conn, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_7dpi_conn, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_conn)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_conn)
sc.tl.pca(HR_Rv_7dpi_conn)
#sc.pp.neighbors(HR_Rv_7dpi_conn, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_7dpi_conn, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_conn)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_conn, color='batch',
          title = '7dpi epicardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_conn, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epicardial niche')


# In[ ]:


#sc.tl.diffmap(HR_Rv_7dpi_conn)
#sc.pp.neighbors(HR_Rv_7dpi_conn, n_neighbors=20, use_rep='X_diffmap')
#sc.tl.draw_graph(HR_Rv_7dpi_conn)


# In[ ]:


sc.tl.leiden(HR_Rv_7dpi_conn, resolution=2)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_conn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_7dpi_conn, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_7dpi_conn, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[ ]:


sc.tl.umap(HR_Rv_7dpi_conn, init_pos = 'paga')


# In[ ]:


scv.pp.moments(HR_Rv_7dpi_conn)


# In[ ]:


scv.tl.velocity(HR_Rv_7dpi_conn, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_7dpi_conn)


# In[ ]:


#scv.pl.velocity_embedding(HR_Rv_7dpi_conn, basis='draw_graph_fa')


# In[ ]:


#scv.pl.velocity_embedding_grid(HR_Rv_7dpi_conn, basis='draw_graph_fa', color = 'Cell_type')


# In[ ]:


#scv.pl.velocity_embedding_stream(HR_Rv_7dpi_conn, basis='draw_graph_fa', title = '', 
#                                 color = 'Cell_type', legend_loc = 'none',
#                                save = '_scvelo_connected_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_conn, title = '', 
                                 color = 'Cell_type', legend_loc = 'none',
                                 save = '_scvelo_epicardial_niche_7dpi_20n_cell_types_umap.png')


# In[ ]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'stra6', 'fn1a', 'nrg1', 'postnb']
for gene in epi_velo_genes:
    scv.pl.velocity(HR_Rv_7dpi_conn, gene, basis='umap', ncols=1, fontsize=16,
                    save = gene+'_velocity_epicardial_niche_7dpi_20n_umap.png')


# In[ ]:


scv.pl.velocity(HR_Rv_7dpi_conn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300,
                save = 'gene_velocity_epicardial_niche_7dpi_umap.png')


# In[ ]:


#epi_velo_genes = ['col11a1a', 'col12a1a', 'stra6', 'fn1a', 'nrg1', 'postnb']
scv.pl.velocity(HR_Rv_7dpi_conn, epi_velo_genes, basis='umap', ncols=1, fontsize=16, layers = '')#,
#                    save = gene+'_velocity_epicardial_niche_7dpi_20n_umap.png')
scv.pl.velocity(HR_Rv_7dpi_conn, ['col11a1a', 'col12a1a', 'stra6', 'fn1a', 'nrg1'], basis='umap', ncols=1, fontsize=12,
               save = 'gene_velocity_epicardial_niche_3dpi_20n_umap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_conn, basis='draw_graph_fa', title = '7dpi connected niche', 
                                 color = 'leiden', legend_loc = 'right margin')


# In[ ]:


HR_Rv_7dpi_conn.write('./write/HR_Rv_7dpi_conn.h5ad')


# # Epicardial niche at all timepoints

# In[ ]:


HR_Rv_epi = HR_Rv_filter[HR_Rv_filter.obs['Cell_type'].isin(connected_3dpi)]
all_genes_but_RFP = [name for name in HR_Rv_epi.var_names if not name == 'RFP']
HR_Rv_epi = HR_Rv_epi[:, all_genes_but_RFP]
#HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
#HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
#HR_Rv_7dpi_conn = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_epi, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_epi, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_epi)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_epi)
sc.tl.pca(HR_Rv_epi)
#sc.pp.neighbors(HR_Rv_epi, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_epi, batch_key='batch')
sc.tl.umap(HR_Rv_epi)


# In[ ]:


sc.pl.umap(HR_Rv_epi, color='batch',
          title = '7dpi epicardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_epi, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_epi.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Epicardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_epi[HR_Rv_epi.obs['dpi'] == '0'], color='Cell_type', #palette = cell_type_colors.loc[HR_Rv_epi.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Epicardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_epi[HR_Rv_epi.obs['dpi'] == '0'], color = ['tcf21'], size = 5)


# In[ ]:


sc.pl.heatmap(HR_Rv_epi[(HR_Rv_epi.obs['dpi'] == '0') & HR_Rv_epi.obs['Cell_type'].isin(['Epicardium (Atrium)', 'Epicardium (Ventricle)', 'Fibroblasts (const.)'])],
              var_names = ['tcf21'], groupby='Cell_type',
             save='Epicardium_fibroblast_tcf21_control.png')


# In[ ]:


sc.pl.dotplot(HR_Rv_filter,
              var_names = ['pdgfrb', 'notch3'], groupby='Cell_type')


# # Endocardial niche at 7dpi

# In[ ]:


endo_7dpi = ['Endocardium (Atrium)', 'Endocardium (frzb)', 'Endocardium (Ventricle)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)']
HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
HR_Rv_7dpi_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_7dpi_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_7dpi_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_endo)
sc.tl.pca(HR_Rv_7dpi_endo)
#sc.pp.neighbors(HR_Rv_7dpi_endo, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_7dpi_endo, batch_key='batch', n_pcs = 30)
sc.tl.umap(HR_Rv_7dpi_endo)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_endo, color='batch',
          title = '7dpi endocardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi endocardial niche')


# In[ ]:


#sc.tl.diffmap(HR_Rv_7dpi_endo)
#sc.pp.neighbors(HR_Rv_7dpi_endo, n_neighbors=20, use_rep='X_diffmap')
#sc.tl.draw_graph(HR_Rv_7dpi_endo)


# In[ ]:


sc.tl.leiden(HR_Rv_7dpi_endo, resolution=2)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_7dpi_endo, groups='Cell_type')


# In[ ]:


sc.pl.paga(HR_Rv_7dpi_endo, show=True, colors = 'Cell_type')
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    '',''], 
#           node_size_scale = 2, save = '_scvelo_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[ ]:


sc.tl.umap(HR_Rv_7dpi_endo, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_endo, color='Cell_type')


# In[ ]:


sc.tl.draw_graph(HR_Rv_7dpi_endo, init_pos='paga')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_7dpi_endo, color = 'leiden', title = '7dpi endocardial niche', save = 'scvelo_endocardial_niche_7dpi_20n_leiden_diffmap.png')


# In[ ]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, color = ['Cell_type'], title = '7dpi endocardial niche')#, save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_7dpi_endo, ncols = 3,
#                 color = ['col1a1a', 'fn1b', 'spock3',
#                          'nppc', 'cyp26b1', 'frzb', 'acta2'],
#                save = '_scvelo_endocardial_niche_7dpi_20n_marker_genes_diffmap.png')


# In[ ]:


#scv.pp.remove_duplicate_cells(HR_Rv_7dpi_endo)
scv.pp.moments(HR_Rv_7dpi_endo)


# In[ ]:


scv.tl.velocity(HR_Rv_7dpi_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_7dpi_endo)


# In[ ]:


#scv.pl.velocity_embedding(HR_Rv_7dpi_endo, basis='draw_graph_fa')


# In[ ]:


#scv.pl.velocity_embedding_grid(HR_Rv_7dpi_endo, basis='draw_graph_fa', color = 'Cell_type')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endo, basis='umap', title = '', 
                                 color = 'Cell_type',
                                save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity(HR_Rv_7dpi_endo, ['aldh1a2', 'krt18', 'kdrl', 'cdh5', 'fli1a', 'nppc', 'spock3'], basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
#                save = 'gene_velocity_epicardial_niche_3dpi_umap.png')


# In[ ]:


scv.tl.rank_velocity_genes(HR_Rv_7dpi_endo, groupby='Cell_type', min_corr=.3)

endo_celltype_velogenes = scv.DataFrame(HR_Rv_7dpi_endo.uns['rank_velocity_genes']['names'])
endo_celltype_velogenes.head()


# In[ ]:


scv.pl.velocity(HR_Rv_7dpi_endo, ['pdgfrb', 'clcf1', 'tspan35'], basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,


# In[ ]:


sc.tl.rank_genes_groups(HR_Rv_7dpi_endo, groupby='Cell_type')


# In[ ]:


HR_Rv_7dpi_endo.uns['rank_genes_groups']['names']['Fibroblasts (nppc)']


# In[ ]:


scv.pl.velocity(HR_Rv_7dpi_endo, ['nppc', 'vmp1', 'serpine1'], basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,


# How are the speed and coherence?

# In[ ]:


scv.tl.velocity_confidence(HR_Rv_7dpi_endo)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_endo, c=keys, cmap='coolwarm', basis='umap', perc=[5, 95])#,
#                                save = 'strength_coherence_scvelo_endo_niche_7dpi_20n_cell_types_diffmap.png')


# # Nppc fibroblasts surrounded by other cells

# In[ ]:


annotations.Cell_type.value_counts()


# In[ ]:


epithelium_fibro = ['Endocardium (Atrium)', 'Endocardium (Ventricle)', 'Fibroblasts (const.)', 
                   'Endocardium (frzb)', 'Bl.ves.EC (apnln)', 'Bl.ves.EC (plvapb)', 'Epicardium (Atrium)',
                   'Epicardium (Ventricle)', 'Valve fibroblasts', 'Fibroblasts (col12a1a)', 'Fibroblasts (cxcl12a)',
                   'Fibroblasts (col11a1a)', 'Fibroblasts (cfd)', 'Bl.ves.EC (lyve1)', 'Perivascular cells',
                   'Fibroblasts (nppc)', 'Fibroblasts (spock3)', 'Fibroblasts (mpeg1.1)', 'Fibroblasts (proliferating)']
HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
HR_Rv_7dpi_ef = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(epithelium_fibro)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_7dpi_ef, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_7dpi_ef, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_ef)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_ef)
sc.tl.pca(HR_Rv_7dpi_ef)
#sc.pp.neighbors(HR_Rv_7dpi_ef, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_7dpi_ef, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_ef)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_ef, color='batch',
          title = '7dpi epithelial and fibroblasts')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_ef, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_ef.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epithelial and fibroblasts')


# In[ ]:


#sc.tl.leiden(HR_Rv_7dpi_ef, resolution=2)


# In[ ]:


#sc.pl.umap(HR_Rv_7dpi_ef, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_7dpi_ef, groups='Cell_type')


# In[ ]:


sc.pl.paga(HR_Rv_7dpi_ef, show=True, color = 'Cell_type',
          node_size_scale = 2, threshold = 0.5)
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    '',''], 
#           node_size_scale = 2, save = '_scvelo_ef_niche_7dpi_20n_leiden_diffmap.png')


# In[ ]:


sc.tl.umap(HR_Rv_7dpi_ef, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_ef, color='Cell_type')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_ef[HR_Rv_7dpi_ef.obs['Cell_type'].isin(['Fibroblasts (const.)', 'Fibroblasts (spock3)',
                                                             'Endocardium (frzb)'])], ncols = 3,
                 color = ['col1a1a', 'spock3',
                          'nppc', 'wif1', 'notum1b'])#,
#                save = '_scvelo_endocardial_niche_7dpi_20n_marker_genes_diffmap.png')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_ef[HR_Rv_7dpi_ef.obs['Cell_type'].isin(['Fibroblasts (const.)', 'Fibroblasts (spock3)',
                                                             'Endocardium (frzb)'])], color = 'Cell_type')


# # Endocardial niche for deep injuries at 7dpi

# In[ ]:


#annotations[annotations.Cell_type == 'Fibroblasts (nppc)']['orig.ident'].value_counts()


# In[ ]:


#annotations['orig.ident'].value_counts()


# In[ ]:


endo_7dpi = ['Endocardium (frzb)', 'Endocardium (Ventricle)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)']
HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
#HR_Rv_7dpi_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)]
deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9'] # All samples with > 50 nppc fibroblasts.
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)]
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['heart'].isin(deep_injury_libraries)]
scv.pp.remove_duplicate_cells(HR_Rv_7dpi_deep_endo)


# In[ ]:


#HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['heart'].isin(['Hr2b'])].obs.Cell_type.value_counts()


# In[ ]:


#HR_Rv_7dpi_deep_endo.obs['Cell_type'].value_counts()


# In[ ]:


sc.pp.filter_genes(HR_Rv_7dpi_deep_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_7dpi_deep_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_deep_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_deep_endo)
sc.tl.pca(HR_Rv_7dpi_deep_endo)
#sc.pp.neighbors(HR_Rv_7dpi_deep_endo, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_7dpi_deep_endo, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_deep_endo)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi deep injury endocardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='batch',
          title = '7dpi deep injury endocardial niche')


# In[ ]:


#sc.tl.diffmap(HR_Rv_7dpi_deep_endo)
#sc.pp.neighbors(HR_Rv_7dpi_deep_endo, n_neighbors=30, use_rep='X_diffmap')
#sc.tl.draw_graph(HR_Rv_7dpi_deep_endo) # Do we need this?
sc.tl.leiden(HR_Rv_7dpi_deep_endo, resolution=2)
#sc.tl.umap(HR_Rv_7dpi_deep_endo)


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


#sc.pl.diffmap(HR_Rv_7dpi_deep_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


#sc.pl.diffmap(HR_Rv_7dpi_deep_endo, color='Cell_type', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


#sc.pl.umap(HR_Rv_7dpi_deep_endo, color='batch',
 #         title = '7dpi deep injury endocardial niche')


# In[ ]:


#sc.pl.umap(HR_Rv_7dpi_deep_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
#          title = '7dpi deep injury endocardial niche')


# In[ ]:


#sc.tl.diffmap(HR_Rv_7dpi_deep_endo)
#sc.pp.neighbors(HR_Rv_7dpi_deep_endo, n_neighbors=20, use_rep='X_diffmap')
#sc.tl.draw_graph(HR_Rv_7dpi_deep_endo)


# In[ ]:


#sc.tl.leiden(HR_Rv_7dpi_deep_endo, resolution=0.4)


# In[ ]:


#sc.pl.umap(HR_Rv_7dpi_deep_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


#sc.tl.umap(HR_Rv_7dpi_deep_endo)


# In[ ]:


#sc.pl.umap(HR_Rv_7dpi_deep_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


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


#sc.pl.umap(HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['Cell_type'].isin(['Endocardium (Ventricle)', 'Fibroblasts (nppc)'])], color = 'Cell_type')


# In[ ]:


#HR_Rv_7dpi_deep_endo_vent = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['Cell_type'].isin(['Endocardium (Ventricle)', 'Fibroblasts (nppc)'])]


# In[ ]:


#sc.pl.umap(HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['Cell_type'].isin(['Endocardium (Atrium)', 'Fibroblasts (nppc)'])], color = 'Cell_type')


# In[ ]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color = ['spock3', 'nppc', 'ptgs2a', 'ptgs2b', 'fabp11a', 'aldh1a2', 'epas1b', 'pcna',
                                         'notum1b', 'wif1'])


# In[ ]:


#sc.tl.draw_graph(HR_Rv_7dpi_deep_endo, init_pos='paga')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['Cell_type'].isin(['Endocardium (Atrium)', 'Fibroblasts (nppc)'])], color = ['Cell_type'], title = '7dpi deep endocardial niche')#, save = '_scvelo_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_7dpi_deep_endo, color = ['Cell_type'], title = '7dpi deep endocardial niche')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_7dpi_deep_endo, color = 'leiden', title = '7dpi deep endocardial niche', save = 'scvelo_deep_endocardial_niche_7dpi_20n_leiden_diffmap.png')


# In[ ]:


#sc.pl.draw_graph(HR_Rv_7dpi_deep_endo, ncols = 3,
#                 color = ['col1a1a', 'fn1b', 'spock3',
#                          'nppc', 'cyp26b1', 'frzb', 'acta2'],
#                save = '_scvelo_deep_endocardial_niche_7dpi_20n_marker_genes_diffmap.png')


# In[ ]:


scv.pp.moments(HR_Rv_7dpi_deep_endo)


# In[ ]:


#scv.tl.recover_dynamics(HR_Rv_7dpi_deep_endo) # For dynamical model, takes ~ 45min
#scv.tl.velocity(HR_Rv_7dpi_deep_endo, mode='dynamical')
scv.tl.velocity(HR_Rv_7dpi_deep_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_7dpi_deep_endo)


# In[ ]:


#HR_Rv_7dpi_deep_endo.write('./write/HR_Rv_7dpi_deep_endo.h5ad')
#HR_Rv_7dpi_deep_endo = sc.read('./write/HR_Rv_7dpi_deep_endo.h5ad')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, basis='umap', title = '', 
                                 color = 'Cell_type', legend_loc = 'none')


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_7dpi_deep_endo, basis='umap', title = '', scale = 0.3,
                                 color = 'Cell_type', legend_loc = 'none',
                         save = 'gene_velocity_endoventricular_deep_niche_7dpi_umap.png')


# In[ ]:


scv.pl.velocity(HR_Rv_7dpi_deep_endo, ['aldh1a2', 'krt18', 'kdrl', 'cdh5', 'fli1a', 'spock3', 'nppc'], basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
#                save = 'gene_velocity_epicardial_niche_3dpi_umap.png')


# In[ ]:


scv.tl.velocity_confidence(HR_Rv_7dpi_deep_endo)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_deep_endo, c=keys, cmap='coolwarm', perc=[5, 95])#,
#                                save = 'strength_coherence_scvelo_deep_endo_niche_7dpi_20n_cell_types_umap.png')


# In[ ]:


# Dynamical
#scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, basis='umap', title = '', 
#                                 color = 'Cell_type', legend_loc = 'none')


# In[ ]:


deep_endo_velofit = HR_Rv_7dpi_deep_endo.var
deep_endo_velofit = deep_endo_velofit[(deep_endo_velofit['fit_likelihood'] > .1) & deep_endo_velofit['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(deep_endo_velofit['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(deep_endo_velofit['fit_beta'] * deep_endo_velofit['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(deep_endo_velofit['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(HR_Rv_7dpi_deep_endo, 'fit*', dropna=True).head()


# A high likelihood means it's an interesting gene, i.e. values close to 1 are better than values close to 0.

# In[ ]:


scv.tl.latent_time(HR_Rv_7dpi_deep_endo)
scv.pl.scatter(HR_Rv_7dpi_deep_endo, color='latent_time', color_map='gnuplot', size=80)


# In[ ]:


sc.tl.rank_genes_groups(HR_Rv_7dpi_deep_endo, groupby='Cell_type')


# In[ ]:


#HR_Rv_7dpi_deep_endo.uns['rank_genes_groups']


# In[ ]:


HR_Rv_7dpi_deep_endo.uns['rank_genes_groups']['names']['Fibroblasts (nppc)']


# In[ ]:


HR_Rv_7dpi_deep_endo.uns['rank_genes_groups']['names']['Endocardium (Ventricle)']


# In[ ]:


scv.pl.velocity(HR_Rv_7dpi_deep_endo_vent, ['aldh1a2', 'krt18', 'kdrl', 'cdh5', 'fli1a'], basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,


# In[ ]:


#top_genes = HR_Rv_7dpi_deep_endo.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(HR_Rv_7dpi_deep_endo, basis=top_genes[:15], ncols=5, frameon=False, color = 'Cell_type')


# In[ ]:


var_names = ['frzb', 'nppc']
scv.pl.scatter(HR_Rv_7dpi_deep_endo, var_names, frameon=False, color = 'Cell_type')
scv.pl.scatter(HR_Rv_7dpi_deep_endo, x='latent_time', y=var_names, frameon=False, color = 'Cell_type')


# In[ ]:


scv.tl.rank_dynamical_genes(HR_Rv_7dpi_deep_endo, groupby='Cell_type')
HR_Rv_7dpi_deep_endo_dynamic_genes = scv.get_df(HR_Rv_7dpi_deep_endo, 'rank_dynamical_genes/names')


# In[ ]:


HR_Rv_7dpi_deep_endo_dynamic_genes


# In[ ]:


#restriction_celltypes = ['Endocardium (Atrium)', 'Endocardium (Ventricle)', 'Endocardium (frzb)', 
#                         'Fibroblasts (nppc)', 'Fibroblasts (spock3)']
#HR_Rv_7dpi_deep_endo_rs = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['Cell_type'].isin(restriction_celltypes)]
scv.pl.scatter(HR_Rv_7dpi_deep_endo_rs, ['brdt', 'vmp1', 'zfand2a', 'ubb', 'sqstm1'], frameon=False, color = 'Cell_type')
scv.pl.scatter(HR_Rv_7dpi_deep_endo_rs, x='latent_time', y=['brdt', 'vmp1', 'zfand2a', 'ubb', 'sqstm1'], frameon=False, color = 'Cell_type')
#HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)


# In[ ]:


s1 = HR_Rv_7dpi_deep_endo_dynamic_genes['Fibroblasts (nppc)']
s2 = HR_Rv_7dpi_deep_endo_dynamic_genes['Endocardium (frzb)']
frzb_nppc = list(set(s1) & set(s2))
#HR_Rv_7dpi_deep_endo_dynamic_genes['Fibroblasts (nppc)'].intersection(HR_Rv_7dpi_deep_endo_dynamic_genes['Endocardium (frzb)'])


# In[ ]:


HR_Rv_7dpi_deep_endo.var.loc[frzb_nppc]#.['fit_likelihood']#[top_genes]


# In[ ]:


#import inspect
#print(inspect.getsource(scv.pl.proportions))


# In[ ]:


#HR_Rv_7dpi_deep_endo


# In[ ]:


#HR_Rv_7dpi_deep_endo.layers['unspliced']


# In[ ]:


#s_genes, g2m_genes = scv.utils.get_phase_marker_genes(HR_Rv_7dpi_deep_endo)


# In[ ]:


#g2m_genes


# In[ ]:


#scv.tl.score_genes_cell_cycle(HR_Rv_7dpi_deep_endo)
#scv.pl.scatter(HR_Rv_7dpi_deep_endo, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])


# In[ ]:


#scv.pl.velocity_embedding(HR_Rv_7dpi_deep_endo, basis='draw_graph_fa')


# In[ ]:


#scv.pl.velocity_embedding_grid(HR_Rv_7dpi_deep_endo, basis='draw_graph_fa', color = 'Cell_type')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, basis='draw_graph_fa', title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
                                #save = '_scvelo_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


endo_velo_genes = ['col1a1a', 'nppc', 'fli1a', 'aldh1a2']
for gene in endo_velo_genes:
    scv.pl.velocity(HR_Rv_7dpi_deep_endo, gene, basis='draw_graph_fa', ncols=1, fontsize=16)


# How are the speed and coherence?

# In[ ]:


scv.tl.velocity_confidence(HR_Rv_7dpi_deep_endo)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_deep_endo, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_deep_endo_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, title = '', 
                                 color = 'Cell_type', legend_loc = 'none')#,
                                 #save = '_scvelo_endocardial_niche_7dpi_deep_20n_cell_types_umap.png')


# In[ ]:


scv.tl.paga(HR_Rv_7dpi_deep_endo, groups='leiden')
df = scv.get_df(HR_Rv_7dpi_deep_endo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[ ]:


scv.pl.paga(HR_Rv_7dpi_deep_endo, color = 'Cell_type', basis = 'umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, legend_loc = 'none')


# In[ ]:


#scv.tl.rank_velocity_genes(HR_Rv_7dpi_deep_endo, groupby='leiden', min_corr=.3)

#dfvg = scv.DataFrame(HR_Rv_7dpi_deep_endo.uns['rank_velocity_genes']['names'])
dfvg['9']


# In[ ]:


endo_velo_genes = ['gfm1', 'il34', 'fhl3a'] #['spock3', 'nppc', 'ptgs2a', 'ptgs2b', 'fabp11a', 'aldh1a2', 'epas1b']
for gene in endo_velo_genes:
    scv.pl.velocity(HR_Rv_7dpi_deep_endo, gene, basis='umap', ncols=1, fontsize=16)#,
                   # save = gene+'_velocity_epicardial_niche_3dpi_20n_umap.png')


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_deep_endo, c=keys, cmap='coolwarm', perc=[5, 95])#,
#                                save = 'strength_coherence_scvelo_deep_endo_niche_7dpi_20n_cell_types_umap.png')


# In[ ]:


#HR_Rv_7dpi_deep_endo.obs


# Velocity length for nppc fibroblasts in Hr9:

# In[ ]:


#HR_Rv_7dpi_deep_endo.obs[(HR_Rv_7dpi_deep_endo.obs['Cell_type'] == "Fibroblast (nppc)") & (HR_Rv_7dpi_deep_endo.obs['heart'] == 'Hr9')]


# In[ ]:


#HR_Rv_Hr9_endo.obs[(HR_Rv_Hr9_endo.obs['Cell_type'] == "Fibroblast (nppc)") & (HR_Rv_Hr9_endo.obs['heart'] == 'Hr9')]


# # Endocardial niche for Hr9 (a deep injury library at 7dpi)

# In[ ]:


HR_Rv_Hr9_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['heart'].isin(['Hr9'])]
HR_Rv_Hr9_endo = HR_Rv_Hr9_endo[HR_Rv_Hr9_endo.obs['Cell_type'].isin(endo_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_Hr9_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_Hr9_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_Hr9_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_Hr9_endo)
sc.tl.pca(HR_Rv_Hr9_endo)
sc.pp.neighbors(HR_Rv_Hr9_endo, n_neighbors=30)
sc.tl.umap(HR_Rv_Hr9_endo)


# In[ ]:


sc.pl.umap(HR_Rv_Hr9_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_Hr9_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Hr9: 7dpi deep injury endocardial niche')


# In[ ]:


sc.tl.leiden(HR_Rv_Hr9_endo, resolution=1)


# In[ ]:


sc.pl.umap(HR_Rv_Hr9_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_Hr9_endo, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_Hr9_endo, show=True, color = 'Cell_type')


# In[ ]:


sc.tl.umap(HR_Rv_Hr9_endo, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_Hr9_endo, color = 'Cell_type')


# In[ ]:


scv.pp.moments(HR_Rv_Hr9_endo)


# In[ ]:


scv.tl.velocity(HR_Rv_Hr9_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_Hr9_endo)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_Hr9_endo, basis='umap', #basis='draw_graph_fa', title = 'Hr9: 7dpi deep endocardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin')#,
                                #save = '_scvelo_Hr9_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_Hr9_endo, basis='umap',
                                 color = 'leiden', legend_loc = 'right margin')#,


# In[ ]:


sc.pl.paga(HR_Rv_Hr9_endo, show=True, color = 'Cell_type', node_size_scale = 2)


# In[ ]:


HR_Rv_Hr9_endo.uns['neighbors']['distances'] = HR_Rv_Hr9_endo.obsp['distances']
HR_Rv_Hr9_endo.uns['neighbors']['connectivities'] = HR_Rv_Hr9_endo.obsp['connectivities']

#scv.tl.paga(adata, groups='clusters')
#df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
#df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.tl.paga(HR_Rv_Hr9_endo, groups='leiden')
df = scv.get_df(HR_Rv_Hr9_endo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[ ]:


scv.pl.paga(HR_Rv_Hr9_endo, color = 'Cell_type', basis = 'umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, legend_loc = 'none')


# # Endocardial niche for Hr2a (a deep injury library at 7dpi)

# In[ ]:


HR_Rv_Hr2a_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['heart'].isin(['Hr2a'])]
HR_Rv_Hr2a_endo = HR_Rv_Hr2a_endo[HR_Rv_Hr2a_endo.obs['Cell_type'].isin(endo_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_Hr2a_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_Hr2a_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_Hr2a_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_Hr2a_endo)
sc.tl.pca(HR_Rv_Hr2a_endo)
sc.pp.neighbors(HR_Rv_Hr2a_endo, n_neighbors=30)
sc.tl.umap(HR_Rv_Hr2a_endo)


# In[ ]:


sc.pl.umap(HR_Rv_Hr2a_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_Hr2a_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Hr2a: 7dpi deep injury endocardial niche')


# In[ ]:


sc.tl.leiden(HR_Rv_Hr2a_endo, resolution=1)


# In[ ]:


sc.pl.umap(HR_Rv_Hr2a_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_Hr2a_endo, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_Hr2a_endo, show=True, color = 'Cell_type')


# In[ ]:


sc.tl.umap(HR_Rv_Hr2a_endo, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_Hr2a_endo, color = 'Cell_type')


# In[ ]:


scv.pp.moments(HR_Rv_Hr2a_endo)


# In[ ]:


scv.tl.velocity(HR_Rv_Hr2a_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_Hr2a_endo)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_Hr2a_endo, basis='umap', #basis='draw_graph_fa', title = 'Hr2a: 7dpi deep endocardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin')#,
                                #save = '_scvelo_Hr2a_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


sc.pl.paga(HR_Rv_Hr2a_endo, show=True, color = 'Cell_type', node_size_scale = 2)


# In[ ]:


HR_Rv_Hr2a_endo.uns['neighbors']['distances'] = HR_Rv_Hr2a_endo.obsp['distances']
HR_Rv_Hr2a_endo.uns['neighbors']['connectivities'] = HR_Rv_Hr2a_endo.obsp['connectivities']

#scv.tl.paga(adata, groups='clusters')
#df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
#df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.tl.paga(HR_Rv_Hr2a_endo, groups='leiden')
df = scv.get_df(HR_Rv_Hr2a_endo, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[ ]:


scv.pl.paga(HR_Rv_Hr2a_endo, color = 'Cell_type', basis = 'umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, legend_loc = 'none')


# # Endocardial niche for Hr2b (a deep injury library at 7dpi)

# In[ ]:


HR_Rv_Hr2b_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['heart'].isin(['Hr2b'])]
HR_Rv_Hr2b_endo = HR_Rv_Hr2b_endo[HR_Rv_Hr2b_endo.obs['Cell_type'].isin(endo_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_Hr2b_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_Hr2b_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_Hr2b_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_Hr2b_endo)
sc.tl.pca(HR_Rv_Hr2b_endo)
sc.pp.neighbors(HR_Rv_Hr2b_endo, n_neighbors=30)
sc.tl.umap(HR_Rv_Hr2b_endo)


# In[ ]:


sc.pl.umap(HR_Rv_Hr2b_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_Hr2b_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Hr2b: 7dpi deep injury endocardial niche')


# In[ ]:


sc.tl.leiden(HR_Rv_Hr2b_endo, resolution=1)


# In[ ]:


sc.pl.umap(HR_Rv_Hr2b_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_Hr2b_endo, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_Hr2b_endo, show=True, color = 'Cell_type')


# In[ ]:


sc.tl.umap(HR_Rv_Hr2b_endo, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_Hr2b_endo, color = 'Cell_type')


# In[ ]:


scv.pp.moments(HR_Rv_Hr2b_endo)


# In[ ]:


scv.tl.velocity(HR_Rv_Hr2b_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_Hr2b_endo)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_Hr2b_endo, basis='umap', 
                                 color = 'Cell_type', legend_loc = 'right margin')


# # Endocardial niche for Hr6v (a deep injury library at 7dpi)

# In[ ]:


HR_Rv_Hr6v_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['heart'].isin(['Hr6v'])]
HR_Rv_Hr6v_endo = HR_Rv_Hr6v_endo[HR_Rv_Hr6v_endo.obs['Cell_type'].isin(endo_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_Hr6v_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_Hr6v_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_Hr6v_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_Hr6v_endo)
sc.tl.pca(HR_Rv_Hr6v_endo)
sc.pp.neighbors(HR_Rv_Hr6v_endo, n_neighbors=30)
sc.tl.umap(HR_Rv_Hr6v_endo)


# In[ ]:


sc.pl.umap(HR_Rv_Hr6v_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_Hr6v_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Hr6v: 7dpi deep injury endocardial niche')


# In[ ]:


sc.tl.leiden(HR_Rv_Hr6v_endo, resolution=1)


# In[ ]:


sc.pl.umap(HR_Rv_Hr6v_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_Hr6v_endo, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_Hr6v_endo, show=True, color = 'Cell_type')


# In[ ]:


sc.tl.umap(HR_Rv_Hr6v_endo, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_Hr6v_endo, color = 'Cell_type')


# In[ ]:


scv.pp.moments(HR_Rv_Hr6v_endo)


# In[ ]:


scv.tl.velocity(HR_Rv_Hr6v_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_Hr6v_endo)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_Hr6v_endo, basis='umap', 
                                 color = 'Cell_type', legend_loc = 'right margin')


# # Endocardial niche for Hr1 (a deep injury library at 7dpi)

# In[ ]:


HR_Rv_Hr1_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['heart'].isin(['Hr1'])]
HR_Rv_Hr1_endo = HR_Rv_Hr1_endo[HR_Rv_Hr1_endo.obs['Cell_type'].isin(endo_7dpi)]


# In[ ]:


sc.pp.filter_genes(HR_Rv_Hr1_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_Hr1_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_Hr1_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_Hr1_endo)
sc.tl.pca(HR_Rv_Hr1_endo)
sc.pp.neighbors(HR_Rv_Hr1_endo, n_neighbors=30)
sc.tl.umap(HR_Rv_Hr1_endo)


# In[ ]:


sc.pl.umap(HR_Rv_Hr1_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_Hr1_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Hr1: 7dpi deep injury endocardial niche')


# In[ ]:


sc.tl.leiden(HR_Rv_Hr1_endo, resolution=1)


# In[ ]:


sc.pl.umap(HR_Rv_Hr1_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_Hr1_endo, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_Hr1_endo, show=True, color = 'Cell_type')


# In[ ]:


sc.tl.umap(HR_Rv_Hr1_endo, init_pos='paga')


# In[ ]:


sc.pl.umap(HR_Rv_Hr1_endo, color = 'Cell_type')


# In[ ]:


scv.pp.moments(HR_Rv_Hr1_endo)


# In[ ]:


scv.tl.velocity(HR_Rv_Hr1_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_Hr1_endo)


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_Hr1_endo, basis='umap', 
                                 color = 'Cell_type', legend_loc = 'right margin')


# # Endocardial niche at 3dpi

# In[ ]:


HR_Rv_3dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '3']
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
HR_Rv_3dpi_endo = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(['Endocardium (V)', 'Fibroblast (nppc)', 'Fibroblast (spock3)', 'Fibroblast-like cells'])]


# In[ ]:


sc.pp.filter_genes(HR_Rv_3dpi_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_3dpi_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_3dpi_endo)


# In[ ]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_endo)
sc.tl.pca(HR_Rv_3dpi_endo)
sc.pp.neighbors(HR_Rv_3dpi_endo, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_3dpi_endo, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_endo)


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_endo, color='batch',
          title = '3dpi endocardial niche')


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi endocardial niche')


# In[ ]:


HR_Rv_3dpi_endo.obs['Cell_type'].value_counts()


# In[ ]:


sc.tl.diffmap(HR_Rv_3dpi_endo)
sc.pp.neighbors(HR_Rv_3dpi_endo, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_Rv_3dpi_endo)


# In[ ]:


sc.tl.leiden(HR_Rv_3dpi_endo, resolution=0.4)


# In[ ]:


sc.pl.umap(HR_Rv_3dpi_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[ ]:


sc.tl.paga(HR_Rv_7dpi_endo, groups='leiden')


# In[ ]:


sc.pl.paga(HR_Rv_7dpi_endo, show=True, 
           labels = ['', '', '', '', '', 
                     '', '', '', '', '', 
                     '', '', '', '', '',
                    '',''], 
           node_size_scale = 2, save = 'scvelo_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[ ]:


sc.tl.draw_graph(HR_Rv_7dpi_endo, init_pos='paga')


# In[ ]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, color = 'leiden', title = '7dpi endocardial niche', save = 'scvelo_endocardial_niche_7dpi_20n_leiden_diffmap.png')


# In[ ]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, color = ['Cell_type'], title = '7dpi endocardial niche', save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, ncols = 3,
                 color = ['col1a1a', 'fn1b', 'spock3',
                          'nppc', 'cyp26b1', 'frzb', 'acta2'],
                save = '_scvelo_endocardial_niche_7dpi_20n_marker_genes_diffmap.png')


# In[ ]:


scv.pp.moments(HR_Rv_7dpi_endo)


# In[ ]:


scv.tl.velocity(HR_Rv_7dpi_endo, mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(HR_Rv_7dpi_endo)


# In[ ]:


scv.pl.velocity_embedding(HR_Rv_7dpi_endo, basis='draw_graph_fa')


# In[ ]:


scv.pl.velocity_embedding_grid(HR_Rv_7dpi_endo, basis='draw_graph_fa', color = 'Cell_type')


# In[ ]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endo, basis='draw_graph_fa', title = '7dpi endocardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[ ]:


scv.pl.proportions(HR_Rv_7dpi_endo, groupby = 'batch')


# How are the speed and coherence?

# In[ ]:


scv.tl.velocity_confidence(HR_Rv_7dpi_endo)


# In[ ]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_endo, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_endo_niche_7dpi_20n_cell_types_diffmap.png')


# # Full niche

# In[ ]:


HR_niche = HR_filter[HR_filter.obs['Cell_type'].isin(trajectory_subset)]


# In[ ]:


sc.pp.filter_genes(HR_niche, min_cells=3)
sc.pp.normalize_per_cell(HR_niche, counts_per_cell_after=1e4,copy=False)
sc.pp.log1p(HR_niche, copy=False)
sc.pp.highly_variable_genes(HR_niche)
sc.tl.pca(HR_niche)
sc.pp.neighbors(HR_niche, n_neighbors=30)
sc.external.pp.bbknn(HR_niche, batch_key='batch')
sc.tl.umap(HR_niche)


# In[ ]:


sc.pl.umap(HR_niche, color='batch',
          title = 'Niche')#, save = 'subset1_3dpi_umap_batch.png')


# In[ ]:


sc.pl.umap(HR_niche, color='dpi',
          title = 'Niche')


# In[ ]:


HR_niche.obs


# In[ ]:


sc.pl.umap(HR_niche, color='Cell_type', palette = cell_type_colors.loc[HR_niche.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[ ]:


HR_niche.obs.Cell_type.cat.categories.tolist()

