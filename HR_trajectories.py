
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
scv.settings.set_figure_params('scvelo')


# In[2]:


annotations = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv', index_col = 0)


# In[3]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors_2.csv', index_col = 0)


# In[4]:


trajectory_subset = ['Fibroblast', 'Fibroblast (cfd)', 'Fibroblast (col11a1a)', 'Fibroblast (col12a1a)', 
                     'Fibroblast (cxcl12a)', 'Fibroblast (mpeg1.1)', 'Fibroblast (nppc)', 'Fibroblast (spock3)',
                    'Fibroblast-like cells', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblast (proliferating)', 'Perivascular cells']


# In[5]:


connected_3dpi = ['Fibroblast', 'Fibroblast (cfd)', 'Fibroblast (col11a1a)', 'Fibroblast (col12a1a)', 
                     'Fibroblast (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblast (proliferating)', 'Perivascular cells']


# In[6]:


connected_7dpi = ['Fibroblast', 'Fibroblast (cfd)', 'Fibroblast (col11a1a)', 'Fibroblast (col12a1a)', 
                     'Fibroblast (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblast (proliferating)', 'Perivascular cells', 'Fibroblast (cxcl12a)']


# In[7]:


#fibro_colors = cell_type_colors[['color', 'setFibro']]
#fibro_colors = fibro_colors[fibro_colors.notna().setFibro]
#fibro_colors = fibro_colors.set_index('setFibro')


# In[8]:


HR_setnames = pd.DataFrame({'batch': ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                                      '11', '12','13', '14', '15', '16', '17', '18', '19', '20',
                                     '21', '22','23', '24', '25', '26', '27', '28', '29', '30',
                                     '31', '32','33', '34', '35', '36', '37', '38', '39', '40',
                                     '41', '42'],
                            'heart': ['H5', 'H6', 'H7', 'H8a', 'H8v', 'Hr1', 'Hr2a', 'Hr2b', 'Hr3', 
                                      'Hr4', 'Hr5', 'Hr6a', 'Hr6v', 'Hr7a', 'Hr7v', 'Hr8', 'Hr9', 'Hr10',
                                     'Hr11', 'Hr12', 'Hr13', 'Hr14', 'Hr15',
                                     'Hr16', 'Hr17', 'Hr18', 'Hr19', 'Hr20',
                                     'Hr21', 'Hr22', 'Hr23', 'Hr24', 'Hr25',
                                     'Hr26', 'Hr27', 'Hr28', 'Hr29', 'Hr30',
                                     'Hr31', 'Hr32', 'Hr33', 'Hr34', 'Hr35'],
                             'dpi': ['0', '0', '0', '0', '0', '7', '7', '7', '30', 
                                    '30', '60', '7', '7', '7', '7', '7', '7', '3',
                                   '3', '3', '7', '7', '7',
                                   '15', '15', '15', '30', '30',
                                   '30', '3', '3', '3', '3',
                                   '3', '3', '3', '3', '7',
                                   '7', '7', '7', '3', '3'],
                           'inhib': ['no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 
                                    'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no',
                                   'no', 'no', 'no', 'no', 'no',
                                   'no', 'no', 'no', 'no', 'no',
                                   'no', 'no', 'no', 'no', 'no',
                                   'no', 'no', 'DMSO', 'IWR1', 'DMSO',
                                   'IWR1', 'IWR1', 'IWR1', 'IWR1', 'IWR1']})


# # Load and annotate single-cell data

# In[25]:


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


# In[26]:


HR =    H5_data.concatenate(H6_data, H7_data, H8a_data, H8v_data,
                        Hr1_data, Hr2a_data, Hr2b_data, Hr3_data, Hr4_data, Hr5_data, 
                        Hr6a_data, Hr6v_data, Hr7a_data, Hr7v_data, Hr8_data, Hr9_data, Hr10_data,
                        Hr11_data, Hr12_data, Hr13_data, Hr14_data, Hr15_data,
                        Hr16_data, Hr17_data, Hr18_data, Hr19_data, Hr20_data,
                        Hr21_data, Hr22_data, Hr23_data, Hr24_data, Hr25_data,
                        Hr26_data, Hr27_data, Hr28_data, Hr29_data, Hr30_data,
                        Hr31_data, Hr32_data, Hr33_data, Hr34_data, Hr35_data)
HR.shape


# In[27]:


HR_obs = HR.obs
HR.obs = HR.obs.reset_index().merge(HR_setnames, how="inner").set_index('index')


# In[28]:


# Rename cells so the cell names correspond to the ones in the annotation file
HR.obs_names = [str(HR.obs.loc[x,'heart'])+'_'+str(x.split('-', 1)[0]) for x in HR.obs_names]


# In[30]:


# Drop annotations that are not in the single-cell object
anno_drop = annotations.index.difference(HR.obs_names)
annotations = annotations.drop(anno_drop)


# In[31]:


HR_filter = HR[annotations.index]
HR_filter


# In[32]:


HR_filter.obs['Cell_type'] = annotations['Cell_type'].tolist()


# In[4]:


#HR_filter.write('./write/HR_filter.h5ad')


# # Trajectory analysis on the 3dpi niche

# In[ ]:


HR_filter = sc.read('./write/HR_filter.h5ad')


# In[34]:


HR_ps_1 = HR_filter[HR_filter.obs['dpi'] == '3']
all_genes_but_RFP = [name for name in HR_ps_1.var_names if not name == 'RFP']
HR_ps_1 = HR_ps_1[:, all_genes_but_RFP]


# In[35]:


sc.pp.filter_genes(HR_ps_1, min_cells=3)
HR_ps_1_norm = sc.pp.normalize_per_cell(HR_ps_1, counts_per_cell_after=1e4,copy=True)
HR_ps_1 = sc.pp.log1p(HR_ps_1_norm, copy=True)


# In[21]:


HR_3d_sub1_cbr = HR_ps_1[HR_ps_1.obs['Cell_type'].isin(trajectory_subset)]
#HR_3d_sub1_cbr = sc.read('./write/HR_3dpi_subset1.h5ad')


# In[22]:


sc.pp.highly_variable_genes(HR_3d_sub1_cbr)
sc.tl.pca(HR_3d_sub1_cbr)
sc.pp.neighbors(HR_3d_sub1_cbr, n_neighbors=30)
sc.external.pp.bbknn(HR_3d_sub1_cbr, batch_key='batch')
sc.tl.umap(HR_3d_sub1_cbr)
# The order here is very important, to my surprise. If I run bbknn before neighbors, I do not get batch correction at all.


# In[28]:


sc.pl.umap(HR_3d_sub1_cbr, color='batch',
          title = 'Subset 1, 3dpi', save = 'subset1_3dpi_umap_batch.png')


# In[112]:


sc.pl.umap(HR_3d_sub1_cbr, color='Cell_type', palette = fibro_colors.loc[HR_3d_sub1_cbr.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Subset 1, 3dpi')#, save = 'subset1_3dpi_umap.png')


# In[31]:


sc.tl.diffmap(HR_3d_sub1_cbr)


# In[32]:


sc.pp.neighbors(HR_3d_sub1_cbr, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_3d_sub1_cbr)


# In[33]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color='Cell_type', title = 'Subset 1, 3dpi', save = 'subset1_3dpi_20n_diffmap.png')


# In[34]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color = 'batch', title = 'Subset 1, 3dpi', save = 'subset1_3dpi_20n_diffmap_batches.png')


# In[35]:


sc.tl.paga(HR_3d_sub1_cbr, groups='Cell_type')


# In[36]:


sc.pl.paga(HR_3d_sub1_cbr, threshold=0.3, show=True, labels = ['', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'subset1_3dpi_20n_diffmap.png')


# In[37]:


sc.tl.draw_graph(HR_3d_sub1_cbr, init_pos='paga')


# In[38]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color = 'Cell_type', title = 'Subset 1, 3dpi', save = 'subset1_3dpi_20n_diffmap_pagainit.png')


# We can see in the PAGA that everything is connected, even cell types such as fibroblast-like cells and spock3 and nppc-fibroblasts that we know are not connected.  
# Maybe this has to do with the coarse grained clustering?

# In[97]:


sc.tl.leiden(HR_3d_sub1_cbr, resolution=0.5)


# In[98]:


sc.pl.umap(HR_3d_sub1_cbr, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[102]:


sc.pl.umap(HR_3d_sub1_cbr, color='Cell_type')


# In[99]:


sc.tl.paga(HR_3d_sub1_cbr, groups='leiden')


# In[107]:


sc.tl.draw_graph(HR_3d_sub1_cbr, init_pos='paga')


# In[106]:


sc.pl.paga(HR_3d_sub1_cbr, threshold = 0.2, color=['leiden'])


# In[110]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color = ['Cell_type'])


# In[111]:


sc.pl.draw_graph(HR_3d_sub1_cbr, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# After some fiddling with parameters, we manage to disconnect the fibroblast-like cells although they are still connected to another cluster of fibroblasts. However, the spock3 and nppc-fibroblasts remain connected, the perivascular cells are also disconnected, and a part of the atrial epicardium is now also disconnected.

# # RNA velocity

# In[128]:


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


# In[129]:


HR_Rv =    H5_Rv_data.concatenate(H6_Rv_data, H7_Rv_data, H8a_Rv_data, H8v_Rv_data,
                        Hr1_Rv_data, Hr2a_Rv_data, Hr2b_Rv_data, Hr3_Rv_data, Hr4_Rv_data, Hr5_Rv_data, 
                        Hr6a_Rv_data, Hr6v_Rv_data, Hr7a_Rv_data, Hr7v_Rv_data, Hr8_Rv_data, Hr9_Rv_data, Hr10_Rv_data,
                        Hr11_Rv_data, Hr12_Rv_data, Hr13_Rv_data, Hr14_Rv_data, Hr15_Rv_data,
                        Hr16_Rv_data, Hr17_Rv_data, Hr18_Rv_data, Hr19_Rv_data, Hr20_Rv_data,
                        Hr21_Rv_data, Hr22_Rv_data, Hr23_Rv_data, Hr24_Rv_data, Hr25_Rv_data,
                        Hr26_Rv_data, Hr27_Rv_data, Hr28_Rv_data, Hr29_Rv_data, Hr30_Rv_data,
                        Hr31_Rv_data, Hr32_Rv_data, Hr33_Rv_data, Hr34_Rv_data, Hr35_Rv_data)
HR_Rv.shape


# In[130]:


#HR_Rv.write('./write/HR_Rv.h5ad')


# In[54]:


HR_Rv = sc.read('./write/HR_Rv.h5ad')


# In[55]:


HR_Rv.obs = HR_Rv.obs.reset_index().merge(HR_setnames, how="inner").set_index('index')


# In[56]:


# Rename cells to match cell names in annotation file
HR_Rv.obs_names = [str(HR_Rv.obs.loc[x,'heart'])+'_'+str(x.split(':', 1)[1])[0:16] for x in HR_Rv.obs_names]
# Drop annotations that are not in the single-cell object
anno_drop_Rv = annotations.index.difference(HR_Rv.obs_names)
annotations_Rv = annotations.drop(anno_drop_Rv)


# In[57]:


HR_Rv_filter = HR_Rv[annotations_Rv.index]
HR_Rv_filter


# In[58]:


HR_Rv_filter.obs['Cell_type'] = annotations_Rv['Cell_type'].tolist()


# In[63]:


HR_Rv_filter.obs['Cell_type'].value_counts()
#HR_Rv_7dpi_deep_endo.obs['Cell_type'].unique()


# In[64]:


#HR_Rv_filter.write('./write/HR_Rv_filter.h5ad')
HR_Rv_filter = sc.read('./write/HR_Rv_filter.h5ad')


# # Trajectories in 3dpi epicardial connected niche

# In[13]:


HR_Rv_3dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '3']
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
HR_Rv_3dpi_conn = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(connected_3dpi)]


# In[14]:


sc.pp.filter_genes(HR_Rv_3dpi_conn, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_3dpi_conn, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_3dpi_conn)


# In[15]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_conn)
sc.tl.pca(HR_Rv_3dpi_conn)
sc.pp.neighbors(HR_Rv_3dpi_conn, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_3dpi_conn, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_conn)


# In[16]:


sc.pl.umap(HR_Rv_3dpi_conn, color='batch',
          title = '3dpi connected niche')


# In[17]:


sc.pl.umap(HR_Rv_3dpi_conn, color='Cell_type', palette = fibro_colors.loc[HR_Rv_3dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi connected niche')


# In[18]:


sc.tl.diffmap(HR_Rv_3dpi_conn)
sc.pp.neighbors(HR_Rv_3dpi_conn, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_Rv_3dpi_conn)


# In[19]:


sc.tl.leiden(HR_Rv_3dpi_conn, resolution=0.4)


# In[20]:


sc.pl.umap(HR_Rv_3dpi_conn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[21]:


sc.tl.paga(HR_Rv_3dpi_conn, groups='leiden')


# In[22]:


sc.pl.paga(HR_Rv_3dpi_conn, show=True, labels = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'scvelo_connected_niche_3dpi_20n_leiden_diffmap.png')


# In[23]:


sc.tl.draw_graph(HR_Rv_3dpi_conn, init_pos='paga')


# In[24]:


sc.pl.draw_graph(HR_Rv_3dpi_conn, color = 'leiden', title = '3dpi connected niche', save = 'scvelo_connected_niche_3dpi_20n_leiden_diffmap.png')


# In[25]:


sc.pl.draw_graph(HR_Rv_3dpi_conn, color = ['Cell_type'], title = '3dpi connected niche', save = '_scvelo_connected_niche_3dpi_20n_cell_types_diffmap.png')


# In[26]:


sc.pl.draw_graph(HR_Rv_3dpi_conn, ncols = 4,
                 color = ['col1a1a', 'tbx18',
                          'notch3', 'pdgfrb', 'mpeg1.1', 'cfd', 
                          'col11a1a', 'col12a1a', 
                          'postnb', 'pcna', 'aldh1a2', 'stra6'],
                save = '_scvelo_connected_niche_3dpi_20n_marker_genes_diffmap.png')


# Slightly different but very similar. I like the scanpy colors better for the heatmap.

# # RNA velocity on 3dpi connected

# In[43]:


#scv.pp.filter_and_normalize(HR_Rv_3dpi_conn)
scv.pp.moments(HR_Rv_3dpi_conn)


# In[44]:


scv.tl.velocity(HR_Rv_3dpi_conn, mode='stochastic')


# In[45]:


scv.tl.velocity_graph(HR_Rv_3dpi_conn)


# In[46]:


scv.pl.velocity_embedding(HR_Rv_3dpi_conn, basis='draw_graph_fa')


# In[47]:


scv.pl.velocity_embedding_grid(HR_Rv_3dpi_conn, basis='draw_graph_fa', color = 'Cell_type')


# In[62]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, basis='draw_graph_fa', title = '3dpi connected niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_connected_niche_3dpi_20n_cell_types_diffmap.png')


# Rename leiden clusters

# In[52]:


sc.pl.umap(HR_Rv_3dpi_conn, color='Cell_type')


# In[50]:


sc.pl.umap(HR_Rv_3dpi_conn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# 0: FB (col11a1a) 1
# 1: EpiC A 1
# 2: FB 1
# 3: EpiC V
# 4: FB (col12a1a)
# 5: EpiC A 2
# 6: FB 2
# 7: FB 3
# 8: FB 4
# 9: FB (col11a1a) 2
# 10: FB 5
# 11: PC
# 12: FB (mpeg1.1)
# 13: FB (col11a1a) 3
# 14: FB 6  
# 
# Used to be:
# 0: EpiC V
# 1: FB 1
# 2: FB (col11a1a) 1
# 3: FB 2
# 4: FB (col11a1a) 2
# 5: FB 3
# 6: FB 4
# 7: FB 5
# 8: FB (col11a1a) 3
# 9: EpiC A 1
# 10: Epic A 2
# 11: FB (col12a1a) 1
# 12: FB 6
# 13: PC
# 14: FB (mpeg1.1)
# 15: FB (col12a1a) 2
# 16: FB 7

# In[13]:


HR_Rv_3dpi_conn.obs['leiden'].cat.categories


# In[53]:


HR_Rv_3dpi_conn.obs['leiden_anno'] = HR_Rv_3dpi_conn.obs['leiden']


# In[54]:


HR_Rv_3dpi_conn.obs['leiden_anno'].cat.categories = ['FB (col11a1a) 1', 'EpiC A 1', 'FB 1', 'EpiC V', 
                                                     'FB (col12a1a)', 'EpiC A 2', 'FB 2', 'FB 3', 
                                                     'FB 4', 'FB (col11a1a) 2', 'FB 5', 'PC', 
                                                     'FB (mpeg1.1)', 'FB (col11a1a)', 'FB 6']


# In[28]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_conn, basis='draw_graph_fa', color = 'leiden_anno', legend_loc = 'on data', title = '3dpi connected niche')


# Intron/exon proportions per cluster - this is very strange; why would it be 9% across the board? It's also 9% for each sample. I hope this is just the plotting calculation going wrong and not an outcome of, say, the batch integration.

# In[29]:


scv.pl.proportions(HR_Rv_3dpi_conn, groupby = 'leiden_anno')


# How are the speed and coherence?

# In[63]:


scv.tl.velocity_confidence(HR_Rv_3dpi_conn)


# In[64]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_3dpi_conn, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_connected_niche_3dpi_20n_cell_types_diffmap.png')


# In[58]:


scv.pl.velocity_graph(HR_Rv_3dpi_conn, threshold=.2, color = 'leiden_anno', basis='draw_graph_fa')
# This doesn't really show me a whole lot; maybe integrating this into PAGA works better?


# Which genes drive the transitions?

# In[68]:


scv.pl.velocity(HR_Rv_3dpi_conn, ['col1a1a', 'postnb', 'col11a1a', 'col12a1a', 
                          'aldh1a2', 'stra6', 'frzb', 'dkk3b'], basis='draw_graph_fa', ncols=2,
                                save = 'marker_genes_scvelo_connected_niche_3dpi_20n_cell_types_diffmap.png')


# In[27]:


#HR_Rv_3dpi_conn.write('./write/HR_Rv_3dpi_conn.h5ad')
HR_Rv_3dpi_conn = sc.read('./write/HR_Rv_3dpi_conn.h5ad')


# In[24]:


scv.tl.rank_velocity_genes(HR_Rv_3dpi_conn, groupby='leiden_anno', min_corr=.3)


# In[27]:


HR_3dpi_conn_velogenes = scv.DataFrame(HR_Rv_3dpi_conn.uns['rank_velocity_genes']['names'])
HR_3dpi_conn_velogenes.head()


# In[33]:


scv.pl.scatter(HR_Rv_3dpi_conn, HR_3dpi_conn_velogenes['FB 1'][:5], ylabel='FB 1', color = 'leiden_anno', add_outline='FB 1, FB 3, FB 5, EpiC V')


# # Same niche at 7dpi

# In[31]:


HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
HR_Rv_7dpi_conn = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_7dpi)]


# In[32]:


sc.pp.filter_genes(HR_Rv_7dpi_conn, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_7dpi_conn, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_conn)


# In[33]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_conn)
sc.tl.pca(HR_Rv_7dpi_conn)
sc.pp.neighbors(HR_Rv_7dpi_conn, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_7dpi_conn, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_conn)


# In[34]:


sc.pl.umap(HR_Rv_7dpi_conn, color='batch',
          title = '7dpi connected niche')


# In[35]:


sc.pl.umap(HR_Rv_7dpi_conn, color='Cell_type', palette = fibro_colors.loc[HR_Rv_7dpi_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi connected niche')


# In[36]:


sc.tl.diffmap(HR_Rv_7dpi_conn)
sc.pp.neighbors(HR_Rv_7dpi_conn, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_Rv_7dpi_conn)


# In[37]:


sc.tl.leiden(HR_Rv_7dpi_conn, resolution=0.4)


# In[38]:


sc.pl.umap(HR_Rv_7dpi_conn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[40]:


scv.pp.moments(HR_Rv_7dpi_conn)


# In[41]:


scv.tl.velocity(HR_Rv_7dpi_conn, mode='stochastic')


# In[42]:


scv.tl.velocity_graph(HR_Rv_7dpi_conn)


# In[43]:


scv.pl.velocity_embedding(HR_Rv_7dpi_conn, basis='draw_graph_fa')


# In[44]:


scv.pl.velocity_embedding_grid(HR_Rv_7dpi_conn, basis='draw_graph_fa', color = 'Cell_type')


# In[45]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_conn, basis='draw_graph_fa', title = '7dpi connected niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_connected_niche_7dpi_20n_cell_types_diffmap.png')


# In[46]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_conn, basis='draw_graph_fa', title = '7dpi connected niche', 
                                 color = 'leiden', legend_loc = 'right margin')


# In[47]:


HR_Rv_7dpi_conn.write('./write/HR_Rv_7dpi_conn.h5ad')


# # Endocardial niche at 7dpi

# In[10]:


endo_7dpi = ['Endocardium (A)', 'Endocardium (frzb)', 'Endocardium (V)', 'Fibroblast', 'Fibroblast (nppc)', 'Fibroblast (spock3)', 'Fibroblast-like cells', 'Smooth muscle cells']


# In[11]:


HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
HR_Rv_7dpi_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)]


# In[50]:


sc.pp.filter_genes(HR_Rv_7dpi_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_7dpi_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_endo)


# In[51]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_endo)
sc.tl.pca(HR_Rv_7dpi_endo)
sc.pp.neighbors(HR_Rv_7dpi_endo, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_7dpi_endo, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_endo)


# In[52]:


sc.pl.umap(HR_Rv_7dpi_endo, color='batch',
          title = '7dpi endocardial niche')


# In[56]:


sc.pl.umap(HR_Rv_7dpi_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi endocardial niche')


# In[74]:


sc.pl.umap(HR_Rv_7dpi_endo, color = ['col1a1a'])


# In[57]:


sc.tl.diffmap(HR_Rv_7dpi_endo)
sc.pp.neighbors(HR_Rv_7dpi_endo, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_Rv_7dpi_endo)


# In[58]:


sc.tl.leiden(HR_Rv_7dpi_endo, resolution=0.4)


# In[59]:


sc.pl.umap(HR_Rv_7dpi_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[60]:


sc.tl.paga(HR_Rv_7dpi_endo, groups='leiden')


# In[61]:


sc.pl.paga(HR_Rv_7dpi_endo, show=True, 
           labels = ['', '', '', '', '', 
                     '', '', '', '', '', 
                     '', '', '', '', '',
                    '',''], 
           node_size_scale = 2, save = 'scvelo_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[62]:


sc.tl.draw_graph(HR_Rv_7dpi_endo, init_pos='paga')


# In[63]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, color = 'leiden', title = '7dpi endocardial niche', save = 'scvelo_endocardial_niche_7dpi_20n_leiden_diffmap.png')


# In[64]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, color = ['Cell_type'], title = '7dpi endocardial niche', save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[86]:


HR_Rv_7dpi_endo.obs['cluster_dummy'] = HR_Rv_7dpi_endo.obs['Cell_type'] == 'Fibroblast (spock3)'
sc.pl.umap(HR_Rv_7dpi_endo, color = 'cluster_dummy')


# In[67]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, ncols = 3,
                 color = ['col1a1a', 'fn1b', 'spock3',
                          'nppc', 'cyp26b1', 'frzb', 'acta2'],
                save = '_scvelo_endocardial_niche_7dpi_20n_marker_genes_diffmap.png')


# In[68]:


scv.pp.moments(HR_Rv_7dpi_endo)


# In[69]:


scv.tl.velocity(HR_Rv_7dpi_endo, mode='stochastic')


# In[70]:


scv.tl.velocity_graph(HR_Rv_7dpi_endo)


# In[71]:


scv.pl.velocity_embedding(HR_Rv_7dpi_endo, basis='draw_graph_fa')


# In[72]:


scv.pl.velocity_embedding_grid(HR_Rv_7dpi_endo, basis='draw_graph_fa', color = 'Cell_type')


# In[73]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endo, basis='draw_graph_fa', title = '7dpi endocardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[79]:


scv.pl.proportions(HR_Rv_7dpi_endo, groupby = 'batch')


# How are the speed and coherence?

# In[76]:


scv.tl.velocity_confidence(HR_Rv_7dpi_endo)


# In[78]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_endo, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_endo_niche_7dpi_20n_cell_types_diffmap.png')


# # Endocardial niche for deep injuries at 7dpi

# In[66]:


#HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
#all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
#HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
#HR_Rv_7dpi_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)]
deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6a', 'Hr6v', 'Hr7a', 'Hr7v', 'Hr9']
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)]
HR_Rv_7dpi_deep_endo = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['heart'].isin(deep_injury_libraries)]


# In[67]:


HR_Rv_7dpi_deep_endo.obs['Cell_type'].value_counts()


# In[68]:


sc.pp.filter_genes(HR_Rv_7dpi_deep_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_7dpi_deep_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_deep_endo)


# In[69]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_deep_endo)
sc.tl.pca(HR_Rv_7dpi_deep_endo)
sc.pp.neighbors(HR_Rv_7dpi_deep_endo, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_7dpi_deep_endo, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_deep_endo)


# In[70]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='batch',
          title = '7dpi deep injury endocardial niche')


# In[71]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi deep injury endocardial niche')


# In[72]:


sc.tl.diffmap(HR_Rv_7dpi_deep_endo)
sc.pp.neighbors(HR_Rv_7dpi_deep_endo, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_Rv_7dpi_deep_endo)


# In[73]:


sc.tl.leiden(HR_Rv_7dpi_deep_endo, resolution=0.4)


# In[74]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[75]:


sc.tl.paga(HR_Rv_7dpi_deep_endo, groups='leiden')


# In[77]:


sc.pl.paga(HR_Rv_7dpi_deep_endo, show=True, 
           labels = ['', '', '', '', '', 
                     '', '', '', '', '', 
                     '', '', '', '', '',
                    '','', ''], 
           node_size_scale = 2, save = 'scvelo_deep_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[78]:


sc.tl.draw_graph(HR_Rv_7dpi_deep_endo, init_pos='paga')


# In[79]:


sc.pl.draw_graph(HR_Rv_7dpi_deep_endo, color = 'leiden', title = '7dpi deep endocardial niche', save = 'scvelo_deep_endocardial_niche_7dpi_20n_leiden_diffmap.png')


# In[80]:


sc.pl.draw_graph(HR_Rv_7dpi_deep_endo, color = ['Cell_type'], title = '7dpi deep endocardial niche', save = '_scvelo_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[81]:


sc.pl.draw_graph(HR_Rv_7dpi_deep_endo, ncols = 3,
                 color = ['col1a1a', 'fn1b', 'spock3',
                          'nppc', 'cyp26b1', 'frzb', 'acta2'],
                save = '_scvelo_deep_endocardial_niche_7dpi_20n_marker_genes_diffmap.png')


# In[82]:


scv.pp.moments(HR_Rv_7dpi_deep_endo)


# In[83]:


scv.tl.velocity(HR_Rv_7dpi_deep_endo, mode='stochastic')


# In[84]:


scv.tl.velocity_graph(HR_Rv_7dpi_deep_endo)


# In[85]:


scv.pl.velocity_embedding(HR_Rv_7dpi_deep_endo, basis='draw_graph_fa')


# In[86]:


scv.pl.velocity_embedding_grid(HR_Rv_7dpi_deep_endo, basis='draw_graph_fa', color = 'Cell_type')


# In[87]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_deep_endo, basis='draw_graph_fa', title = '7dpi deep endocardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[106]:


scv.pl.proportions(HR_Rv_7dpi_deep_endo, groupby = 'batch')


# How are the speed and coherence?

# In[107]:


scv.tl.velocity_confidence(HR_Rv_7dpi_deep_endo)


# In[108]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_deep_endo, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_deep_endo_niche_7dpi_20n_cell_types_diffmap.png')


# In[112]:


HR_Rv_7dpi_deep_endo.obs


# Velocity length for nppc fibroblasts in Hr9:

# In[116]:


HR_Rv_7dpi_deep_endo.obs[(HR_Rv_7dpi_deep_endo.obs['Cell_type'] == "Fibroblast (nppc)") & (HR_Rv_7dpi_deep_endo.obs['heart'] == 'Hr9')]


# In[117]:


HR_Rv_Hr9_endo.obs[(HR_Rv_Hr9_endo.obs['Cell_type'] == "Fibroblast (nppc)") & (HR_Rv_Hr9_endo.obs['heart'] == 'Hr9')]


# # Endocardial niche for Hr9 (a deep injury library at 7dpi)

# In[88]:


#HR_Rv_7dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '7']
#all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
#HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
#HR_Rv_7dpi_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(endo_7dpi)]
#deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6a', 'Hr6v', 'Hr7a', 'Hr7v', 'Hr9']
HR_Rv_Hr9_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['heart'].isin(['Hr9'])]
HR_Rv_Hr9_endo = HR_Rv_Hr9_endo[HR_Rv_Hr9_endo.obs['Cell_type'].isin(endo_7dpi)]


# In[89]:


sc.pp.filter_genes(HR_Rv_Hr9_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_Hr9_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_Hr9_endo)


# In[90]:


sc.pp.highly_variable_genes(HR_Rv_Hr9_endo)
sc.tl.pca(HR_Rv_Hr9_endo)
sc.pp.neighbors(HR_Rv_Hr9_endo, n_neighbors=30)
#sc.external.pp.bbknn(HR_Rv_7dpi_deep_endo, batch_key='batch')
sc.tl.umap(HR_Rv_Hr9_endo)


# In[91]:


sc.pl.umap(HR_Rv_Hr9_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_Hr9_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Hr9: 7dpi deep injury endocardial niche')


# In[92]:


sc.tl.diffmap(HR_Rv_Hr9_endo)
sc.pp.neighbors(HR_Rv_Hr9_endo, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_Rv_Hr9_endo)


# In[93]:


sc.tl.leiden(HR_Rv_Hr9_endo, resolution=0.4)


# In[94]:


sc.pl.umap(HR_Rv_Hr9_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[95]:


sc.tl.paga(HR_Rv_Hr9_endo, groups='leiden')


# In[96]:


sc.pl.paga(HR_Rv_Hr9_endo, show=True, 
           labels = ['', '', '', '', '', 
                     '', '', '', '', '', 
                     '', '', '', '', '',
                    '', ''], 
           node_size_scale = 2, save = 'scvelo_Hr9_deep_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[97]:


sc.tl.draw_graph(HR_Rv_Hr9_endo, init_pos='paga')


# In[98]:


sc.pl.draw_graph(HR_Rv_Hr9_endo, color = 'leiden', title = 'Hr9: 7dpi deep endocardial niche', save = 'scvelo_Hr9_deep_endocardial_niche_7dpi_20n_leiden_diffmap.png')


# In[99]:


sc.pl.draw_graph(HR_Rv_Hr9_endo, color = ['Cell_type'], title = 'Hr9: 7dpi deep endocardial niche', save = '_scvelo_Hr9_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[100]:


scv.pp.moments(HR_Rv_Hr9_endo)


# In[101]:


scv.tl.velocity(HR_Rv_Hr9_endo, mode='stochastic')


# In[102]:


scv.tl.velocity_graph(HR_Rv_Hr9_endo)


# In[103]:


scv.pl.velocity_embedding(HR_Rv_Hr9_endo, basis='draw_graph_fa')


# In[104]:


scv.pl.velocity_embedding_grid(HR_Rv_Hr9_endo, basis='draw_graph_fa', color = 'Cell_type')


# In[105]:


scv.pl.velocity_embedding_stream(HR_Rv_Hr9_endo, basis='draw_graph_fa', title = 'Hr9: 7dpi deep endocardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_Hr9_deep_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[109]:


scv.pl.proportions(HR_Rv_Hr9_endo, groupby = 'Cell_type')


# How are the speed and coherence?

# In[110]:


scv.tl.velocity_confidence(HR_Rv_Hr9_endo)


# In[111]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_Hr9_endo, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_Hr9_endo_niche_7dpi_20n_cell_types_diffmap.png')


# # Endocardial niche at 3dpi

# In[111]:


HR_Rv_3dpi = HR_Rv_filter[HR_Rv_filter.obs['dpi'] == '3']
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
HR_Rv_3dpi_endo = HR_Rv_3dpi[HR_Rv_3dpi.obs['Cell_type'].isin(['Endocardium (V)', 'Fibroblast (nppc)', 'Fibroblast (spock3)', 'Fibroblast-like cells'])]


# In[112]:


sc.pp.filter_genes(HR_Rv_3dpi_endo, min_cells=3)
sc.pp.normalize_per_cell(HR_Rv_3dpi_endo, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_3dpi_endo)


# In[113]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_endo)
sc.tl.pca(HR_Rv_3dpi_endo)
sc.pp.neighbors(HR_Rv_3dpi_endo, n_neighbors=30)
sc.external.pp.bbknn(HR_Rv_3dpi_endo, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_endo)


# In[114]:


sc.pl.umap(HR_Rv_3dpi_endo, color='batch',
          title = '3dpi endocardial niche')


# In[115]:


sc.pl.umap(HR_Rv_3dpi_endo, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi endocardial niche')


# In[116]:


HR_Rv_3dpi_endo.obs['Cell_type'].value_counts()


# In[96]:


sc.tl.diffmap(HR_Rv_3dpi_endo)
sc.pp.neighbors(HR_Rv_3dpi_endo, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_Rv_3dpi_endo)


# In[97]:


sc.tl.leiden(HR_Rv_3dpi_endo, resolution=0.4)


# In[98]:


sc.pl.umap(HR_Rv_3dpi_endo, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[60]:


sc.tl.paga(HR_Rv_7dpi_endo, groups='leiden')


# In[61]:


sc.pl.paga(HR_Rv_7dpi_endo, show=True, 
           labels = ['', '', '', '', '', 
                     '', '', '', '', '', 
                     '', '', '', '', '',
                    '',''], 
           node_size_scale = 2, save = 'scvelo_endo_niche_7dpi_20n_leiden_diffmap.png')


# In[62]:


sc.tl.draw_graph(HR_Rv_7dpi_endo, init_pos='paga')


# In[63]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, color = 'leiden', title = '7dpi endocardial niche', save = 'scvelo_endocardial_niche_7dpi_20n_leiden_diffmap.png')


# In[64]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, color = ['Cell_type'], title = '7dpi endocardial niche', save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[67]:


sc.pl.draw_graph(HR_Rv_7dpi_endo, ncols = 3,
                 color = ['col1a1a', 'fn1b', 'spock3',
                          'nppc', 'cyp26b1', 'frzb', 'acta2'],
                save = '_scvelo_endocardial_niche_7dpi_20n_marker_genes_diffmap.png')


# In[68]:


scv.pp.moments(HR_Rv_7dpi_endo)


# In[69]:


scv.tl.velocity(HR_Rv_7dpi_endo, mode='stochastic')


# In[70]:


scv.tl.velocity_graph(HR_Rv_7dpi_endo)


# In[71]:


scv.pl.velocity_embedding(HR_Rv_7dpi_endo, basis='draw_graph_fa')


# In[72]:


scv.pl.velocity_embedding_grid(HR_Rv_7dpi_endo, basis='draw_graph_fa', color = 'Cell_type')


# In[73]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endo, basis='draw_graph_fa', title = '7dpi endocardial niche', 
                                 color = 'Cell_type', legend_loc = 'right margin',
                                save = '_scvelo_endocardial_niche_7dpi_20n_cell_types_diffmap.png')


# In[79]:


scv.pl.proportions(HR_Rv_7dpi_endo, groupby = 'batch')


# How are the speed and coherence?

# In[76]:


scv.tl.velocity_confidence(HR_Rv_7dpi_endo)


# In[78]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(HR_Rv_7dpi_endo, c=keys, cmap='coolwarm', basis='draw_graph_fa', perc=[5, 95],
                                save = 'strength_coherence_scvelo_endo_niche_7dpi_20n_cell_types_diffmap.png')


# # Full niche

# In[5]:


HR_niche = HR_filter[HR_filter.obs['Cell_type'].isin(trajectory_subset)]


# In[7]:


sc.pp.filter_genes(HR_niche, min_cells=3)
sc.pp.normalize_per_cell(HR_niche, counts_per_cell_after=1e4,copy=False)
sc.pp.log1p(HR_niche, copy=False)
sc.pp.highly_variable_genes(HR_niche)
sc.tl.pca(HR_niche)
sc.pp.neighbors(HR_niche, n_neighbors=30)
sc.external.pp.bbknn(HR_niche, batch_key='batch')
sc.tl.umap(HR_niche)


# In[8]:


sc.pl.umap(HR_niche, color='batch',
          title = 'Niche')#, save = 'subset1_3dpi_umap_batch.png')


# In[9]:


sc.pl.umap(HR_niche, color='dpi',
          title = 'Niche')


# In[13]:


HR_niche.obs


# In[15]:


sc.pl.umap(HR_niche, color='Cell_type', palette = cell_type_colors.loc[HR_niche.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[16]:


HR_niche.obs.Cell_type.cat.categories.tolist()

