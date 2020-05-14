
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


# In[2]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors_unmerged.csv', index_col = 0)


# In[19]:


trajectory_subset = ['Fibroblast', 'Fibroblast (cfd)', 'Fibroblast (col11a1a)', 'Fibroblast (col12a1a)', 
                     'Fibroblast (cxcl12a)', 'Fibroblast (mpeg1.1)', 'Fibroblast (nppc)', 'Fibroblast (spock3)',
                    'Fibroblast-like cells', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblast (proliferating)', 'Perivascular cells']


# In[ ]:


connected_3dpi = ['Fibroblast', 'Fibroblast (cfd)', 'Fibroblast (col11a1a)', 'Fibroblast (col12a1a)', 
                     'Fibroblast (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblast (proliferating)', 'Perivascular cells']


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
HR_obs = HR.obs
HR.obs = HR.obs.reset_index().merge(HR_setnames, how="inner").set_index('index')


# In[28]:


# Rename cells so the cell names correspond to the ones in the annotation file
HR.obs_names = [str(HR.obs.loc[x,'heart'])+'_'+str(x.split('-', 1)[0]) for x in HR.obs_names]


# In[29]:


annotations = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged.csv', index_col = 0)


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

# # Pseudotime on connected cell types 3dpi

# In[40]:


HR_3d_conn = HR_ps_1[HR_ps_1.obs['Cell_type'].isin(connected_3dpi)]


# In[41]:


sc.pp.highly_variable_genes(HR_3d_conn)
sc.tl.pca(HR_3d_conn)
sc.pp.neighbors(HR_3d_conn, n_neighbors=30)
sc.external.pp.bbknn(HR_3d_conn, batch_key='batch')
sc.tl.umap(HR_3d_conn)


# In[94]:


sc.pl.umap(HR_3d_conn, color='batch',
          title = '3dpi connected niche', save = 'connected_niche_3dpi_umap_batch.png')


# In[95]:


sc.pl.umap(HR_3d_conn, color='Cell_type', palette = fibro_colors.loc[HR_3d_conn.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '3dpi connected niche', save = 'connected_niche_3dpi_umap.png')


# In[45]:


sc.tl.diffmap(HR_3d_conn)
sc.pp.neighbors(HR_3d_conn, n_neighbors=20, use_rep='X_diffmap')
sc.tl.draw_graph(HR_3d_conn)


# In[53]:


sc.tl.leiden(HR_3d_conn, resolution=0.5)


# In[79]:


sc.pl.umap(HR_3d_conn, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# 0: shared col12, col11, proliferating
# 1: part of fibroblast
# 2: part of col11
# 3: part of fibroblast
# 4: part of col11
# 5: part of epi A
# 6: part of epi A
# 7: part of fibroblast
# 8: part of fibroblast
# 9: epi V
# 10: part of fibroblast
# 11: part of fibroblast
# 12: part of col12
# 13: perivascular
# 14: mpeg1.1
# 15: shared fibroblast/prolif

# In[55]:


sc.tl.paga(HR_3d_conn, groups='leiden')


# In[89]:


sc.pl.paga(HR_3d_conn, show=True, labels = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''], node_size_scale = 2, save = 'connected_niche_3dpi_20n_leiden_diffmap.png')


# In[83]:


sc.tl.draw_graph(HR_3d_conn, init_pos='paga')


# In[92]:


sc.pl.draw_graph(HR_3d_conn, color = 'leiden', title = '3dpi connected niche', save = 'connected_niche_3dpi_20n_diffmap_leidenshow_pagainit.png')


# In[90]:


sc.pl.draw_graph(HR_3d_conn, color = ['Cell_type'], title = '3dpi connected niche', save = 'connected_niche_3dpi_20n_leiden_diffmap_pagainit.png')


# In[124]:


sc.pl.draw_graph(HR_3d_conn, ncols = 3,
                 color = ['col1a1a', 'tcf21', 'tbx18',
                          'notch3', 'pdgfrb', 'mpeg1.1', 'cfd', 
                          'postnb', 'col11a1a', 'col12a1a', 
                          'pcna', 'aldh1a2', 'stra6'])


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

