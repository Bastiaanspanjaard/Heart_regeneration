
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


# # Load and annotate single-cell data

# In[2]:


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


# In[48]:


HR =    H5_data.concatenate(H6_data, H7_data, H8a_data, H8v_data,
                        Hr1_data, Hr2a_data, Hr2b_data, Hr3_data, Hr4_data, Hr5_data, 
                        Hr6a_data, Hr6v_data, Hr7a_data, Hr7v_data, Hr8_data, Hr9_data, Hr10_data,
                        Hr11_data, Hr12_data, Hr13_data, Hr14_data, Hr15_data,
                        Hr16_data, Hr17_data, Hr18_data, Hr19_data, Hr20_data,
                        Hr21_data, Hr22_data, Hr23_data, Hr24_data, Hr25_data,
                        Hr26_data, Hr27_data, Hr28_data, Hr29_data, Hr30_data,
                        Hr31_data, Hr32_data, Hr33_data, Hr34_data, Hr35_data)
HR.shape


# In[49]:


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


# In[50]:


HR.obs_names = [str(HR.obs.loc[x,'heart'])+'_'+str(x.split('-', 1)[0]) for x in HR.obs_names]


# In[51]:


annotations = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged.csv', index_col = 0)


# In[52]:


# Drop annotations that are not in the single-cell object
#set_A.difference(set_B) for (A - B)
anno_drop = annotations.index.difference(HR.obs_names)
annotations = annotations.drop(anno_drop)


# In[53]:


HR_filter = HR[annotations.index]
#annotations_filter = annotations[annotations.index.isin(HR.obs_names.intersection(annotations.index))]
HR_filter


# In[54]:


HR_filter.obs['Cell_type'] = annotations_filter['Cell_type'].tolist()


# # Pseudotime on fibroblasts, fibroblast-like cells and epicardium at 3dpi

# In[55]:


HR_ps_1 = HR_filter[HR_filter.obs['dpi'] == '3']


# In[56]:


sc.pp.filter_genes(HR_ps_1, min_cells=3)
HR_ps_1_norm = sc.pp.normalize_per_cell(HR_ps_1, counts_per_cell_after=1e4,copy=True)
HR_ps_1 = sc.pp.log1p(HR_ps_1_norm, copy=True)
sc.pp.highly_variable_genes(HR_ps_1)
sc.tl.pca(HR_ps_1)
sc.pp.neighbors(HR_ps_1, n_neighbors=30)
sc.tl.umap(HR_ps_1)


# In[57]:


sc.pl.umap(HR_ps_1, color='Cell_type')


# In[44]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors.csv', index_col = 0)


# In[36]:


cell_type_colors.index


# In[47]:


np.setdiff1d(np.unique(HR_ps_1.obs.Cell_type), cell_type_colors.index)


# In[32]:


(np.unique(HR_ps_1.obs.Cell_type))


# In[ ]:


# Add correct colors
# Subset celltypes
# Do trajectory

