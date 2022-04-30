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
import cellrank as cr
scv.settings.set_figure_params('scvelo', dpi_save = 300)


# In[2]:


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


# In[3]:


annotations = pd.read_csv('~/Documents/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv', index_col = 0)


# In[4]:


cell_type_colors = pd.read_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Cell_type_colors_2.csv', index_col = 0)


# In[5]:


epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (cfd)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (cxcl12a)', 'Fibroblasts (mpeg1.1)', 'Fibroblasts (nppc)', 'Fibroblasts (spock3)',
                    'Valve fibroblasts', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[6]:


connected_epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Fibroblasts (mpeg1.1)','Epicardium (Atrium)', 'Epicardium (Ventricle)','Fibroblasts (cxcl12a)',
                    'Fibroblasts (proliferating)', 'Perivascular cells']


# In[7]:


local_connected_epifibro_types = ['Fibroblasts (const.)', 'Fibroblasts (col11a1a)', 'Fibroblasts (col12a1a)', 
                     'Epicardium (Ventricle)',
                    'Fibroblasts (proliferating)', 'Perivascular cells', 'Fibroblasts (cxcl12a)']


# In[8]:


connected_endofibro_types = ['Endocardium (frzb)', 'Endocardium (Ventricle)', 
                            'Fibroblasts (nppc)', 'Fibroblasts (spock3)']


# In[9]:


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


# In[11]:


#HR_Rv_filter.write('./write/HR_Rv_filter.h5ad')
HR_Rv_filter = sc.read('./write/HR_Rv_filter.h5ad')


# In[12]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]


# In[13]:


HR_Rv_3dpi_winh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi_winh.var_names if not name == 'RFP']
HR_Rv_3dpi_winh = HR_Rv_3dpi_winh[:, all_genes_but_RFP]


# In[83]:


HR_Rv_7dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]


# In[15]:


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
#HR_Rv_3dpi_epiconn.obs['n_counts'] = HR_Rv_3dpi_epiconn.X.sum(axis=1)
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


sc.pl.umap(HR_Rv_3dpi_epiconn, color='tcf21')


# In[33]:


sc.pl.violin(HR_Rv_3dpi_epiconn, keys = 'tcf21', groupby = 'Cell_type')


# In[34]:


scv.pp.moments(HR_Rv_3dpi_epiconn)


# In[35]:


# For dynamic model
#scv.tl.recover_dynamics(HR_Rv_3dpi_epiconn, n_jobs=8)


# In[36]:


#scv.tl.velocity(HR_Rv_3dpi_epiconn, mode="dynamical")
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


# In[51]:


# Dynamic
scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', 
                                 color = 'Cell_type')


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


# In[45]:


scv.tl.rank_velocity_genes(HR_Rv_3dpi_epiconn, groupby='leiden')

df = scv.DataFrame(HR_Rv_3dpi_epiconn.uns['rank_velocity_genes']['names'])
df.head()
# Note the difference between col12a1A and col12a1B!


# In[50]:


df['6'][0:5].tolist()


# In[51]:


scv.pl.velocity(HR_Rv_3dpi_epiconn, df['6'][0:5].tolist(), basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)


# ## Determining transitions with CellRank

# In[52]:


cr.tl.terminal_states(HR_Rv_3dpi_epiconn, show_plots=True) # This does not really work


# In[42]:


from cellrank.tl.kernels import VelocityKernel
vk = VelocityKernel(HR_Rv_3dpi_epiconn)


# In[43]:


vk.compute_transition_matrix()


# In[44]:


from cellrank.tl.kernels import ConnectivityKernel
ck = ConnectivityKernel(HR_Rv_3dpi_epiconn).compute_transition_matrix()


# In[45]:


combined_kernel = 0.8 * vk + 0.2 * ck


# In[46]:


from cellrank.tl.estimators import GPCCA
g = GPCCA(combined_kernel)
print(g)


# In[47]:


g.compute_schur(n_components=20)
g.plot_spectrum()


# In[54]:


g.compute_macrostates(n_states=7, cluster_key="leiden")
g.plot_macrostates()


# In[50]:


g.plot_macrostates(same_plot=False)


# In[50]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', color = 'Cell_type')


# In[51]:


scv.pl.velocity_embedding_stream(HR_Rv_3dpi_epiconn, title = '', color = 'leiden')


# In[51]:


g.set_terminal_states_from_macrostates()


# In[52]:


g.compute_absorption_probabilities()


# In[53]:


cr.pl.circular_projection(HR_Rv_3dpi_epiconn, keys="Cell_type", legend_loc="right")


# In[53]:


scv.tl.paga(
    HR_Rv_3dpi_epiconn,
    groups="leiden")


# In[54]:


cr.pl.cluster_fates(
    HR_Rv_3dpi_epiconn,
    mode="paga_pie",
    cluster_key="Cell_type")
#,
#    basis="umap",
#    legend_kwargs={"loc": "top right out"},
#    legend_loc="top left out",
#    node_size_scale=5,
#    edge_width_scale=1,
#    max_edge_width=4,
#    title="directed PAGA",
#)


# In[41]:


HR_Rv_3dpi_epiconn.obs


# In[40]:


cr.pl.terminal_states(HR_Rv_3dpi_epiconn)


# In[42]:


cr.tl.initial_states(HR_Rv_3dpi_epiconn, cluster_key="Cell_type")
cr.pl.initial_states(HR_Rv_3dpi_epiconn, discrete=True)


# In[43]:


cr.tl.lineages(HR_Rv_3dpi_epiconn)
cr.pl.lineages(HR_Rv_3dpi_epiconn, same_plot=False)


# ## Which genes drive the transitions?

# In[ ]:





# In[37]:


sc.pl.paga(HR_Rv_3dpi_epiconn, show=True, color = 'Cell_type', node_size_scale = 2)


# In[38]:


sc.pl.umap(HR_Rv_3dpi_epiconn, color=['leiden', 'Cell_type'], legend_loc='on data', legend_fontsize='x-large')


# For each cell of type 1, extract the maximum transition probability to type 2. The mean and distribution of these probabilities tell us how likely a source type 1 is for type 2. We compare these maximum transition probabilities to transition probabilities of type 1 to 'anything but type 1+2' to estimate whether the transition probabilities from type 1 to type 2 are significantly higher than a null hypothesis.
# The questions we want to answer with this are:  
# 1) Is there a difference between transitions from Epi V to col11/col12 at 3dpi and at 7dpi?  
# 2) Is the seeming transition from constitutional fibroblasts to epicardium as meaningful as that from constitutional fibroblasts to col11/col12 at 3dpi and 7dpi? At 3dpi, the constitutional fibroblast clusters involved are 2 and 7; col11/col12 are 4,5,11,13; epicardium is 0,3,10
# 
# Slightly different plan: for each cell, which cell does it best 'point' to? Then aggregate over cell types or clusters to understand either where cell types want to go.

# In[59]:


HR_Rv_3dpi_epiconn.uns['velocity_graph_neg']


# In[60]:


HR_Rv_3dpi_epiconn.uns['velocity_graph']


# In[57]:


HR_Rv_3dpi_epiconn.uns['velocity_graph'][0,]


# In[58]:


HR_Rv_3dpi_epiconn.uns['neighbors']['connectivities'][0, ]


# In[63]:


ref_ix = 7
C = HR_Rv_3dpi_epiconn.obsp['connectivities']
neighbors_ix = C[ref_ix, :].indices


# In[72]:


# get the computed velocity graph values
vgraph = HR_Rv_3dpi_epiconn.uns['velocity_graph']
vgraph_neg = HR_Rv_3dpi_epiconn.uns['velocity_graph_neg']
vgraph_tot = HR_Rv_3dpi_epiconn.uns['velocity_graph'] + HR_Rv_3dpi_epiconn.uns['velocity_graph_neg']


# In[76]:


vgraph[ref_ix, :].data[0:10]


# In[77]:


vgraph_neg[ref_ix, :].data[0:10]


# In[78]:


vgraph_tot[ref_ix, :].data[0:10]


# In[39]:


cell_number = HR_Rv_3dpi_epiconn.obs.shape[0]
# Extract transition probabilities
transition_probabilities = (HR_Rv_3dpi_epiconn.uns['velocity_graph']) # Note there's also a graph that contains negative cosine similarities. 
# For our purposes, it makes sense to only consider the positive ones.
# Find which cell maximum probability points to and make linking matrix (maximum_probability_index_matrix)
maximum_probability_index_array = [transition_probabilities.getrow(i).argmax() for i in range(cell_number)]
maximum_probability_index_matrix = np.zeros([cell_number, cell_number])
maximum_probability_index_matrix[range(cell_number), maximum_probability_index_array] = 1
# M is cell * cell that links every row cell to its maximum transition probability cell (column)
# To aggregate this, we need a cell * Cell_type matrix C, and then C^T ** M ** C counts how many cells from each cluster go where.
# Create matrix using pivot.
df_for_pivot = HR_Rv_3dpi_epiconn.obs[['Cell_type']]
df_for_pivot = df_for_pivot.assign(Presence = pd.Series([1] * df_for_pivot.shape[0]).values)
cells_types_matrix = df_for_pivot.pivot(columns = 'Cell_type', values = 'Presence')
cells_types_matrix = cells_types_matrix.fillna(0)
# Aggregate over cell types
inter_result = pd.DataFrame(data = maximum_probability_index_matrix.dot(cells_types_matrix), index = cells_types_matrix.index, columns = cells_types_matrix.columns)
max_prob_trans = (cells_types_matrix.transpose()).dot(inter_result)
max_prob_trans


# In[40]:


transition_fractions = (max_prob_trans.transpose()/max_prob_trans.sum(axis = 1)).transpose()


# In[124]:


# This normalization doesn't work:
ax = sns.heatmap(transition_fractions, linewidth=0.5)
pl.show()


# In[134]:


max_prob_trans.sum(axis = 1)


# In[41]:


expected_trans = np.outer(max_prob_trans.sum(axis = 1), max_prob_trans.sum(axis = 1)/cell_number)
max_prob_trans - expected_trans


# In[42]:


type_fractions = max_prob_trans.sum(axis = 1)/cell_number
standard_devs = np.sqrt(np.outer(max_prob_trans.sum(axis = 1), type_fractions * (1 - type_fractions)))


# In[43]:


z_scaled_transitions = (max_prob_trans - expected_trans)/standard_devs
z_scaled_transitions


# In[44]:


ax = sns.heatmap(z_scaled_transitions, linewidth=0.5, center = 0, cbar_kws={'label': 'Transition z-score'})
ax.set(xlabel='Transition to', ylabel='Transition from')
#pl.savefig('./figures/Z_scale_transitions_3dpi_epiniche.png', bbox_inches = 'tight')
pl.show()


# In[45]:


p_scaled_transitions = max_prob_trans/expected_trans
p_scaled_transitions


# In[46]:


ax = sns.heatmap(p_scaled_transitions, linewidth=0.5, center = 1, cbar_kws={'label': 'Transition p-score'})
ax.set(xlabel='Transition to', ylabel='Transition from')
#pl.savefig('./figures/Z_scale_transitions_3dpi_epiniche.png', bbox_inches = 'tight')
pl.show()


# In[123]:


transition_fractions


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


# In[45]:


# Epicardium to col11/col12
leiden_1 = ['0', '3', '10']
leiden_2 = ['4', '5', '11', '13']
# Find (numeric) indices for both types
leiden_1_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_1] 
leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if HR_Rv_3dpi_epiconn.obs.leiden[i] in leiden_2]
# Subset velocity graph for both types
transitions_1to2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, leiden_2_index]
# Take row maxima
#max_transitions_1to2 = transitions_1to2.max(axis=1).toarray()
# Subset velocity graph for type 1 and not-type-2
#not_leiden_2_index = [i for i in range(len(HR_Rv_3dpi_epiconn.obs.leiden)) if not HR_Rv_3dpi_epiconn.obs.leiden[i] in (leiden_1 + leiden_2)]
#transitions_1tonot2 = (HR_Rv_3dpi_epiconn.uns['velocity_graph'])[leiden_1_index,:][:, not_leiden_2_index]
# Take row maxima
#max_transitions_1tonot2 = transitions_1tonot2.max(axis=1).toarray()
# Plot distributions
#fig, (ax1, ax2) = pl.subplots(1, 2, sharex=True)
#fig.suptitle('3dpi', y = 1.05)
#ax1.set_title('Epicardium to col11/col12')
#ax2.set_title('Epicardium to rest')
#ax1.hist(max_transitions_1to2[max_transitions_1to2 != 0])
#ax2.hist(max_transitions_1tonot2)
#pl.savefig('figures/Epi_col1112_transitions_3dpi.png')
#print(np.mean(max_transitions_1to2[max_transitions_1to2 != 0]))
#print(np.mean(max_transitions_1tonot2[max_transitions_1tonot2 != 0]))


# In[55]:


sum(HR_Rv_3dpi_epiconn.obs.leiden.isin(leiden_2))


# In[50]:


transitions_1to2


# In[64]:


#transitions_1to2.todense()[0, :]
(transitions_1to2.todense()[1, :] == 0).sum(axis = 1)


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


# # Trajectories in 3dpi epicardial connected WNT-inhibited niche

# In[16]:


HR_Rv_3dpi_epiconn_winh = HR_Rv_3dpi_winh[HR_Rv_3dpi_winh.obs['Cell_type'].isin(local_connected_epifibro_types)] #connected_epifibro_types)]
sc.pp.filter_genes(HR_Rv_3dpi_epiconn_winh, min_cells=3)


# In[17]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_3dpi_epiconn_winh.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_3dpi_epiconn_winh.obs['percent_mito'] = np.sum(
    HR_Rv_3dpi_epiconn_winh[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_3dpi_epiconn_winh.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_3dpi_epiconn_winh, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_3dpi_epiconn_winh, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[18]:


sc.pp.normalize_per_cell(HR_Rv_3dpi_epiconn_winh, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_3dpi_epiconn_winh)


# In[19]:


sc.pp.regress_out(HR_Rv_3dpi_epiconn_winh, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[20]:


sc.pp.highly_variable_genes(HR_Rv_3dpi_epiconn_winh)
sc.tl.pca(HR_Rv_3dpi_epiconn_winh)
sc.external.pp.bbknn(HR_Rv_3dpi_epiconn_winh, batch_key='batch')
sc.tl.umap(HR_Rv_3dpi_epiconn_winh)


# In[21]:


sc.pl.umap(HR_Rv_3dpi_epiconn_winh, color='batch',
          title = '7dpi epicardial niche')


# In[22]:


sc.pl.umap(HR_Rv_3dpi_epiconn_winh, color=['total_counts', 'n_genes_by_counts', 'percent_mito'],
          title = '7dpi epicardial niche')


# In[23]:


sc.pl.umap(HR_Rv_3dpi_epiconn_winh, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_3dpi_epiconn_winh.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi epicardial niche')


# In[24]:


sc.tl.leiden(HR_Rv_3dpi_epiconn_winh, resolution=2)


# In[25]:


sc.pl.umap(HR_Rv_3dpi_epiconn_winh, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[26]:


sc.tl.paga(HR_Rv_3dpi_epiconn_winh, groups='leiden')


# In[27]:


sc.pl.paga(HR_Rv_3dpi_epiconn_winh, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[28]:


sc.tl.umap(HR_Rv_3dpi_epiconn_winh, init_pos = 'paga')


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


# In[110]:


epi_velo_genes = ['col11a1a', 'col12a1a', 'fn1a', 'postnb', 'stra6','nrg1']
scv.pl.velocity(HR_Rv_3dpi_epiconn_winh, epi_velo_genes, basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)#,
                #save = 'gene_velocity_local_epicardial_niche_7dpi_noinhib_umap.png')


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


# In[73]:


scv.pl.velocity(HR_Rv_7dpi_epiconn, ['serpine1', 'col12a1b', 'map1ab'], basis='umap', ncols=1, fontsize=16, figsize=(9,6), dpi = 300)


# In[72]:


scv.tl.rank_velocity_genes(HR_Rv_7dpi_epiconn, groupby='leiden')

df = scv.DataFrame(HR_Rv_7dpi_epiconn.uns['rank_velocity_genes']['names'])
df.head()


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


# ## Transitions

# In[175]:


sc.pl.paga(HR_Rv_7dpi_deep_endo, show=True, color = 'Cell_type', node_size_scale = 2)


# In[195]:


sc.pl.umap(HR_Rv_7dpi_deep_endo, color=['leiden', 'Cell_type'], legend_loc='on data', legend_fontsize='x-large', ncols = 1, save = '_leiden_celltypes_7dpi_endo.png')


# For each cell of type 1, extract the maximum transition probability to type 2. The mean and distribution of these probabilities tell us how likely a source type 1 is for type 2. We compare these maximum transition probabilities to transition probabilities of type 1 to 'anything but type 1+2' to estimate whether the transition probabilities from type 1 to type 2 are significantly higher than a null hypothesis.
# The questions we want to answer with this are:  
# 1) Is there a difference between transitions from Epi V to col11/col12 at 3dpi and at 7dpi?  
# 2) Is the seeming transition from constitutional fibroblasts to epicardium as meaningful as that from constitutional fibroblasts to col11/col12 at 3dpi and 7dpi? At 3dpi, the constitutional fibroblast clusters involved are 2 and 7; col11/col12 are 4,5,11,13; epicardium is 0,3,10
# 
# Slightly different plan: for each cell, which cell does it best 'point' to? Then aggregate over cell types or clusters to understand either where cell types want to go.

# In[184]:


cell_number = HR_Rv_7dpi_deep_endo.obs.shape[0]
# Extract transition probabilities
transition_probabilities = (HR_Rv_7dpi_deep_endo.uns['velocity_graph'])
# Find which cell maximum probability points to and make linking matrix (maximum_probability_index_matrix)
maximum_probability_index_array = [transition_probabilities.getrow(i).argmax() for i in range(cell_number)]
maximum_probability_index_matrix = np.zeros([cell_number, cell_number])
maximum_probability_index_matrix[range(cell_number), maximum_probability_index_array] = 1
# M is cell * cell that links every row cell to its maximum transition probability cell (column)
# To aggregate this, we need a cell * Cell_type matrix C, and then C^T ** M ** C counts how many cells from each cluster go where.
# Create matrix using pivot.
df_for_pivot = HR_Rv_7dpi_deep_endo.obs[['leiden']]
df_for_pivot = df_for_pivot.assign(Presence = pd.Series([1] * df_for_pivot.shape[0]).values)
cells_types_matrix = df_for_pivot.pivot(columns = 'leiden', values = 'Presence')
cells_types_matrix = cells_types_matrix.fillna(0)
# Aggregate over cell types
inter_result = pd.DataFrame(data = maximum_probability_index_matrix.dot(cells_types_matrix), index = cells_types_matrix.index, columns = cells_types_matrix.columns)
max_prob_trans = (cells_types_matrix.transpose()).dot(inter_result)
max_prob_trans


# In[185]:


transition_fractions = (max_prob_trans.transpose()/max_prob_trans.sum(axis = 1)).transpose()


# In[186]:


# This normalization doesn't work:
ax = sns.heatmap(transition_fractions, linewidth=0.5)
pl.show()


# In[187]:


max_prob_trans.sum(axis = 1)


# In[188]:


expected_trans = np.outer(max_prob_trans.sum(axis = 1), max_prob_trans.sum(axis = 1)/cell_number)
max_prob_trans - expected_trans


# In[189]:


type_fractions = max_prob_trans.sum(axis = 1)/cell_number
standard_devs = np.sqrt(np.outer(max_prob_trans.sum(axis = 1), type_fractions * (1 - type_fractions)))


# In[190]:


z_scaled_transitions = (max_prob_trans - expected_trans)/standard_devs
z_scaled_transitions


# In[191]:


ax = sns.heatmap(z_scaled_transitions, linewidth=0.5, center = 0, cbar_kws={'label': 'Transition z-score'})
ax.set(xlabel='Transition to', ylabel='Transition from')
pl.savefig('./figures/Z_scale_transitions_7dpi_endoniche.png', bbox_inches = 'tight')
pl.show()


# # Nppc and ventricular endocardium

# In[12]:


nppc_endo = ['Endocardium (Ventricle)', 'Fibroblasts (nppc)']


# In[30]:


adata = HR_Rv_filter[(HR_Rv_filter.obs['Cell_type'].isin(nppc_endo)) & (HR_Rv_filter.obs['inhib'] != 'IWR1') & HR_Rv_filter.obs['dpi'].isin(['0','7'])]
#deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9'] # All samples with > 50 nppc fibroblasts.
#HR_Rv_7dpi_deep_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(all_endocardial_lineage)]#connected_endofibro_types)]
#HR_Rv_7dpi_deep_endo = HR_Rv_7dpi_deep_endo[HR_Rv_7dpi_deep_endo.obs['heart'].isin(deep_injury_libraries)]
scv.pp.remove_duplicate_cells(adata)
sc.pp.filter_genes(adata, min_cells=3)


# In[31]:


# Find mito_genes in dataset.
mito_in_index = list(set(adata.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_in_index].X, axis=1) / np.sum(adata.X, axis=1)
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[32]:


adata


# In[33]:


sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)


# In[34]:


sc.pp.regress_out(adata, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[35]:


sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
sc.tl.umap(adata)


# In[36]:


sc.pl.umap(adata, color='Cell_type', palette = cell_type_colors.loc[adata.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Ventricular endocardium and nppc fibroblasts')


# In[37]:


sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[38]:


sc.pl.umap(adata, color=['batch','dpi'],
          title = 'Ventricular endocardium and nppc fibroblasts')


# In[80]:


endmt_genes = ["nppc", "cdh5", "vegfaa"]#"vwf", "snai2", "twist1b", "sparc", "acta2", "icn"]


# In[81]:


sc.pl.umap(adata, color = endmt_genes, ncols = 2)


# In[82]:


sc.pl.violin(adata[adata.obs['Cell_type']=='Endocardium (Ventricle)'], keys = ['nppc', 'cdh5'], groupby = 'dpi')


# In[73]:


sc.pl.violin(adata[adata.obs['Cell_type']=='Fibroblasts (nppc)'], keys = ['nppc', 'cdh5'], groupby = 'dpi')


# In[101]:


HR_Rv_7dpi = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(['Fibroblasts (nppc)', 'Endocardium (Ventricle)', 
                                                          'Endocardium (Atrium)', 'Endocardium (frzb)'])]
scv.pp.remove_duplicate_cells(HR_Rv_7dpi)
sc.pp.filter_genes(HR_Rv_7dpi, min_cells=3)
mito_in_index = list(set(HR_Rv_7dpi.var.index.values) & set(mito_genes))
HR_Rv_7dpi.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi, percent_top=None, log1p=True, inplace=True)
sc.pp.normalize_per_cell(HR_Rv_7dpi, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi)


# In[102]:


sc.pl.violin(HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(['Fibroblasts (nppc)', 'Endocardium (Ventricle)', 
                                                          'Endocardium (Atrium)', 'Endocardium (frzb)'])], 
                        keys = 'nppc', groupby = 'Cell_type')


# In[103]:


sc.pp.highly_variable_genes(HR_Rv_7dpi)
sc.tl.pca(HR_Rv_7dpi)
sc.pp.neighbors(HR_Rv_7dpi)
sc.tl.umap(HR_Rv_7dpi)


# In[105]:


sc.pl.umap(HR_Rv_7dpi, color = ['batch'])


# In[110]:


sc.pl.umap(HR_Rv_7dpi, color = ['Cell_type'],
          palette = cell_type_colors.loc[HR_Rv_7dpi.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[106]:


endmt_genes = ["nppc", "cdh5", "vegfaa", "vwf", "snai2", "twist1b", "sparc", "acta2", "icn"]


# In[108]:


sc.pl.umap(HR_Rv_7dpi, color = endmt_genes, ncols = 2)


# In[118]:


sc.pl.violin(HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(['Endocardium (Ventricle)', 'Fibroblasts (nppc)'])], 
             keys = 'cdh5', groupby = 'Cell_type',
            save = 'cdh5_nppc_endo_7dpi')


# In[ ]:


['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9']


# In[116]:


sc.pl.violin(HR_Rv_7dpi[HR_Rv_7dpi.obs['heart'] == 'Hr6v'], keys = 'cdh5', groupby = 'Cell_type')


# In[88]:


sc.pl.violin(HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type']=='Endocardium (Ventricle)'], keys = ['nppc', 'cdh5'])


# In[ ]:





# In[52]:


adata_2 = adata[adata.obs['Cell_type'] == "Fibroblasts (nppc)"]


# In[53]:


sc.tl.pca(adata_2)


# In[55]:


sc.pl.pca(adata_2, components = ['1,2', '3,4', '5,6', '7,8'], color = ['nppc', 'cdh5', 'vwf'])


# In[63]:


sc.pl.pca_loadings(adata_2, components = [1, 2])


# In[61]:


sc.pl.pca(adata, components = ['1,2', '3,4', '5,6', '7,8'], color = ['nppc', 'cdh5', 'vwf', 'Cell_type'])


# # Endocardial niche for WNT inhibition at 7dpi

# In[97]:


all_endocardial_lineage = ['Endocardium (frzb)', 'Endocardium (Ventricle)', 
                            'Fibroblasts (nppc)']


# In[153]:


#deep_injury_libraries = ['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9'] # All samples with > 50 nppc fibroblasts.
HR_Rv_7dpi_endo_winh = HR_Rv_7dpi_winh[HR_Rv_7dpi_winh.obs['Cell_type'].isin(all_endocardial_lineage)]#connected_endofibro_types)]
#HR_Rv_7dpi_endo_winh = HR_Rv_7dpi_endo_winh[HR_Rv_7dpi_endo_winh.obs['heart'].isin(deep_injury_libraries)]
scv.pp.remove_duplicate_cells(HR_Rv_7dpi_endo_winh)
sc.pp.filter_genes(HR_Rv_7dpi_endo_winh, min_cells=3)


# In[134]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_7dpi_endo_winh.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_7dpi_endo_winh.obs['percent_mito'] = np.sum(
    HR_Rv_7dpi_endo_winh[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_7dpi_endo_winh.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_7dpi_endo_winh, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_7dpi_endo_winh, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[135]:


sc.pp.normalize_per_cell(HR_Rv_7dpi_endo_winh, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_7dpi_endo_winh)


# In[136]:


sc.pp.regress_out(HR_Rv_7dpi_endo_winh, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[137]:


sc.pp.highly_variable_genes(HR_Rv_7dpi_endo_winh)
sc.tl.pca(HR_Rv_7dpi_endo_winh)
sc.external.pp.bbknn(HR_Rv_7dpi_endo_winh, batch_key='batch')
sc.tl.umap(HR_Rv_7dpi_endo_winh)


# In[138]:


sc.pl.umap(HR_Rv_7dpi_endo_winh, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_7dpi_endo_winh.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi deep injury endocardial niche')


# In[139]:


sc.pl.umap(HR_Rv_7dpi_endo_winh, color=['total_counts', 'n_genes_by_counts', 'percent_mito'], 
          title = '7dpi deep injury endocardial niche')


# In[140]:


sc.pl.umap(HR_Rv_7dpi_endo_winh, color='batch',
          title = '7dpi deep injury endocardial niche')


# In[141]:


sc.tl.leiden(HR_Rv_7dpi_endo_winh, resolution=2)


# In[142]:


sc.tl.paga(HR_Rv_7dpi_endo_winh, groups='leiden')


# In[143]:


sc.pl.paga(HR_Rv_7dpi_endo_winh, show=True, color = 'Cell_type', node_size_scale = 2)#, 
    #       labels = ['', '', '', '', '', 
   #                  '', '', '', '', '', 
   #                  '', '', '', '', '', '', ''], 
  #         node_size_scale = 2)#, save = 'scvelo_endo_winh_niche_7dpi_20n_leiden_diffmap.png')


# In[144]:


sc.tl.umap(HR_Rv_7dpi_endo_winh, init_pos='paga')


# In[145]:


sc.pl.umap(HR_Rv_7dpi_endo_winh, color = 'leiden')


# In[146]:


sc.pl.umap(HR_Rv_7dpi_endo_winh, color = 'batch')


# In[147]:


sc.pl.umap(HR_Rv_7dpi_endo_winh, color = 'Cell_type')


# In[148]:


sc.pl.umap(HR_Rv_7dpi_endo_winh, color = ['pecam1', 'cdh5', 'nppc'])


# In[149]:


scv.pp.moments(HR_Rv_7dpi_endo_winh)


# In[150]:


scv.tl.velocity(HR_Rv_7dpi_endo_winh, mode='stochastic')


# In[151]:


scv.tl.velocity_graph(HR_Rv_7dpi_endo_winh)


# In[152]:


scv.pl.velocity_embedding_stream(HR_Rv_7dpi_endo_winh, basis='umap', title = '', 
                                 color = 'Cell_type',
                                save = 'gene_velocity_local_endocardial_winhib_niche_7dpi_umap.png')


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


# # All timepoints in epicardial niche on same UMAP & clusters

# In[73]:


HR_Rv_epiconn_all = HR_Rv_filter[HR_Rv_filter.obs['Cell_type'].isin(local_connected_epifibro_types)]
all_genes_but_RFP = [name for name in HR_Rv_epiconn_all.var_names if not name == 'RFP']
HR_Rv_epiconn_all = HR_Rv_epiconn_all[:, all_genes_but_RFP]
sc.pp.filter_genes(HR_Rv_epiconn_all, min_cells=3)


# In[74]:


HR_Rv_epiconn_all


# In[75]:


# Find mito_genes in dataset.
mito_in_index = list(set(HR_Rv_epiconn_all.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
HR_Rv_epiconn_all.obs['percent_mito'] = np.sum(
    HR_Rv_epiconn_all[:, mito_in_index].X, axis=1) / np.sum(HR_Rv_epiconn_all.X, axis=1)
sc.pp.calculate_qc_metrics(HR_Rv_epiconn_all, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(HR_Rv_epiconn_all, ['total_counts', 'n_genes_by_counts', 'percent_mito'], 
             jitter=0.4, multi_panel=True, size = 0.1)


# In[76]:


sc.pp.normalize_per_cell(HR_Rv_epiconn_all, counts_per_cell_after=1e4)
sc.pp.log1p(HR_Rv_epiconn_all)


# In[77]:


sc.pp.regress_out(HR_Rv_epiconn_all, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[78]:


sc.pp.highly_variable_genes(HR_Rv_epiconn_all)
sc.tl.pca(HR_Rv_epiconn_all)
sc.external.pp.bbknn(HR_Rv_epiconn_all, batch_key='batch')
sc.tl.umap(HR_Rv_epiconn_all)


# In[79]:


sc.pl.umap(HR_Rv_epiconn_all, color='batch',
          title = 'Epicardial niche')


# In[83]:


sc.pl.umap(HR_Rv_epiconn_all, color='dpi',
          title = 'Epicardial niche')


# In[81]:


sc.pl.umap(HR_Rv_epiconn_all, color='Cell_type', palette = cell_type_colors.loc[HR_Rv_epiconn_all.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Epicardial niche')


# In[86]:


sc.pl.umap(HR_Rv_epiconn_all[HR_Rv_epiconn_all.obs['dpi'] == '0'], color='Cell_type', 
          title = 'Epicardial niche control')


# In[84]:


sc.pl.umap(HR_Rv_epiconn_all[HR_Rv_epiconn_all.obs['dpi'] == '3'], color='Cell_type', 
          title = 'Epicardial niche 3 dpi')


# In[87]:


sc.pl.umap(HR_Rv_epiconn_all[HR_Rv_epiconn_all.obs['inhib'] == 'IWR1'], color='Cell_type', 
          title = 'Epicardial niche Wnt inhibited')


# In[82]:


sc.pl.umap(HR_Rv_epiconn_all, color='inhib',
          title = 'Epicardial niche')


# In[55]:


sc.tl.leiden(HR_Rv_epiconn_all, resolution=2)


# In[56]:


sc.pl.umap(HR_Rv_epiconn_all, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[57]:


sc.tl.paga(HR_Rv_epiconn_all, groups='leiden')


# In[58]:


sc.pl.paga(HR_Rv_epiconn_all, show=True, node_size_scale = 2)#,
#           labels = ['', '', '', '', '', 
#                     '', '', '', '', '', 
#                     '', '', '', '', '',
#                    ''], 
#           save = 'scvelo_epicardial_niche_3dpi_20n_leiden_diffmap_recolor.png')


# In[59]:


sc.tl.umap(HR_Rv_epiconn_all, init_pos = 'paga')


# # Single dataset endocardium transitions

# In[79]:


deep_injury_libraries = ['Hr2a'] #['Hr1', 'Hr2a', 'Hr2b', 'Hr6v', 'Hr9'] # All samples with > 50 nppc fibroblasts.
Single_lib_7dpi_deep_endo = HR_Rv_7dpi[HR_Rv_7dpi.obs['Cell_type'].isin(connected_endofibro_types)]
Single_lib_7dpi_deep_endo = Single_lib_7dpi_deep_endo[Single_lib_7dpi_deep_endo.obs['heart'].isin(deep_injury_libraries)]
scv.pp.remove_duplicate_cells(Single_lib_7dpi_deep_endo)
#sc.pp.filter_genes(Single_lib_7dpi_deep_endo, min_cells=3)


# In[80]:


scv.pp.filter_and_normalize(Single_lib_7dpi_deep_endo, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(Single_lib_7dpi_deep_endo, n_pcs=30, n_neighbors=30)


# In[81]:


scv.tl.recover_dynamics(Single_lib_7dpi_deep_endo)


# In[82]:


scv.tl.velocity(Single_lib_7dpi_deep_endo, mode='dynamical')


# In[83]:


scv.tl.velocity_graph(Single_lib_7dpi_deep_endo)


# In[84]:


scv.tl.umap(Single_lib_7dpi_deep_endo)


# In[85]:


sc.tl.leiden(Single_lib_7dpi_deep_endo)


# In[77]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap',
                                 color='Cell_type', title = deep_injury_libraries,
                                 palette = cell_type_colors.loc[Single_lib_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[78]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap', color = 'leiden')


# In[40]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap',
                                 color='Cell_type', title = deep_injury_libraries,
                                 palette = cell_type_colors.loc[Single_lib_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[86]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap',
                                 color='Cell_type', title = deep_injury_libraries,
                                 palette = cell_type_colors.loc[Single_lib_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[87]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap', color = 'leiden', title = deep_injury_libraries)


# In[46]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap',
                                 color='Cell_type', title = deep_injury_libraries,
                                 palette = cell_type_colors.loc[Single_lib_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[58]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap',
                                 color='Cell_type', title = deep_injury_libraries,
                                 palette = cell_type_colors.loc[Single_lib_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[64]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap',
                                 color='Cell_type', title = deep_injury_libraries,
                                 palette = cell_type_colors.loc[Single_lib_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist())


# In[67]:


sc.tl.leiden(Single_lib_7dpi_deep_endo)


# In[69]:


scv.pl.velocity_embedding_stream(Single_lib_7dpi_deep_endo, basis='umap', color = 'leiden')


# In[22]:


# Find mito_genes in dataset.
mito_in_index = list(set(Single_lib_7dpi_deep_endo.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
Single_lib_7dpi_deep_endo.obs['percent_mito'] = np.sum(
    Single_lib_7dpi_deep_endo[:, mito_in_index].X, axis=1) / np.sum(Single_lib_7dpi_deep_endo.X, axis=1)
sc.pp.calculate_qc_metrics(Single_lib_7dpi_deep_endo, percent_top=None, log1p=True, inplace=True)
#sc.pl.violin(Single_lib_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
#             jitter=0.4, multi_panel=True, size = 0.1)


# In[23]:


sc.pp.normalize_per_cell(Single_lib_7dpi_deep_endo, counts_per_cell_after=1e4)
sc.pp.log1p(Single_lib_7dpi_deep_endo)


# In[24]:


sc.pp.regress_out(Single_lib_7dpi_deep_endo, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[ ]:


sc.pp.highly_variable_genes(Single_lib_7dpi_deep_endo)
sc.tl.pca(Single_lib_7dpi_deep_endo)
sc.external.pp.bbknn(Single_lib_7dpi_deep_endo, batch_key='batch')
sc.tl.umap(Single_lib_7dpi_deep_endo)


# In[27]:


sc.pl.umap(Single_lib_7dpi_deep_endo, color='Cell_type', palette = cell_type_colors.loc[Single_lib_7dpi_deep_endo.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = '7dpi deep injury endocardial niche')


# In[ ]:





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


# In[14]:


secretome = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Alliance_secretome_gene_names_noDRduplicates.scsv', sep = ';')
#secretome = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Alliance_secretome_gene_names.scsv', sep = '\t')
secretome = secretome.rename(columns={'DR_name': 'external_gene_name'})
secretome = secretome[secretome['external_gene_name'].notna()]
secretome = (secretome[['external_gene_name']]).drop_duplicates()


# In[17]:


secretome_old = pd.read_csv('~/Documents/Projects/heart_Bo/Data/Alliance_secretome_gene_names.scsv', sep = '\t')
secretome_old = secretome_old.rename(columns={'DR_name': 'external_gene_name'})
secretome_old = secretome_old[secretome_old['external_gene_name'].notna()]
secretome_old = (secretome_old[['external_gene_name']]).drop_duplicates()


# In[18]:


secretome_old


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

# In[15]:


def CountSecretome(adata, secretome):
    #secretome_ind = np.unique([adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index])
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


# In[16]:


def ExpressedSecretome(adata, secretome):
    secretome_ind = [adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index]
    #secretome_ind = np.unique([adata.var.index.get_loc(x) for x in np.array(secretome.external_gene_name) if x in adata.var.index])

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


# In[17]:


HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '0') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)
secretome_ctrl = CountSecretome(HR_Rv_ctrl, secretome)
secretome_averages_ctrl, secretome_z_averages_ctrl = ExpressedSecretome(HR_Rv_ctrl, secretome)


# In[21]:


#secretome_old_ctrl = CountSecretome(HR_Rv_ctrl, secretome_old)
#secretome_old_averages_ctrl, secretome_old_z_averages_ctrl = ExpressedSecretome(HR_Rv_ctrl, secretome_old)


# In[108]:


#secretome_z_averages_ctrl


# In[109]:


#secretome_old_z_averages_ctrl


# In[107]:


secretome_averages_ctrl_save = secretome_averages_ctrl[secretome_averages_ctrl.max(axis = 1) > 10]
secretome_averages_ctrl_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_averages_alliance_conversion_nodup_ctrl.csv')
secretome_z_averages_ctrl_save = secretome_z_averages_ctrl[secretome_z_averages_ctrl.max(axis = 1) > 2]
secretome_z_averages_ctrl_save.to_csv('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Data/Secretome_z_averages_alliance_conversion_nodup_ctrl.csv')


# In[27]:


#secretome_old_averages_ctrl_save = secretome_old_averages_ctrl[secretome_old_averages_ctrl.max(axis = 1) > 10]
#secretome_old_z_averages_ctrl_save = secretome_old_z_averages_ctrl[secretome_old_z_averages_ctrl.max(axis = 1) > 2]


# In[110]:


#secretome_z_averages_ctrl_save


# In[111]:


#secretome_old_z_averages_ctrl_save


# In[112]:


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


# # Multiplied transcript counts

# In[33]:


genes_of_interest = ['nrg1', 'fn1a', 'aldh1a2', 'cxcl12b']


# In[112]:


def CellTypeExpression(adata, genes):
    gene_ind = [adata.var.index.get_loc(x) for x in genes if x in adata.var.index]

    gene_slice = adata.X[:,gene_ind]

    ctype_sums = pd.DataFrame(index = [x for x in genes if x in adata.var.index])
    for ctype in adata.obs['Cell_type'].cat.categories:
        ctype_index = [x for x in range(0, len(adata.obs) - 1) if adata.obs.Cell_type[x] == ctype]
        if(len(ctype_index) == 0):
            # Add zeroes in case there are no cells of this type
            ctype_sum = pd.DataFrame({ctype : np.zeros(len(genes))},
                                      index = [x for x in np.array(genes) if x in adata.var.index])
        else:
            ctype_sum = pd.DataFrame({ctype : 10000 * np.squeeze(np.asarray(gene_slice[ctype_index, :].sum(axis = 0)))/len(adata.obs)},
                                      index = [x for x in np.array(genes) if x in adata.var.index])
        ctype_sums = ctype_sums.join(ctype_sum)
        
    return(ctype_sums)


# ## Control

# In[115]:


# Normalized transcripts
HR_Rv_ctrl = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '0') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_ctrl.var_names if not name == 'RFP']
HR_Rv_ctrl = HR_Rv_ctrl[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_ctrl, target_sum=1e4)

# Sum specific transcripts over cells from all celltypes
expression_sum_ctrl = CellTypeExpression(HR_Rv_ctrl, genes_of_interest)


# In[116]:


HR_Rv_3dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi.var_names if not name == 'RFP']
HR_Rv_3dpi = HR_Rv_3dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi, target_sum=1e4)

expression_sum_3dpi = CellTypeExpression(HR_Rv_3dpi, genes_of_interest)


# In[117]:


HR_Rv_7dpi = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] != 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi.var_names if not name == 'RFP']
HR_Rv_7dpi = HR_Rv_7dpi[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi, target_sum=1e4)

expression_sum_7dpi = CellTypeExpression(HR_Rv_7dpi, genes_of_interest)


# In[118]:


HR_Rv_3dpi_inh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '3') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_3dpi_inh.var_names if not name == 'RFP']
HR_Rv_3dpi_inh = HR_Rv_3dpi_inh[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_3dpi_inh, target_sum=1e4)

expression_sum_3dpi_inh = CellTypeExpression(HR_Rv_3dpi_inh, genes_of_interest)


# In[119]:


HR_Rv_7dpi_inh = HR_Rv_filter[(HR_Rv_filter.obs['dpi'] == '7') & (HR_Rv_filter.obs['inhib'] == 'IWR1')]
all_genes_but_RFP = [name for name in HR_Rv_7dpi_inh.var_names if not name == 'RFP']
HR_Rv_7dpi_inh = HR_Rv_7dpi_inh[:, all_genes_but_RFP]
sc.pp.normalize_total(HR_Rv_7dpi_inh, target_sum=1e4)

expression_sum_7dpi_inh = CellTypeExpression(HR_Rv_7dpi_inh, genes_of_interest)


# In[120]:


# Plot preparation
# Make sure all cell types are in all dataframes.
import matplotlib.pyplot as plt

expression_sum_ctrl_plot = expression_sum_ctrl.transpose()
expression_sum_ctrl_plot = expression_sum_ctrl_plot.drop(['Neuronal cells', 'Myelin cells'])
expression_sum_ctrl_plot = expression_sum_ctrl_plot.join(cell_type_colors, how = 'inner')

expression_sum_3dpi_plot = expression_sum_3dpi.transpose()
expression_sum_3dpi_plot = expression_sum_3dpi_plot.drop(['Neuronal cells', 'Myelin cells'])
expression_sum_3dpi_plot = expression_sum_3dpi_plot.join(cell_type_colors, how = 'inner')

expression_sum_7dpi_plot = expression_sum_7dpi.transpose()
expression_sum_7dpi_plot = expression_sum_7dpi_plot.drop(['Neuronal cells', 'Myelin cells'])
expression_sum_7dpi_plot = expression_sum_7dpi_plot.join(cell_type_colors, how = 'inner')

expression_sum_3dpi_inh_plot = expression_sum_3dpi_inh.transpose()
expression_sum_3dpi_inh_plot = expression_sum_3dpi_inh_plot.drop(['Neuronal cells', 'Myelin cells'])
expression_sum_3dpi_inh_plot = expression_sum_3dpi_inh_plot.join(cell_type_colors, how = 'inner')

expression_sum_7dpi_inh_plot = expression_sum_7dpi_inh.transpose()
expression_sum_7dpi_inh_plot = expression_sum_7dpi_inh_plot.drop(['Neuronal cells', 'Myelin cells'])
expression_sum_7dpi_inh_plot = expression_sum_7dpi_inh_plot.join(cell_type_colors, how = 'inner')

expression_array = [expression_sum_ctrl_plot, expression_sum_3dpi_plot, expression_sum_7dpi_plot, 
                   expression_sum_3dpi_inh_plot, expression_sum_7dpi_inh_plot]
annotation_array = ['ctrl', '3dpi', '7dpi', '3dpi IWR1', '7dpi IWR1']


# In[123]:


# Plot to-do
fig, axs = plt.subplots(len(genes_of_interest), len(expression_array), sharex='col', sharey='row', figsize=(25,10), gridspec_kw={'hspace': 0.1, 'wspace': 0.1})
for i in np.arange(len(expression_array)):
    plot_data = expression_array[i]
    condition = annotation_array[i]
    for gene_number in np.arange(len(genes_of_interest)):
        gene = genes_of_interest[gene_number]
        axs[gene_number, i].bar(x = np.arange(len(plot_data)),
                             height = plot_data[gene],
                             color = plot_data.color)
        axs[gene_number, i].set_xticks(np.arange(len(plot_data))) 
        axs[gene_number, i].set_xticklabels(plot_data.index, rotation=270)

for ax, col in zip(axs[0], annotation_array):
    ax.set_title(col)

for ax, row in zip(axs[:,0], genes_of_interest):
    ax.set_ylabel(row)

for ax in axs.flat:
    ax.label_outer()
    
plt.savefig('/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Images/Expressed_genes_times_cells_norm10k_cells.pdf', bbox_inches = 'tight', transparent = True)


# In[97]:


len(expression_sum_ctrl_plot), len(expression_sum_3dpi_plot)


# In[99]:


expression_sum_ctrl_plot.index


# In[111]:


expression_sum_3dpi_plot


# # Subclusters in endocardium

# In[83]:


adata = HR_Rv_filter[HR_Rv_filter.obs['inhib'] != 'IWR1']
all_genes_but_RFP = [name for name in adata.var_names if not name == 'RFP']
adata = adata[:, all_genes_but_RFP]


# In[84]:


HR_Rv_filter


# In[85]:


adata


# In[86]:


adata = adata[adata.obs['Cell_type'].isin(['Endocardium (Ventricle)', 'Endocardium (Atrium)', 'Endocardium (frzb)'])] #connected_epifibro_types)]
sc.pp.filter_genes(adata, min_cells=3)


# In[87]:


# Find mito_genes in dataset.
mito_in_index = list(set(adata.var.index.values) & set(mito_genes))
# Get mitochondrial percentages and read counts; plot violin plots
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_in_index].X, axis=1) / np.sum(adata.X, axis=1)
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, size = 0.1)


# In[88]:


sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)


# In[89]:


sc.pp.regress_out(adata, ['total_counts', 'n_genes_by_counts', 'percent_mito'])


# In[90]:


sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
#sc.pp.neighbors(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
sc.tl.umap(adata)


# In[91]:


adata


# In[92]:


sc.pl.umap(adata, color='Cell_type', palette = cell_type_colors.loc[adata.obs.Cell_type.cat.categories.tolist()].color.tolist(),
          title = 'Endocardial niche')


# In[93]:


sc.pl.umap(adata, color='batch')


# In[94]:


sc.tl.leiden(adata, resolution=2)


# In[95]:


sc.pl.umap(adata, color='leiden', legend_loc='on data', legend_fontsize='x-large')


# In[96]:


EV_I_H5 = ['H5_AAACGGGAGCCTTGAT', 'H5_AACTCTTAGTCGAGTG', 'H5_ACTATCTCACTACAGT', 'H5_ACTGAGTAGCGATGAC', 'H5_ATCCGAATCCTCAATT', 
           'H5_ATTATCCGTGAGCGAT', 'H5_CAACCAATCTTGCAAG', 'H5_CAAGAAAGTCACACGC', 'H5_CACTCCAAGACTGGGT', 'H5_CATCCACAGGAGTCTG', 
           'H5_CCCATACGTTGGTGGA', 'H5_CGATGGCAGAGTGACC', 'H5_CGCGGTACACGTCAGC', 'H5_CGGAGTCGTAGCGATG', 'H5_CGGTTAACACCGATAT', 
           'H5_CTCGTCAAGACGCAAC', 'H5_CTTAACTCACGAGGTA', 'H5_GAAATGACATCACGTA', 'H5_GACGCGTAGAAGGCCT', 'H5_GACGGCTCAGCCTTGG', 
           'H5_GCATACAAGTGGGATC', 'H5_GCGCGATGTCCAACTA', 'H5_GGCTGGTAGTGGAGAA', 'H5_GGCTGGTTCCAAGTAC', 'H5_GGTATTGCACCGGAAA', 
           'H5_GTACTCCCAGTTCATG', 'H5_GTATCTTCAAAGGCGT', 'H5_GTCTTCGCAAGCCGCT', 'H5_GTTCTCGAGCTCCCAG', 'H5_TACGGATTCACCGTAA', 
           'H5_TGGTTCCAGCGTGAAC', 'H5_TTATGCTCATAACCTG', 'H5_TTCGAAGAGGAGCGAG', 'H5_TTCTCAAGTAAGTAGT']
EV_I_Hr27 = ['Hr27_AACAAGATCGAGCCAC', 'Hr27_AAGGAATCAACACGTT', 'Hr27_AATGCCACAGTCTTCC', 'Hr27_AATTCCTCAAATCAAG', 'Hr27_ACCACAACATAATCCG', 
             'Hr27_ACGTACAGTCTCGACG', 'Hr27_ACGTTCCGTTCCACAA', 'Hr27_AGAGCAGGTTAGAGTA', 'Hr27_AGGTGTTTCAACGTGT', 'Hr27_AGTACCAGTAACGCGA', 
             'Hr27_AGTACCAGTGTGATGG', 'Hr27_ATACTTCAGATCGCTT', 'Hr27_ATCAGGTAGCTTAGTC', 'Hr27_ATGATCGTCAAGATAG', 'Hr27_ATGCATGGTATGCTTG', 
             'Hr27_ATGCCTCCAGAGTTCT', 'Hr27_ATGGAGGCAGTCAGAG', 'Hr27_ATTACCTTCGACGAGA', 'Hr27_CAAGCTATCATAGGCT', 'Hr27_CAATGACGTGCCGAAA', 
             'Hr27_CACGGGTTCCGGCTTT', 'Hr27_CATCCACTCGAGAGAC', 'Hr27_CCAAGCGTCCTCACGT', 'Hr27_CCAATGATCTCACTCG', 'Hr27_CCAATTTAGTGATGGC', 
             'Hr27_CCACAAACATGACAAA', 'Hr27_CCACGTTGTCGCGTCA', 'Hr27_CCACGTTTCTACGCGG', 'Hr27_CCGAACGAGAACGCGT', 'Hr27_CCTGCATTCCATAAGC', 
             'Hr27_CCTTGTGAGACCAGAC', 'Hr27_CGAGGCTAGACACACG', 'Hr27_CGAGTTATCGATTGAC', 'Hr27_CGGCAGTAGGTTCACT', 'Hr27_CGGGCATAGCCATTTG', 
             'Hr27_CGTGTCTGTTCGGGTC', 'Hr27_CTAAGTGAGCCGTAAG', 'Hr27_CTACATTAGGATTTAG', 'Hr27_CTATCCGTCTGTAACG', 'Hr27_CTCATCGCAGTGGGTA', 
             'Hr27_CTGCCATGTCCGGACT', 'Hr27_CTTAGGACACACTGGC', 'Hr27_CTTCTAAGTTGTGTTG', 'Hr27_CTTCTCTAGCGGGTTA', 'Hr27_GAAGCGAGTTGACTGT', 
             'Hr27_GACCAATCACAAGCAG', 'Hr27_GATCACACAAATGATG', 'Hr27_GATCATGAGGGTTAAT', 'Hr27_GATGGAGAGTGAATAC', 'Hr27_GATTGGTAGATACATG', 
             'Hr27_GCACGTGTCAGGGATG', 'Hr27_GCAGGCTAGGTTATAG', 'Hr27_GCATCTCCATATTCGG', 'Hr27_GCCCGAAAGTACCGGA', 'Hr27_GGAGGATAGTCTCGTA', 
             'Hr27_GGATCTAAGATCCAAA', 'Hr27_GGTAACTGTACGTGAG', 'Hr27_GGTGTCGTCTTCACGC', 'Hr27_GGTGTTATCTATGTGG', 'Hr27_GTACAACCAGCTTTCC', 
             'Hr27_GTAGAAACACACTTAG', 'Hr27_GTCACTCGTTTATGCG', 'Hr27_GTCAGCGTCCTGTACC', 'Hr27_GTCCACTAGTACCGGA', 'Hr27_GTCGAATAGCCATTCA', 
             'Hr27_GTCTAGAAGTCTGTAC', 'Hr27_TAAGTCGTCATTCACT', 'Hr27_TAGTGCAAGCACTCAT', 'Hr27_TCACATTGTATGATCC', 'Hr27_TCACTCGAGGTTACCT', 
             'Hr27_TCATGTTGTGAGGCAT', 'Hr27_TCATTGTAGACCAAGC', 'Hr27_TCCGGGAGTGGTCTCG', 'Hr27_TGACTCCAGGAGACCT', 'Hr27_TGAGCGCTCCTTCTGG', 
             'Hr27_TGGGATTCAGACAATA', 'Hr27_TGTGAGTTCAATCGGT', 'Hr27_TTACCGCAGCATCTTG', 'Hr27_TTCATGTCACCAAAGG', 'Hr27_TTCATTGTCGTACCTC', 
             'Hr27_TTGAGTGGTCAACCTA', 'Hr27_TTGGGATAGCCTGAAG', 'Hr27_TTGGTTTGTTATGGTC', 'Hr27_TTTCCTCAGCGCCTTG', 'Hr27_TTTGGTTTCGAGAAGC']
EV_II_H5 = ['H5_CAGAGAGCACTGTGTA', 'H5_CAGGTGCAGTCCTCCT', 'H5_CATGGCGGTCCGTCAG', 'H5_CCTTCCCGTCAGAAGC', 'H5_CGCTTCAGTGTGGCTC', 
            'H5_CTCGAAACACTGTGTA', 'H5_CTCGGAGCATCGATGT', 'H5_GACCTGGAGCTGCAAG', 'H5_GACCTGGCAGCCAGAA', 'H5_GACGGCTTCATGCATG', 
            'H5_GCAGTTAAGGAGTACC', 'H5_GCTCTGTCATTGCGGC', 'H5_GGGTCTGAGAGATGAG', 'H5_GTATCTTTCTCCAACC', 'H5_TACTCATAGAGGTAGA', 
            'H5_TACTTGTAGTTAAGTG']
EV_II_Hr27 = ['Hr27_AAAGTCCCAGTTCTAG', 'Hr27_AAAGTGAAGGCTCAAG', 'Hr27_AACCTTTCAGGTGACA', 'Hr27_AACCTTTGTAAGGTCG', 'Hr27_AACGTCAGTCTGTAGT', 
              'Hr27_AACTTCTAGCAGATAT', 'Hr27_AACTTCTAGGGTCTTT', 'Hr27_AACTTCTTCCCTTTGG', 'Hr27_AAGCGTTTCAGCTTGA', 'Hr27_AAGGAATTCGCTGCGA', 
              'Hr27_AAGGTAACACACCTGG', 'Hr27_AAGGTAAGTGCATACT', 'Hr27_AATGACCAGCTGAGTG', 'Hr27_AATGGAATCAGCATTG', 'Hr27_ACAAAGATCACCTCGT', 
              'Hr27_ACAAGCTGTATGTCAC', 'Hr27_ACATCGACATACAGAA', 'Hr27_ACATTTCCAAGACAAT', 'Hr27_ACATTTCTCGCTATTT', 'Hr27_ACCATTTGTCGCTTAA', 
              'Hr27_ACCGTTCGTCTCGGGT', 'Hr27_ACGTAACTCTAGACAC', 'Hr27_ACGTAGTTCCATCTCG', 'Hr27_AGAGAGCGTTAAGGAT', 'Hr27_AGAGCAGTCCGCCTAT', 
              'Hr27_AGATCCATCAACTTTC', 'Hr27_AGATCCATCCTACGAA', 'Hr27_AGATCGTCACTGGACC', 'Hr27_AGATCGTTCTCGTGAA', 'Hr27_AGCCACGGTCGAATGG', 
              'Hr27_AGCGCTGCACGGTCTG', 'Hr27_AGCGTATTCCACCTGT', 'Hr27_AGGCATTCAAATCAGA', 'Hr27_AGGCTGCTCTTACTGT', 'Hr27_AGGGCCTAGACTCTAC', 
              'Hr27_AGGTAGGGTTACCGTA', 'Hr27_AGGTCATGTCCAGCGT', 'Hr27_AGGTTGTCATTGTAGC', 'Hr27_AGTCAACCATGCGTGC', 'Hr27_AGTCTCCAGCGTTGTT', 
              'Hr27_ATCCCTGGTTGTTGTG', 'Hr27_ATCCGTCGTTAGGACG', 'Hr27_ATCGCCTCACCTTCGT', 'Hr27_ATGCATGCAGGTGGAT', 'Hr27_ATGCGATGTGCGCTCA', 
              'Hr27_ATTGTTCCAATGTTGC', 'Hr27_ATTTCACTCCTATTGT', 'Hr27_ATTTCTGTCAAGAATG', 'Hr27_CAACGATTCGGAATGG', 'Hr27_CAATCGAGTAGATCGG', 
              'Hr27_CACAACACACCTGAAT', 'Hr27_CACACAATCTCCGAGG', 'Hr27_CACAGGCTCGGACAAG', 'Hr27_CACGTTCCATGAATAG', 'Hr27_CACTGAATCCAAGCCG', 
              'Hr27_CAGGCCATCGTCAACA', 'Hr27_CAGTTAGGTTCAGGTT', 'Hr27_CAGTTAGTCCACGAAT', 'Hr27_CATACCCCACATAACC', 'Hr27_CATGCCTGTTGGCCGT', 
              'Hr27_CATTCATAGGCCACTC', 'Hr27_CATTGTTTCAGTGTTG', 'Hr27_CCAAGCGCACTGGCCA', 'Hr27_CCGGTGATCGCAGTGC', 'Hr27_CGAAGGATCGGACGTC', 
              'Hr27_CGAGGAAAGATTAGAC', 'Hr27_CGAGTTACATTGACCA', 'Hr27_CGATCGGAGTCCCAGC', 'Hr27_CGATGCGAGCATCCTA', 'Hr27_CGATGCGTCGTAACAC', 
              'Hr27_CGTAAGTCACGTATAC', 'Hr27_CGTGAATGTCATCGCG', 'Hr27_CGTGTCTCATAATGCC', 'Hr27_CGTTAGAAGTAGGAAG', 'Hr27_CGTTAGAGTCATCGGC', 
              'Hr27_CTAGACATCGTGGCGT', 'Hr27_CTATAGGGTAGGACTG', 'Hr27_CTATAGGTCCTTGACC', 'Hr27_CTATCCGAGTTCATCG', 'Hr27_CTATCTAAGAGCATAT', 
              'Hr27_CTCAAGACACTTGAGT', 'Hr27_CTCCGATCATGACAGG', 'Hr27_CTCCGATGTTTACACG', 'Hr27_CTGCGAGCAGACTCTA', 'Hr27_CTGTCGTTCAACTGAC', 
              'Hr27_CTGTGAACATCATCTT', 'Hr27_CTGTGGGAGAGCATAT', 'Hr27_CTTCAATTCTTGGTGA', 'Hr27_CTTCCTTTCCTGTTAT', 'Hr27_GAACGTTGTGAACTAA', 
              'Hr27_GAATAGAGTCCAAGAG', 'Hr27_GACACGCAGGAGATAG', 'Hr27_GACTTCCGTGAGACGT', 'Hr27_GAGGCCTTCAAGGTGG', 'Hr27_GAGTTGTTCCGCTAGG', 
              'Hr27_GATCGTATCTGGTGGC', 'Hr27_GCAGGCTTCAAATAGG', 'Hr27_GCTGCAGGTACCGCGT', 'Hr27_GCTTCACAGAGCAGAA', 'Hr27_GGAGGTAGTATTTCTC', 
              'Hr27_GGCAGTCAGGCCTGAA', 'Hr27_GGGATCCCACCTGCAG', 'Hr27_GGGATCCTCGTACCTC', 'Hr27_GGGTGAAGTCGTACTA', 'Hr27_GGTAATCCACTCCTGT', 
              'Hr27_GGTCACGGTTTGGGAG', 'Hr27_GGTGATTTCTTGCAGA', 'Hr27_GGTGGCTTCTCTGCCA', 'Hr27_GTAAGTCTCCGCGAGT', 'Hr27_GTAGTACCAATAACCC', 
              'Hr27_GTCAAACCATACCGTA', 'Hr27_GTGAGCCGTGCCCTTT', 'Hr27_GTGAGCCGTGGTTCTA', 'Hr27_GTGAGGATCAAACGTC', 'Hr27_GTGCACGTCGATCCCT', 
              'Hr27_GTGCAGCTCCGTGTCT', 'Hr27_GTGCTGGCACTGCGTG', 'Hr27_GTGGGAATCAAGGTGG', 'Hr27_GTGTCCTCAGCACCCA', 'Hr27_GTGTGGCAGCGTCGAA', 
              'Hr27_GTTGAACCAATAGTAG', 'Hr27_TAACTTCAGGTGCTTT', 'Hr27_TAAGCGTGTTATAGCC', 'Hr27_TACTGCCCATAGAAAC', 'Hr27_TAGACCATCTCCTGTG', 
              'Hr27_TATCTTGAGCCTGAAG', 'Hr27_TATTGCTTCCGTATGA', 'Hr27_TCACACCCAATTCACG', 'Hr27_TCAGGGCTCAGCGCAC', 'Hr27_TCAGGTATCTATCACT', 
              'Hr27_TCATATCAGCGAAACC', 'Hr27_TCATCCGCAACCAATC', 'Hr27_TCATGAGAGACGCATG', 'Hr27_TCATGCCCAATCACGT', 'Hr27_TCATGGATCTCGTCGT', 
              'Hr27_TCATTTGCAGCGATTT', 'Hr27_TCCCACAAGCCTCTTC', 'Hr27_TCGACGGTCGAGATGG', 'Hr27_TCGATTTTCACGGACC', 'Hr27_TCGCAGGAGACAGCGT', 
              'Hr27_TCGCTTGCATCCTTCG', 'Hr27_TCGGATACACATAACC', 'Hr27_TCGGGACCAAGCGATG', 'Hr27_TCGGGCAGTTTGCAGT', 'Hr27_TCGGTCTCACAGCCTG', 
              'Hr27_TCGTAGACACTCATAG', 'Hr27_TCGTAGAGTTGGTAGG', 'Hr27_TGACTCCAGGTATTGA', 'Hr27_TGAGCGCAGATAACGT', 'Hr27_TGAGGGAAGTGCCCGT', 
              'Hr27_TGATCTTCATAATCCG', 'Hr27_TGCAGATAGAAGCCAC', 'Hr27_TGCAGATAGAGCCATG', 'Hr27_TGCCGAGTCTGCTCTG', 'Hr27_TGCGACGAGGACATCG', 
              'Hr27_TGCGACGGTTGCACGC', 'Hr27_TGCTTCGCAACGTTAC', 'Hr27_TGCTTCGTCGAGAGAC', 'Hr27_TGGGAAGCACAATGTC', 'Hr27_TGTACAGCAAATACAG', 
              'Hr27_TGTACAGCACGAAGAC', 'Hr27_TTACGTTTCACTTATC', 'Hr27_TTCAATCGTCAGGTAG', 'Hr27_TTCATGTTCTAGCAAC', 'Hr27_TTCGATTGTCGAACGA', 
              'Hr27_TTCGATTGTCGATTTG', 'Hr27_TTCGATTTCTTCCACG', 'Hr27_TTCTTGACAAGGGTCA', 'Hr27_TTGATGGGTCTTCATT', 'Hr27_TTGCTGCTCCTACGGG', 
              'Hr27_TTGTGTTAGGTCGTAG', 'Hr27_TTGTGTTGTGGTTTGT', 'Hr27_TTTACTGGTTAGGGTG', 'Hr27_TTTATGCGTGCATCTA', 'Hr27_AAACGAATCAGCTGAT', 
              'Hr27_AAAGAACTCTAGGCAT', 'Hr27_AAAGGGCGTGCGACAA', 'Hr27_AACCAACTCTGCTTAT', 'Hr27_AACCCAAGTATACCTG', 'Hr27_AACCTTTAGACATAGT', 
              'Hr27_AACGTCATCTACGCAA', 'Hr27_AAGCCATTCTTTGATC', 'Hr27_AAGTTCGTCTTCCGTG', 'Hr27_AATCACGTCTCTATAC', 'Hr27_AATCGACTCACTCTTA', 
              'Hr27_ACGATGTAGCCATATC', 'Hr27_ACTGTCCAGGGAGTTC', 'Hr27_AGCGCCAGTCCGATCG', 'Hr27_AGGAAATCAACTAGAA', 'Hr27_AGGCATTTCTGAATGC', 
              'Hr27_AGTACTGAGCAACAAT', 'Hr27_AGTGCCGCAACTTCTT', 'Hr27_ATACTTCAGTCCCGAC', 'Hr27_ATCGTCCGTCGAGTTT', 'Hr27_ATTCCTATCAGGAGAC', 
              'Hr27_ATTCGTTTCCGTGTCT', 'Hr27_ATTTCTGAGTGTTGAA', 'Hr27_CACGGGTGTCATCGCG', 'Hr27_CACTGGGTCCGCTTAC', 'Hr27_CACTTCGCATAGGAGC', 
              'Hr27_CAGCAATAGCCACTCG', 'Hr27_CATCCACCAGGGATAC', 'Hr27_CATCCCACACGTAGAG', 'Hr27_CATCGCTCAAGAGTGC', 'Hr27_CATTCTAGTACGTGTT', 
              'Hr27_CCACGAGCACGCGCAT', 'Hr27_CCCATTGTCCCGAGGT', 'Hr27_CCGGGTACAGCTATTG', 'Hr27_CCTCTAGCAGGATCTT', 'Hr27_CGCAGGTAGATTACCC', 
              'Hr27_CGGACACTCTGCTTTA', 'Hr27_CGTCAAAGTATCTCGA', 'Hr27_CTACAGACAGCTAACT', 'Hr27_CTACTATAGAGAGGTA', 'Hr27_CTAGACATCTTACCAT', 
              'Hr27_CTCAGGGAGAACTTCC', 'Hr27_CTCATGCAGGCACAAC', 'Hr27_CTGAATGCAATCCTTT', 'Hr27_CTGAGGCTCATCTGTT', 'Hr27_CTGCATCGTTGGAGAC', 
              'Hr27_CTGCCTACAATCAGCT', 'Hr27_CTGTCGTCATATTCGG', 'Hr27_GACTCAATCGATGCTA', 'Hr27_GAGGGATAGACCATAA', 'Hr27_GATGCTACATCCTTCG', 
              'Hr27_GATTCGACATCGTGGC', 'Hr27_GCAACCGCACCTCGTT', 'Hr27_GCACGTGTCATGCCCT', 'Hr27_GCAGTTATCTGAGTCA', 'Hr27_GCATCGGCACGCTGTG', 
              'Hr27_GCATCTCTCCCTATTA', 'Hr27_GCGTGCACAGCTAACT', 'Hr27_GCTTTCGGTCTGCGCA', 'Hr27_GGAACCCTCCCAACTC', 'Hr27_GGAGCAAGTCAATGGG', 
              'Hr27_GGCGTCACAGCCGTTG', 'Hr27_GGGAAGTGTTCATCGA', 'Hr27_GGGCTCAGTAGAGTTA', 'Hr27_GGGTGAATCCTGGCTT', 'Hr27_GGTAACTCACCCAAGC', 
              'Hr27_GGTGATTAGCACCGAA', 'Hr27_GTACAGTAGATGGCGT', 'Hr27_GTAGCTAAGGAAAGAC', 'Hr27_GTCACTCCATAGCACT', 'Hr27_GTCCTCACATTCACAG', 
              'Hr27_GTCTCACCATTGCAAC', 'Hr27_GTGACGCGTCTACTGA', 'Hr27_GTGACGCTCATAGAGA', 'Hr27_GTGACGCTCTCTCTAA', 'Hr27_GTGATGTAGAAGAACG', 
              'Hr27_GTGCTTCAGTAAACTG', 'Hr27_GTTCATTAGTTCCGGC', 'Hr27_TAAGCCAGTGACTGAG', 'Hr27_TAATTCCCAAGAGTAT', 'Hr27_TACCGAACAATTTCTC', 
              'Hr27_TATCAGGGTGTAGGAC', 'Hr27_TATCTGTAGGAGGTTC', 'Hr27_TCACTCGTCTCCGATC', 'Hr27_TCCACCAGTTTAGACC', 'Hr27_TCGCTCACACAGAGAC', 
              'Hr27_TCGGGTGCACCTCAGG', 'Hr27_TCGTGCTAGGCATCTT', 'Hr27_TGAGGTTGTATGTCAC', 'Hr27_TGCAGGCGTTCCACAA', 'Hr27_TGCATGATCTCTGGTC', 
              'Hr27_TGCGATACAGAAACCG', 'Hr27_TGGAACTTCTACTGCC', 'Hr27_TGGCGTGAGATAACGT', 'Hr27_TGGGAAGGTGTTCAGT', 'Hr27_TGGTACACAACGACAG', 
              'Hr27_TGGTTAGGTCAACACT', 'Hr27_TGTAAGCAGAGTCCGA', 'Hr27_TGTAAGCTCATCGTAG', 'Hr27_TGTCAGAAGACCAGCA', 'Hr27_TGTGAGTAGATGCGAC', 
              'Hr27_TGTGAGTGTAGGCAAC', 'Hr27_TGTGAGTGTTGCATAC', 'Hr27_TTAGTCTGTTGGAGGT', 'Hr27_TTCCTTCCAGCAGACA', 'Hr27_TTTACCAGTCCACACG']


# In[38]:


#adata = sc.datasets.pbmc68k_reduced()
# create new categorical column called `selection`
#adata.obs['selection'] = pd.Categorical((adata.obs_vector('CD3G') > 2) & (adata.obs_vector('CD4') < 3))
# adjust colors
#adata.uns['selection_colors'] = ['blue', 'yellow']
#sc.pl.umap(adata, color='selection', add_outline=True, s=20)


# In[79]:


ininputput = 'H5_AAACGGGAGCCTTGAT'# adata.obs.index[0]
if ininputput in EV_I_H5:
    print("I_H5")
elif ininputput in EV_I_Hr27:
    print("I_Hr27")
elif ininputput in EV_II_H5:
    print("II_H5")
elif ininputput in EV_II_Hr27:
    print("II_Hr27")
else:
    print("No")


# In[ ]:


def cellchecker(x):
    if ininputput in EV_I_H5:
        print("I_H5")
    elif ininputput in EV_I_Hr27:
        print("I_Hr27")
    elif ininputput in EV_II_H5:
        print("II_H5")
    elif ininputput in EV_II_Hr27:
        print("II_Hr27")
    else:
        print("No")


# In[82]:


adata.obs['selection'] = pd.Categorical(np.where(adata.obs.index.isin(EV_I_H5), 'I_H5', adata.obs.index.isin(EV_II_H5), 'II_H5', 'No'))#adata.obs.index.isin(EV_I_H5) b)
adata.uns['selection_colors'] = ['red', 'black', 'grey']
sc.pl.umap(adata, color = ['selection'], groups = ['I_H5', 'II_H5'])


# In[46]:


sum(adata.obs.index.isin(EV_I_H5))


# In[ ]:




