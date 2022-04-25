# Heart regeneration
Data analysis scripts for single cell heart regeneration in zebrafish. External packages used: [Seurat](http://satijalab.org/seurat/) for single-cell transcriptome analysis, [AutoGeneS](https://autogenes.readthedocs.io/en/latest/) for bulk deconvolution, [LINNAEUS](https://github.com/Bastiaanspanjaard/LINNAEUS) for lineage tree building, [CellPhoneDB](https://github.com/Teichlab/cellphonedb) for ligand-receptor analysis, and [scanpy/PAGA](https://scanpy.readthedocs.io/en/stable/) and [scvelo](https://scvelo.readthedocs.io) for single-cell trajectory analysis.

# Data analysis workflow
1. __Determine cell types and marker gene expression.__ (Seurat, add Bo's scripts)
2. __Determine cell type dynamics.__ (Seurat, add Bo's scripts)
3. __Determine cell type locations.__ (tomoDeconv.py)
4. __Calculate secretome.__ (Translate_secretome.R, HR_trajectories.py)
5. __Ligand-receptor analysis.__ (Prep_CellPhoneDB.R, CPDB_alliance_ctrl_zoom.sh, CellPhoneDB_outcome_analysis.R)
6. __Build lineage trees.__ (./Library_scar_filtering/ + Iterative_tree_building.R + ./collapsibleTree/ + scar_helper_functions.R, refer to LINNAEUS)
7. __Lineage tree analysis.__ (Cell_type_tree_relations.R + HR_library.R) Include technical validations.
8. __Single-cell trajectory analysis.__ (HR_trajectories.py)

# Data description
We provide a single-cell annotation file /Data/final_celltypes.tsv for cell type annotation and to extract scars on valid transcriptome barcodes. Other datasets can be downloaded from GEO. 
