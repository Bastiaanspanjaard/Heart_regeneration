# Heart regeneration
Data analysis scripts for single cell heart regeneration in zebrafish. External packages used: [Seurat](http://satijalab.org/seurat/) for single-cell transcriptome analysis, [AutoGeneS](https://autogenes.readthedocs.io/en/latest/) for bulk deconvolution, [LINNAEUS](https://github.com/Bastiaanspanjaard/LINNAEUS) for lineage tree building, [CellPhoneDB](https://github.com/Teichlab/cellphonedb) for ligand-receptor analysis, and [scanpy/PAGA](https://scanpy.readthedocs.io/en/stable/) and [scvelo](https://scvelo.readthedocs.io) for single-cell trajectory analysis.

# Data analysis workflow
1. __Determine cell types and marker gene expression.__ This is done in the script heartregen_clustering.R; our results were obtained using Seurat 3.2.0.
2. __Determine marker gene expression and cell type dynamics.__ This builds further on the cell type determination and is done in Heart_regen_Seurat_analysis.R.
3. __Determine cell type locations.__ In tomoDeconv.py, single-cell transcriptomes are combined with spatial RNA-seq to determine the atrium-ventricle location of cell types. 
4. __Calculate secretome.__ Cellular secretomes are calculated in HR_trajectories.py after an orthologue converion in Translate_secretome.R.
5. __Ligand-receptor analysis.__ After preparing data using Prep_CellPhoneDB.R, a ligand-receptor analysis is performed using CellPhoneDB through CPDB_alliance_ctrl_zoom.sh. The resulting interactions are analysed using CellPhoneDB_outcome_analysis.R.
6. __Build lineage trees.__ We use LINNAEUS to create single-cell lineage trees and here provide all necessary scripts (initial filtering scripts in ./Library_scar_filtering/, tree building in Iterative_tree_building.R using the package provided in ./collapsibleTree/, and the function library scar_helper_functions.R)
7. __Lineage tree analysis.__ Lineage trees are analyzed in Cell_type_tree_relations.R, together with the function library HR_library.R.
8. __Single-cell trajectory analysis.__ Single-cell trajectories based on PAGA and intron velocity are created in HR_trajectories.py.

# Data description
Transcriptome and scar sequencing data can be downloaded from GEO (accession numbers GSE159032 and GSE158919). In folder /Data in this repository, we provide a single-cell annotation file (final_celltypes.tsv) for cell type annotation, cell type color data files (color_scheme_seurat.rds and Cell_type_colors_2.csv - same underlying color scheme), a list of mitochondrial zebrafish genes (mito.genes.vs.txt), a list of zebrafish secretome genes (Alliance_secretome_gene_names_noDRduplicates.scsv), and the object containing all single-cell lineage trees built for this dataset (Tree_list_oneEndo.rds, scripts for visualization and analysis are in Cell_type_tree_relations.R).
