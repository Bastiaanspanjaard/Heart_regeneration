# Heart regeneration
Data analysis scripts for single cell heart regeneration in zebrafish. Includes cell type determination using Seurat [Seurat](http://satijalab.org/seurat/), bulk RNA deconvolution using [AutoGeneS](https://autogenes.readthedocs.io/en/latest/), lineage tree building using [LINNAEUS](https://github.com/Bastiaanspanjaard/LINNAEUS), lineage tree analysis using scripts in this repository and trajectory analysis using [scanpy/PAGA](https://scanpy.readthedocs.io/en/stable/) and RNA velocity using [scvelo](https://scvelo.readthedocs.io).

# Data analysis workflow
1. __Determine cell types and marker gene expression.__ (Seurat, Bo's scripts)
2. __Determine cell type dynamics.__ (Seurat, Bo's scripts)
3. __Determine cell type locations.__ (python, AutoGeneS)
4. __Calculate secretome.__ (Seurat, my script)
5. __Ligand-receptor analysis.__ (CellPhoneDB, analysis script)
6. __Build lineage trees.__ (Existing code, refer to LINNAEUS)
7. __Lineage tree analysis.__ (R, my script) Include technical validations.
8. __Single-cell trajectory analysis.__ (python, my script)

# Data description
We provide a single-cell annotation file /Data/final_celltypes.tsv for cell type annotation and to extract scars on valid transcriptome barcodes. Other datasets can be downloaded from GEO.
