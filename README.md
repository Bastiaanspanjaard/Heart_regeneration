# Heart regeneration
Data analysis scripts for single cell heart regeneration in zebrafish. Includes cell type determination using Seurat [Seurat](http://satijalab.org/seurat/), bulk RNA deconvolution using [AutoGeneS](https://autogenes.readthedocs.io/en/latest/), lineage tree building using [LINNAEUS](https://github.com/Bastiaanspanjaard/LINNAEUS), lineage tree analysis using scripts in this repository and trajectory analysis using [scanpy/PAGA](https://scanpy.readthedocs.io/en/stable/) and RNA velocity using [scvelo](https://scvelo.readthedocs.io).

# Tree-building workflow
A description of the full tree-building workflow, starting from sequenced scars and sequenced and (Cell Ranger) mapped mRNA data.

1. __Determine cell types using Seurat.__ OLD - Run the Larvae_Seurat_mt.R script in the R_scripts folder. This script combines the Cell Ranger-generated count tables of all sequenced larval 10X libraries, filters, clusters and determines differentially expressed genes for the cells sequenced. Uses [Seurat](http://satijalab.org/seurat/).
1. __Extract scars.__ Run the pipe_scar_H5.sh UNIX shell script in the scar_extraction folder. This script uses predetermined cellular barcodes from step 1) (see /Data/final_celltypes.tsv as example) and calls scar_CIGAR_sc_10X_v2.pl and scar_filter.py, present in the same folder. The script scar_CIGAR_sc_10X_v2.pl requires [bwa](http://bio-bwa.sourceforge.net) to be installed, and the (Python 2.7) script scar_filter.py assumes the existence of a virtual environment to run in, and packages optparse, pandas and distance to be installed in that environment.
3. __Filter scars within the library and determine creation probability of scars.__ Run the Preprocess_droplet_scar_allreads.R scripts and the Scar_comparison_between_fish.R script in the library_scar_filtering folder to combine scar libraries with mRNA libraries and perform further filtering steps.
4. __Run the tree-building algorithm.__ Run Iterative_tree_building.R to build and plot tree using cell type and filtered scar data.

# Data description
We provide a barcode annotation file /Data/final_celltypes.tsv that is used to extract scars on valid transcriptome barcodes. Other datasets can be downloaded from GEO.
