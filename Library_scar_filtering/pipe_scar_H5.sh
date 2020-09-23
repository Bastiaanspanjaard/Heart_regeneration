#!/bin/bash
# Copy and concatenate fastq-files
cp /data/junker/users/bspanja/sc_scars_wt/2017_10X_7_CR/reads/outs/fastq_path/Heart/H5_scar*R* ./
gunzip H5_scar*R*.gz
cat H5_scar_part*R1* > H5_scar_R1.fastq
cat H5_scar_part*R2* > H5_scar_R2.fastq
rm H5_scar_part*
# Select Seurat barcodes as cellular barcodes
awk 'BEGIN{bc = 1}{if($2=="H5"){printf("%s\t%s\n", bc, $3);bc++}}' /local/Bastiaan/Projects/heart_Bo/Data/final_celltypes.tsv > H5_scar_barcodes.csv
# Create UMIBC-file
awk '(NR % 4 == 1) {print $1, $2} (NR % 4 == 3) {print $1} (NR % 2 == 0) {printf("%s%s\n", substr($1, 17, 10), substr($1, 1, 16))}' H5_scar_R1.fastq > H5_scar_UMIBC.fastq
# Map and extract scars
/local/Bastiaan/Scripts/scar_CIGAR_sc_10X_v2.pl -R1=H5_scar_UMIBC.fastq -R2=H5_scar_R2.fastq -op=H5_scar -t=1 -r=/local/gene_models/lintrace_hGFP_RFP_ERCC92.fa -bc=H5_scar_barcodes.csv -g=RFP -k=75 -mU=1 -l=75 -ps=GAGTTCAAGACCATCTACATGGCC
# Count how often every scar is sequenced, remove those that only occur once, filter the rest
tail -n +2 H5_scar_scars.txt | sort -k2,2 -k3,3 -k1,1 -k7,7 | uniq -c | cut -f1-4,7 > H5_scar_reads.txt
awk '$1>1' H5_scar_reads.txt > H5_scar_reads_over1.txt
source /local/Bastiaan/Scripts/scar_filter/bin/activate
/local/Bastiaan/bitbucket_scripts/sc_scar/scar_filter.py -i H5_scar_reads_over1.txt -o H5_scar_filtered_scars.csv
