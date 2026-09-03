#!/bin/bash
tfs=../data/work/SCENIC/mm_mgi_tfs.txt
#Downloaded from https://github.com/aertslab/pySCENIC/blob/master/resources/mm_mgi_tfs.txt
feather=../data/work/SCENIC/*.genes_vs_motifs.rankings.feather
# Two cisTarget reference files, "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather" and "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather", 
# were adopted for the analysis and downloaded from https://resources.aertslab.org/cistarget/databases/.
tbl=../data/work/SCENIC/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
#Downloaded from https://resources.aertslab.org/cistarget/motif2tf/.

input_loom1=../data/work/pySCENIC/HSC.loom
ls $tfs $feather $tbl $input_loom1

# GRN
pyscenic grn \
--num_workers 16 \
--output ../data/work/pySCENIC/HSC_adj.samplen.tsv \
--method grnboost2 \
$input_loom1 \
$tfs

# cisTarget
pyscenic ctx \
../data/work/pySCENIC/HSC_adj.samplen.tsv \
$feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom1 \
--mode dask_multiprocessing \
--output ../data/work/pySCENIC/HSC_regn.csv \
--num_workers 16 \
--mask_dropouts

# AUCell
pyscenic aucell \
$input_loom1 \
../data/work/pySCENIC/HSC_regn.csv \
--output ../data/work/pySCENIC/HSC_out_SCENICn.loom \
--num_workers 16

input_loom1=../data/work/pySCENIC/hep_50k.loom
ls $tfs $feather $tbl $input_loom1

# GRN
pyscenic grn \
--num_workers 16 \
--output ../data/work/pySCENIC/hep_adj.samplen.tsv \
--method grnboost2 \
$input_loom1 \
$tfs

# cisTarget
pyscenic ctx \
../data/work/pySCENIC/hep_adj.samplen.tsv \
$feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom1 \
--mode dask_multiprocessing \
--output ../data/work/pySCENIC/hep_regn.csv \
--num_workers 16 \
--mask_dropouts

# AUCell
pyscenic aucell \
$input_loom1 \
../data/work/pySCENIC/hep_regn.csv \
--output ../data/work/pySCENIC/hep_out_SCENICn.loom \
--num_workers 16

input_loom1=../data/work/pySCENIC/LSECs.loom
ls $tfs $feather $tbl $input_loom1

# GRN
pyscenic grn \
--num_workers 16 \
--output ../data/work/pySCENIC/LSEC_adj.samplen.tsv \
--method grnboost2 \
$input_loom1 \
$tfs

# cisTarget
pyscenic ctx \
../data/work/pySCENIC/LSEC_adj.samplen.tsv \
$feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom1 \
--mode dask_multiprocessing \
--output ../data/work/pySCENIC/LSEC_regn.csv \
--num_workers 16 \
--mask_dropouts

# AUCell
pyscenic aucell \
$input_loom1 \
../data/work/pySCENIC/LSEC_regn.csv \
--output ../data/work/pySCENIC/LSEC_out_SCENICn.loom \
--num_workers 16
