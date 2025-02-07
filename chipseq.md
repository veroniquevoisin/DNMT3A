#!/bin/bash
#BATCH -t 5:0:0
#SBATCH --mem=10G 
#SBATCH -J FLT3bam
#SBATCH -p himem 
#SBATCH -c 1
#SBATCH -N 1 
#SBATCH -o %x-%j.out

# Load Deeptools
module load deeptools/3.5.1

computeMatrix reference-point --referencePoint TSS \
-b 5000 -a 5000 \
-R allgenesmm10.bed \
-S 1_0EGB_020IPMCC_WT-VEH-1_H3K27me3_mm10_i76_dmnorm_signal.bw  3_0EGD_020IPMCC_WT-MET-1_H3K27me3_mm10_i78_dmnorm_si
gnal.bw 5_0EFL_020IPMCC_RH-VEH-1_H3K27me3_mm10_i52_dmnorm_signal.bw 7_0EGG_020IPMCC_RH-MET-1_H3K27me3_mm10_i82_dmnor
m_signal.bw         \
--skipZeros \
-o output.gz \
-p 6 \
--binSize 20 \
--outFileSortedRegions mergedintervals_center.bed
