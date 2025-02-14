# CHIP-seq
## code from data included in manuscript :
***Metformin reduces the competitive advantage of Dnmt3a R878H HSPCs*** <br>
Mohsen Hosseini , Veronique Voisin , Ali Chegini , Angelica Varesi , Severine Cathelin ,
Dhanoop Manikoth Ayyathan , Alex C.H. Liu , Yitong Yang , Vivian Wang , Abdula Maher,
Eric Grignano , Julie A. Reisz , Angelo D’Alessandro , Kira Young , Yiyan Wu , Martina
Fiumara , Samuele Ferrari , Luigi Naldini , Federico Gaiti , Shraddha Pai , Grace Egan ,
Aaron D. Schimmer , Gary D. Bader , John E. Dick , Stephanie Z. Xie , Jennifer J.
Trowbridge , and Steven M. Chan 


### Merging the bigwig files
```Ruby
#!/bin/bash

#conda install -c bioconda ucsc-bigwigmerge
#conda install -c bioconda ucsc-bedGraphToBigWig
#download https://github.com/igvteam/igv/blob/master/genomes/sizes/mm10.chrom.sizes


##WT VEH
bigWigMerge /Users/veroniquevoisin/bigwig/1_0EGB_020IPMCC_WT-VEH-1_H3K27me3_mm10_i76_dmnorm_signal.bw /Users/veroniquevoisin/bigwig/2_0EGC_020IPMCC_WT-VEH-2_H3K27me3_mm10_i77_dmnorm_signal.bw  WT_VEH.bedgraph;

sort -k1,1 -k2,2n WT_VEH.bedgraph | \
awk 'BEGIN{OFS = "\t"}($5 = $4/3){print $1,$2,$3,$5}' > WT_VEH_adjusted.bedgraph;

bedGraphToBigWig WT_VEH_adjusted.bedgraph mm10.chrom.sizes WT_VEH_aggregated.bw;


##WT MET
bigWigMerge /Users/veroniquevoisin/bigwig/3_0EGD_020IPMCC_WT-MET-1_H3K27me3_mm10_i78_dmnorm_signal.bw /Users/veroniquevoisin/bigwig/4_0EGE_020IPMCC_WT-MET-2_H3K27me3_mm10_i80_dmnorm_signal.bw  WT_MET.bedgraph;

sort -k1,1 -k2,2n WT_MET.bedgraph | \
awk 'BEGIN{OFS = "\t"}($5 = $4/3){print $1,$2,$3,$5}' > WT_MET_adjusted.bedgraph;

bedGraphToBigWig WT_MET_adjusted.bedgraph mm10.chrom.sizes WT_MET_aggregated.bw;

##RH VEH
bigWigMerge /Users/veroniquevoisin/bigwig/5_0EFL_020IPMCC_RH-VEH-1_H3K27me3_mm10_i52_dmnorm_signal.bw /Users/veroniquevoisin/bigwig/6_0EGF_020IPMCC_RH-VEH-2_H3K27me3_mm10_i81_dmnorm_signal.bw RH_VEH.bedgraph;

sort -k1,1 -k2,2n RH_VEH.bedgraph | \
awk 'BEGIN{OFS = "\t"}($5 = $4/3){print $1,$2,$3,$5}' > RH_VEH_adjusted.bedgraph;

bedGraphToBigWig RH_VEH_adjusted.bedgraph mm10.chrom.sizes RH_VEH_aggregated.bw;

##RH MET
bigWigMerge /Users/veroniquevoisin/bigwig/8_0EGH_020IPMCC_RH-MET-2_H3K27me3_mm10_i84_dmnorm_signal.bw /Users/veroniquevoisin/bigwig/7_0EGG_020IPMCC_RH-MET-1_H3K27me3_mm10_i82_dmnorm_signal.bw RH_MET.bedgraph;

sort -k1,1 -k2,2n RH_MET.bedgraph | \
awk 'BEGIN{OFS = "\t"}($5 = $4/3){print $1,$2,$3,$5}' > RH_MET_adjusted.bedgraph;

bedGraphToBigWig RH_MET_adjusted.bedgraph mm10.chrom.sizes RH_MET_aggregated.bw;
```

## Calculate scores per genome regions
```Ruby
# Load Deeptools
module load deeptools/3.5.1

computeMatrix reference-point --referencePoint TSS \
-b 2000 -a 2000 \
-R allgenesmm10.bed \
-S 0_0EH4_020IPMCC_BM-Pooled_Input_mm10_i85_uniqnorm_signal.bw WT_VEH_aggregated.bw WT_MET_aggregated.bw RH_VEH_aggregated.bw RH_MET_aggregated.bw  \
--skipZeros \
-o aggregated_promoter3.gz \
-p 6 \
--binSize 20 \
--outFileSortedRegions aggregatedpromoter_center2.bed  \
--missingDataAsZero \
--numberOfProcessors 4 


```
## Draw the heatmap
```Ruby
plotHeatmap -m aggregated_promoter3.gz \
-out aggregated_promoter_heatmap3v2.pdf \
--outFileSortedRegions sorted_aggregated_promoter.bed \
--zMin 0 --zMax 0.8 \
--boxAroundHeatmaps no \
--kmeans 2
```

## Optional: select regions from one cluster
```Ruby
##unzip aggregated_promoter3.gz
data <- read.delim("~/aggregated_promoter3", header=FALSE)
sorted = read.delim("~/sorted_aggregated_promoter.bed")
cluster1= sorted[ which(sorted$deepTools_group == "cluster_1"),]
data_selected = data[which(data$V4 %in% cluster1$name) ,]
write.table(data_selected, "data_selected_aggregated.txt" ,sep="\t", quote=FALSE, row.names=FALSE, col.names=F)
gzip -c data_selected_aggregated.txt > data_selected_aggregated2.gz
```

## Re-draw the heatmap
```Ruby
plotHeatmap -m data_selected_aggregated2.gz \
-out aggregated_promoter_heatmap3v3.pdf \
-min -0.5 -max 2.5 \
--boxAroundHeatmaps no 
--heatmapHeight 10
--heatmapWidth 2
```


