# CHIP-seq
## code from data included in manuscript :
***Metformin reduces the competitive advantage of Dnmt3a R878H HSPCs*** <br>
Mohsen Hosseini , Veronique Voisin , Ali Chegini , Angelica Varesi , Severine Cathelin ,
Dhanoop Manikoth Ayyathan , Alex C.H. Liu , Yitong Yang , Vivian Wang , Abdula Maher,
Eric Grignano , Julie A. Reisz , Angelo D’Alessandro , Kira Young , Yiyan Wu , Martina
Fiumara , Samuele Ferrari , Luigi Naldini , Federico Gaiti , Shraddha Pai , Grace Egan ,
Aaron D. Schimmer , Gary D. Bader , John E. Dick , Stephanie Z. Xie , Jennifer J.
Trowbridge , and Steven M. Chan 

## Calculate scores per genome regions
```Ruby
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
```

## Optional: select regions from cluster 1
```Ruby
data <- read.delim("~/nonaggregated_promoter2", header=FALSE)
sorted = read.delim("~/Dropbox (Bader Lab)/CSCTeam/CHAN/SC02/chipseq/sorted_aggregated_promoter.bed")
cluster1= sorted[ which(sorted$deepTools_group == "cluster_1"),]
dim(cluster1)
data_selected = data[which(data$V4 %in% cluster1$name) ,]
write.table(data_selected, "data_selected_aggregated.txt" ,sep="\t", quote=FALSE, row.names=FALSE, col.names=F)
gzip -c /Users/veroniquevoisin/Dropbox\ \(Bader\ Lab\)/CSCTeam/CHAN/SC02/chipseq/data_selected_aggregated.txt > data_selected_aggregated2.gz
```

## Draw the heatmap
```Ruby
plotHeatmap -m data_selected_aggregated2.gz \
-out aggregated_promoter_heatmap3v3.pdf \
-min -0.5 -max 2.5 \
--boxAroundHeatmaps no 
--heatmapHeight 10
--heatmapWidth 2
```

```Ruby
plotHeatmap -m data_selected_aggregated.gz \
-out aggregated_promoter_heatmap3v2.pdf \
-min -0.5 -max 2.5 \
--boxAroundHeatmaps no
```
