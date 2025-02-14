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
data <- read.delim("~/nonaggregated_promoter2", header=FALSE)
sorted = read.delim("~/sorted_aggregated_promoter.bed")
cluster1= sorted[ which(sorted$deepTools_group == "cluster_1"),]
data_selected = data[which(data$V4 %in% cluster1$name) ,]
write.table(data_selected, "data_selected_aggregated.txt" ,sep="\t", quote=FALSE, row.names=FALSE, col.names=F)
gzip -c data_selected_aggregated.txt > data_selected_aggregated2.gz
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
