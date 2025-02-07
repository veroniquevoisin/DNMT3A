# CITE-seq
## code from data included in manuscript :
***Metformin reduces the competitive advantage of Dnmt3a R878H HSPCs*** <br>
Mohsen Hosseini 1 , Veronique Voisin 2# , Ali Chegini 1,3# , Angelica Varesi 1,4# , Severine Cathelin 1# , Dhanoop Manikoth Ayyathan 1 , Alex C.H. Liu 1,3 , Yitong Yang 1,3 , Vivian Wang 1,3 , Abdula Maher 1 , Eric Grignano 1 , Julie A. Reisz 5 , Angelo D’Alessandro 5 , Kira Young 6 , Yiyan Wu 3 , Martina Fiumara 7,8 , Samuele Ferrari 7,8 , Luigi Naldini 7,8 , Federico Gaiti 1,3 , Shraddha Pai 2,9 , Grace Egan 1,10 , Aaron D. Schimmer 1,3 , Gary D. Bader 2 , John E. Dick 1,4 , Stephanie Z. Xie 1 , Jennifer J. Trowbridge 6 , and Steven M. Chan 1,3 *
Steven M. Chan1,3*

### R libraries
```Ruby
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
```
### RNA data (GEX) and CITseq data processing: 

Pipeline adapted from Seurat (V4): [tutorial](https://satijalab.org/seurat/articles/get_started.html)

Description of the steps:
 - Step1: for each sample (2 MET and 2 VEH), a Seurat object was initialized with the raw non-normalized data. 
 - Step2: data exploration, quality control plots were performed by exploring the data distribution of number of count per gene, number of feature per cell and and the percentage of mitochondria per cell
 - Step3: Unwanted cells were removed by applying filters to retain cells with nFeature_RNA > 500 & nFeature_RNA <8000 & percent.mt < 15
 - Step4: Normalization and Scaling of each individual sample using SCTransform; followed by dimension reduction and clustering
 - Step5: Integration of the 4 samples together using an anchor method (RPCA) followed by dimension reduction and clustering and UMAP on all cells of all samples
 - Step6: Detection of the  hematpoietic cell types on the integrated UMAP using AUCell.
 - Step7: Quantification of the mutant and wild type cells in each sample and in each cell population using the tag quantification of CD45.1 and CD45.2 

### Reading the Cell ranger output
```Ruby
data <- Read10X(data.dir = inputdir)
```
### Initializing the Seurat object with the RNA raw (non-normalized data).
```Ruby
sob <- CreateSeuratObject(counts = data$`Gene Expression`, project = "DNMT3a")
```
### Quality control metrics
```Ruby
sob[["percent.mt"]] <- PercentageFeatureSet(sob, pattern = "^MT-")
print(VlnPlot(sob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
plot1 <- FeatureScatter(sob, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```
### Adding CITEseq data to the RNA data Seurat object
```Ruby
cite <- CreateAssayObject(counts = data$`Antibody Capture`)
sob[["ADT"]] = cite
```
### Removing unwanted cells
```Ruby
DefaultAssay(sob) <- "RNA"
sob <- subset(sob, subset = nFeature_RNA > 500 & nFeature_RNA <8000 & percent.mt < 15)
```
### Normalizing and scaling the data using SCTransform
```Ruby
sob <- SCTransform(sob)
```
### Dimension reduction using the top 30 principal components
```Ruby
sob <- RunPCA(sob, features = VariableFeatures(object = sob))
```
### Clustering using the Louvain algorithm
```Ruby
sob <- FindNeighbors(sob, dims = 1:30)
sob <- FindClusters(sob, resolution = 0.5)
```
### Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
```Ruby
sob <- RunUMAP(sob, dims = 1:30)
```
### Integration of the 4 samples (2 MET and 2 VEH) using Seurat reciprocal PCA (‘RPCA’)[https://satijalab.org/seurat/articles/integration_rpca.html]
When determining anchors between any two datasets using RPCA, each dataset is projected into the others PCA space and constrain the anchors by the same mutual neighborhood requirement.
Anchors are identified using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and these anchors are used to integrate the two datasets together with IntegrateData().
```Ruby
sob12 = list(MET_1, MET_2,VFH_1, VFH_2)
features <- SelectIntegrationFeatures(object.list = sob12)
sob12 <- PrepSCTIntegration(sob12,anchor.features = features)
sob.anchors <- FindIntegrationAnchors(object.list = sob12, anchor.features = features,normalization.method = "SCT", dims = 1:50, reductio
n = "rpca", k.anchor = 5)
sob.combined <- IntegrateData(anchorset = sob.anchors, normalization.method = "SCT", dims = 1:50)
DefaultAssay(sob.combined) <- "integrated"
```
### Running the standard workflow for visualization and clustering
```Ruby
#sob.combined <- ScaleData(sob.combined, verbose = FALSE)
#sob.combined <- SCTransform(sob.combined, verbose = FALSE)
sob.combined <- RunPCA(sob.combined, npcs = 30, verbose = FALSE)
sob.combined <- RunUMAP(sob.combined, reduction = "pca", dims = 1:30)
sob.combined <- FindNeighbors(sob.combined, reduction = "pca", dims = 1:30)
sob.combined <- FindClusters(sob.combined, resolution = 0.5)
```
### CITE-seq: normalization of CD45.1 and CD45.2 and ratio calculation
```Ruby
sob.combined = readRDS(sob.combined)
DefaultAssay(sob.combined) <- "ADT"
sob.combined <- subset(x = sob.combined, subset = `B0178-CD45-1-TotalSeqB` <600 & `B0157-CD45-2-TotalSeqB` <600)
sob.combined <- NormalizeData(sob.combined, normalization.method = "LogNormalize",margin = 2)
sob.combined <- ScaleData(sob.combined, assay="ADT", display.progress = FALSE)
result = FetchData(object = sob.combined, vars = c("B0157-CD45-2-TotalSeqB", "B0178-CD45-1-TotalSeqB"))
result$ratio = result[,1] - result[,2]
```

###Cell type annotation (code: Federico Gaiti)
```Ruby
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(AUCell)
library(RColorBrewer)


#input data
sob.combined = readRDS("sob.combined.rds")

inp_data = GetAssayData(sob.combined, assay = 'RNA', slot = 'counts') 

##genesets

genesets = read.csv("genesets.csv", stringsAsFactors=FALSE) ##genesets from Landau's paper containing genes specific to each cell type

mylist = list()
for (i in c(1:ncol(genesets))){
  print(i)
  genesets2 = unlist(genesets[,i])
  genelist2 = rownames(inp_data)[which((rownames(inp_data) %in% genesets2))]
  mylist[[i]] =  genelist2
  print( genelist2)
}

names(mylist) = colnames(genesets)
length(mylist)

##parameters (running in batches to limit memory consumption)
num_batches = 100
num_cells <- ncol(inp_data)
batch_size <- ceiling(num_cells/num_batches)
score_mat <- c()
print('Running AUCell scoring')


##running AUCell
for (i in 1:num_batches) {
      print(paste('batch', i, Sys.time()))
      ind1 <- (i-1)*batch_size + 1
      ind2 <- i*batch_size
      if (ind2 > num_cells) {
        ind2 <- num_cells
      }
      gene_rankings <- AUCell::AUCell_buildRankings(inp_data[,ind1:ind2], plotStats = FALSE)
      score_mat_i <- AUCell::AUCell_calcAUC(geneSets = mylist, rankings = gene_rankings, aucMaxRank = ceiling(0.25 * nrow(gene_rankings)) )
      score_mat_i <- t(SummarizedExperiment::assay(score_mat_i, 'AUC'))
      score_mat <- rbind(score_mat, score_mat_i)
      gc(full = TRUE, verbose = TRUE)
}

##plots
FeaturePlot(object = sob.combined_AUC, features = names(mylist)[i],order=TRUE, raster=TRUE, cols = c("blue", "red"),min.cutoff = "q10
", max.cutoff = "q90"))

FeaturePlot(object = sob.combined_AUC, features = names(mylist)[i],order=TRUE,split.by = "protocol", raster=TRUE, cols = c("blue", "red"),min.cutoff = "q10", max.cutoff = "q90"))
```



###  UMAP (code: Alex C.H. Liu)
```Ruby
d %>% 
  group_by(CellType2) %>% 
  dplyr::summarise(xpos = median(UMAP1),
                   ypos = median(UMAP2)) -> umap.lab.pos
umap.lab.pos$lab <- umap.lab.pos$CellType2 %>% gsub("_|Cd3d\\.", "", .)

p = umap.gene.expression.celltype2 + scale_color_viridis(
  discrete = TRUE, 
  option="H",
  alpha = 1,
  begin = 0,
  end = 1,
  direction = 1,) 


```
