# CITE-seq
## code from data included in manuscript titled '''Metformin reduces the clonal fitness of Dnmt3aR878H hematopoietic stem and progenitor cells by reversing their aberrant metabolic and epigenetic state''' 
Mohsen Hosseini1, Veronique Voisin2#, Ali Chegini1,3#, Angelica Varesi1,4#, Severine Cathelin1#,
Dhanoop Manikoth Ayyathan1, Alex C.H. Liu1,3, Yitong Yang1,3, Vivian Wang1,3, Abdula Maher1,
Eric Grignano1, Julie A. Reisz5, Angelo D’Alessandro5, Kira Young6, Yiyan Wu3, Martina
Fiumara7,8, Samuele Ferrari7,8, Luigi Naldini7,8, Federico Gaiti1,3, Shraddha Pai2,9, Aaron D.
Schimmer1,3, Gary D. Bader2, John E. Dick1,4, Stephanie Z. Xie1, Jennifer J. Trowbridge6, and
Steven M. Chan1,3*

### R libraries
```Ruby
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
```
### RNA data (GEX) and CITseq data processing: pipeline adapted from Seurat (V4): link to tutorial, Seurat V4)
Description of the steps:
 - Step1: for each sample (2 MET and 2 VEH), a Seurat object was initialized with the raw non-normalized data. 
 - Step2: data exploration, quality control plots were performed by exploring the data distribution of number of count per gene, number of feature per cell and and the percentage of mitochondria per cell
 - Step3: Unwanted cells were removed by applying filters to retain cells with nFeature_RNA > 500 & nFeature_RNA <8000 & percent.mt < 15
- Step3: Normalization and Scaling of each individual sample using SCTransform; followed by dimension reduction and clustering
- Step4: Integration of the 4 samples together using an anchor method (RPCA) followed by dimension reduction and clustering and UMAP on all cells of all samples
- Step5: Detection of the  hematpoietic cell types on the integrated UMAP using AUCell.
- Step6: Quantification of the mutant and wild type cells in each sample and in each cell population using the tag quantification of CD45.1 and CD45.2 

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
