# RRBS
## code [Ali Chegini] from data included in manuscript :
***Metformin reduces the competitive advantage of Dnmt3a R878H HSPCs*** <br>
Mohsen Hosseini , Veronique Voisin , Ali Chegini , Angelica Varesi , Severine Cathelin ,
Dhanoop Manikoth Ayyathan , Alex C.H. Liu , Yitong Yang , Vivian Wang , Abdula Maher,
Eric Grignano , Julie A. Reisz , Angelo D’Alessandro , Kira Young , Yiyan Wu , Martina
Fiumara , Samuele Ferrari , Luigi Naldini , Federico Gaiti , Shraddha Pai , Grace Egan ,
Aaron D. Schimmer , Gary D. Bader , John E. Dick , Stephanie Z. Xie , Jennifer J.
Trowbridge , and Steven M. Chan 


### Loading libraries
```Ruby
library(methylKit)
library(ggvenn)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(viridis)
library(BiocParallel)
register(MulticoreParam(20))

source('Settings.R')
source(file.path(METHYLOMEANALYSIS_PATH, 'local_settings.R'))
```


### Functions
```Ruby


MyAnnotateRegions <- function (regions, TxDb, annoDb, tssRegion = c(-3000, 3000))
{
  ### The original function is DMRichR::annotateRegions
  ### I modified it to be able to change tssRegion window
  ### The original function only consider 3kb up and down stream
  
  genome <- TxDb %>% GenomeInfoDb::genome() %>% unique()
  print(glue::glue("Annotating {tidyRegions} regions from {genome} with gene symbols", 
    tidyRegions = length(regions)))
  if (is(TxDb, "EnsDb")) {
    genome <- dplyr::case_when(GenomeInfoDb::genome(TxDb) == 
       "GRCh38" ~ "hg38", GenomeInfoDb::genome(TxDb) == 
       "GRCm38" ~ "mm10", GenomeInfoDb::genome(TxDb) == 
       "Mmul_10" ~ "rheMac10", GenomeInfoDb::genome(TxDb) == 
       "Mmul_8.0.1" ~ "rheMac8", GenomeInfoDb::genome(TxDb) == 
       "Rnor_6.0" ~ "rn6", GenomeInfoDb::genome(TxDb) == 
       "GRCz11" ~ "danRer11", GenomeInfoDb::genome(TxDb) == 
       "GRCg6a" ~ "galGal6", GenomeInfoDb::genome(TxDb) == 
       "ARS-UCD1.2" ~ "bosTau9", GenomeInfoDb::genome(TxDb) == 
       "BDGP6.28" ~ "dm6", GenomeInfoDb::genome(TxDb) == 
       "Sscrofa11.1" ~ "susScr11", GenomeInfoDb::genome(TxDb) == 
       "CanFam3.1" ~ "canFam3") %>% unique()
  }
  CpGs <- DMRichR::getCpGs(genome)
  regionsCpG <- regions %>% plyranges::join_overlap_left(CpGs %>% 
    plyranges::filter(type == "islands") %>% plyranges::select(CpG.Island = type)) %>% 
    unique() %>% plyranges::join_overlap_left(CpGs %>% plyranges::filter(type == 
    "shores") %>% plyranges::select(CpG.Shore = type)) %>% 
    unique() %>% plyranges::join_overlap_left(CpGs %>% plyranges::filter(type == 
    "shelves") %>% plyranges::select(CpG.Shelf = type)) %>% 
    unique() %>% plyranges::join_overlap_left(CpGs %>% plyranges::filter(type == 
    "inter") %>% plyranges::select(Open.Sea = type)) %>% 
    unique() %>% plyranges::mutate(CpG.Island = dplyr::case_when(CpG.Island == 
    "islands" ~ "Yes", TRUE ~ "No"), CpG.Shore = dplyr::case_when(CpG.Shore == 
    "shores" ~ "Yes", TRUE ~ "No"), CpG.Shelf = dplyr::case_when(CpG.Shelf == 
    "shelves" ~ "Yes", TRUE ~ "No"), Open.Sea = dplyr::case_when(Open.Sea == 
    "inter" ~ "Yes", TRUE ~ "No"))
  if (is(TxDb, "EnsDb")) {
    GenomeInfoDb::seqlevelsStyle(regionsCpG) <- "Ensembl"
  }
  
  
  regionsCpG %>% ChIPseeker::annotatePeak(tssRegion = tssRegion, TxDb = TxDb, annoDb = annoDb, 
    overlap = "all", verbose = FALSE) %>% dplyr::as_tibble() %>% 
    dplyr::select(-tidyselect::any_of(c("strand", "index.start", 
      "index.end", "index.width", "area", "geneChr", "geneStart", 
      "geneEnd", "geneLength", "geneStrand", "transcriptId", 
      "transcriptBiotype", "ENTREZID"))) %>% dplyr::mutate(annotation = gsub(" \\(.*", 
    "", annotation)) %>% dplyr::rename_with(~dplyr::case_when(. == 
    "seqnames" ~ "chr", . == "L" ~ "CpGs", . == "beta" ~ 
    "betaCoefficient", . == "stat" ~ "statistic", . == "pval" ~ 
    "p.value", . == "qval" ~ "q.value", . == "SYMBOL" ~ 
    "geneSymbol", . == "GENENAME" ~ "gene", TRUE ~ .)) %>% 
    return()
}

MakeFileName <- function(
    region_,
    min_cov_,
    min_per_group_,
    comparison_,
    dataset_name)
{
  
  dataset_name <- glue::glue('{dataset_name}_Comp:{comparison_}_Region:{region_}_MinCov:{min_cov_}_MinPerGroup:{min_per_group_}')
  
  return(dataset_name)
}

AnnoRegionFilter <- list(
  'ALL' = list(
    'cond' = function(x) rep(T, nrow(x)),
    'slug' = 'ALL'),
  'PRO' = list(
    'cond' = function(x) (x$annotation %in% c('Promoter')),
    'slug' = 'PRO'),
  'ISL' = list(
    'cond' = function(x) (x$CpG.Island %in% c('Yes')), 
    'slug' = 'ISL'),
  'Shelf' = list(
    'cond' = function(x) (x$CpG.Shelf %in% c('Yes')), 
    'slug' = 'Shelf'),
  'Shore' = list(
    'cond' = function(x) (x$CpG.Shore %in% c('Yes')), 
    'slug' = 'Shore'),
  'Open.Sea' = list(
    'cond' = function(x) (x$Open.Sea == 'Yes'), 
    'slug' = 'Open.Sea')
)

MyGalaxyPlot <- function(
  data_set_,
  region_cond,
  col_name_map,
  comparison_slug_map,
  lrf,
  results_slug = 'Galaxy'
  )
{
  
  myplot <- function(dtt, condition_slug = 'All'){
    
    dtt <- dtt[dtt$V1 <=0, ]
    maxx <- max(abs(dtt))

    title_ <- condition_slug
    y_label <- str_replace_all(comparison_slug_map[col_name_map['V1']], '\n', '    ')
    x_label <- str_replace_all(comparison_slug_map[col_name_map['V2']], '\n', '    ')
    
    
    p1 <- ggscatter(
      data = dtt,
      x = 'V2',
      y = 'V1',
      size = 1
    ) + 
      # geom_hline(yintercept=0, linetype="dashed", color = "red") +
      geom_vline(xintercept=0, linetype="dashed", color = "red") +
      scale_color_viridis(option="viridis") + theme(legend.position="none") +
      scale_y_continuous(n.breaks = 10, 
                         limits = c(0, -maxx),
                         trans = "reverse"
      ) +
      scale_x_continuous(n.breaks = 10, limits = c(-maxx, maxx)) +
      xlab(x_label) + 
      ylab(y_label) +
      ggtitle(title_)
    
    file_name <- glue::glue('{results_slug}_{condition_slug}.png')
    file_dir <- file.path(lrf, file_name)
    ggsave(file_dir, p1, units="in", width=7, height=5, dpi=300)
    
    return(p1)
  }
  
  
  
  
  ## results dir preparation
  lrf <- file.path(
    lrf,
    'Plots'
  )
  if(!(dir.exists(lrf))){
    dir.create(lrf)
  }
  
  
  anno_cols_temp <- c(
    "CpG.Island",
    "CpG.Shelf", "Open.Sea", "annotation",
    "geneId", "distanceToTSS", "geneSymbol",
    "gene")

  data_set_ <- data_set_[,c('V1', 'V2', anno_cols_temp)]
  
  
  cond_slug <- glue::glue('{region_cond$slug}')
  data_set_ <- data_set_[region_cond$cond(data_set_),]
  
  p <- myplot(
    dtt = data_set_[,c('V1','V2')],
    condition_slug = cond_slug
  )
  
  return(p)
}


LoadDMRs <- function(
  file_name,
  lrf)
{
  
  lrf <- file.path(
    lrf, 
    'Comparisons',
    file_name)
  
  dtt <- read_tsv(lrf) %>%
    as.data.frame(dtt)
  
  cov_columns <- grepl(
    pattern = 'cov.',
    x = colnames(dtt),
    fixed = T)
  cov_ <- rowSums(dtt[,cov_columns], na.rm = T)
  dtt$cov.Total<- cov_
  
  
  return(dtt)
}

StoreDMRResults <- function(
    region_meth,
    region_diff,
    comparison_,
    comparison_dir,
    organism_,
    txdb_,
    results_slug = 'Annotated'){
  
  
  
  ### Preparing coverage data
  dtt <- methylKit::getData(region_meth)

  cov_index <- paste('cov', region_meth@sample.ids, sep = '.')
  names(cov_index) <- region_meth@coverage.index
  
  Cs_index <- paste('numCs', region_meth@sample.ids, sep = '.')
  names(Cs_index) <- region_meth@numCs.index
  
  Ts_index <- paste('numTs', region_meth@sample.ids, sep = '.')
  names(Ts_index) <- region_meth@numTs.index
  
  anno_col <- c(cov_index,Cs_index,Ts_index)
  colnames(dtt)[as.numeric(names(anno_col))] <- anno_col
  
  ### Preparing comparison data
  region_diff <- methylKit::getData(region_diff)
  region_diff$dataset_name <- results_slug
  region_diff$comparison_name <- comparison_
  
  ### Merging coverage and comparison
  merge_columns <- c("chr", "start", "end", "strand")
  
  dtt <- merge(
    x = dtt,
    y = region_diff,
    by = merge_columns
  )
  
  region_diff_granges <- as(dtt,"GRanges")
  seqlevelsStyle(region_diff_granges) <- "UCSC"
  
  region_anno <- region_diff_granges %>% MyAnnotateRegions(
    TxDb = txdb_, 
    annoDb = organism_, 
    tssRegion = tssRegion
    ) %>%
    as.data.frame()

  file_name <- glue::glue('{results_slug}.tsv')
  write_tsv(region_anno, file.path(comparison_dir,file_name))

  return(1)
}

Filter4Comparison <- function(
    myobj,
    comparison_
){
  groups_ <- strsplit(
    x = as.character(comparison_),
    split = '_vs_',
    fixed = T
  )[[1]]
  
  sample_ids <- names(filtered_obj@treatment)
  treated_slug <- groups_[1]
  control_slug <- groups_[2]
  cond1 <- grepl(pattern = control_slug, x = sample_ids, fixed = T)
  cond2 <- grepl(pattern = treated_slug, x = sample_ids, fixed = T)
  sub_sample_ids <- sample_ids[cond1 | cond2]
  names(sub_sample_ids) <- sub_sample_ids
  sub_treatments <- rep(0, length(sub_sample_ids))
  names(sub_treatments) <- sub_sample_ids
  sub_treatments[grepl(treated_slug, sub_sample_ids, fixed = T)] <- 1
  sub_myobj  <- reorganize(
    myobj,
    sample.ids=sub_sample_ids,
    treatment=sub_treatments
  )
  
  return(sub_myobj)
  
}
```
### Loading Data
```Ruby
### Load hyper parameters
N_CPU = 1

tssRegion <- c(-1000,150)

ComparisonSets <- c(
  'RH_MET_vs_RH_VEH', 
  'RH_VEH_vs_WT_VEH'
)
ComparisonDirMap <- c(
  RH_MET_vs_RH_VEH = +1,
  RH_VEH_vs_WT_VEH = -1
)
ComparisonSlugMap <- c(
  RH_MET_vs_RH_VEH = 'R878/+ Metformin\nvs\nR878/+ Ctrl',
  RH_VEH_vs_WT_VEH = 'R878/+ Ctrl\nvs\n+/+ Ctrl')
ComparisonSlugMap <- ComparisonSlugMap[ComparisonSets]

PvalueCutoff <- 0.01
CX <- 'CpG'

WindowSize <- 1000
RegionSlug <- '1kbBin'

RegionConds <- c('ALL', 'ISL', 'PRO')
DeduplicationApp <- 'COV'

MinCov <- 3 
MinPerGroup <- 0L
Strands <- c('*')
OverLap <- 0
### Loading Bismark cytosine reports
files_list <- list.files(APP_DB_PATH)
cond1 <- grepl('.txt.gz', files_list)
files_list <- files_list[cond1]

treated_slug <- 'MET'
control_slug <- 'VEH'

dataset_name_temp <- 'BismarkCpG'

cond1 <- grepl(pattern = control_slug, x = files_list, fixed = T)
cond2 <- grepl(pattern = treated_slug, x = files_list, fixed = T)
files_list <- files_list[cond1 | cond2]
files_dir <- sapply(files_list, function(x) file.path(APP_DB_PATH,x))
sample_ids <- str_replace_all(files_list, '_CX_report.txt.gz' , '')
names(sample_ids) <- names(files_dir)

treatments <- rep(0, length(sample_ids))
treatments[grepl(pattern = treated_slug, x = sample_ids, fixed = T)] <- 1
names(treatments) <- sample_ids

sample_ids_bk <- sample_ids
treatments_bk <- treatments

if(!('myobj_bk' %in% ls())){
  ## be careful with this! you are gonna need lots of memory!
  myobj_bk <- methRead(
    location = as.list(files_dir),
    sample.id=as.list(sample_ids),
    assembly="mm10",
    treatment=treatments,
    context=CX,
    mincov = 3,
    pipeline = 'bismarkCytosineReport',
    header = F
    # dbtype = 'tabix',
    # dbdir = dirname(APP_DB_PATH)
  )
}

```

### Excluding low quality samples.
```Ruby
### removing low quality samples
TotalBases <- do.call(rbind,lapply(myobj_bk,dim)) %>%
  as.data.frame()
colnames(TotalBases) <- c('NumBases', 'NumCols')
rownames(TotalBases) <- sample_ids
TotalBases$NumBases <- TotalBases$NumBases / mean(TotalBases$NumBases) * 100
cond1 <- TotalBases$NumBases < 50
print(TotalBases[cond1,])
# RH_MET4 and WT_VEH1 should be excluded due to insufficient read counts.

excluded_samples <- c('RH_MET4', 'WT_VEH1')
cond1 <- which(!(sample_ids %in% excluded_samples))
filtered_sample_ids <- sample_ids[cond1]
filtered_treatments <- treatments[cond1]
filtered_obj <- reorganize(
  methylObj = myobj_bk,
  sample.ids = filtered_sample_ids,
  treatment = filtered_treatments
)
```

### Estimating Differentially Methylated regions (DMRs) on 'RH_VEH_vs_WT_VEH'
```Ruby
selected_comparison <- 'RH_VEH_vs_WT_VEH'

### Preparing results directory
comparison_dir <- file.path(results_dir, 'Comparisons')
if(!(dir.exists(comparison_dir))){
  dir.create(comparison_dir)
}

### Selecting samples involve in comparison
sub_obj <- Filter4Comparison(
  myobj = filtered_obj,
  comparison_ = selected_comparison
)

### Aggregate bases to regions
region_methylation <- tileMethylCounts(
  sub_obj,
  win.size=WindowSize,
  step.size=WindowSize*(1 - OverLap),
  cov.bases = 3,
  mc.cores = 10)

### Filtering low quality reads
region_methylation <- filterByCoverage(
  region_methylation,
  lo.count=MinCov,
  lo.perc=NULL,
  hi.count=NULL,
  hi.perc = NULL
)

### Uniting read counts
min_per_group_ <- MinPerGroup
if(min_per_group_ != 0L){
  conditions_freq <- table(region_methylation@treatment)
  max_freq <- max(conditions_freq)
  min_per_group_ <- as.integer(max_freq - min_per_group_)
}else{
  min_per_group_ <- NULL
}
region_meth <- methylKit::unite(
  region_methylation,
  destrand = T,
  min.per.group = min_per_group_
)

### Finding differentially methylated regions
region_diff <- calculateDiffMeth(
  region_meth,
  mc.cores=10,
  adjust = 'fdr'
)

results_slug <- MakeFileName(
  region_ = RegionSlug,
  min_cov_ = MinCov,
  min_per_group_ = MinPerGroup,
  comparison_ = selected_comparison,
  dataset_name = dataset_name_temp
)

StoreDMRResults(
  region_meth = region_meth,
  region_diff = region_diff,
  comparison_ = selected_comparison,
  comparison_dir = comparison_dir,
  organism_ = ORGANISM,
  txdb_ = TXDB_mm10,
  results_slug = results_slug
  )

```

### Estimating Differentially Methylated regions (DMRs) on 'RH_MET_vs_RH_VEH'
```Ruby
selected_comparison <- 'RH_MET_vs_RH_VEH'

### Preparing results directory
comparison_dir <- file.path(results_dir, 'Comparisons')
if(!(dir.exists(comparison_dir))){
  dir.create(comparison_dir)
}

### Selecting samples involve in comparison
sub_obj <- Filter4Comparison(
  myobj = filtered_obj,
  comparison_ = selected_comparison
)

### Aggregate bases to regions
region_methylation <- tileMethylCounts(
  sub_obj,
  win.size=WindowSize,
  step.size=WindowSize*(1 - OverLap),
  cov.bases = 3,
  mc.cores = 10)

### Filtering low quality reads
region_methylation <- filterByCoverage(
  region_methylation,
  lo.count=MinCov,
  lo.perc=NULL,
  hi.count=NULL,
  hi.perc = NULL
)

### Uniting read counts
min_per_group_ <- MinPerGroup
if(min_per_group_ != 0L){
  conditions_freq <- table(region_methylation@treatment)
  max_freq <- max(conditions_freq)
  min_per_group_ <- as.integer(max_freq - min_per_group_)
}else{
  min_per_group_ <- NULL
}
region_meth <- methylKit::unite(
  region_methylation,
  destrand = T,
  min.per.group = min_per_group_
)

#### Finding differentially methylated regions
region_diff <- calculateDiffMeth(
  region_meth,
  mc.cores=10,
  adjust = 'fdr'
)

results_slug <- MakeFileName(
  region_ = RegionSlug,
  min_cov_ = MinCov,
  min_per_group_ = MinPerGroup,
  comparison_ = selected_comparison,
  dataset_name = dataset_name_temp
)

StoreDMRResults(
  region_meth = region_meth,
  region_diff = region_diff,
  comparison_ = selected_comparison,
  comparison_dir = comparison_dir,
  organism_ = ORGANISM,
  txdb_ = TXDB_mm10,
  results_slug = results_slug
  )

```

### Visualization using scatter plots
```Ruby
selected_comps <- ComparisonSets
selected_region_conds <- c('ALL', 'PRO', 'ISL')
selected_columns <- c(
  'chr', 'start', 'end', 'width',
  'pvalue', 'qvalue', 'meth.diff',
  "CpG.Island", "CpG.Shore", "CpG.Shelf", "Open.Sea", 
  "annotation", "geneId", "distanceToTSS", "geneSymbol", "gene",
  "cov.Total"
)


dmrs <- list()
for(comp_ in selected_comps){
  print(comp_)
  results_slug <- MakeFileName(
    region_ = RegionSlug,
    min_cov_ = MinCov,
    min_per_group_ = MinPerGroup,
    comparison_ = comp_,
    dataset_name = dataset_name_temp
  )
  file_name <- glue::glue('{results_slug}.tsv')
  dmr_df <- LoadDMRs(
    file_name = file_name,
    lrf = results_dir
  )
  dmr_df <- dmr_df[,selected_columns]
  new_col_names <- colnames(dmr_df)
  names(new_col_names) <- new_col_names
  new_col_names[c('pvalue', 'qvalue', 'meth.diff', 'cov.Total')] <- paste(
    c('pvalue', 'qvalue', 'meth.diff', 'cov.Total'),
    comp_,
    sep = '.'
  )
  colnames(dmr_df) <- unname(new_col_names)
  dmrs[[comp_]] <- dmr_df
}

dmrs <- merge(
  x = dmrs[['RH_VEH_vs_WT_VEH']],
  y = dmrs[['RH_MET_vs_RH_VEH']],
  # by = c('chr', 'start', 'end', 'width', 'strand')
)

selected_columns <- c(
  'meth.diff.RH_VEH_vs_WT_VEH','meth.diff.RH_MET_vs_RH_VEH',
  "CpG.Island", "CpG.Shore", "CpG.Shelf", "Open.Sea", 
  "annotation", "geneId", "distanceToTSS", "geneSymbol", "gene"
)

cond1 <- dmrs$pvalue.RH_MET_vs_RH_VEH <= PvalueCutoff
cond2 <- dmrs$pvalue.RH_VEH_vs_WT_VEH <= PvalueCutoff
cond <- cond1 & cond2
data_set <- dmrs[cond ,selected_columns]
colnames(data_set)[1:2] <- c('V1', 'V2')

col_name_map <- c(V1 = 'RH_VEH_vs_WT_VEH', V2 = 'RH_MET_vs_RH_VEH')

plots_list <- list()
for(selected_region_cond in selected_region_conds){
  region_cond <- AnnoRegionFilter[[selected_region_cond]]
  p <- MyGalaxyPlot(
    data_set_ = data_set,
    region_cond = region_cond,
    col_name_map = col_name_map,
    comparison_slug_map = ComparisonSlugMap,
    lrf = results_dir,
    results_slug = 'Galaxy'
    )
  
  plots_list[[length(plots_list) + 1]] <- p
}

plots_list

```

### Visualizaton using violin plots

```Ruby
selected_comps <- ComparisonSets
selected_region_conds <- c('ALL', 'PRO', 'ISL')
selected_columns <- c(
  'chr', 'start', 'end', 'width',
  'pvalue', 'qvalue', 'meth.diff', 'comparison_name',
  "CpG.Island", "CpG.Shore", "CpG.Shelf", "Open.Sea", 
  "annotation", "geneId", "distanceToTSS", "geneSymbol", "gene",
  "cov.Total"
)


dmrs <- list()
for(comp_ in selected_comps){
  print(comp_)
  results_slug <- MakeFileName(
    region_ = RegionSlug,
    min_cov_ = MinCov,
    min_per_group_ = MinPerGroup,
    comparison_ = comp_,
    dataset_name = dataset_name_temp
  )
  file_name <- glue::glue('{results_slug}.tsv')
  dmr_df <- LoadDMRs(
    file_name = file_name,
    lrf = results_dir
  )
  dmr_df <- dmr_df[,selected_columns]
  dmrs[[comp_]] <- dmr_df
}
dmrs <- do.call(rbind, dmrs) %>%
  as.data.frame()

cond1 <- dmrs$pvalue < PvalueCutoff
dmrs <- dmrs[cond1,]
dmrs$comparison_name <- factor(
  dmrs$comparison_name,
  levels = c('RH_VEH_vs_WT_VEH', 'RH_MET_vs_RH_VEH')
)

plots_list <- list()
for(selected_region_cond in selected_region_conds){
  region_cond <- AnnoRegionFilter[[selected_region_cond]]
  dtt <- dmrs[region_cond$cond(dmrs),]
  
  p <- ggviolin(
    data = dtt,
    x = 'comparison_name',
    y = 'meth.diff'
    ) +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') +
    geom_boxplot(
      aes(x = comparison_name, y = meth.diff),
      width = 0.1,
      outlier.shape = NA
    ) +
    labs(
      y = 'Change in beta value',
      x = 'Comparisons',
      title = region_cond$slug
    )
  
  
  plots_list[[length(plots_list) + 1]] <- p
}

plots_list

```


### Calculation of Methylation Percentage

DMRs were restricted to the ones that were hypomethylated in vehicle-treated Dnmt3aR878H/+ 655
samples relative to Dnmt3aR878H/+ 656 samples. Each dot represents the average value of the DMRs in
657 one sample. n=3 or 4 biological replicates for each condition.

```Ruby
selected_comps <- ComparisonSets
selected_region_conds <- c('ALL', 'PRO', 'ISL')
selected_columns <- c(
  'chr', 'start', 'end', 'width',
  'pvalue', 'meth.diff',
  "CpG.Island", "CpG.Shore", "CpG.Shelf", "annotation"
)



dmrs <- list()
for(comp_ in selected_comps){
  print(comp_)
  results_slug <- MakeFileName(
    region_ = RegionSlug,
    min_cov_ = MinCov,
    min_per_group_ = MinPerGroup,
    comparison_ = comp_,
    dataset_name = dataset_name_temp
  )
  file_name <- glue::glue('{results_slug}.tsv')
  dmr_df <- LoadDMRs(
    file_name = file_name,
    lrf = results_dir
  )
  dmr_df <- dmr_df[,selected_columns]

  dmrs[[comp_]] <- dmr_df
}

dmrs <- merge(
  x = dmrs[[1]],
  y = dmrs[[2]],
  by = c(
    'chr', 'start', 'end', 'width',
    "CpG.Island", "CpG.Shore", "CpG.Shelf", "annotation"
  ),
  suffixes = paste('.', selected_comps, sep = '')
)

cond1 <- dmrs$meth.diff.RH_VEH_vs_WT_VEH < 0
cond2 <- dmrs$pvalue.RH_VEH_vs_WT_VEH < PvalueCutoff
cond3 <- dmrs$pvalue.RH_MET_vs_RH_VEH < PvalueCutoff
cond <- cond1 & cond2 & cond3
dmrs <- dmrs[cond,]

id_cols <- c(
  'chr', 'start', 'end', 'width',
  "CpG.Island", "CpG.Shore", "CpG.Shelf", "annotation"
  )
selected_ids <- dmrs[,id_cols]
selected_ids$id <- paste(
    selected_ids$chr,
    selected_ids$start,
    selected_ids$end,
    sep = '.'
  )
selected_ids$id <- str_replace_all(
  selected_ids$id,
  'chr',
  ''
)


percentage_dt_list <- NA
for(Cond in c('WT_VEH', 'WT_MET', 'RH_VEH', 'RH_MET')){
  selected_samples <- grepl(Cond, filtered_sample_ids)
  selected_samples <- filtered_sample_ids[selected_samples]
  sub_obj <- reorganize(
    filtered_obj,
    sample.ids=selected_samples,
    treatment=filtered_treatments[selected_samples]
  )
  
  ### Aggregate bases to regions
  region_methylation <- tileMethylCounts(
    sub_obj,
    win.size=WindowSize,
    step.size=WindowSize*(1 - OverLap),
    cov.bases = 3,
    mc.cores = 10)
  
  ### Filtering low quality reads
  region_methylation <- filterByCoverage(
    region_methylation,
    lo.count=MinCov,
    lo.perc=NULL,
    hi.count=NULL,
    hi.perc = NULL
  )
  
  ### Uniting read counts
  min_per_group_ <- MinPerGroup
  if(min_per_group_ != 0L){
    conditions_freq <- table(region_methylation@treatment)
    max_freq <- max(conditions_freq)
    min_per_group_ <- as.integer(max_freq - min_per_group_)
  }else{
    min_per_group_ <- NULL
  }
  region_meth <- methylKit::unite(
    region_methylation,
    destrand = T,
    min.per.group = min_per_group_
  )
  
  percentage_dt <- percMethylation(
    region_meth,
    rowids = T) %>%
    as.data.frame()
  
  percentage_dt$id <- rownames(percentage_dt)
  
  percentage_dt <- merge(
    x = selected_ids,
    y = percentage_dt,
    by = 'id'
  )
  
  percentage_dt_list <- merge(
    x = percentage_dt_list,
    y = percentage_dt,
    all = T
  )
}

plots_list <- list()
for(selected_region_cond in selected_region_conds){
  
  region_cond <- AnnoRegionFilter[[selected_region_cond]]
  dtt <- percentage_dt_list[region_cond$cond(percentage_dt_list),filtered_sample_ids]
  dtt <- reshape2::melt(dtt)
  dtt$Cond <- sapply(
    dtt$variable,
    function(x) substr(x, 1,6)
  )
  
  dtt <- dtt %>%
    group_by(Cond, variable) %>%
    summarise(Mean = mean(value, na.rm = T), N = n()) %>%
    as.data.frame()
  
  dtt$Cond <- factor(
    dtt$Cond,
    levels = c('WT_VEH', 'WT_MET', 'RH_VEH', 'RH_MET')
    )
  
  p <- ggscatter(
    data = dtt,
    x = 'Cond',
    y = 'Mean'
  ) +
    labs(
      title = region_cond$slug,
      y = 'Mean Methylation %' 
    )
  
  plots_list[[length(plots_list) + 1]] <- p
  
}

plots_list

```
























































