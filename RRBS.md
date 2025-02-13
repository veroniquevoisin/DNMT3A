# CITE-seq
## code from data included in manuscript - Ali Chegini :
***Metformin reduces the competitive advantage of Dnmt3a R878H HSPCs*** <br>
Mohsen Hosseini , Veronique Voisin , Ali Chegini , Angelica Varesi , Severine Cathelin ,
Dhanoop Manikoth Ayyathan , Alex C.H. Liu , Yitong Yang , Vivian Wang , Abdula Maher,
Eric Grignano , Julie A. Reisz , Angelo D’Alessandro , Kira Young , Yiyan Wu , Martina
Fiumara , Samuele Ferrari , Luigi Naldini , Federico Gaiti , Shraddha Pai , Grace Egan ,
Aaron D. Schimmer , Gary D. Bader , John E. Dick , Stephanie Z. Xie , Jennifer J.
Trowbridge , and Steven M. Chan 



```{r}
library(methylKit)
library(ggvenn)

source('Settings.R')
source(file.path(METHYLOMEANALYSIS_PATH, 'local_settings.R'))

source(file.path(modules_dir, 'Filter4Comparison.R'))
source(file.path(modules_dir, 'MakeFileName.R'))
source(file.path(modules_dir, 'MyAnnotateRegions.R'))
source(file.path(modules_dir, 'StoreDMRResults.R'))
source(file.path(modules_dir, 'LoadDMRs.R'))
source(file.path(modules_dir, 'AnnotatedRegionFilter.R'))
source(file.path(modules_dir, 'MyGalaxyPlot.R'))
```

# Loading Data

```{r}
### Load hyper parameters
source(file.path(modules_dir, 'GlobalMethylationHyperParameters.R'))


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

# excluding low quality samples.
```{r}
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

# DMRs
```{r}
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

```{r}
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

Galaxy plots
```{r fig.width= 7, fig.height= 5}
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

# Violin plot

```{r fig.width= 5, fig.height= 6}
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

# Venn plot

```{r fig.width= 5, fig.height= 6}
selected_comps <- ComparisonSets
selected_region_conds <- c('ALL', 'PRO', 'ISL')
selected_columns <- c(
  'chr', 'start', 'end', 'width',
  'pvalue', 'comparison_name',
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
  cond1 <- dmr_df$pvalue < PvalueCutoff
  dmr_df <- dmr_df[cond1, ]

  dmrs[[comp_]] <- dmr_df
}

dmrs <- merge(
  x = dmrs[[1]],
  y = dmrs[[2]],
  by = c(
    'chr', 'start', 'end', 'width',
    "CpG.Island", "CpG.Shore", "CpG.Shelf", "annotation"
  ),
  all = T
)

dmrs$id <- paste(
    dmrs$chr,
    dmrs$start,
    dmrs$end,
    dmrs$width,
    sep = '.'
  )

plots_list <- list()
for(selected_region_cond in selected_region_conds){
  
  region_cond <- AnnoRegionFilter[[selected_region_cond]]
  dtt <- dmrs[region_cond$cond(dmrs),]
  
  cond1 <- !is.na(dtt$comparison_name.x)
  set1 <- dtt$id[cond1]
  
  cond1 <- !is.na(dtt$comparison_name.y)
  set2 <- dtt$id[cond1]
  
  sets <- list(set1, set2)
  names(sets) <- selected_comps
  
  p <- ggvenn(
    sets, 
    fill_color = c("#0073C2FF", "#EFC000FF"),
    stroke_size = 0.5, set_name_size = 4
    ) + ggtitle(region_cond$slug)
  
  plots_list[[length(plots_list) + 1]] <- p
}

plots_list


```

# Methylation Percentage

DMRs were restricted to the ones that were hypomethylated in vehicle-treated Dnmt3aR878H/+ 655
samples relative to Dnmt3aR878H/+ 656 samples. Each dot represents the average value of the DMRs in
657 one sample. n=3 or 4 biological replicates for each condition.

```{r fig.width= 5, fig.height=5}
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

```{r}
sessionInfo()
```






















































