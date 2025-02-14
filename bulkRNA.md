---
title: "Bulk-RNAseq Analysis"
author: "Ali Chegini"
date: "2025-01-01"
output: html_document
---


```{r setup, include=FALSE}
rm(list = ls())
gc()
knitr::opts_chunk$set(echo = TRUE)
root_dir <- dirname(getwd())
knitr::opts_knit$set(root.dir = root_dir)
```

```{r}
library(biomaRt)
library(enrichplot)

source('Settings.R')
source(file.path(BULKRNASEQ_PATH, 'local_settings.R'))
```

# Functions
```{r}
MyGSEAplot <- function(gsea_obj, pw_id){

  row_number <- which(gsea_obj$ID == pw_id)
  if(length(row_number) != 0){
    
    pw_comp <- gsea_obj@result$comparison_name[row_number]
    pw_des <- gsea_obj@result$Description[row_number]  
    pw_pval <- gsea_obj@result$pvalue[row_number]
    pw_pval <- format(pw_pval, scientific = T, digits = 4)
    pw_qval <- gsea_obj@result$qvalue[row_number]
    pw_qval <- format(pw_qval, scientific = T, digits = 4)
    pw_nes <- gsea_obj@result$NES[row_number]
    pw_nes <- format(pw_nes, scientific = T, digits = 4)
    pw_es <- gsea_obj@result$enrichmentScore[row_number]
    pw_es <- format(pw_es, scientific = T, digits = 4)
    
    title_ <- glue::glue('{pw_comp} - {pw_des}\n(pval: {pw_pval}, ES: {pw_es}, NES: {pw_nes})')
    figure <- gseaplot2(gsea_obj, title = title_, geneSetID = pw_id)
    figure[[1]] <- figure[[1]] + geom_hline(yintercept = 0, linetype = 'dashed')
    return(figure)
  }else{
    return(NA)
  }

}

mapIt <- function(x) {
    require("Orthology.eg.db", character.only = TRUE)
    require("org.Mm.eg.db", character.only = TRUE)
    require("org.Hs.eg.db", character.only = TRUE)
    human <- mapIds(org.Hs.eg.db, x, "ENTREZID", "ALIAS")
    mapped <- select(Orthology.eg.db, human, "Mus.musculus","Homo.sapiens")
    mouse <- mapIds(org.Mm.eg.db, as.character(mapped[,2]), "SYMBOL","ENTREZID")
    mouse <- do.call(c, lapply(mouse, function(x) if(is.null(x)) return(NA) else return(x)))
    cbind(x, mapped, mouse)
}

LoadGeneSets <- function(){
  
  go_gs_file_name <- 'GO_geneSets_2024_Oct_16.tsv'
  kegg_gs_file_name <- 'KEGG_geneSets_2024_Oct_16.tsv'
  reactome_gs_file_name <- 'Reactome_geneSets_2024_Oct_16.tsv'
  wp_gs_file_name <- 'WikiPathway_geneSets_2024_Oct_16.tsv'
  
  
  ## Loading KEGG, Reactome, and WikiPathways gene sets
  entrez_gs_dbs <- c(
    kegg_gs_file_name,
    wp_gs_file_name,
    reactome_gs_file_name
  )
  
  entrez_gs <- list()
  
  for(gs_file_name in entrez_gs_dbs){
  
    file_path <- file.path(MOUSE_GENE_SETS_PATH, gs_file_name)
    term2gene_collapsed <- read_tsv(file_path) %>%
      as.data.frame()
    
    cond1 <- is.na(term2gene_collapsed$Description)
    term2gene_collapsed[cond1, 'Description'] <- term2gene_collapsed[cond1, 'ID']
    
    dtt <- list()
    for(i in c(1:nrow(term2gene_collapsed))){
      gs_id <- term2gene_collapsed[i, 'ID']
      genes_ <- term2gene_collapsed[i, 'gene_set']
      genes_ <- strsplit(genes_, ', ')[[1]]
      gs_description <- term2gene_collapsed[i, 'Description']
      
      dtt[[length(dtt) + 1]] <- data.frame(
        gs_id = rep(gs_id, length(genes_)),
        gene_name = genes_,
        gs_description = rep(gs_description, length(genes_))
      )
    }
    
    dtt <- do.call(rbind, dtt) %>%
      as.data.frame()
    
    entrez_gs[[length(entrez_gs) + 1]] <- dtt
  }
  
  entrez_gs <- do.call(rbind, entrez_gs) %>%
    as.data.frame()
  
  
  ## Loading GO gene sets
  file_path <- file.path(MOUSE_GENE_SETS_PATH, go_gs_file_name)
  term2gene_collapsed <- read_tsv(file_path) %>%
    as.data.frame()
  
  cond1 <- is.na(term2gene_collapsed$Description)
  term2gene_collapsed[cond1, 'Description'] <- term2gene_collapsed[cond1, 'ID']
  
  symbol_gs <- list()
  for(i in c(1:nrow(term2gene_collapsed))){
    gs_id <- term2gene_collapsed[i, 'ID']
    genes_ <- term2gene_collapsed[i, 'gene_set']
    genes_ <- strsplit(genes_, ', ')[[1]]
    gs_description <- term2gene_collapsed[i, 'Description']
    
    symbol_gs[[length(symbol_gs) + 1]] <- data.frame(
      gs_id = rep(gs_id, length(genes_)),
      gene_name = genes_,
      gs_description = rep(gs_description, length(genes_))
    )
  }
  
  symbol_gs[[length(symbol_gs) + 1]] <- data.frame(
        gs_id = rep(gs_id, length(genes_)),
        gene_name = genes_,
        gs_description = rep(gs_description, length(genes_))
      )
  
  symbol_gs <- do.call(rbind, symbol_gs)
  
  return(list(entrez_gs,symbol_gs))
}

```


# Loading Data

```{r}
# load gene sets

all_gs <- LoadGeneSets()
entrez_gs <- all_gs[[1]]
symbol_gs <- all_gs[[2]]

entrez_term2name <- entrez_gs %>%
  group_by(gs_id, gs_description) %>%
  summarise(N = n()) %>%
  as.data.frame()
entrez_term2name <- entrez_term2name[,c('gs_id', 'gs_description')]
entrez_gs <- entrez_gs[,c('gs_id','gene_name')]

symbol_term2name <- symbol_gs %>%
  group_by(gs_id, gs_description) %>%
  summarise(N = n()) %>%
  as.data.frame()
symbol_term2name <- symbol_term2name[,c('gs_id', 'gs_description')]
symbol_gs <- symbol_gs[,c('gs_id','gene_name')]

aa <- c('GO:0006119', 'GO:0006120', 'GO:0046653', 'GO:0035999', 'GO:0006760', 'R-MMU-1592230')
cond1 <- symbol_term2name$gs_id %in% aa
symbol_term2name <- symbol_term2name[cond1,]

cond2 <- symbol_gs$gs_id %in% aa
symbol_gs <- symbol_gs[cond2,]

rm(all_gs)
```

```{r}
file_name <- 'gene_count.csv'
file_dir <- file.path(APP_DB_PATH, file_name)
org_dataset <- read_csv(file_dir) %>%
  as.data.frame()

head(org_dataset)
```

# Compering RH_MET vs RH_VEH (treatment effect on mutant samples)

## Separating ReadCounts and Metadata matrixes
```{r}
comparison_name <- 'RH_MET vs RH_VEH'

# Making expression matrix

selected_cols <- c(
  'RH_MET1',
  'RH_MET2',
  'RH_MET3',
  'RH_VEH1',
  'RH_VEH2',
  'RH_VEH3'
)

expression_dt <- org_dataset[,selected_cols] %>%
  as.data.frame()

# Making genes metadata matrix
selected_cols <- c(
  'gene_id',
  "gene_name",
  "gene_chr",
  "gene_start",
  "gene_end",
  "gene_strand",
  "gene_length",
  "gene_biotype"
)
gene_metadata <- org_dataset[,selected_cols] %>%
  as.data.frame()

rownames(expression_dt) <- gene_metadata$gene_id
rownames(gene_metadata) <- gene_metadata$gene_id

head(expression_dt)
```


## Cleaning expression matrix
```{r}
# Removing probs with zero read counts across all the samples
prob_sum <- rowSums(expression_dt)

cond1 <- prob_sum == 0
selected_gene_ids <- names(prob_sum)[!cond1]

filtered_exp <- expression_dt[selected_gene_ids,]
filtered_metadata <- gene_metadata[selected_gene_ids,]

# Selecting probs with maximum readcount to represent the gene
prob_sum <- rowSums(filtered_exp)
prob_sum <- prob_sum[filtered_metadata$gene_id]
filtered_metadata$ProbSum <- prob_sum 
prob_max <- filtered_metadata %>%
  group_by(gene_name) %>%
  summarise(ProbSum = max(ProbSum))
prob_max <- merge(
  x = filtered_metadata,
  y = prob_max
)

# The gene-set I am going to use for the GSEA used different gene symbols for the 
# mitochondrial protein coding genes for instance it used COX1 instead of mt-co1.
# Therefore, I am going to convert the gene symbols of expression matrix.
cond1 <- grepl('mt-', prob_max$gene_name)
mt_genes <- prob_max[cond1, 'gene_name']

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_info <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                   filters = "external_gene_name",
                   values = mt_genes,
                   mart = ensembl)
gene_info <- gene_info[!is.na(gene_info$entrezgene_id),] # keeping mitochondrial protein coding genes
ids<-bitr(
  gene_info$entrezgene_id, 
  fromType = "ENTREZID", 
  toType = "SYMBOL", 
  OrgDb=ORGANISM)

ids <- merge(
  x = ids,
  y = gene_info,
  by.x = 'ENTREZID',
  by.y = 'entrezgene_id'
)
ids <- as.data.frame(ids)
rownames(ids) <- ids$external_gene_name

cond1 <- prob_max$gene_name %in% ids$external_gene_name
prob_max[cond1, 'gene_name'] <- ids[prob_max[cond1, 'gene_name']  ,'SYMBOL']

print(ids)

# # Selecting genes with Entrezid
ids<-bitr(
  prob_max$gene_name,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb=ORGANISM)

## remove duplicate IDS
c1 <- duplicated(ids[c("ENTREZID")])
c2 <- duplicated(ids[c("SYMBOL")])
dedup_ids = ids[!(c1 | c2), ]
prob_max <- merge(
  x = prob_max,
  y = dedup_ids,
  by.x = 'gene_name',
  by.y = 'SYMBOL'
  )
rownames(prob_max) <- prob_max$gene_id

cond1 <- duplicated(prob_max$ENTREZID)
cond2 <- duplicated(prob_max$gene_name)
cond <- cond1 | cond2

prob_max <- prob_max[!cond,]


# keeping expression data of max probs
selected_gene_ids <- intersect(rownames(filtered_exp), prob_max$gene_id)
filtered_exp <- filtered_exp[selected_gene_ids ,]
filtered_exp <- filtered_exp[complete.cases(filtered_exp), ]
prob_max <- prob_max[selected_gene_ids,]
```

## Making Comparison

```{r}
# Creating Condition Matrix
colData_ <- data.frame(
  row.names = colnames(filtered_exp), 
  condition = sapply(
    colnames(filtered_exp),
    function(x) substr(x, 1, nchar(x) - 1)
  ),
  batch = c(
    '1', '2', '2',
    '1', '2', '2'
  )
)

# Removing zero std genes
std_ex <- apply(filtered_exp,1,sd)
cond1 <- std_ex > 0
# cond1 <- std_ex >= quantile(std_ex, 0.1)
# cond2 <- std_ex <= quantile(std_ex, 0.9)
ex <- filtered_exp[cond1,]

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = ex,
                              colData = colData_, 
                              design = ~ 1 + condition + batch)

dds <- DESeq(
  dds, 
  test = "LRT",
  reduced = ~ 1,
  # fitType = "mean",
  # sfType = "iterate"
  )
```

```{r}
res <- results(
  dds,
  contrast=c("condition", 'RH_MET', 'RH_VEH'),
  parallel = T
  ) %>%
  as.data.frame()

res$gene_id <- rownames(res)
res <- merge(
  x = res,
  y = prob_max,
  by = 'gene_id'
)

head(res)

# Creating Ranked gene list
## filtering for protein coding genes
cond1 <- res$gene_biotype == 'protein_coding'
f_res <- res[cond1,]

## Ranked gene list for kegg and wikipathway
entrez_ranked_gene_list <- -log10(f_res$pvalue) * sign(f_res$log2FoldChange)
names(entrez_ranked_gene_list) <- f_res$ENTREZID
entrez_ranked_gene_list <- entrez_ranked_gene_list[order(entrez_ranked_gene_list, decreasing = T)]

## Ranked gene list for GO
symbol_ranked_gene_list <- -log10(f_res$pvalue) * sign(f_res$log2FoldChange)
names(symbol_ranked_gene_list) <- f_res$gene_name
symbol_ranked_gene_list <- symbol_ranked_gene_list[order(symbol_ranked_gene_list, decreasing = T)]


```

```{r}

cond1 <- entrez_ranked_gene_list != 0
cond2 <- !is.na(entrez_ranked_gene_list)

set.seed(333)
entrez_gsea_ <- GSEA(
  geneList = entrez_ranked_gene_list[cond1 & cond2],
  TERM2GENE = entrez_gs,
  TERM2NAME = entrez_term2name,
  minGSSize = 15,
  maxGSSize = 10000,
  pvalueCutoff = 10,
  verbose = T,
  pAdjustMethod = "fdr",
  seed = T
)
entrez_gsea_@result$comparison_name <- comparison_name

pw_id <- 'WP435'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'WP1248'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'mmu00190'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'mmu00670'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'R-MMU-1428517'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)


```

```{r}

cond1 <- symbol_ranked_gene_list != 0
cond2 <- !is.na(symbol_ranked_gene_list)

set.seed(333)
symbol_gsea_ <- GSEA(
  geneList = symbol_ranked_gene_list[cond1 & cond2],
  TERM2GENE = symbol_gs,
  TERM2NAME = symbol_term2name,
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 10,
  verbose = T,
  pAdjustMethod = "fdr",
  seed = T
)
symbol_gsea_@result$comparison_name <- comparison_name

pw_id <- 'GO:0006119'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0006120'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0046653'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0035999'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0006760'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'R-MMU-1592230'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

```

```{r}
# EPPERT_CE_HSC_LSC
gs_genes <- 'Abcb1, Adgrg6, Alcam, Baalc, Bcl11a, Cacnb2, Crhbp, Dapk1, Dram1, Elk3, Erg, Fam30a, Flt3, Frmd4b, Gucy1a1, Hla-drb4, Hlf, Hoxa5, Hoxb2, Hoxb3, Htr1f, Inpp4b, Kat6a, Mecom, Meis1, Myo5c, Plscr4, Ppp1r16b, Prkch, Rbpms, Rnf125, Slc25a36, Smarca1, Socs2, Spink2, Sptbn1, Tceal9, Tfpi, Tmem38b, Yes1, Znf165'
gs_genes <- str_split(gs_genes, ', ')[[1]]
# gs_genes <- mapIt(gs_genes)
# gs_genes <- gs_genes$mouse
# gs_genes <- gs_genes[!is.na(gs_genes)]

term2gene <- data.frame(
  gs_name = rep('EPPERT_CE_HSC_LSC', length(gs_genes)),
  SYMBOL = gs_genes
)


cond1 <- symbol_ranked_gene_list != 0
cond2 <- !is.na(symbol_ranked_gene_list)

set.seed(333)
gsea_ <- GSEA(
  geneList = symbol_ranked_gene_list[cond1 & cond2],
  TERM2GENE = term2gene,
  nPerm =  500,
  minGSSize = 15,
  maxGSSize = 10000,
  pvalueCutoff = 10,
  verbose = T,
  pAdjustMethod = "fdr",
  seed = T,
  # scoreType = 'neg'
)
gsea_@result$comparison_name <- comparison_name

pw_id <- 'EPPERT_CE_HSC_LSC'
MyGSEAplot(gsea_obj = gsea_, pw_id = pw_id)
```

```{r}
# EPPERT_CE_HSC_LSC
gs_genes <- 'Abcc12, Agxt, Arhgap32, Arhgap42, Arhgef28, Arl4c, Ash1l, Baalc, Bcl2l11, Camta1, Ccar1, Ccdc136, Ccdc185, Ccl19, Cd164l2, Cdip1, Cdk14, Cers4, Chn1, Chrna6, Chrnb1, Clps, Cobll1, Col18a1, Col4a1, Cpsf2, Cpxm1, Ctnnal1, Cxcl14, Cyp2d6, Cyp2e1, Ddit4, Def8, Dnajb13, Dnajc21, Dnmt3a, Dusp5, Dydc1, Efna3, Egr1, Emcn, Emp1, Eng, Esr1, Etv1, Etv3, Fam110b, Fam168b, Fam181b, Fbxw12, Foxb2, Frmd8, Fundc2, Gapdhs, Gata3, Gatad2a, Gem, Glb1l, Gnl1, Hnrnpr, Hook1, Hoxb2, Hyal1, Ifi44, Ifih1, Ift81, Igdcc4, Il36a, Inppl1, Iqgap2, Itga3, Itga9, Itih5, Jade1, Jun, Kcnj6, Kcnk15, Kif5a, Klhl26, Laptm4b, Lpar4, Lrrc49, Lrrtm4, Lypd1, Maff, Marchf9, Mef2d, Mex3a, Mirlet7a1, Misfa, Msi2, Mycbp2, Myo10, Ncoa2, Nfatc2, Nfix, Nfkbia, Nfkbie, Nkain3, Nkx2-8, Nr4a1, Nr4a2, Nrxn1, Nrxn3, Nsd3, Ntrk3, Paip1, Palmd, Pard6g, Pck1, Pdc, Pde4b, Pdf, Pdzd2, Peak1, Pex11a, Pglyrp2, Pitpnc1, Pkd2, Plac8l1, Polr2a, Ppp1r9a, Prcd, Prkacb, Prkce, Prkg1, Ptgr1, Ptprk, Pygm, Rasd1, Rbbp6, Rbm28, Rbx1, Rest, Rfpl4b, Rilpl1, Robo4, Rps6ka3, S100a5, Samd10, Sash1, Scaf11, Sec14l3, Selenbp1, Sema3d, Sema4c, Septin3, Serpinb8, Sgsm1, Siah2, Skil, Slc16a5, Slc23a2, Slc25a30, Slc35f1, Slc41a1, Slit2, Smad7, Smarca2, Smc5, Socs5, Sos1, Srek1, Sult4a1, Sv2a, Syn2, Syt11, Tbxa2r, Tcf25, Tcf4, Tctn2, Tgoln2, Tha1p, Trim61, Trim9, Tspan13, Tspan6, Usp34, Usp38, Vamp2, Vldlr, Vps37a, Wnt6, Wwc2, Zbtb20, Zbtb37, Zc2hc1a, Zc3h12c, Zcchc7, Zdbf2, Zmynd8, Zrsr2, Zxdb'
gs_genes <- str_split(gs_genes, ', ')[[1]]
term2gene <- data.frame(
  gs_name = rep('IVANOVA_HEMATOPOIESIS_STEM_CELL', length(gs_genes)),
  SYMBOL = gs_genes
)

cond1 <- symbol_ranked_gene_list != 0
cond2 <- !is.na(symbol_ranked_gene_list)

set.seed(1234)
gsea_ <- GSEA(
  geneList = symbol_ranked_gene_list[cond1 & cond2],
  TERM2GENE = term2gene,
  # nPerm =  500,
  minGSSize = 15,
  maxGSSize = 10000,
  pvalueCutoff = 10,
  verbose = T,
  pAdjustMethod = "fdr",
  seed = T
)
gsea_@result$comparison_name <- comparison_name

pw_id <- 'IVANOVA_HEMATOPOIESIS_STEM_CELL'
MyGSEAplot(gsea_obj = gsea_, pw_id = pw_id)
```

# Compering RH_VEH vs WT_VEH

## Separating ReadCounts and Metadata matrixes
```{r}

comparison_name <- 'RH_VEH vs WT_VEH'

# Making expression matrix

selected_cols <- c(
  'RH_VEH1',
  'RH_VEH2',
  'RH_VEH3',
  'WT_VEH1',
  'WT_VEH2',
  'WT_VEH3'
)

expression_dt <- org_dataset[,selected_cols] %>%
  as.data.frame()

# Making genes metadata matrix
selected_cols <- c(
  'gene_id',
  "gene_name",
  "gene_chr",
  "gene_start",
  "gene_end",
  "gene_strand",
  "gene_length",
  "gene_biotype"
)
gene_metadata <- org_dataset[,selected_cols] %>%
  as.data.frame()

rownames(expression_dt) <- gene_metadata$gene_id
rownames(gene_metadata) <- gene_metadata$gene_id

head(expression_dt)
```


## Cleaning expression matrix
```{r}
# Removing probs with zero read counts across all the samples
prob_sum <- rowSums(expression_dt)

cond1 <- prob_sum == 0
selected_gene_ids <- names(prob_sum)[!cond1]

filtered_exp <- expression_dt[selected_gene_ids,]
filtered_metadata <- gene_metadata[selected_gene_ids,]

# Selecting probs with maximum readcount to represent the gene
prob_sum <- rowSums(filtered_exp)
prob_sum <- prob_sum[filtered_metadata$gene_id]
filtered_metadata$ProbSum <- prob_sum 
prob_max <- filtered_metadata %>%
  group_by(gene_name) %>%
  summarise(ProbSum = max(ProbSum))
prob_max <- merge(
  x = filtered_metadata,
  y = prob_max
)

# The gene-set I am going to use for the GSEA used different gene symbols for the 
# mitochondrial protein coding genes for instance it used COX1 instead of mt-co1.
# Therefore, I am going to convert the gene symbols of expression matrix.
cond1 <- grepl('mt-', prob_max$gene_name)
mt_genes <- prob_max[cond1, 'gene_name']

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_info <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                   filters = "external_gene_name",
                   values = mt_genes,
                   mart = ensembl)
gene_info <- gene_info[!is.na(gene_info$entrezgene_id),] # keeping mitochondrial protein coding genes
ids<-bitr(
  gene_info$entrezgene_id, 
  fromType = "ENTREZID", 
  toType = "SYMBOL", 
  OrgDb=ORGANISM)

ids <- merge(
  x = ids,
  y = gene_info,
  by.x = 'ENTREZID',
  by.y = 'entrezgene_id'
)
ids <- as.data.frame(ids)
rownames(ids) <- ids$external_gene_name

cond1 <- prob_max$gene_name %in% ids$external_gene_name
prob_max[cond1, 'gene_name'] <- ids[prob_max[cond1, 'gene_name']  ,'SYMBOL']

# # Selecting genes with Entrezid
ids<-bitr(
  prob_max$gene_name,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb=ORGANISM)

## remove duplicate IDS
c1 <- duplicated(ids[c("ENTREZID")])
c2 <- duplicated(ids[c("SYMBOL")])
dedup_ids = ids[!(c1 | c2), ]
prob_max <- merge(
  x = prob_max,
  y = dedup_ids,
  by.x = 'gene_name',
  by.y = 'SYMBOL'
  )
rownames(prob_max) <- prob_max$gene_id

cond1 <- duplicated(prob_max$ENTREZID)
cond2 <- duplicated(prob_max$gene_name)
cond <- cond1 | cond2

prob_max <- prob_max[!cond,]


# keeping expression data of max probs
selected_gene_ids <- intersect(rownames(filtered_exp), prob_max$gene_id)
filtered_exp <- filtered_exp[selected_gene_ids ,]
filtered_exp <- filtered_exp[complete.cases(filtered_exp), ]
prob_max <- prob_max[selected_gene_ids,]
```

## Making Comparison
```{r}
# Creating Condition Matrix
colData_ <- data.frame(
  row.names = colnames(filtered_exp), 
  condition = sapply(
    colnames(filtered_exp),
    function(x) substr(x, 1, nchar(x) - 1)
  ),
  batch = c(
    '1', '2', '2',
    '1', '2', '2'
  )
)

# Removing zero std genes
std_ex <- apply(filtered_exp,1,sd)
cond1 <- std_ex > 0
# cond1 <- std_ex >= quantile(std_ex, 0.1)
# cond2 <- std_ex <= quantile(std_ex, 0.9)
ex <- filtered_exp[cond1,]

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = ex,
                              colData = colData_, 
                              design = ~ 1 + condition + batch)

dds <- DESeq(
  dds, 
  test = "LRT",
  reduced = ~ 1,
  # fitType = "parametric",
  # sfType = "poscounts"
  )
```

```{r}
res <- results(
  dds,
  contrast=c("condition", 'RH_VEH', 'WT_VEH'),
  parallel = T
  ) %>%
  as.data.frame()

res$gene_id <- rownames(res)
res <- merge(
  x = res,
  y = prob_max,
  by = 'gene_id'
)

head(res)

# Creating Ranked gene list
## filtering for protein coding genes
cond1 <- res$gene_biotype == 'protein_coding'
f_res <- res[cond1,]

## Ranked gene list for kegg and wikipathway
entrez_ranked_gene_list <- -log10(f_res$pvalue) * sign(f_res$log2FoldChange)
names(entrez_ranked_gene_list) <- f_res$ENTREZID
entrez_ranked_gene_list <- entrez_ranked_gene_list[order(entrez_ranked_gene_list, decreasing = T)]

## Ranked gene list for GO
symbol_ranked_gene_list <- -log10(f_res$pvalue) * sign(f_res$log2FoldChange)
names(symbol_ranked_gene_list) <- f_res$gene_name
symbol_ranked_gene_list <- symbol_ranked_gene_list[order(symbol_ranked_gene_list, decreasing = T)]


```


```{r}

cond1 <- entrez_ranked_gene_list != 0
cond2 <- !is.na(entrez_ranked_gene_list)

set.seed(333)
entrez_gsea_ <- GSEA(
  geneList = entrez_ranked_gene_list[cond1 & cond2],
  TERM2GENE = entrez_gs,
  TERM2NAME = entrez_term2name,
  minGSSize = 15,
  maxGSSize = 10000,
  pvalueCutoff = 10,
  verbose = T,
  pAdjustMethod = "fdr",
  seed = T
)
entrez_gsea_@result$comparison_name <- comparison_name

pw_id <- 'WP435'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'WP1248'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'mmu00190'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'mmu00670'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

pw_id <- 'R-MMU-1428517'
MyGSEAplot(gsea_obj = entrez_gsea_, pw_id = pw_id)

```

```{r}

cond1 <- symbol_ranked_gene_list != 0
cond2 <- !is.na(symbol_ranked_gene_list)

set.seed(333)
symbol_gsea_ <- GSEA(
  geneList = symbol_ranked_gene_list[cond1 & cond2],
  TERM2GENE = symbol_gs,
  TERM2NAME = symbol_term2name,
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 10,
  verbose = T,
  pAdjustMethod = "fdr",
  seed = T
)
symbol_gsea_@result$comparison_name <- comparison_name

pw_id <- 'GO:0006119'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0006120'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0046653'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0035999'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'GO:0006760'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

pw_id <- 'R-MMU-1592230'
MyGSEAplot(gsea_obj = symbol_gsea_, pw_id = pw_id)

```































