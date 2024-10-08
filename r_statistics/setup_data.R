library(tidyverse)
library(matrixStats)
library(biomaRt)
library(GO.db)

# For t-test ----

set.seed(1)
n_sample_groupA <- 5
n_sample_groupB <- 5
mean_groupA <- 2
sd_groupA <- 0.2
mean_groupB <- 3
sd_groupB <- 1
gene_groupA <- rnorm(n_sample_groupA, mean_groupA, sd_groupA)
gene_groupB <- rnorm(n_sample_groupB, mean_groupB, sd_groupB)
df <- data.frame(
    group = rep(c("groupA", "groupB"), c(n_sample_groupA, n_sample_groupB)),
    gene_exprs = c(gene_groupA, gene_groupB)
)
write.csv(df, "data/gene_exprs.csv", row.names = FALSE, quote = FALSE)

# Import source data ----

## Read full normalised matrix ----

logcounts_matrix <- read.csv("../data/aulicino/logcounts_matrix.txt", header = TRUE, sep = "\t")
logcounts_matrix <- as.matrix(logcounts_matrix)
logcounts_matrix <- round(logcounts_matrix, digits = 2)
dim(logcounts_matrix)

## Re-export gene information ----

gene_metadata <- read.csv("../data/aulicino/gene_metadata.txt", header = TRUE, sep = "\t")
dim(gene_metadata)

## Read cell metadata ----

cell_metadata <- read.csv("../data/aulicino/cell_metadata.txt", header = TRUE, sep = "\t")
dim(cell_metadata)

## Gene sets ----

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# subset to biological processes associated with genes in the universe
out <- getBM(
  attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
  filters = c("ensembl_gene_id", "go_parent_term"),
  values = list(rownames(logcounts_matrix), "GO:0008150"),
  mart = mart)

# Subset and reorder data ----

## Gene information ----

gene_metadata_subset <- gene_metadata[rownames(logcounts_matrix), c("gene_id", "gene_name")]
dim(gene_metadata_subset)

## Cell information ----

cell_metadata_subset <- subset(cell_metadata,
  Time == "6h" &
    Infection %in% c("Mock", "STM-LT2") &
    Status %in% c("Uninfected", "Infected") &
    Sample %in% colnames(logcounts_matrix),
  select = c("Sample", "Infection"))
dim(cell_metadata_subset)

logcounts_matrix_subset <- logcounts_matrix[gene_metadata_subset$gene_id, cell_metadata_subset$Sample]

## subset to gene sets with 10 to 100 genes in the universe ----

keep_go_ids <- out %>% 
  group_by(go_id) %>% 
  summarise(n = n()) %>% 
  filter(n >= 10 & n <= 100) %>% 
  pull(go_id)
gene_go_bp_subset <- subset(out, go_id %in% keep_go_ids, c("ensembl_gene_id", "go_id"))
dim(gene_go_bp_subset)

## Gene set information ----

# columns(GO.db)
go_subset <- AnnotationDbi::select(
  GO.db,
  keys = unique(gene_go_bp_subset$go_id),
  columns = c("GOID", "TERM"))

## Export files ----

write.csv(logcounts_matrix_subset, "data/logcounts_matrix.csv", quote = FALSE)
write.csv(cell_metadata_subset, "data/cell_metadata.csv", row.names = FALSE, quote = FALSE)
write.csv(gene_metadata_subset, "data/gene_metadata.csv", row.names = FALSE, quote = FALSE)
write.csv(gene_go_bp_subset, "data/human_go_bp.csv", quote = FALSE, row.names = FALSE)
write.csv(go_subset, "data/go_info.csv", quote = TRUE, row.names = FALSE)
