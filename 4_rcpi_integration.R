library(data.table)
library(ggplot2)
library(patchwork)
library(RISC)
library(Seurat)

tc_norm <- readRDS('timecourse_norm_seu.RDS')
sp_norm <- readRDS('sigpert_norm_seu.RDS')

# get raw counts, metadata, and genes as input for integration
tc_mat <- tc_norm@assays$RNA@counts # count matrix
tc_cmt <- tc_norm[[]] # metadata
tc_gmt <- data.frame('gene' = row.names(tc_mat), row.names = row.names(tc_mat)) # genes

sp_mat <- sp_norm@assays$RNA@counts
sp_cmt <- sp_norm[[]]
sp_gmt <- data.frame('gene' = row.names(sp_mat), row.names = row.names(sp_mat))

# read raw data (each experiment separately) into risc objects, filter and normalise using risc
process <- function(mat, cmt, gmt){
	# read in raw data, filter and normalise
	risc <- readscdata(mat, cmt, gmt)
	risc <- scFilter(risc, min.UMI = 10, max.UMI = Inf, min.gene = 1, min.cell = 1)
	risc <- scNormalize(risc)
	# detect variable genes and plot
	risc <- scDisperse(risc)
	risc <- scPCA(risc)
	risc <- scUMAP(risc, npc = 10)
} 
tc_risc <- process(tc_mat, tc_cmt, tc_gmt)
sp_risc <- process(sp_mat, sp_cmt, sp_gmt)

data0 <- list(tc_risc, sp_risc)

# make gene list and integrate
genelist <- Reduce(intersect, list(sp_risc@vargene, tc_risc@vargene))
data0 <- scMultiIntegrate(objects = data0, eigens = 15, var.gene = genelist, ncore = 4, npc = 50) 

# UMAP reduction of integrated data
data0 <- scUMAP(data0, npc = 15, use = 'PLS')

saveRDS(data0, 'int_risc.RDS')) 
saveRDS(as.matrix(riscvis@assay$logcount), 'int_matrix.RDS') # save matrix of the normalised data for ease


