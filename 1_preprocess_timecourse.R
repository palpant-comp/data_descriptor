tc_gex_dir <- '/path/to/cellranger/count/sample/outs/filtered_feature_bc_matrix/'
tc_hto_dir <- '/path/to/cite-count/sample/umi_count/'

library(Seurat)
library(ggplot2)
library(ggforce)
library(patchwork)
library(data.table)
library(Matrix)

# load in data and perform preprocesing and dimensionality reduction for QC ###########
dat <- Read10X(data.dir = tc_gex_dir)
seu <- CreateSeuratObject(counts = dat)

seu[['percent.mt']] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[['percent.rb']] <- PercentageFeatureSet(seu, pattern = "^RPS|^RPL")

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)
seu <- ScaleData(seu, features = row.names(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu))
seu <- RunUMAP(seu, dims = 1:10)
DimPlot(seu, reduction = 'umap')

# scds doublet detection ##############################################################
library(scds)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(seu)

# run three doublet prediction methods, save in seurat object metadata
sce <- cxds(sce, retRes = TRUE, verb = TRUE, estNdbl = TRUE)
sce <- bcds(sce, retRes = TRUE, verb = TRUE, estNdbl = TRUE)
sce <- cxds_bcds_hybrid(sce, estNdbl = TRUE)

calls <- as.data.frame(colData(sce)[,c('cxds_call', 'bcds_call', 'hybrid_call')])
seu@meta.data <- cbind(seu@meta.data, calls)

# demultiplexing HTO counts ###########################################################
read10X <- function(matrix_dir){
	features.path <- paste0(matrix_dir, 'features.tsv.gz')
	barcodes.path <- paste0(matrix_dir, 'barcodes.tsv.gz')
	matrix.path <- paste0(matrix_dir, 'matrix.mtx.gz')
	features.names <- read.delim(features.path, header=F, stringsAsFactors=F)
	barcodes.names <- read.delim(barcodes.path, header=F, stringsAsFactors=F)
	ab <- readMM(matrix.path)
	colnames(ab) <- barcodes.names$V1
	row.names(ab) <- features.names$V1
	return(ab)
}

# load in and format hto counts
hto <- read10X(tc_hto_dir)
colnames(hto) <- paste0(colnames(hto), '-1')
row.names(hto) <- gsub('-.*', '', rownames(hto))

# get overlapping cells from gene expression and hto matrix, run hto demux
sharedCells <- intersect(colnames(seu), colnames(hto)) # 19997 cells
seu_hto <- seu[,which(colnames(seu) %in% sharedCells)]
seu_hto[['BC']] <- CreateAssayObject(counts = hto[,which(colnames(hto) %in% sharedCells)])
seu_hto <- NormalizeData(seu_hto, assay = 'BC', normalization.method = 'CLR')
seu_hto <- HTODemux(seu_hto, assay = 'BC', positive.quantile = 0.99)

# generate tsnes based on hto reads
dist.mtx <- as.matrix(dist(t(GetAssayData(seu_hto, assay = 'BC'))))
seu_hto <- RunTSNE(seu_hto, distance.matrix = dist.mtx, perplexity = 100)

# "unmapped" is treated as a HTO, so reassign cells labelled as "unmapped" as "Negative" and doublets with "unmapped" are singlets for the top barcode.
correctedID <- gsub('.*_.*', 'Doublet', gsub('unmapped', 'Negative', gsub('_unmapped', '', seu_hto[[]]$BC_classification)))
correctedClass <- seu_hto[[]]$BC_classification.global
correctedClass[grep('_unmapped', seu_hto[[]]$BC_classification)] <- 'Singlet'

seu_hto[['correctedID']] <- correctedID
seu_hto[['correctedClass']] <- correctedClass

# summarise QC data and remove low quality cells #####################################
doubletCount <- apply(seu[[]][c('correctedClass', 'cxds_call', 'bcds_call', 'hybrid_call')], 1, function(x){length(x[x %in% c(TRUE, 'Doublet')])})
seu_hto[['doubletCount']] <- doubletCount

df <- data.frame(seu@reductions$umap@cell.embeddings, seu@reductions$tsne@cell.embeddings, seu[[]])
seu_hto[['filter']] <- ifelse(df$nCount_RNA > 1e4 & df$nCount_RNA < 1e5 & df$nFeature_RNA > 2500 & df$nFeature_RNA < 1e4 & df$percent.mt < 20 & df$percent.rb < 45 & doubletCount < 2 & df$correctedClass != 'Negative', 'pass', 'fail')

seu_filt <- subset(seu_hto, subset = filter == 'pass') # 13682 cells remain 

# seurat scaling and normalisation on filtered dataset ###########################
seu_norm <- NormalizeData(seu_filt)
seu_norm <- FindVariableFeatures(seu_norm)
seu_norm <- ScaleData(seu_norm, features = row.names(seu_norm))
seu_norm <- RunPCA(seu_norm, dims = 1:50)
seu_norm <- RunUMAP(seu_norm, dims = 1:50)

# assign timepoints to barcode and save in metadata
seu_norm[['day']] <- factor(c(9,2:8)[match(seu[[]]$BC_maxID, c(paste0('A025', 1:8)))], levels=c(2:9))
seu_norm[['experiment']] <- 'timecourse'

# save output
saveRDS(seu_hto, 'timecourse_raw_seu.RDS') 
saveRDS(seu_filt, 'timecourse_filt_seu.RDS') 
saveRDS(seu_norm, 'timecourse_norm_seu.RDS') 

