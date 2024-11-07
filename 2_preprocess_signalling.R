sp_dir <- '/path/to/cellranger/count/'

library(Matrix)
library(Seurat)

# FUNCTIONS ##################################################
# function for reading in 10X libraries; returns sparse matrix in gene x cell format with ensembl gene IDs
read10X <- function(libname, type){
	matrix_dir = paste0(dir, libname, 'filtered_feature_bc_matrix/')
	barcode.path <- paste0(matrix_dir, 'barcodes.tsv.gz') 
	features.path <- paste0(matrix_dir, 'features.tsv.gz') 
	matrix.path <- paste0(matrix_dir, 'matrix.mtx.gz') 

	mat <- readMM(matrix.path)
	barcode.names <- read.delim(barcode.path, header=F, stringsAsFactors=F)
	feature.names <- read.delim(features.path, header=F, stringsAsFactors=F)

	colnames(mat) <- barcode.names$V1
	feature.names <- as.data.frame(lapply(feature.names, function(y) gsub('GRCh38_|addSeq_|_', '', y)))
	if (type == 'gex'){
		feature.names$comb <- paste0(feature.names$V2, '--', feature.names$V1)
		bcidx <- (nrow(feature.names)-19):nrow(feature.names)
		feature.names$comb[bcidx] <- gsub('--.*', '', feature.names$comb[bcidx])
		rownames(mat) <- feature.names$comb
	} else if (type == 'bc') {rownames(mat) <- feature.names$V1}
	return(mat)
}

# get same cells in lib and bc
getSameCells <- function(lib_dat, bc_dat){
	lib_sub <- lib_dat[,which(colnames(lib_dat) %in% colnames(bc_dat))]
	bc_sub <- bc_dat[,which(colnames(bc_dat) %in% colnames(lib_dat))]
	return(list(lib_sub, bc_sub))
}

# create seurat objects. goodbcs is the column numbers of the barcodes that 'worked'
# i.e. we know that BC01, 2, and 5 (in lib2) were excluded from sequencing, so all reads mapping to those barcodes are incorrect. 
mkSeu <- function(both, projectname, goodbcs){
	lib_seu <- CreateSeuratObject(counts = both[[1]], project = projectname)
	lib_seu[['BC']] <- CreateAssayObject(counts = both[[2]][goodbcs,])
	return(lib_seu)
}

# process barcoding libraries and print distribution of barcodes, doublets and negatives
assignBCs <- function(seu){
	seu <- NormalizeData(seu, assay = 'BC', normalization.method = 'CLR')
	seu <- HTODemux(seu, assay = 'BC', positive.quantile = 0.99)
	print(table(seu$BC_classification.global))
	print(table(seu$hash.ID))
	return(seu)
}

# run tSNE on the barcodes from each object
mkTSNEs <- function(seu){
	dist.mtx <- as.matrix(dist(t(GetAssayData(seu, assay = 'BC'))))
	seu <- RunTSNE(seu, distance.matrix = dist.mtx, perplexity = 100)
	return(seu)
}

#############################################################
# load gene expression and barcoding matrices, get overlapping cells and put into seurat object
libnames <- c(paste0('2_signalling_lib', 1:3, '_mRNA_scRNAseq'))
barnames <- c(paste0('2_signalling_lib', 1:3, '_bc_scRNAseq'))

libs <- lapply(libnames, read10X, type='gex')
bars <- lapply(barnames, read10X, type='bc')
bots <- mapply(getSameCells, libs, bars)

seu1 <- mkSeu(both1, 'lib1', c(3:20))
seu2 <- mkSeu(both2, 'lib2', c(3:4,6:20))
seu3 <- mkSeu(both3, 'lib3', c(3:20))
seus <- list(seu1, seu2, seu3)

seus <- lapply(seus, assignBCs)
seus <- lapply(seus, mkTSNEs)

seu <- merge(seu1, c(seu2, seu3), add.cell.ids = c('rxn1', 'rxn2', 'rxn3'))

# doublet annotation and filtering ############################
library(scds)
library(SingleCellExperiment)

# doublet annotation by scds
sce <- as.SingleCellExperiment(seu)
sce <- cxds(sce, retRes = TRUE, verb = TRUE, estNdbl = TRUE)
sce <- bcds(sce, retRes = TRUE, verb = TRUE, estNdbl = TRUE)
sce <- cxds_bcds_hybrid(sce, estNdbl = TRUE)

calls <- as.data.frame(colData(sce)[,c('BC_classification.global', 'cxds_call', 'bcds_call', 'hybrid_call')])
calls$BC_classification.global[which(calls$BC_classification.global == 'Doublet')] <- TRUE
calls$BC_classification.global[which(calls$BC_classification.global != 'TRUE')] <- FALSE
calls$BC_classification.global <- as.logical(calls$BC_classification.global)

# number of algorithms that call each barcode a doublet (including barcoding doublets)
calls$count <- apply(calls, 1, function(x){length(x[x == TRUE])})

# cell filtering based on number of genes, library size and percent mitochondrial genes (as per seurat vignette)
seu[['percent.mt']] <- PercentageFeatureSet(seu, features = row.names(seu)[grep('^MT-', row.names(seu))])
filt <- (seu$nFeature_RNA > 2000 & seu$nFeature_RNA < 7500 & seu$nCount_RNA > 5000 & seu$nCount_RNA < 50000 & seu$percent.mt < 25) # vector of true = pass filter threshold
calls$filt <- filt # 4475 are FALSE (don't pass filtering)

seu$callcount <- calls$count
seu$filt <- calls$filt

seu_filt <- subset(seu, subset = callcount < 2 & filt == TRUE) # 33714 48526

# seurat scaling and normalisation on filtered dataset ###########################
seu_norm <- NormalizeData(seu_filt)
seu_norm <- FindVariableFeatures(seu_norm)
seu_norm <- ScaleData(seu_norm, features = row.names(seu_norm))
seu_norm <- RunPCA(seu_norm, dims = 1:50)
seu_norm <- RunUMAP(seu_norm, dims = 1:50)

# assign treatment group and timepoint for each barcode ###########################
oldbc <- c('bcbg03', 'bcbg12', 'bcbg05', 'bcbg14', 'bcbg06', 'bcbg15', 'bcbg08', 'bcbg17', 'bcbg09', 'bcbg18', 'bcbg10', 'bcbg19', 'bcbg13', 'bcbg20', 'bcbg04', 'bcbg11', 'bcbg07', 'bcbg16')
allcond <- rep(c('no_treatment', 'XAV', 'CHIR', 'Dorso', 'K02288', 'BMP4', 'VEGF', 'lowXAV', 'lowDorso'), each = 2)
oldbc_ord <- c(paste0('bcbg0', c(3:9)), paste0('bcbg', c(10, 13, 11, 12, 14:20)))
oldrxnbc <- paste0(rep(c('rxn1', 'rxn2', 'rxn3'), each = 18), oldbc_ord)
newtp <- rep(c(5, 2, 2, 9, 9, 5), each = 9)

b2t <- function(barcodes, newLabel, oldLabel){newLabel[match(barcodes, oldLabel)]}

meta <- seu_norm[[]]
seu_norm[['condition']] <- b2t(meta$BC_maxID, newcond, oldbc)
seu_norm[['timepoint']] <- b2t(paste0(gsub('_.*', '', row.names(meta)), meta$BC_maxID), newtp, oldrxnbc)
seu_norm[['allconditions']] <- b2t(meta$BC_maxID, allcond, oldbc)
seu_norm[['experiment']] <- 'SigPert'

# save output
saveRDS(seu, 'sigpert_raw_seu.RDS') 
saveRDS(seu_filt, 'sigpert_filt_seu.RDS') 
saveRDS(seu_norm, 'sigpert_norm_seu.RDS'))

