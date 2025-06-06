# label transfer analysis with three in vivo datasets

library(Seurat)
library(data.table)
library(Matrix)
library(ggrastr)
library(TRIAGE)
set.seed(7)

# load in query (our) data, as integrated gene expression matrix and metadata with clusters: 
mat <- readRDS('int_matrix.RDS')
seu <- readRDS('int_clust0.3_seu.RDS')
met <- seu[[]]

# load in reference data -- example code provided
ref_norm <- readRDS('ref_matrix.RDS') # normalised expression matrix of reference dataset
ref_meta <- readRDS('ref_metadata.RDS') # metadata

# TRIAGE-transformation of both datasets
process <- function(in_mat, species, in_cells, in_meta){
	disc <- TRIAGEgene(in_mat, species=species, log=FALSE)
	colnames(disc) <- in_cells
	dseu <- CreateSeuratObject(counts=disc)
	dseu@meta.data <- in_meta
	dseu <- FindVariableFeatures(dseu)
	dseu <- ScaleData(dseu)
	dseu <- RunPCA(dseu, npcs = 30)
	return(dseu)}

seu_query <- process(in_mat=mat, species='Human', in_cells=gsub('--.*', '', row.names(mat)), in_meta=met)
seu_refer <- process(in_mat=ref_norm, species='Mouse', in_cells=colnames(ref_norm), in_meta=ref_meta)

# label transfer
# edit refdata based on metadata column containing the reference data's cell type annotations
anc <- FindTransferAnchors(reference=seu_refer, query=seu_query, dims=1:30)
prd <- TransferData(anchorset=anc, refdata=seu_refer$celltype, dims=1:30) 

# save results; repeat for all datasets
met_LT <- data.table(met, lt1.id=prd$predicted.id, lt1.maxscore=pred$prediction.score.max) 

# look at histogram of prediction scores to select thresholds
ggplot(met_LT, aes(lt1.maxscore, colour=lt1.id)) + geom_density() + theme_classic()

# apply selected thresholds
lt1_th <- 0.4
met_LT$lt1.id.filt <- ifelse(met_LT$lt1.maxscore > lt1_th, met_LT$lt1.id, NA)

# plot output in bar plot (figure 4d)
ggplot(met_LT[,.N,by=.(RNA_snn_res.0.3, lt1.id.filt)], aes(RNA_snn_res.0.3, N, fill=lt1.id.filt)) + geom_col(position='fill', colour='black') + theme_classic() + theme(legend.position='bottom') + coord_flip()

# look at top labels per group
met_LT[,.N,by=lt1.id.filt][order(-N)]
