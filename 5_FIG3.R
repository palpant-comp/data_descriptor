library(data.table)
library(ggplot2)
library(patchwork)
library(Seurat)
library(ggrastr)

# colours
treatment_lvls <- c('XAV', 'lowXAV', 'CHIR', 'BMP4', 'Dorso', 'lowDorso', 'K02288', 'VEGF', 'no_treatment') 
treatment_cols <- setNames(c('red', '#ff6969', 'dodgerblue', '#00c79c', 'darkviolet', '#9A6EAF', 'darkmagenta', '#dbba00', 'grey60'), treatment_lvls)
day_lvls <- c(2:9)
day_cols <- setNames(c('#000080', '#24006d', '#48005b', '#6d0049', '#910036', '#b60024', '#da0012', '#ff0000'), day_lvls)

########################################################
# Figure 3a: pre- and post-integration PCA and umaps

# combine the two pre-integration datasets and generate pca and umaps
tc_norm <- readRDS('timecourse_norm_seu.RDS')
sp_norm <- readRDS('sigpert_norm_seu.RDS')

both_raw <- cbind(as.matrix(tc_norm@assays$RNA@counts,sp_norm@assays$RNA@counts))
both_seu <- CreateSeuratObject(both_raw)
both_seu@meta.data <- rbind(tc_norm[[]], sp_norm[[]])

both_seu <- NormalizeData(both_seu)
both_seu <- FindVariableFeatures(both_seu)
both_seu <- ScaleData(both_seu, features = row.names(both_seu))
both_seu <- RunPCA(both_seu, dims = 1:50)
both_seu <- RunUMAP(both_seu, dims = 1:50)

# extract pre-integration coordinates
rawdim <- do.call(cbind, list(both_seu[[]], both_seu@reductions$umap@cell.embeddings, both_seu@reductions$pca@cell.embeddings[,1:2]))
colnames(rawdim)[grep('PC|UMAP', colnames(rawdim))] <- paste0('raw_', grep('PC|UMAP', colnames(rawdim),value=T))

# get post-integration coordinates #####
int_risc <- readRDS('int_risc.RDS')
int_meta <- as.data.table(riscvis@coldata)

int_pcs <- data.table('cell_barcode'=gsub('Set1_|Set2_', '', row.names(riscvis@DimReduction$cell.pls)), riscvis@DimReduction$cell.pls[,1:2])
int_ums <- data.table('cell_barcode'=gsub('Set1_|Set2_', '', row.names(riscvis@DimReduction$cell.umap)), riscvis@DimReduction$cell.umap)
intdim <- Reduce(merge, list(int_meta, int_pcs, int_ums))
colnames(intdim)[grep('PC|UMAP', colnames(intdim))] <- paste0('int_', grep('PC|UMAP', colnames(intdim),value=T))

# combine pre and post and all pca umap into one table
alldim <- as.data.table(merge(rawdim, intdim[,.(cell_barcode, int_PC1, int_PC2, int_UMAP1, int_UMAP2)]))

# define colour groups for plots
alldim$d0 <- ifelse(alldim$experiment=='timecourse', alldim$day, 'other')
alldim$db <- ifelse(alldim$experiment=='SigPert', as.character(alldim$condition), 'other')
alldim$db2 <- ifelse(alldim$experiment=='SigPert', as.character(alldim$day), 'other')

# reorder so that grey points are always at the bottom (they have to come first in the table) and ensure time/conditionsa are shuffled
alldim0 <- rbind(alldim[experiment=='SigPert'], alldim[experiment=='timecourse']) 
alldimb <- rbind(alldim[experiment=='timecourse'], alldim[experiment=='SigPert'][sample(length(which(alldim$experiment=='SigPert')))]) 

# set other experiment as grey
dcols <- c(day_cols, 'other'='grey90')
tcols <- c(treatment_cols, 'other'='grey90')

# plot bfore integrate
p1 <- ggplot(alldim0,aes(raw_PC_1,raw_PC_2,colour=d0))+ggplot(alldim0,aes(raw_UMAP_1,raw_UMAP_2,colour=d0))&scale_colour_manual(values=dcols)
p2<-ggplot(alldimb,aes(raw_PC_1,raw_PC_2,colour=db))+ggplot(alldimb,aes(raw_UMAP_1,raw_UMAP_2,colour=db))&scale_colour_manual(values=tcols)
p3<-ggplot(alldimb,aes(raw_PC_1,raw_PC_2,colour=db2))+ggplot(alldimb,aes(raw_UMAP_1,raw_UMAP_2,colour=db2))&scale_colour_manual(values=dcols)
pp1<-p1/p3/p2&geom_point_rast(size=0.5)&theme_classic()&theme(legend.position='bottom');pp1

# plot after integrate
p4<-ggplot(alldim0,aes(int_PC1,int_PC2,colour=d0))+ggplot(alldim0,aes(int_UMAP1,int_UMAP2,colour=d0))&scale_colour_manual(values=dcols)
p5<-ggplot(alldimb,aes(int_PC1,int_PC2,colour=db))+ggplot(alldimb,aes(int_UMAP1,int_UMAP2,colour=db))&scale_colour_manual(values=tcols)
p6<-ggplot(alldimb,aes(int_PC1,int_PC2,colour=db2))+ggplot(alldimb,aes(int_UMAP1,int_UMAP2,colour=db2))&scale_colour_manual(values=dcols)
pp2<-p4/p6/p5&geom_point_rast(size=0.5)&theme_classic()&theme(legend.position='bottom');pp2

##############################################################
# Figure 3b: integrated umaps with each treatment group coloured separately

hl2 <- function(focus, df){
	'%!in%' <- function(x,y){!('%in%'(x,y))}
	colours <- treatment_cols
	colours[which(names(colours) %!in% focus)] <- 'grey90'
	d <- rbind(df[df$condition %!in% focus,], df[df$condition %in% focus,])
	p <- ggplot(d, aes(UMAP1, UMAP2, colour=condition)) + geom_point_rast(size=0.1) + theme_void() + scale_colour_manual(values=colours) + ggtitle(focus) + theme(panel.border=element_rect(fill=NA), legend.position='none')
	return(p)}

f3b <- wrap_plots(lapply(names(treatment_cols), hl2, df=int_meta))

##############################################################
# Figure 3c: umap plots of different marker genes in integrated data

pltex <- data.table(riscvis@DimReduction$cell.umap, riscvis@coldata, t(as.matrix(riscvis@assay$logcount))) # rcpi-corrected values

pltg <- function(ingene, df=pltex, fac='.~experiment'){
	gene <- grep(paste0('^', ingene, '--'), colnames(df), value = T)[1]
	if (is.na(gene)==TRUE){gene <- ingene; df[[gene]] <- NaN}
	df <- df[,c('UMAP1', 'UMAP2', 'experiment', gene), with=F]; df2 <- df
	df2$experiment <- ifelse(df2$experiment=='0Xav', 'SigPert', '0Xav')
	df2[[gene]] <- NA
	df3 <- rbind(df2, df)
	p <- ggplot(df3, aes(UMAP1, UMAP2, colour=get(gene))) + geom_point_rast(size=0.05) + theme_void() + facet_grid(fac) + scale_colour_viridis_c(option='mako', na.value='grey90') + theme(panel.border=element_rect(fill=NA), legend.position='none') + ggtitle(ingene)
	return(p)
}


genes <- c('MIXL1', 'NOG', 'CDH5', 'PAX3', 'PRRX1', 'MESP1', 'FOXA2', 'SOX2', 'TTR', 'HAND1', 'MYH6', 'T', 'GSC', 'FOXD3', 'EOMES', 'DNAAF3', 'SOX11', 'TTYH1', 'NR2F2', 'ACTC1', 'AFP')

for (gene in genes){
	print(gene)
	plt <- pltg(gene)
	pdf(paste0(wdir, 'outs/F3c_', gene, '.pdf'), height = 7, width = 13)
	print(plt)
	dev.off()
}
