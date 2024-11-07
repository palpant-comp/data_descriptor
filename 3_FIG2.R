library(Seurat)
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrastr)
library(scales)

# commonly used colours
treatment_lvls <- c('XAV', 'lowXAV', 'CHIR', 'BMP4', 'Dorso', 'lowDorso', 'K02288', 'VEGF', 'no_treatment') 
treatment_cols <- setNames(c('red', '#ff6969', 'dodgerblue', '#00c79c', 'darkviolet', '#9A6EAF', 'darkmagenta', '#dbba00', 'grey60'), treatment_lvls)
day_lvls <- c(2:9)
day_cols <- setNames(c('#000080', '#24006d', '#48005b', '#6d0049', '#910036', '#b60024', '#da0012', '#ff0000'), day_lvls)

############################################################################################
# figure 2a: tsnes of hto reads (time course data) #########################################
seu_hto <- readRDS('timecourse_raw_seu.RDS')
tsn_hto <- data.table(cell_barcode = colnames(seu_hto), cbind(seu_hto@reductions$tsne@cell.embeddings, seu_hto[[]]))

hto_cols <- setNames(c('#003f5c', '#58508d', '#8a508f', '#bc5090', '#de5a79', '#ff6361', '#ff8531', '#ffa600', 'grey50', 'grey80'), c(paste0('A025', 1:8), 'Doublet', 'Negative'))

f2a <- ggplot(tsn_hto, aes(tSNE_1, tSNE_2, colour = correctedID)) + geom_point_rast(size=0.5) + theme_classic() + scale_colour_manual(values=hto_cols)

############################################################################################
# figure 2b: tsnes of barcoding reads (signalling perturbation data) #######################
seu_bcs <- readRDS('sigpert_raw_seu.RDS')
tsn_bcs <- data.table(cell_barcode = colnames(seu_bcs), cbind(seu_bcs@reductions$tsne@cell.embeddings, seu_bcs[[]]))

bcs_cols <- setNames(c(hue_pal()(18), 'grey50', 'grey90'), levels(as.factor(tsn_bcs$hash.ID)))

f2b <- ggplot(tsn_bcs, aes(tSNE_1, tSNE_2, colour=hash.ID)) + geom_point_rast(size=0.5) + theme_classic() + facet_grid('orig.ident') + scale_colour_manual(values=bcs_cols)

############################################################################################
# figure 2c: umap coloured by each time point (time course) ################################
tc_norm <- readRDS('timecourse_norm_seu.RDS')
tc_meta <- data.table(cell_barcode = colnames(tc_norm), cbind(tc_norm@reductions$umap@cell.embeddings, tc_norm[[]]))

timeplot <- function(d){
	cols <- day_cols
	cols[which(names(cols) != d)] <- 'grey90'
	dt <- rbind(tc_meta[day!=d], tc_meta[day==d])
	p <- ggplot(dt, aes(UMAP_1, UMAP_2, colour=day)) + geom_point_rast(size=0.5) + theme_void() + scale_colour_manual(values=cols) + ggtitle(paste0('day ', d)) + theme(legend.position='none', panel.border=element_rect(fill=NA)) 
	return(p)
}
f2c <- wrap_plots(lapply(names(day_cols), timeplot), ncol=2)

############################################################################################
# figure 2d: umaps coloured by treatment group split by library
sp_norm <- readRDS('sigpert_norm_seu.RDS')
sp_meta <- data.table(cell_barcode = colnames(sp_norm), cbind(sp_norm@reductions$umap@cell.embeddings, sp_norm[[]]))

# focus on one library at a time, colour by condition, and make the non library cells in grey90
rplot <- function(lib, time){
	sp_meta$hl <- sp_meta$allconditions
	sp_meta[get('orig.ident') != lib][['hl']] <- 'other'
	sp_meta[get('timepoint') != time][['hl']] <- 'other'
	treatment_cols[['other']] <- 'grey90'
	sp_meta <- rbind(sp_meta[hl == 'other'], sp_meta[hl != 'other'])
	p <- ggplot(sp_meta, aes(UMAP_1, UMAP_2, colour=hl)) + geom_point_rast(size=0.5) + theme_void()  + scale_colour_manual(values=treatment_cols) + ggtitle(paste0(lib,' day', time)) + theme(legend.position='none', panel.border=element_rect(fill=NA)) 
	return(p)
}
# each valid combination of library and timepoint:
sets <- sp_meta[,.N,by=c('orig.ident', 'timepoint')]; setkey(sets, timepoint) 
f2d <- wrap_plots(mapply(rplot, sets$orig.ident, sets$timepoint, SIMPLIFY=F), ncol=2)

