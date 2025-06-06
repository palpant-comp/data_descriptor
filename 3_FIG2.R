library(Seurat)
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrastr)
library(scales)
libarary(pheatmap)

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
# figure 2d: umaps coloured by treatment group, split by timepoint and barcode
sp_norm <- readRDS('sigpert_norm_seu.RDS')
sp_meta <- data.table(cell_barcode = colnames(sp_norm), cbind(sp_norm@reductions$umap@cell.embeddings, sp_norm[[]]))

# create background cells (light grey) by repeating coordinates across all facet levels
facet_combos <- unique(sp_meta[, .(BC_maxID, timepoint)])
background_data <- facet_combos[, cbind(.SD, sp_meta[, .(UMAP_1, UMAP_2)]), by = .(BC_maxID, timepoint)]

# defne replicate order for visualisation in both fore- and background points
sp_meta$allconditions <- factor(sp_meta$allconditions, levels=treatment_lvls)
sp_meta$BC_maxID <- factor(sp_meta$BC_maxID, levels=sp_meta[,.N,by=c('BC_maxID', 'allconditions')][order(allconditions)]$BC_maxID)
bc_order <- sp_meta[, .N, by = .(allconditions, BC_maxID)][order(allconditions, -N)  ]$BC_maxID # group by treatment, sort largest on top
sp_meta$BC_maxID <- factor(sp_meta$BC_maxID, levels = bc_order) # Apply factor levels to BC_maxID
background_data$BC_maxID <- factor(background_data$BC_maxID, levels = bc_order)

f2d <- ggplot() + geom_point_rast(data = background_data, aes(UMAP_1, UMAP_2), colour = "grey90", size = 0.1) + geom_point_rast(data = sp_meta, aes(UMAP_1, UMAP_2, colour = allconditions), size = 0.1) + facet_grid(BC_maxID ~ timepoint) + theme_void() + theme(legend.position = "none", panel.border = element_rect(fill = NA),strip.text.y = element_text(angle = -90)) + scale_colour_manual(values = treatment_cols)


############################################################################################
# figure 2e: heatmap of correlation between mean expression in each sample
sp_norm <- readRDS('sigpert_norm_seu.RDS')
sp_meta <- data.table(cell_barcode = colnames(sp_norm), cbind(sp_norm@reductions$umap@cell.embeddings, sp_norm[[]]))
expr_mat <- GetAssayData(sp_norm)

# update barcode labels and append to time point and treamtent information 
sp_meta[,newbc:=paste0('BC', c(paste0(0, 1:9), 10:18))[match(sp_meta$BC_maxID, paste0('bcbg', c(paste0(0, 3:9), 10:20)))]]
sp_meta[,newconds:=gsub('no_treatment', 'Control', gsub('low', 'Low ', allconditions))]
sp_meta[,group := paste(newbc, paste0('Day ', timepoint), newconds, sep = ", ")]
sp_meta$group <- factor(sp_meta$group, levels=unique(sp_meta[,c('BC_maxID', 'timepoint', 'group')])[order(timepoint, BC_maxID)]$group)
group_labels <- setNames(sp_meta$group, sp_meta$rn)
group_levels <- unique(sp_meta$group)

# calculate mean expresion in each sample (group)
avg_expr_list <- lapply(group_levels, function(g) {
  cells <- names(group_labels[group_labels == g])
  rowMeans(expr_mat[, cells, drop = FALSE])})

# Convert to a matrix: genes x groups and calculate pairwise pearson correlation 
avg_expr_mat <- do.call(cbind, avg_expr_list)
cor_mat <- cor(avg_expr_mat, method = "pearson")

# plot heatmap of correlations, with and without hierarchical clustering on the rows
pheatmap(cor_mat[levels(sp_meta$group), levels(sp_meta$group)], cluster_rows=T, cluster_cols=T, color=colorRampPalette(c('black', 'red'))(100), breaks=seq(0.9, 1, length.out=100))
pheatmap(cor_mat[levels(sp_meta$group), levels(sp_meta$group)], cluster_rows=F, cluster_cols=T, color=colorRampPalette(c('black', 'red'))(100), breaks=seq(0.9, 1, length.out=100))

# note on the selected heatmap scale: the lowest value in cor_mat is 0.9052171
min(cor_mat)

