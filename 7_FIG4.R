library(data.table)
library(ggplot2)
library(patchwork)
library(ggrastr)
library(Seurat)
library(scales)
set.seed(7)

seu <- readRDS('int_clust0.3_seu.RDS')

# prepare metadata data.table for ease of plotting
met <- as.data.table(seu[[c('UMAP1', 'UMAP2', 'cell_barcode', 'run', 'day', 'condition', 'experiment', 'Set', 'Barcode', grep('RNA_snn', colnames(seu[[]]), value =T))]])
met$day <- as.factor(met$day)
met$condition <- factor(ifelse(met$experiment=='0Xav', 'no_treatment_tc', as.character(met$condition)), levels=treatment_lvls)

##########################################
# Figure 4a: plot clustering output umaps (different resolutions and each cluster for chosen 0.3 resolution)
pltres <- function(res, df){
	resname <- paste0('RNA_snn_res.', res)
	p1 <- ggplot(df, aes(UMAP1, UMAP2, colour=get(resname))) + ggtitle(resname)
	p2 <- p1 + facet_wrap(resname)
	p1 + p2 & geom_point_rast(size=0.1) & theme_void() & theme(legend.position='none', panel.border=element_rect(fill=NA))}

hl <- function(focus, df){
	'%!in%' <- function(x,y){!('%in%'(x,y))}
	colours <- setNames(hue_pal()(length(levels(df$RNA_snn_res.0.3))), levels(df$RNA_snn_res.0.3))
	colours[which(names(colours) %!in% focus)] <- 'grey90'
	d <- rbind(df[df$RNA_snn_res.0.3 %!in% focus,], df[df$RNA_snn_res.0.3 %in% focus,])
	p <- ggplot(d, aes(UMAP1, UMAP2, colour=RNA_snn_res.0.3)) + geom_point_rast(size=0.2) + theme_void() + scale_colour_manual(values=colours) + ggtitle(focus) + theme(panel.border=element_rect(fill=NA), legend.position='none')
	return(p)}

f4a_top <- wrap_plots(lapply(seq(0.1, 0.4, 0.1), pltres, df=met))
f4a_bot <- wrap_plots(lapply(levels(pltdf$RNA_snn_res.0.3), hl, df = met))

##########################################
# Figure 4b: gene expression bubble plots
mat <- readRDS('int_matrix.RDS')
markers <- fread('T4_markergenes.txt')

bplt <- function(genes, ensembl=TRUE){
	pltex <- as.data.table(cbind(met, t(mat[genes,])))
	long <- melt(pltex, measure.vars=grep('--', colnames(pltex), value=T), id.vars=c('RNA_snn_res.0.3'))
	summ <- long[,.(mean_expr=mean(value), pct_pos=length(which(value>0))/length(value)),by=.(RNA_snn_res.0.3, variable)]
	# gene order so that the clusters highly expressing that gene are together
	geneord <- summ[summ[,.I[which.max(mean_expr)],by=variable]$V1][order(RNA_snn_res.0.3)]$variable
	if (ensembl == FALSE){summ$variable <- factor(gsub('--.*', '', summ$variable), levels=gsub('--.*', '', geneord))
	} else {summ$variable <- factor(summ$variable, levels=geneord)}
	p_out <- ggplot(summ, aes(RNA_snn_res.0.3, variable, size=pct_pos, colour=mean_expr)) + geom_point() + theme_minimal() + theme(axis.line = element_line(), axis.ticks=element_line()) + scale_colour_gradientn(colours=c('grey90', 'cornflowerblue'))
	return(p_out)
}

# select known genes from the marker genes lists
selectgenes_names <- c('NANOG', 'NKX1-2', 'LEFTY1', 'GSC', 'MIXL1', 'LHX1', 'GATA3', 'PAX1', 'NKX2-3', 'NOG', 'SHH', 'T', 'AFP', 'APOB', 'TTR', 'MESP2', 'MESP1', 'HAND1', 'HOXB3', 'FOXC2', 'SNAI2', 'KRT7' , 'IGF1', 'NKX2-1', 'PRRX1', 'TTYH1', 'MYL3', 'MYH6', 'MYH7', 'CCNA1', 'CDH5', 'PECAM1', 'FLI1', 'ACTC1', 'DNAAF3', 'EOMES', 'FOXA2', 'NR2F2', 'SOX2')
selectgenes <- grep(paste(paste0('^', selectgenes_names, '--'), collapse='|'), row.names(mat), value = T)

# top 3 genes: select marker genes and top DE genes in each cluster by fold change and significance
top3_genes <- markers[p_val_adj==0&avg_log2FC>2,head(.SD,3),by='cluster']$gene
allgenes <- unique(c(selectgenes, top3_genes))

f4b <- bplt(allgenes)

##########################################
# F4c: kegg and GO output bubble plots
kgs <- readRDS(paste0(ddir, 'clust0.3_degokegg.RDS'))

clustcols <- setNames(hue_pal()(length(unique(kgs[[1]]$cluster))), unique(kgs[[1]]$cluster))

plt <- function(dt, n=2){
	topn <- dt[order(p.adjust), head(.SD,n),by=cluster][order(cluster, p.adjust)]
	topn$Description <- factor(topn$Description, levels=unique(topn$Description))
	topn$filtclust <- ifelse(topn$p.adjust < 0.05, as.character(topn$cluster), NA)
	topn$filtclust <- factor(topn$filtclust, levels=names(clustcols))
	p <- ggplot(topn, aes(cluster, Description, size=-log10(p.adjust), colour=filtclust)) + geom_point() + theme_minimal() + theme(axis.line = element_line()) + scale_colour_manual(values=clustcols, drop=FALSE, na.value='grey90')
	return(p)
}

f4c <- plt(kgs[[1]]) + plt(kgs[[2]])

##########################################
# Fig4d: proportions of each time point and treatment group in each cluster, with normalisation for number of cells captured at each timepoint
prop1 <- met[,.N,by=c('RNA_snn_res.0.3', 'day', 'condition')]
prop2 <- met[,.('totalday'=.N),by=day]
prop3 <- met[,.('totaldaycond'=.N),by=c('condition', 'day')]

props1 <- merge(prop1, prop2, by='day')
props2 <- merge(prop1, prop3, by=c('condition', 'day'))

props1 <- props1[,.(RNA_snn_res.0.3, day, 'pr_clust_per_day'=N/totalday)]
props2 <- props2[,.(RNA_snn_res.0.3, day, condition, 'pr_clust_per_daycond'=N/totaldaycond)]

f4d <- ggplot(props1, aes(RNA_snn_res.0.3, pr_clust_per_day, fill=day)) + geom_col(position='fill') + theme_classic() + scale_fill_manual(values=day_cols) + labs(y='Normalised proportion of cells by day') +
	   ggplot(props2, aes(RNA_snn_res.0.3, pr_clust_per_daycond, fill=condition)) + geom_col(position='fill') + theme_classic() + scale_fill_manual(values=treatment_cols) + labs(y='Normalised proportion of cells by day and treatment group')

