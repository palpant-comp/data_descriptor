library(data.table)
library(ggplot2)
library(patchwork)
library(ggrastr)
library(Seurat)
library(scales)
set.seed(7)

int_risc <- readRDS('int_risc.RDS')

# extract data from risc object (esp. PLS dimensionality reduction) to use as input for seurat clustering
int_meta <- cbind(as.data.table(int_risc@DimReduction$cell.umap), int_risc@coldata)
row.names(int_meta) <- int_meta$Barcode

pcs <- int_risc@DimReduction$cell.pls
colnames(pcs) <- gsub('PC', 'PC_', colnames(pcs))
ums <- int_risc@DimReduction$cell.umap
colnames(ums) <- gsub('UMAP', 'UMAP_', colnames(ums))

# put into seurat object
seu <- CreateSeuratObject(counts = int_risc@assay$logcount, meta.data=int_meta)
seu[['pca']] <- CreateDimReducObject(embeddings=pcs, stdev=apply(pcs, 2, sd), key='PC_', assay='RNA')
seu[['umap']] <- CreateDimReducObject(embeddings=ums, key='UMAP_', assay='RNA')

# run processing and clustering
seu <- FindVariableFeatures(seu)
seu <- FindNeighbors(seu, dims = 1:10) # because umap coords were generated using 10 pcs
for (res in seq(0.1, 0.4, 0.1)){seu <- FindClusters(seu, resolution=res)}

saveRDS(seu, paste0(ddir, 'int_clust0.3_seu.RDS')) 

#########################################################
# get marker genes

Idents(seu) <- 'RNA_snn_res.0.3'
markers <- FindAllMarkers(seu, only.pos=TRUE)
write.table(markers, 'T4_markergenes.txt'), quote = F, sep = '\t', col.names=T, row.names=F)

#########################################################
# do kegg and go with clusterprofiler
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(scales)

mat <- readRDS('int_matrix.RDS')

allensembl <- gsub('.*--', '', row.names(mat))
ids <- bitr(gsub('.*--', '', allensembl), fromType = 'ENSEMBL', toType = 'UNIPROT', OrgDb = 'org.Hs.eg.db') # convert ensembl ids to uniprot ids for kegg

dokg <- function(clust, method, n=100){
	mar <- markers[cluster==clust][order(p_val_adj,-avg_log2FC),head(.SD,n)]
	print(head(mar))
	genes <- gsub('.*--', '', mar$gene)
	if (method == 'go'){
		kg <- enrichGO(gene=genes, universe=allensembl, keyType='ENSEMBL', OrgDb=org.Hs.eg.db, ont='BP', pAdjustMethod='BH', pvalueCutoff=0.01, qvalueCutoff=0.05)
	} else if (method == 'kegg'){
		useids <- ids[grep(paste(genes,collapse = '|'), ids$ENSEMBL),]$UNIPROT
		print(paste('number of uniprot ids:', length(useids)))
		kg <- enrichKEGG(gene = useids, organism = 'hsa', pvalueCutoff = 0.05, keyType = 'uniprot')
	}
	out <- as.data.table(kg@result[,c('ID', 'Description', 'p.adjust', 'qvalue')])
	out$cluster <- clust
	return(out)
}

gos <- do.call(rbind, lapply(0:12, dokg, method='go'))
keggs <- do.call(rbind, lapply(0:12, dokg, method='kegg'))
gos$cluster <- factor(gos$cluster, levels = 0:12)
keggs$cluster <- factor(keggs$cluster, levels = 0:12)

saveRDS(list('go' = gos, 'kegg'=keggs), paste0(ddir, 'clust0.3_degokegg.RDS'))

write.table(gos, 'T4_clust0.3_gos.txt'), quote = F, sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(keggs, 'T4_clust0.3_keggs.txt'), quote = F, sep = '\t', col.names = TRUE, row.names = FALSE)

