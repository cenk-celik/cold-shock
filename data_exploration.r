 ##----Install libraries----
libraries <- rlang::quos(ggplot2, DESeq2, AnnotationDbi, org.Ce.eg.db, pheatmap,
                         genefilter, rafalib, viridis,
                         EnhancedVolcano, MetBrewer)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Load saved featureCounts----
load(paste0(rnaLocal, "Rdata/featurecounts2.RData"))

countMatrix <- featurecounts$counts
colnames(countMatrix) <- gsub(".bam$", "", colnames(countMatrix))
coldata <- read.csv("coldata.csv", header = TRUE)
coldata <- coldata[order(coldata$test), ]
countMatrix <- countMatrix[, coldata$X]

##----Differential Expression Analysis----
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = coldata,
                              design = ~ test + batch)

##----Pre-filtering----
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]

dds <- DESeq(dds)
colnames(dds) <- coldata$sample # rename matrix columns
  
##----normalisation----
rld <- rlog(dds)
assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)

#----Heat maps---
topVarGenes <- head(order(-rowVars(assay(rld)), decreasing = T), 10000)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- mapIds(org.Ce.eg.db,
                        keys=row.names(mat),
                        column="SYMBOL",
                        keytype="WORMBASE",
                        multiVals="first")

df <- as.data.frame(colData(rld)[,c("test","mutant")])
colnames(df) <- c("Test", "Mutant")

col <- met.brewer("Pissaro")

cols <- list(Mutant = c(ire.ok799 = col[1], wt = col[6]),
             Test = c(ire.ok799_cs = col[2], ire.ok799 = col[3], wt_cs = col[5], wt = col[4]))

pheatmap(mat, cluster_rows = T, show_rownames = F,# kmeans_k = 999,
         cluster_cols = F, annotation_col = df, annotation_names_col = F, annotation_names_row = F,
         angle_col = 45, color = met.brewer("Hokusai", n = 1000, type = "continuous"), gaps_col = c(6), annotation_colors = cols) # 650 x 700

##----Principle Component Analysis----
pca <- plotPCA(rld, intgroup = c("test","mutant"), ntop = 5000, returnData = TRUE)
jitter <- c(rep(col[2],3), rep(col[3],3), rep(col[5], 6), rep(col[4],6))
percentVar <- round(100 * attr(pca, "percentVar"))

ggplot(pca, aes(PC1, PC2, fill = group)) + 
  geom_point(size=4, color = jitter) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_fill_manual(labels=c("ire1(-) (2°C)", "ire1(-)", "WT (2°C)", "WT"), values = c(col[2],col[3],col[5],col[4])) +
  coord_fixed() + theme_minimal() 

##----Session info----
sessionInfo()