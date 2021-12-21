##----Install libraries----
libraries <- rlang::quos(ggplot2, 
                         DESeq2, 
                         AnnotationDbi, 
                         org.Ce.eg.db, 
                         pheatmap, 
                         genefilter, 
                         RColorBrewer, 
                         rafalib, 
                         MetBrewer, 
                         EnhancedVolcano,
                         vsn)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Results----
res <- results(dds, contrast = c("test", "wt_cs", "wt"), lfcThreshold = 0.5,
               altHypothesis = "greaterAbs")

summary(res)

res$symbol <- mapIds(org.Ce.eg.db,
                        keys=row.names(res),
                        column="SYMBOL",
                        keytype="WORMBASE",
                        multiVals="first")
res$name <- mapIds(org.Ce.eg.db,
                      keys=row.names(res),
                      column="GENENAME",
                      keytype="WORMBASE",
                      multiVals="first")
res$entrez <- mapIds(org.Ce.eg.db,
                        keys=row.names(res),
                        column="ENTREZID",
                        keytype="WORMBASE",
                        multiVals="first")
res$GO_term <- mapIds(org.Ce.eg.db,
                     keys=row.names(res),
                     column="GO",
                     keytype="WORMBASE",
                     multiVals="first")

res <- res[order(res$padj), ]
hist(res$pvalue[res$baseMean > 1],
     col = "grey", border = "white", xlab = "", ylab = "", main = "")

res <- na.omit(res)
res <- subset(res, pvalue < 0.05)
write.csv(as.data.frame(res), file = "results/filtered_treatment1_vs_control.csv")

rownames(res) <- mapIds(org.Ce.eg.db,
                        keys=row.names(res),
                        column="SYMBOL",
                        keytype="WORMBASE",
                        multiVals="first")

col <- met.brewer("Pissaro")
col <- c(col[4], col[4], col[4], col[5])

WT.CSvsNone<-EnhancedVolcano(res, lab = rownames(res),
                #selectLab = c('ire-1','pek-1', "atf-6", "jnk-1", "kgb-1", "kgb-2", "nlp-3", "fat-5", "fat-6","xbp-1"),
                x = 'log2FoldChange', title = "WT", pCutoff = tail(res[res$padj < 0.1, ])$pvalue[1] + 1e-10,
                y = 'pvalue', legendPosition = "right", FCcutoff = 0.5,
                drawConnectors = TRUE, boxedLabels = TRUE,xlim = c(-6,8), ylim = c(0, 21),
                col = col,
                subtitle = "", caption = "") +
  theme(panel.background = element_rect(colour = "black", fill = "transparent", size = 0.6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) #500x400

sessionInfo()