##----Install libraries----
libraries <- rlang::quos(clusterProfiler, DOSE, ggplot2, enrichplot,org.Ce.eg.db,
                         ReactomePA, forcats, dplyr, ggstance, MetBrewer)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Create gene list----
d <- read.csv("results/filtered_treatment1_vs_control.csv")
genelist <- d[,3]
names(genelist) <- as.character(d[,10])
genelist <- na.omit(genelist) 
genelist <- sort(genelist, decreasing = TRUE)

#----Gene Ontology----
gse <- gseGO(geneList = genelist,
             ont = "BP",
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Ce.eg.db, 
             pAdjustMethod = "none",
             eps = 0) 

#write.csv(gse@result, file = "WT_gsea.csv")

categories <- c("response to stress",
                "response to external stimulus",
                "cell communication",
                "signaling",
                "signal transduction",
                "regulation of lipid localization",
                "regulation of biological process",
                "fatty acid metabolic process",
                "lipid localization",
                "lipid storage",
                "posttranscriptional regulation of gene expression",
                "posttranscriptional gene silencing",
                "gene silencing by RNA",
                "regulation of transcription by RNA polymerase II",
                "regulation of nucleic acid-templated transcription",
                "negative regulation of transcription, DNA-templated",
                "RNA biosynthetic process"
                )

col <- met.brewer("Pissaro")

##----Dot plot for enrichment result----
dotplot(gse, showCategory = categories, split = ".sign", orderBy = "x" ) + 
  facet_grid(.~.sign) + labs(title = "WT (2°C)") + 
  scale_color_gradient(low = col[5], high = "grey")

##----Enrichment analysis----
genes <- names(genelist)
ego <- enrichGO(gene = genes,
                OrgDb = org.Ce.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

##----Emap Plots----
pt <- pairwise_termsim(ego,
                       method = "JC",
                       semData = NULL,
                       showCategory = 200)
pt <- simplify(pt, cutoff=0.7, by="p.adjust", select_fun=min)

emapplot(pt, color = "p.adjust", node_label = "category", shadowtext = TRUE,
         clusterFunction = cluster::pam, split = "category", pie = "Count",
         repel = TRUE) + 
  labs(title = "WT (2°C)") + 
  scale_color_gradient(low = col[5], high = "grey")

##----Gene Ontology----
edo <- gseGO(geneList = genelist,
             ont = "BP",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = org.Ce.eg.db, 
             pAdjustMethod = "none",
             eps = 0) 

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Ce.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
#edox2<- simplify(edox2, cutoff=0.7, by="p.adjust", select_fun=min)

treeplot(edox2, showCategory = 50) + 
  scale_color_continuous(low = col[5], high = "grey") + 
  labs(title = "WT (2°C)")

gseaplot2(edo, geneSetID = 61, pvalue_table = TRUE, subplots = 1:2)

sessionInfo()