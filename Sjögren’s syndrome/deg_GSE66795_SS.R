gse <- getGEO("GSE66795")

gse1em <- getGEO("GSE66795", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse <- gse1em[[1]]

exprs_data <- exprs(gse)
exprs_data <- normalizeBetweenArrays(exprs_data)

metadata <- pData(gse)

group <- ifelse(grepl("control", metadata$characteristics_ch1.2, ignore.case = TRUE), "control", "patient")
group <- factor(group)
table(group)  # grup sayısını kontrol et

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)

# Grup karşılaştırması
contrast_matrix <- makeContrasts(patient - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)
head(results)

probe_ids <- rownames(results)

results$GeneSymbol <- mapIds(illuminaHumanv4.db, 
                             keys = probe_ids, 
                             column = "SYMBOL", 
                             keytype = "PROBEID", 
                             multiVals = "first")

deg <- subset(results, abs(logFC) > 0.05 & adj.P.Val < 0.05)
deg <- deg %>%
  filter(!is.na(deg$GeneSymbol))
head(deg)


results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value")

# Heatmap (ilk 50 DEG)
top50 <- deg[order(deg$adj.P.Val), ][1:50, ]
deg_ids <- rownames(top50)

heatmap_data <- exprs_data[deg_ids, ]

heatmap_data <- t(scale(t(heatmap_data)))

pheatmap(heatmap_data, 
         scale = "none", 
         show_rownames = TRUE, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Top 50 DEG Heatmap")

write_xlsx(deg, "deg_GSE66795_SS.xlsx")


logfc_threshold <- 1
pval_threshold <- 0.05

# Upregulated genler (pozitif logFC)
upregulated_genes <- subset(deg, logFC > logfc_threshold & adj.P.Val < pval_threshold)

# Downregulated genler (negatif logFC)
downregulated_genes <- subset(deg, logFC < -logfc_threshold & adj.P.Val < pval_threshold)

genes <- results$GeneSymbol[!is.na(results$GeneSymbol)]



go_results <- enrichGO(gene         = genes,
                       OrgDb        = org.Hs.eg.db,
                       keyType      = "SYMBOL",
                       ont          = "BP", "MF", "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

head(go_enrich)

#Barplot
barplot(go_enrich, showCategory = 20, title = "GO BP Enrichment")

#kegg enrichment

entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

entrez_ids <- entrez_ids[!is.na(entrez_ids)]

kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(kegg_enrich)

# Barplot
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(go_enrich), "deg_GSE66795_SS_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "deg_GSE66795_SS_KEGG_enrichment.xlsx")




