gse <- getGEO("GSE81622")

gse1em <- getGEO("GSE81622", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse <- gse1em[[1]]

exprs_data <- exprs(gse)
exprs_data <- log2(exprs_data + 1)
exprs_data <- normalizeBetweenArrays(exprs_data)

metadata <- pData(gse)

group <- ifelse(grepl("control", metadata$title, ignore.case = TRUE), "control", "SLE")
group <- factor(group)
table(group)  # grup sayısını kontrol et

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)

contrast_matrix <- makeContrasts(SLE - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)

probe_ids <- rownames(results)

results$GeneSymbol <- mapIds(illuminaHumanv4.db, 
                             keys = probe_ids, 
                             column = "SYMBOL", 
                             keytype = "PROBEID", 
                             multiVals = "first")

deg <- subset(results, abs(logFC) > 0.05 & adj.P.Val < 0.05)

deg <- deg %>%
  filter(!is.na(deg$GeneSymbol))


# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

# Heatmap (ilk 50 DEG)
deg_ids <- rownames(deg)[1:50]
pheatmap(exprs_data[deg_ids, ], scale = "row")

write_xlsx(deg, "deg_GSE81622_SLE.xlsx")


up_genes <- deg %>% filter(logFC > 1)
down_genes <- deg %>% filter(logFC < -1)




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

write_xlsx(as.data.frame(go_enrich), "deg_GSE81622_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "deg_GSE81622_KEGG_enrichment.xlsx")





