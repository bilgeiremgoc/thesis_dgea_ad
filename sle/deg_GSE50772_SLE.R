
gse <- getGEO("GSE50772", GSEMatrix = TRUE)
gse_data <- gse[[1]]

exprs_data <- exprs(gse_data)
exprs_data <- log2(exprs_data + 1)
exprs_data <- normalizeBetweenArrays(exprs_data)

metadata <- pData(gse_data)

group <- ifelse(grepl("control", metadata$characteristics_ch1, ignore.case = TRUE), "control", "SLE")
group <- factor(group)
table(group)  # grup sayısını kontrol et

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(SLE - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)

results$GeneSymbol <- mapIds(hgu133plus2.db, 
                             keys = rownames(results), 
                             column = "SYMBOL", 
                             keytype = "PROBEID", 
                             multiVals = "first")

deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)
deg <- na.omit(deg)
head(deg)

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme(legend.position = "none")

# Heatmap (ilk 50 DEG)
deg_ids <- rownames(deg)[1:50]
pheatmap(exprs_data[deg_ids, ], scale = "row")


write_xlsx(deg, "deg_GSE50772_SLE.xlsx")


deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)


deg_filtered$ENTREZID <- mapIds(hgu133a.db,
                                keys = rownames(deg_filtered),
                                column = "ENTREZID",
                                keytype = "PROBEID",
                                multiVals = "first")


entrez_ids <- na.omit(deg_filtered$ENTREZID)

go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

head(go_enrich)

#Barplot
barplot(go_enrich, showCategory = 20, title = "GO BP Enrichment")

#kegg enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Entrez ID'den gen adlarına dönüştür (okunabilirlik)
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# İlk sonuçlara bak
head(kegg_enrich)

# Barplot
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(go_enrich), "GSE50772_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "GSE50772_KEGG_enrichment.xlsx")



