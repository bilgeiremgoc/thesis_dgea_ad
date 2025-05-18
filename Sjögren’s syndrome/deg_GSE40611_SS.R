gse <- getGEO("GSE40611")

gse1em <- getGEO("GSE40611", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse<- gse1em[[1]]
  
exprs_data <- exprs(gse)

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")
exprs_data <- log2(exprs_data + 1) 
exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

metadata <- pData(gse)
head(metadata)

group <- ifelse(grepl("control", metadata$`disease status:ch1`, ignore.case = TRUE), "control", "SS")
group <- factor(group)
table(group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(SS - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)
results$gene_symbol <- getSYMBOL(rownames(results), "hgu133plus2.db")
head(results)

deg <- subset(results, abs(logFC) > 0.5 & adj.P.Val < 0.05)
head(deg)


deg <- deg %>%
  filter(!is.na(deg$gene_symbol))
head(deg)


# Volkan Plotu

results$threshold <- as.factor(abs(results$logFC) > 0.5 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "logFC",
       y = "-log10(adjusted P-value)",
       color = "Significant")



# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg)[1:14], ], scale = "row")


write_xlsx(deg, "deg_GSE40611_SS.xlsx")


genes <- results$gene_symbol[!is.na(results$gene_symbol)]


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

write_xlsx(as.data.frame(go_enrich), "deg_GSE40611_SS_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "deg_GSE40611_SS_KEGG_enrichment.xlsx")

reactome_enrich <- enrichPathway(entrez_ids, organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
barplot(reactome_enrich, showCategory = 20, title = "Reactome Pathway Enrichment")

write_xlsx(as.data.frame(reactome_enrich), "deg_GSE40611_SS_reactome_enrichment.xlsx")


# Bağışıklık Hücresi İnfiltrasyonu (xCell)
exprs_data_with_symbols <- exprs_data[!is.na(gene_symbol), ]

# Gen sembollerini satır isimlerine ekleyelim
rownames(exprs_data_with_symbols) <- gene_symbols[!is.na(gene_symbols)]
head(exprs_data_with_symbols)

xcell_df <- as.data.frame(xcell_results)

# Satır isimlerini ayrı bir kolon olarak ekleme (örn. hücre tipleri)
xcell_df <- cbind(Cell_Type = rownames(xcell_results), xcell_df)

# Excel olarak kaydetme
write_xlsx(xcell_df, "xcell_results_GSE84844.xlsx")

write_xlsx(xcell_results, "xcell_results_GSE84844.xlsx")



# LogFC ve adj.P.Val eşiklerine göre filtreleme
logfc_threshold <- 1
adj_pval_threshold <- 0.05

deg_filtered <- results %>%
  filter(!is.na(gene_symbol)) %>%
  filter(abs(logFC) > logfc_threshold & adj.P.Val < adj_pval_threshold)


# Up ve down gen listelerini ayır
up_genes <- deg_filtered$gene_symbol[deg_filtered$logFC > 0]
down_genes <- deg_filtered$gene_symbol[deg_filtered$logFC < 0]

write.table(up_genes, "up_genes_GSE84844.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(down_genes, "down_genes_GSE84844.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

print(dim(deg_filtered))  # Gen sayısını gör
table(deg_filtered$logFC > 0)  # Up ve down genlerin sayısını gör








