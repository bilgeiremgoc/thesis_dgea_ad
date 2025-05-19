#Data Alımı ve Normalizasyonu
gse <- getGEO("GSE138198")

gse1em <- getGEO("GSE138198", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse<- gse1em[[1]]

exprs_data <- exprs(gse)

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")
exprs_data <- log2(exprs_data + 1) 
exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

#Gruplandırma ve DGEA
metadata <- pData(gse)
head(metadata)

group <- ifelse(metadata$source_name_ch1== "Thyroid, normal histology", "control", 
                ifelse(metadata$source_name_ch1 == "HT histology", "HT", NA))
group <- factor(group)
table(group)

valid_samples <- !is.na(group)
exprs_data <- exprs_data[, valid_samples]
group <- group[valid_samples]

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(HT - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)
if (!require("hugene10sttranscriptcluster.db")) {
  BiocManager::install("hugene10sttranscriptcluster.db")
}
library(hugene10sttranscriptcluster.db)

results$gene_symbol <- mapIds(hugene10sttranscriptcluster.db,
                              keys = rownames(results),
                              column = "SYMBOL",
                              keytype = "PROBEID",
                              multiVals = "first")
head(results)

deg <- subset(results, abs(logFC) > 0.5 & adj.P.Val < 0.05)
deg <- deg %>% filter(!is.na(deg$gene_symbol))
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

write_xlsx(deg, "deg_GSE138198_HT.xlsx")


#GO, KEGG ve Reactome Pathway Enrichment
genes <- results$gene_symbol[!is.na(results$gene_symbol)]


go_results_BP <- enrichGO(gene         = genes,
                       OrgDb        = org.Hs.eg.db,
                       keyType      = "SYMBOL",
                       ont          = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

head(go_results_BP)
barplot(go_results_BP, showCategory = 20, title = "GO BP Enrichment")
write_xlsx(as.data.frame(go_results_BP), "deg_GSE138198_HT_GO_BP_enrichment.xlsx")

go_results_MF <- enrichGO(gene         = genes,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.2)

head(go_results_MF)
barplot(go_results_MF, showCategory = 20, title = "GO MF Enrichment")
write_xlsx(as.data.frame(go_results_MF), "deg_GSE138198_HT_GO_MF_enrichment.xlsx")

go_results_CC <- enrichGO(gene         = genes,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.2)

head(go_results_CC)
barplot(go_results_CC, showCategory = 20, title = "GO CC Enrichment")
write_xlsx(as.data.frame(go_results_CC), "deg_GSE138198_HT_GO_CC_enrichment.xlsx")

#KEGG enrichment
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
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(kegg_enrich), "deg_GSE138198_HT_KEGG_enrichment.xlsx")

reactome_enrich <- enrichPathway(gene = entrez_ids, 
                                 organism = "human", 
                                 pvalueCutoff = 0.1, 
                                 qvalueCutoff = 0.2)
barplot(reactome_enrich, showCategory = 20, title = "Reactome Pathway Enrichment")

write_xlsx(as.data.frame(reactome_enrich), "deg_GSE138198_HT_reactome_enrichment.xlsx")


# Bağışıklık Hücresi İnfiltrasyonu (xCell)
gene_symbols <- fData(gse)$`Gene symbol`

exprs_data_with_symbols <- exprs_data[!is.na(gene_symbols) & gene_symbols != "", ]
rownames(exprs_data_with_symbols) <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
exprs_data_with_symbols <- as.data.frame(exprs_data_with_symbols)
head(exprs_data_with_symbols)

# xCell analizi
xcell_results <- xCellAnalysis(exprs_data_with_symbols)

xcell_df <- as.data.frame(xcell_results)
xcell_df <- cbind(Cell_Type = rownames(xcell_results), xcell_df)

# Excel olarak kaydetme
write_xlsx(xcell_df, "xcell_results_deg_GSE138198_HT.xlsx")


# LogFC ve adj.P.Val eşiklerine göre filtreleme
logfc_threshold <- 1
adj_pval_threshold <- 0.05

deg_filtered <- results %>%
  filter(!is.na(gene_symbol)) %>%
  filter(abs(logFC) > logfc_threshold & adj.P.Val < adj_pval_threshold)


# Up ve down gen listelerini ayırma
up_genes <- deg_filtered$gene_symbol[deg_filtered$logFC > 0]
down_genes <- deg_filtered$gene_symbol[deg_filtered$logFC < 0]

print(dim(deg_filtered))  # Gen sayısını görme
table(deg_filtered$logFC > 0)  # Up ve down genlerin sayısını görme
