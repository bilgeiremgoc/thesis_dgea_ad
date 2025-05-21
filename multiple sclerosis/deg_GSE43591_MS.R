gse <- getGEO("GSE43591")

gse1em <- getGEO("GSE43591", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse <- gse1em[[1]]

exprs_data <- exprs(gse)
boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")
exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

metadata <- pData(gse)

group <- ifelse(grepl("control", metadata$characteristics_ch1, ignore.case = TRUE), "control", "MS")
group <- factor(group)
table(group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(MS - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)
results$GeneSymbol <- getSYMBOL(rownames(results),"hgu133plus2.db")
head(results)

deg <- subset(results, abs(logFC) > 0.5 & adj.P.Val < 0.05)
deg <- deg %>% filter(!is.na(deg$GeneSymbol))
head(deg)

write_xlsx(deg, "deg_GSE43591_MS_05.xlsx")


deg_filtered <- deg %>% filter(!is.na(deg$GeneSymbol))

deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)
nrow(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC> 0.5)
down_genes <- deg_filtered %>% filter(logFC < -0.5)

write_xlsx(up_genes, "deg_GSE43591_MS_up_genes_05.xlsx")
write_xlsx(down_genes, "deg_GSE43591_MS_down_genes_05.xlsx")

# Volkan Plot
results$threshold <- as.factor(abs(results$logFC) > 0.5 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

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

deg$ENTREZID <- mapIds(
  hgu133plus2.db,
  keys = rownames(deg),
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)


entrez_ids <- deg$ENTREZID

#GO_BP
go_result_BP <- enrichGO(gene = entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

head(go_result_BP)
barplot(go_result_BP, showCategory = 20, title = "GO BP Enrichment")
write_xlsx(as.data.frame(go_result_BP), "deg_GSE43591_MS_GO_BP_enrichment.xlsx")

#GO_MF
go_result_MF<- enrichGO(gene = entrez_ids,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

head(go_result_MF)
barplot(go_result_MF, showCategory = 20, title = "GO MF Enrichment")
write_xlsx(as.data.frame(go_result_MF), "deg_GSE43591_MS_GO_MF_enrichment.xlsx")

#GO_CC
go_result_CC <- enrichGO(gene = entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "CC",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

head(go_result_CC)
barplot(go_result_CC, showCategory = 20, title = "GO CC Enrichment")
write_xlsx(as.data.frame(go_result_CC), "deg_GSE43591_MS_GO_CC_enrichment.xlsx")

#KEGG enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(kegg_enrich)
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(kegg_enrich), "deg_GSE43591_MS_KEGG_enrichment.xlsx")

#Reactome Enrichment
reactome_enrich <- enrichPathway(gene = entrez_ids, 
                                 organism = "human", 
                                 pvalueCutoff = 0.1, 
                                 qvalueCutoff = 0.2)
barplot(reactome_enrich, showCategory = 20, title = "Reactome Pathway Enrichment")

write_xlsx(as.data.frame(reactome_enrich), "deg_GSE43591_MS_Reactome_enrichment.xlsx")







#xCell İnfiltrasyonu
gene_symbols <- getSYMBOL(rownames(exprs_data), "hgu133plus2.db")
exprs_data_symbols <- exprs_data
rownames(exprs_data_symbols) <- gene_symbols

exprs_data_symbols <- exprs_data_symbols[!is.na(rownames(exprs_data_symbols)), ]
exprs_data_symbols <- rowsum(exprs_data_symbols, group = rownames(exprs_data_symbols))

xcell_result <- xCellAnalysis(exprs_data_symbols)

xcell_df <- as.data.frame(xcell_result)
xcell_df <- cbind(Cell_Type = rownames(xcell_result), xcell_df)

pheatmap(xcell_result, 
         main = "xCell Immune Infiltration - GSE43591",
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

write_xlsx(xcell_df, "deg_GSE43591_MS_SLE_xCell_Results.xlsx")




xcell_result <- xCellAnalysis(exprs_data_symbols)

xcell_matrix <- as.data.frame(t(xcell_result))
xcell_matrix$sample <- rownames(xcell_matrix)
xcell_matrix$group <- group 
head(xcell_matrix)

long_xcell <- xcell_matrix %>%
  pivot_longer(cols = -c(sample, group), names_to = "cell_type", values_to = "proportion")


stat_tests <- long_xcell %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(proportion ~ group)$p.value,
    mean_SLE = mean(proportion[group == "SLE"]),
    mean_control = mean(proportion[group == "control"])
  ) %>%
  mutate(adj_p = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)

write_xlsx(stat_tests, "deg_GSE43591_MS_xcell_stat_tests.xlsx")


ggplot(long_xcell, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Immune Cell Infiltration (xCell) - GSE43591") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






