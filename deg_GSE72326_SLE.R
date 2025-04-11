gse <- getGEO("GSE72326")

gse1em <- getGEO("GSE72326", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse <- gse1em[[1]]

pheno_data <- pData(gse)

head(pheno_data)

# Grup bilgisinin bulunduğu sütunu kontrol etmek
table(pheno_data$characteristics_ch1)

# Sadece SLE ve ona ait kontrolleri filtreleme
selected_pheno <- pheno_data[pheno_data$`characteristics_ch1` %in% c("group: SLE", "group: Healthy control of SLE"), ]

# Ekspresyon matrisini bu örneklerle sınırlama
expr_matrix <- exprs(gse)
filtered_expr_matrix <- expr_matrix[, rownames(selected_pheno)]

group_labels <- selected_pheno$characteristics_ch1
group_labels <- factor(ifelse(group_labels == "group: SLE", "SLE", "Control"))

design <- model.matrix(~0 + group_labels)
colnames(design) <- levels(group_labels)
design

boxplot(expr_matrix, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(expr_matrix)

boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")


fit <- lmFit(filtered_expr_matrix, design)


# Grup karşılaştırması
contrast_matrix <- makeContrasts(SLEvsControl = SLE - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Sonuçları al

results <- topTable(fit2, adjust = "fdr", number = Inf)

results$GeneSymbol <- getSYMBOL(rownames(results), "hgdb")
head(results)

deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)
head(deg)

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg)[1:50], ], scale = "row")

gpl <- getGEO("GPL10558", AnnotGPL = TRUE)
gpl_table <- Table(gpl)


head(colnames(gpl_table)) 

probe2symbol <- gpl_table[, c("ID", "Gene Symbol")]
colnames(probe2symbol) <- c("ID", "GeneSymbol")

deg_results$ID <- rownames(deg_results)

# Merge işlemi
deg_annotated <- merge(deg_results, probe2symbol, by = "ID", all.x = TRUE)

# GeneSymbol.x sütununu yeniden adlandır
deg_annotated$GeneSymbol <- deg_annotated$GeneSymbol.x

deg_annotated$GeneSymbol.x <- NULL
deg_annotated$GeneSymbol.y <- NULL

deg_annotated

# Artık gen sembolleri deg_annotated içinde
head(deg_annotated[, c("ID", "GeneSymbol", "logFC", "adj.P.Val")])


write.csv(deg, "deg_GSE50772_SLE.csv")

write.csv(deg_annotated, "DEG_results_with_symbols.csv")


