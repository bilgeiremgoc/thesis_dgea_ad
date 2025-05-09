gse <- getGEO("GSE51092")

gse1em <- getGEO("GSE51092", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])

gse_data <- gse[[1]]

colnames(exprs(gse_data))  # örnek adları
pData(gse_data)[, c("title", "source_name_ch1")]  # grup bilgisi var mı?

metadata <- pData(gse_data)
group <- ifelse(grepl("control", metadata$title, ignore.case = TRUE), "control", "case")
group <- factor(group)
table(group)  # grup sayısını kontrol et


design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)

boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

fit <- lmFit(exprs_data, design)

# Grup karşılaştırması
contrast_matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


# Sonuçları al

gpl <- getGEO("GPL6884", AnnotGPL = TRUE)
gpl_table <- Table(gpl)


head(colnames(gpl_table)) 

probe2symbol <- gpl_table[, c("ID", "Gene symbol")]
colnames(probe2symbol) <- c("ID", "Gene symbol")

results$ID <- rownames(results)


library(limma)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)

# DEG seç
deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)

# Merge işlemi
deg_annotated <- merge(results, probe2symbol, by = "ID", all.x = TRUE)


deg_annotated$GeneSymbol <- NULL

deg_annotated

deg <- subset(deg_annotated, abs(logFC) > 1 & adj.P.Val < 0.05)
head(deg)

results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

top_genes <- rownames(deg)[1:50]
pheatmap(exprs_data[top_genes, ], scale = "row", annotation_col = data.frame(Group = group))

write.csv(deg, "deg_GSE51092_SS.csv")
