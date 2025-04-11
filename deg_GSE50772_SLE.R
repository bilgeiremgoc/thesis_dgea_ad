gse <- getGEO("GSE50772")

gse1em <- getGEO("GSE50772", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])

group <- factor(c(rep("SLE", 61), rep("control", 20)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)

boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

fit <- lmFit(exprs_data, design)

# Grup karşılaştırması
contrast_matrix <- makeContrasts(SLE - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Sonuçları al

results <- topTable(fit2, adjust = "fdr", number = Inf)

results$GeneSymbol <- getSYMBOL(rownames(results), "hgu133plus2.db")
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


write.csv(deg, "deg_GSE50772_SLE.csv")

