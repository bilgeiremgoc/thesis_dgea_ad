gse <- getGEO("GSE17449")

gse1em <- getGEO("GSE17449", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse_data <- gse1em[[1]]

exprs_data <- exprs(gse_data)
boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")
exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

metadata <- pData(gse_data)
group <- ifelse(grepl("healthy", metadata$`disease state:ch1`, ignore.case = TRUE), "healthy", "MS")
group <- factor(group)
table(group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(MS - healthy, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)
results$GeneSymbol <- getSYMBOL(rownames(results), "hgu133plus2.db")
head(results)

deg <- subset(results, abs(results$logFC) > 1 & results$adj.P.Val < 0.05)
deg <- deg %>% filter(!is.na(deg$GeneSymbol))
head(deg)

write_xlsx(deg, "deg_61635_SLE.xlsx")


deg_table <- read.delim(file = "GSE17449.top.table.tsv")
deg <- subset(deg_table, abs(logFC) > 0.5 & adj.P.Val < 0.05)
deg <- deg %>% filter(!is.na(deg$GeneSymbol))
head(deg)






























#GSE111972

deg_table_GSE111972 <- read.delim("GSE111972.top.table.tsv")

deg <- subset(deg_table_GSE111972,log2FoldChange > 1 & padj < 0.5)
deg <- deg %>% filter(!is.na(deg$Symbol))
head(deg)


write_xlsx(deg, "deg_GSE111972_MS.xlsx")
