if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "affy", "limma", "annotate", "hgu133plus2.db"))

library(GEOquery)
library(limma)
library(annotate)
library(ggplot2)
library(pheatmap)
library(hgu133a.db)
library(affy)
library(hgu133plus2.db)

gse <- getGEO("GSE61635")

gse1em <- getGEO("GSE61635", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse_data <- gse1em[[1]]

exprs_data <- exprs(gse_data)
exprs_data <- log2(exprs_data + 1)
exprs_data <- normalizeBetweenArrays(exprs_data)

metadata <- pData(gse_data)
group <- ifelse(grepl("healthy", metadata$characteristics_ch1, ignore.case = TRUE), "healthy", "SLE")
group <- factor(group)
table(group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(gse_data, design)
contrast_matrix <- makeContrasts(SLE - healthy, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


results <- topTable(fit2, adjust = "fdr", number = Inf)

results$GeneSymbol <- getSYMBOL(rownames(results), "hgu133a.db")


deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)
deg_clean <- na.omit(deg$GeneSymbol)
head(deg)
deg_clean <- as.data.frame(deg_clean)

# Volkan Plot
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

# Heatmap (ilk 50 DEG)
deg_ids <- rownames(deg)[1:50]
pheatmap(exprs_data[deg_ids, ], scale = "row")

write_xlsx(deg_clean, "deg_GSE61635_SLE.xlsx")

deg_filtered <- deg_clean %>% filter(!is.na(deg_clean))
                                     
deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1))

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter( logFC> 1)
down_genes <- deg_filtered %>% filter(logFC < -1)


deg$ENTREZID <- mapIds(
  hgu133a.db,
  keys = deg$ID,
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)
                               

# NA olmayanları al
entrez_ids <- deg$ENTREZID

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

write_xlsx(as.data.frame(go_enrich), "deg_GSE61635_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "deg_GSE61635_KEGG_enrichment.xlsx")




