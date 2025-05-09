gse <- getGEO("GSE66795")

gse1em <- getGEO("GSE66795", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])

colnames(exprs(gse_data))  # örnek adları
pData(gse_data)[, c("title", "source_name_ch1")]  # grup bilgisi var mı?

metadata <- pData(gse_data)
group <- ifelse(grepl("control", metadata$title, ignore.case = TRUE), "control", "case")
group <- factor(group)
table(group)  # grup sayısını kontrol et

