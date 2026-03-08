# =============================================================================
# Computational Analysis of Gene Expression in Alzheimer's Disease
# Dataset: GSE33000 (NCBI GEO) — Prefrontal Cortex Microarray
# 310 AD + 157 Controls = 467 samples | 25,929 probes
# Author: theaxonaut
# =============================================================================

# --- 1. PACKAGES --------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages("pak", quiet = TRUE)
pak::pkg_install(c(
  "bioc::GEOquery",
  "bioc::limma",
  "bioc::EnhancedVolcano",
  "bioc::clusterProfiler",
  "bioc::org.Hs.eg.db",
  "bioc::enrichplot"
))

library(GEOquery)
library(limma)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# --- 2. LOAD DATA -------------------------------------------------------------

gse <- getGEO("GSE33000", GSEMatrix = TRUE, getGPL = FALSE)

expr <- exprs(gse[[1]])
meta <- pData(gse[[1]])
feat <- fData(gse[[1]])

# --- 3. FILTER SAMPLES --------------------------------------------------------

meta_filt <- meta[meta$`disease state:ch1` %in% c("Alzheimer's disease", "non-demented"), ]

meta_filt$disease_state <- factor(
  meta_filt$`disease state:ch1`,
  levels = c("non-demented", "Alzheimer's disease")
)

table(meta_filt$disease_state)

# --- 4. FILTER PROBES ---------------------------------------------------------

# Note: gene name column is 'ORF' — no Gene Symbol column in this dataset
feat_filt <- feat[!is.na(feat$Gene.ID) & feat$Gene.ID != "", ]

# --- 5. ALIGN EXPRESSION MATRIX -----------------------------------------------

expr_matched <- expr[rownames(expr) %in% rownames(feat_filt),
                     colnames(expr) %in% rownames(meta_filt)]

expr_matched <- expr_matched[, rownames(meta_filt)]

# Verify alignment — must return TRUE
stopifnot(all(colnames(expr_matched) == rownames(meta_filt)))

dim(expr_matched) # 25929 x 467

# --- 6. DIFFERENTIAL EXPRESSION — limma ---------------------------------------

design <- model.matrix(~ disease_state, data = meta_filt)

fit  <- lmFit(expr_matched, design)
fit2 <- eBayes(fit)

# coef=2 = disease_stateAlzheimer's disease coefficient
# positive logFC = higher in AD relative to non-demented
results <- topTable(fit2, coef = 2, number = Inf, sort.by = "logFC")

# --- 7. ANNOTATE RESULTS ------------------------------------------------------

results_annot <- merge(
  results,
  feat_filt[, c("ID", "ORF", "Gene.ID")],
  by.x = "row.names", by.y = "ID",
  all.x = TRUE
)

colnames(results_annot)[colnames(results_annot) == "Row.names"] <- "ID"
colnames(results_annot)[colnames(results_annot) == "Gene.ID"]   <- "EntrezGeneID"

results_annot <- results_annot[order(results_annot$logFC, decreasing = TRUE), ]

# --- 8. VOLCANO PLOT ----------------------------------------------------------

# FCcutoff = 0.5 (not standard 1.0) — brain microarray has small fold changes
dir.create("figures", showWarnings = FALSE)

EnhancedVolcano(
  results_annot,
  lab      = results_annot$ORF,
  x        = "logFC",
  y        = "adj.P.Val",
  pCutoff  = 0.05,
  FCcutoff = 0.5,
  title    = "GSE33000: AD vs Non-Demented",
  subtitle = "Prefrontal Cortex | limma + eBayes"
)

ggsave("figures/volcano_plot.png", width = 10, height = 8, dpi = 200)

# Note: XIST upregulation is a sex imbalance artefact — not a disease signal

# --- 9. EXPORT RESULTS --------------------------------------------------------

write.csv(results_annot, "GSE33000_DEG_results.csv", row.names = FALSE)

# --- 10. PATHWAY ENRICHMENT — GO & KEGG ---------------------------------------

# Filter: adj.P.Val < 0.05 AND |logFC| > 0.3
# logFC threshold lowered from standard 1.0 for brain microarray data
sig <- results_annot[results_annot$adj.P.Val < 0.05 & abs(results_annot$logFC) > 0.3, ]

up_genes   <- as.character(sig$EntrezGeneID[sig$logFC > 0])
down_genes <- as.character(sig$EntrezGeneID[sig$logFC < 0])
background <- as.character(results_annot$EntrezGeneID)

cat("Upregulated:  ", length(up_genes), "\n")
cat("Downregulated:", length(down_genes), "\n")

# GO Biological Process
go_up <- enrichGO(gene = up_genes, universe = background, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

go_down <- enrichGO(gene = down_genes, universe = background, OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

# KEGG
kegg_up <- enrichKEGG(gene = up_genes, universe = background, organism = "hsa",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)

kegg_down <- enrichKEGG(gene = down_genes, universe = background, organism = "hsa",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)

# Dot plots
p1 <- dotplot(go_up,    showCategory = 20, title = "GO BP — Upregulated in AD")
p2 <- dotplot(go_down,  showCategory = 20, title = "GO BP — Downregulated in AD")
p3 <- dotplot(kegg_up,  showCategory = 15, title = "KEGG — Upregulated in AD")
p4 <- dotplot(kegg_down,showCategory = 15, title = "KEGG — Downregulated in AD")

ggsave("figures/GO_BP_upregulated.png",   p1, width = 9, height = 8, dpi = 200)
ggsave("figures/GO_BP_downregulated.png", p2, width = 9, height = 8, dpi = 200)
ggsave("figures/KEGG_upregulated.png",    p3, width = 9, height = 7, dpi = 200)
ggsave("figures/KEGG_downregulated.png",  p4, width = 9, height = 7, dpi = 200)

# Export enrichment tables
write.csv(as.data.frame(go_up),    "GO_BP_upregulated.csv",   row.names = FALSE)
write.csv(as.data.frame(go_down),  "GO_BP_downregulated.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_up),  "KEGG_upregulated.csv",    row.names = FALSE)
write.csv(as.data.frame(kegg_down),"KEGG_downregulated.csv",  row.names = FALSE)
