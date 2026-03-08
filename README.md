# GSE33000 — Differential Expression Analysis in Alzheimer's Disease

## Background

Alzheimer's disease (AD) is a progressive neurodegenerative disorder and the leading cause of dementia worldwide. Its molecular basis remains incompletely understood, with competing hypotheses centred on amyloid-beta accumulation and tau tangle formation. Transcriptomic profiling of post-mortem brain tissue offers a hypothesis-free approach to identifying gene expression changes associated with disease.

This analysis uses the GSE33000 dataset (NCBI GEO), a large microarray study of human prefrontal cortex tissue comprising 310 AD and 157 non-demented control samples.

---

## Dataset

| Parameter | Details |
|-----------|---------|
| Accession | GSE33000 |
| Tissue | Prefrontal cortex (post-mortem) |
| Platform | Microarray |
| Samples | 467 total (310 AD, 157 non-demented) |
| Probes | 25,929 |
| Source | NCBI GEO |

---

## Methods

### Differential Expression
Differential expression analysis was performed using the **limma** package in R. A linear model was fitted with disease state (AD vs non-demented) as the predictor. Empirical Bayes moderation was applied via `eBayes()`. Probes without a valid EntrezGeneID were excluded prior to analysis. Results were corrected for multiple testing using the Benjamini-Hochberg FDR method.

A fold change cutoff of |logFC| > 0.5 was applied for visualisation (volcano plot), reduced from the standard 1.0 to account for the characteristically small fold changes observed in brain microarray data.

### Pathway Enrichment
Over-representation analysis was performed using **clusterProfiler** against Gene Ontology Biological Process (GO BP) and KEGG pathway databases. Significant DEGs were defined as adj.P.Val < 0.05 and |logFC| > 0.3, yielding 262 upregulated and 249 downregulated genes. The background gene set comprised all 25,929 tested probes.

---

## Key Findings

### Downregulated in AD
Gene ontology analysis revealed consistent downregulation of neuronal processes including regulation of membrane potential, axonogenesis, synapse assembly, synapse organisation, and neuropeptide signalling. KEGG analysis corroborated this, with neuroactive ligand-receptor interaction as the most significantly enriched downregulated pathway, alongside GABAergic, serotonergic, cholinergic, and glutamatergic synapse pathways.

### Upregulated in AD
Two distinct themes emerged in upregulated genes. First, a strong adaptive immune signature including humoral immune response, B cell mediated immunity, leukocyte mediated immunity, and complement activation. Second, a metal ion stress response involving copper, zinc, and cadmium response pathways. KEGG enrichment identified complement and coagulation cascades, cytokine-cytokine receptor interaction, and TNF signalling as the most significantly upregulated pathways.

### Interpretation
In the AD prefrontal cortex, we observe concurrent downregulation of neuronal functions — including membrane potential regulation, ligand-receptor signalling, and synapse formation — alongside upregulation of adaptive immune response and metal ion stress pathways. This dataset does not establish causal direction between these processes; both may represent downstream consequences of amyloid or tau pathology, or independent parallel processes.

---

## Important Caveats

- **XIST artefact:** XIST appears upregulated in AD in this dataset due to a sex imbalance between AD and control groups, not a disease signal. It has been excluded from biological interpretation.
- **Gene name column:** This dataset uses `ORF` as the gene name field. There is no Gene Symbol column in the feature table.
- **Fold change direction:** Positive logFC = higher expression in AD relative to non-demented controls.
- **Causality:** This is an observational transcriptomic dataset. Co-occurring expression changes cannot establish causation.

---

## Repository Structure
```
GSE33000-AD-analysis/
├── README.md                    
├── GSE33000_analysis.R          
├── GSE33000_DEG_results.csv     
└── figures/
    ├── volcano_plot.png
    ├── GO_BP_upregulated.png
    ├── GO_BP_downregulated.png
    ├── KEGG_upregulated.png
    └── KEGG_downregulated.png
```

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| GEOquery | Bioconductor | Data loading |
| limma | Bioconductor | Differential expression |
| EnhancedVolcano | Bioconductor | Volcano plot |
| clusterProfiler | Bioconductor | Pathway enrichment |
| org.Hs.eg.db | Bioconductor | Human gene annotation |
| ggplot2 | CRAN | Plotting |

---

## Environment
R / Google Colab
