# Rheumatoid Arthritis RNA-Seq Analysis Using DESeq2


##objectives
Rheumatoid Arthritis (RA) is a chronic autoimmune disease characterized by persistent inflammation and progressive joint destruction.

This project aims to:

Identify differentially expressed genes (DEGs) in RA patients.
Investigate biological processes associated with disease development.
Explore enriched molecular pathways and signaling networks.
Generate biologically meaningful insights from RNA-seq data.
## dataset
Dataset ID: GSE89408

Source: Gene Expression Omnibus (GEO)

Study Design:
Samples include:
.Healthy Controls
.Osteoarthritis (OST)
.Arthralgia (ARTH)
.Undifferentiated Arthritis
.Early Rheumatoid Arthritis
.Established Rheumatoid Arthritis

For differential expression analysis, the comparison performed in this project is:
Early Rheumatoid Arthritis vs Healthy Controls
##Bioinformatics Workflow
1. Data Import
.RNA-seq count matrix loading
.Sample metadata loading
.Metadata matching with count matrix
2. Data Preprocessing
.Disease category cleaning
.Sample annotation organization
.Conversion of variables to factors
.Construction of DESeq2 dataset
3. Gene Filtering
Lowly expressed genes were removed using:
.Minimum count threshold ≥ 10
.Expression required in at least the smallest experimental group
4. Quality Control

.Library size inspection
.Count distribution visualization
.Boxplots before normalization
.Sample-to-sample distance heatmap
.Principal Component Analysis (PCA)
5. Normalization
Variance Stabilizing Transformation (VST) was applied to reduce heteroscedasticity and improve downstream visualization.
6. Differential Expression Analysis

Differential expression analysis was performed using DESeq2.

Comparison:

Early Rheumatoid Arthritis vs Healthy Controls

Statistical criteria:
Adjusted p-value < 0.05
|log2 Fold Change| > 1

7. LFC Shrinkage
Log fold changes were stabilized using the apeglm shrinkage method.

8. Gene Annotation

Gene identifiers were mapped to:
.Gene Symbols
.Gene Names
.Entrez IDs
using org.Hs.eg.db.
9. Visualization
Generated visualizations include:
.PCA plots
.Sample distance heatmaps
.Volcano plots
.Differential expression heatmaps

10. Functional Enrichment Analysis

Gene Ontology (GO):
.Biological Process (BP)
.Molecular Function (MF)
.Cellular Component (CC)
11. Pathway Analysis

KEGG pathway enrichment analysis was performed to identify significantly affected biological pathways.

12. Gene Set Enrichment Analysis (GSEA)

Ranked-gene GSEA was used to identify coordinated pathway-level changes associated with Rheumatoid Arthritis.

13. Reactome Analysis

Reactome pathway enrichment analysis was performed to investigate higher-order biological mechanisms.
##Tools and Packages

Programming Language
R
Main Packages
DESeq2
clusterProfiler
ReactomePA
enrichplot
EnhancedVolcano
org.Hs.eg.db
ggplot2
pheatmap
dplyr
apeglm
IHW
##Output Files
The workflow generates:

Differential Expression Results
desq_allrheumatiod.csv
desq_up_allrheumatiod.csv
desq_down_allrheumatiod.csv
##Annotated Results
desq_annotat_rheumatiod.csv
desq_annotat_uprheumatiod.csv
desq_annotat_downrheumatiod.csv
##GO Enrichment
up_biologicalprocces.csv
down_biologicalprocces.csv
up_molecularfunction.csv
down_molecularfunction.csv
##KEGG Analysis
kggpathway_up.csv
kggpathway_down.csv
