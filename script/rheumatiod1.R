#set up my directory
getwd()
setwd("/home/rehabashraf/R/x86_64-pc-linux-gnu-library/4.5/ rheumatoid arthritis")
list.files()
# instaalling package 
install.packages("tidydr")
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("ggrepel")
BiocManager::install("RColorBrewer")
BiocManager::install("IHW")
BiocManager::install("apeglm")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("ReactomePA")
#calliing  package  
library(tidydr)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(IHW)
library(apeglm)
library(org.Hs.eg.db) 
library(tibble)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
#load data 
rh_count=read.table("GSE89408_GEO_count_matrix_rename.txt.gz",row.names = 1, header =TRUE,check.names =FALSE)
rh_count=as.matrix(rh_count)
rh_meta=read.csv("SraRunTable.csv",row.names = 1)
#  match data 
all(rownames(rh_meta)==colnames(rh_count))
colnames(rh_count)=rownames(rh_meta)
all(rownames(rh_meta)==colnames(rh_count))
## organize data 
table(rh_meta$disease)
rh_meta$disease=as.character(rh_meta$disease)
rh_meta$disease[rh_meta$disease== "Arthralgia"] ="ARTH"
rh_meta$disease[rh_meta$disease== "Osteoarthritis"] ="OST"
rh_meta$disease[rh_meta$disease== "Rheumatoid arthritis (early)"] ="RH_early"
rh_meta$disease[rh_meta$disease== "Rheumatoid arthritis (established)"] ="RH_established"
rh_meta$disease[rh_meta$disease== "Undifferentiated arthritis"] ="UNDIFF_arthritis"
#convert categories data to factor 
rh_meta$disease=as.factor(rh_meta$disease)
rh_meta$disease=factor(rh_meta$disease,levels = c("Normal","OST","ARTH","UNDIFF_arthritis","RH_early","RH_established"))
#deseq data set object
type(rh_count)
rh_count=round(rh_count)
dds=DESeqDataSetFromMatrix(countData = rh_count,colData = rh_meta,design = ~disease)
#pre_filteration
smallsizegroup=min(table(rh_meta$disease))
keep=rowSums(counts(dds)>=10) >=smallsizegroup
dds=dds[keep,]
dds
# inspect library size 
barplot(colSums(rh_count)/1e6, las=2,ylab = "million of reads",main = "librarysize")
# examine log10 count distribution 
hist(log10(rh_count[rh_count>0]),breaks = 50,main = "log10 count distribution")
# per sample 
boxplot(log2(rh_count+1),las=2,outline=FALSE)
#normaliation
vsd=vst(dds,blind = TRUE)
vsd
#sample distance 
sample_distance =dist(t(assay(vsd)))
sample_distance
sampledistance_matrix=as.matrix(sample_distance)
rownames(sampledistance_matrix)=paste(vsd$disease,sep = "_")
sampledistance_matrix
colnames(sampledistance_matrix)=NULL
sampledistance_matrix
add_colors=list(disease=c(Normal="blue",OST="red",ARTH="pink",UNDIFF_arthritis="orange",RH_early="yellow",RH_established="green"))
#heatmap
pheatmap(sampledistance_matrix,clustering_distance_rows = sample_distance,clustering_distance_cols = sample_distance,
         color = colorRampPalette(rev(brewer.pal(9,"Blues")))(255),main = "samplet-to sample diatance matrix",angle_col = 45,frontsize=10)

#pca 

pcaData<- plotPCA(vsd, intgroup = "disease", returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x=PC1, y=PC2, color=disease, label=name)) +geom_point(size=4, alpha=0.85) + geom_text(vjust=-0.8, size=3.2) +xlab(paste0("PC1:" , percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  
  theme_bw(base_size=13) + ggtitle("PCA of VST-normalized Expression") + scale_color_brewer(palette="Set1") 
# library size 
lib_size=data.frame(sample=colnames(dds),lib_size=colSums(counts(dds)/1e6),condition=dds$disease)
ggplot(lib_size, aes(x=reorder(sample, lib_size), y=lib_size, fill=condition)) +geom_col() + coord_flip() + labs(title="Library Sizes", x="Sample"  , y="Library Size (Millions)") +
  theme_bw()
# boxplot before normalization 
log2_counts=log2(counts(dds,normalized=FALSE)+1)
boxplot(log2_counts,las=2,col=as.integer(rh_meta$disease),main="raw counts(log2+1)",ylab="log2(counts+1",cex.axis=.7)
# after normalization 
dds_tmp=estimateSizeFactors(dds)
boxplot(log2(counts(dds_tmp,normalized=TRUE)+1),las=2,col=as.integer(rh_meta$disease),ylab="log2 count+1",main="normalized count(log2+1)")
#desq
ddss=DESeq(dds)
sizeFactors(ddss)
plotDispEsts(ddss)
#extract result
res=results(ddss,contrast = c("disease","RH_early","Normal"),alpha = 0.05,lfcThreshold = 0,pAdjustMethod = "BH")
View(res)
res=as.data.frame(res)
summary(res)
resultsNames(ddss)
#result with IHW
ihw_result=results(ddss,contrast = c("disease","RH_early","Normal"),filterFun = ihw)
View(ihw_result)
ihw_result=as.data.frame(ihw_result)
#lfc shrinkage 
res_shrinkage=lfcShrink(ddss,coef = "disease_RH_early_vs_Normal" ,type = "apeglm",res = res)
res_shrinkage=as.data.frame(res_shrinkage)
#gene annotation 
res_shrinkage$genename=mapIds(org.Hs.eg.db,keys = rownames(res_shrinkage),column = "GENENAME",keytype = "SYMBOL",multiVals = "first")
#filter ensemble_id 
ens=  res_shrinkage [grep("ENSG",rownames(res_shrinkage),value = TRUE),]
ens$ensamble=rownames(ens)
ens
rownames(ens)=sub(".*(ENSG[0-9]+).*","\\1",rownames(ens))
#gene anottation for ensamble 
ens$symbol =mapIds(org.Hs.eg.db,keys = rownames(ens),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ens$genename =mapIds(org.Hs.eg.db,keys = rownames(ens),column = "GENENAME",keytype = "ENSEMBL",multiVals = "first")
ens$entrezid=mapIds(org.Hs.eg.db,keys = rownames(ens),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
#filter for downstream analysis before annotation
res_filter=as.data.frame(res_shrinkage)%>%
  rownames_to_column("gene_id") %>%
  arrange(padj)%>%
  filter(!is.na(padj))
#filter for downstream analysis after annotation
ens_filter=as.data.frame(ens)%>%
  rownames_to_column("gene_id")%>%
  arrange(padj)%>%
  filter(!is.na(padj))
#filter signficant gene after genes anotation 
siges_genes=ens_filter %>% filter(padj < 0.05, abs(log2FoldChange) >1)
up_siges=filter(siges_genes,log2FoldChange>1)
down_siges=filter(siges_genes,log2FoldChange < -1)
#filter signficatnt gene before annotation
signf_genes=res_filter %>% filter(padj < 0.05 ,abs(log2FoldChange)>1)
signf_genes$entrezid=mapIds(org.Hs.eg.db,keys = signf_genes$gene_id,column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
up_signf=filter(signf_genes,log2FoldChange  >1)
down_signf=filter(signf_genes,log2FoldChange < -1)
#save file 
write.csv(signf_genes,"desq_allrheumatiod.csv",row.names = FALSE)
write.csv(up_signf,"desq_up_allrheumatiod.csv",row.names = FALSE)
write.csv(down_signf,"desq_down_allrheumatiod.csv",row.names = FALSE)
write.csv(siges_genes,"desq_annotat_rheumatiod.csv",row.names=FALSE)
write.csv(up_siges,"desq_annotat_uprheumatiod.csv",row.names=FALSE)
write.csv(down_siges,"desq_annotat_downrheumatiod.csv",row.names=FALSE)
#volcano plot
EnhancedVolcano(signf_genes,x = "log2FoldChange",lab = rownames(signf_genes),y = "padj",title = "disease_RH_early_vs_Normal",subtitle = "DESeq2 | apeglm shrinkage",pCutoff = 0.05,FCcutoff = 1, pointSize      = 2.5,
                labSize         = 3,col            = c("grey40", "forestgreen", "royalblue", "red2"),
                legendLabels   = c("NS", "LFC only", "p-value only", "p & LFC"),
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                max.overlaps   = 20)
#volcano plot after annotat
EnhancedVolcano(siges_genes,x = "log2FoldChange",lab =siges_genes$symbol,y = "padj",title = "disease_RH_early_vs_Normal",subtitle = "DESeq2 | apeglm shrinkage",pCutoff = 0.05,FCcutoff = 1, pointSize      = 2.5,
                labSize         = 3,col            = c("grey40", "forestgreen", "royalblue", "red2"),
                legendLabels   = c("NS", "LFC only", "p-value only", "p & LFC"),
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                max.overlaps   = 20)
#heatmap after anotat
top_genes= signf_genes %>%
        arrange(padj)%>%
  head(50)%>%
  pull(gene_id)
mat   <- assay(vsd)[top_genes, ]
mat   <- mat - rowMeans(mat)  


ann_col <- data.frame(
  disease = rh_meta$disease,
  row.names = colnames(vsd)
)
# heat map vst
pheatmap(mat,
         annotation_col   = ann_col,
         show_rownames    = TRUE,
         show_colnames    = TRUE,
         scale            = "row",
         clustering_method = "complete",
         color            = colorRampPalette(c("navy","white","firebrick3"))(100),
         fontsize_row     = 7,
         main             = "Top 50 DE Genes (DESeq2)",
         border_color     = NA
)
view(mat)

#functional enrichment analysis
up_entrez=up_signf$entrezid
up_entrez=na.omit(up_entrez)
down_entrez=down_signf$entrezid
down_entrez=na.omit(down_entrez)
down_entrez
#universe for down and high gene 
up_universe=na.omit(signf_genes$entrezid)
up_universe
down_universe=na.omit(signf_genes$entrezid)
#functional enrichment for up genes 
go_upen_BP=enrichGO(up_entrez,universe = up_universe,OrgDb = org.Hs.eg.db,ont = "Bp",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,readable = TRUE)
go_upen_MF=enrichGO(up_entrez,universe = up_universe,OrgDb = org.Hs.eg.db,ont = "MF",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,readable = TRUE)
go_upen_cc=enrichGO(up_entrez,universe = up_universe,OrgDb = org.Hs.eg.db,ont = "CC",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)
go_upen_all=enrichGO(up_entrez,universe = up_universe,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)
#functional enrichment  analysis for down gene
go_downen_BP=enrichGO(down_entrez,universe = down_universe,OrgDb = org.Hs.eg.db,ont = "Bp",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,readable = TRUE)
go_downen_MF=enrichGO(down_entrez,universe = down_universe,OrgDb = org.Hs.eg.db,ont = "MF",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,readable = TRUE)
go_downen_cc=enrichGO(down_entrez,universe = down_universe,OrgDb = org.Hs.eg.db,ont = "CC",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)
go_dowwnen_all=enrichGO(down_entrez,universe = down_universe,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)
#save functional enrichment for up gene
write.csv(go_upen_BP,"up_biologicalprocces.csv",row.names = FALSE)
write.csv(go_upen_cc,"up_cellularcomponent.csv",row.names = FALSE)
write.csv(go_upen_MF,"up_molecularfunction.csv",row.names = FALSE)
write.csv(go_upen_all,"up_allfunctional_enrichment.csv",row.names = FALSE)
#save  functional enrichment for down gene
write.csv(go_downen_BP,"down_biologicalprocces.csv",row.names = FALSE)
write.csv(go_downen_cc,"down_cellularcomponent.csv",row.names = FALSE)
 write.csv(go_downen_MF,"down_molecularfunction.csv",row.names = FALSE)
 write.csv(go_dowwnen_all,"down_allfunctional_enrichment.csv",row.names = FALSE)
#visualization up enrichment by dotplot
 dotplot(go_upen_BP,showCategory=10,title="up_biologicalproccess")
dotplot(go_upen_MF,showCategory=10,title="up_molecularfunction")
dotplot(go_upen_cc,showCategory=10,title="up_cellularcomponent")
dotplot(go_upen_all,showCategory=10,title="up_allenrichment")
#visualization  down enrichment
dotplot(go_downen_BP,showCategory=10,title="down_biologicalproccess")
dotplot(go_downen_MF,showCategory=10,title="down_molecularfunction")
dotplot(go_downen_cc,showCategory=10,title="down_cellularcomponent")
dotplot(go_dowwnen_all,showCategory=10,title="down_allenrichment")
#visualization up enrichment bybarplot
barplot(go_upen_BP,showCategory=10,title="up_biologicalproccess")
barplot(go_upen_MF,showCategory=10,title="up_molecularfunction")
barplot(go_upen_cc,showCategory=10,title="up_cellularcomponent")
barplot(go_upen_all,showCategory=10,title="up_allenrichment")
#visualization  down enrichment by barplot
barplot(go_downen_BP,showCategory=10,title="down_biologicalproccess")
barplot(go_downen_MF,showCategory=10,title="down_molecularfunction")
barplot(go_downen_cc,showCategory=10,title="down_cellularcomponent")
barplot(go_dowwnen_all,showCategory=10,title="down_allenrichment")
#emaplot up
emapplot(pairwise_termsim(go_upen_BP),showCategory=10)
emapplot(pairwise_termsim(go_upen_MF),showCategory=10)
emapplot(pairwise_termsim(go_upen_cc),showCategory=10)
emapplot(pairwise_termsim(go_upen_all),showCategory=10)
#emaplot down
emapplot(pairwise_termsim(go_downen_BP),showCategory=10)
emapplot(pairwise_termsim(go_downen_MF),showCategory=10)
emapplot(pairwise_termsim(go_downen_cc),showCategory=10)
emapplot(pairwise_termsim(go_dowwnen_all),showCategory=10)
#kggpathway enrichment
kggp_up_entrez=enrichKEGG(gene = up_entrez,organism = "hsa",pvalueCutoff = 0.05,universe = up_universe)
kggp_down_entrez=enrichKEGG(gene = down_entrez,organism = "hsa",pvalueCutoff = 0.05,universe = down_universe)
#save kggpathway
write.csv(kggp_up_entrez,"kggpathway_up.csv",row.names = FALSE)
write.csv(kggp_down_entrez,"kggpathway_down.csv",row.names = FALSE)
#visualization kggpathway
dotplot(kggp_up_entrez,showCategory=10,title="kggp_up")
dotplot(kggp_down_entrez,showCategory=10,title="kggp_down")
#GSEA
ranked_gene=signf_genes %>%
  filter(!is.na(entrezid),!is.na(pvalue))%>%
  mutate(rank=-log10(pvalue)*sign(log2FoldChange))%>%
  arrange(desc(rank))
ranked_gene_list=ranked_gene$rank
names(ranked_gene_list)=ranked_gene$entrezid
gsea_res <- gseKEGG(
  geneList      = ranked_gene_list,
  organism      = "hsa",
  nPermSimple   = 1000,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)
gseaplot(gsea_res,geneSetID = 1,title = gsea_res$Description[1])
dotplot(gsea_res)
ridgeplot(gsea_res)
emapplot(gsea_res,showCategory=10)
#Reactome up__gene
reactom_res_up=enrichPathway(
  gene     = up_entrez,
  organism = "human",
  universe = up_universe,
  pvalueCutoff = 0.05
)
#reactome pathway  for down ggenne
reactome_res_down=enrichPathway(gene = down_entrez,organism = "human",universe = down_universe,pvalueCutoff = 0.05)
#visuualizatiion for  up  
dotplot(reactom_res_up,showCategory=10,title="reactom_res_up")
#visuualizatiion for  down
dotplot(reactome_res_down,showCategory=10,title="reactom_res_up")









