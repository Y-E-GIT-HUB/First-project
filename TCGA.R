rm(list = ls())
if(!require(datTraits))BiocManager::install('datTraits')
if(!require(ggplotify))BiocManager::install('ggplotify')
if(!require(patchwork))install.packages("patchwork")
if(!require(cowplot))install.packages("cowplot")
if(!require(DESeq2))BiocManager::install('DESeq2')
if(!require(edgeR))BiocManager::install('edgeR')
if(!require(datTraits))BiocManager::install('datTraits')
if(!require(stringr))install.packages("stringr")
library(tidyr)
library(dplyr)
# Loading data
a= read.table("C:/TCGA-LIHC_MergeCOUNT.txt"
              ,header=T)
colnames(a)[1]<-'gene_id'
dim(a)
count1<-separate(a,gene_id,into = c("gene_id"),sep="\\.")

#Note: Gene name
load("C:/id_symble_type.Rda")
count2=merge(geneid_df_nopoint,count1,by = "gene_id" )
count2 <- count2[, c(-1,-3)]
# Replace the first column with row names
row.names(count2) <- count2[, 1]
count2 <- count2[, -1]

# Group the samples
library(stringr)
group_list <- ifelse(as.numeric(str_sub(colnames(count2),14,15))<10,"tumor","normal")
group_list <- factor(group_list,levels = c("normal","tumor"))
table(group_list)

# Performing differential analysis using the DESeq2 method
library(DESeq2)
colData <- data.frame(row.names =colnames(count2), 
                      condition=group_list) 
dds <- DESeqDataSetFromMatrix(
  countData = count2,
  colData = colData,
  design = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$padj),] 
DEG <- as.data.frame(resOrdered)
head(DEG)
DEG <- na.omit(DEG)
DEG$change = as.factor(
  ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > 1.5,
         ifelse(DEG$log2FoldChange > 1.5 ,'UP','DOWN'),'NOT')
)
head(DEG)
DESeq2_DEG <- DEG
write.csv(DESeq2_DEG, "C:/TCGA-Differentially expressed genes.csv")

cg = rownames(DESeq2_DEG)[DESeq2_DEG$change !="NOT"]

# Heatmap
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#3F2BFF','white','#A20011'))(100) 
mat=count2[cg,]
n=t(scale(t(mat)))
n[n>1]=1
n[n< -1]= -1
ac=data.frame(group=group_list)
rownames(ac)=colnames(mat)
ht1 <- pheatmap(n,show_rownames = F,show_colnames = F, 
                cluster_rows = F,cluster_cols = T,
                annotation_col = ac,color=color)
# Volcano Chart
library(ggplot2)
library(ggrepel)
library(dplyr)
data <- DESeq2_DEG
data$gene <- rownames(DESeq2_DEG)
ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj)))

data$Significant <- ifelse(data$padj < 0.05 & abs(data$log2FoldChang) >= 1.5, 
                           ifelse(data$log2FoldChang > 1.5, "Up", "Down"), "Stable")
ggplot(
  data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(position="jitter",aes(color = Significant), size=1) +
  scale_color_manual(values = c("#7815FF","#BBBBBB", "#EC004E")) +
  geom_vline(xintercept=c(-1.5,1.5),lty=2,col=c("black","black"),lwd=0.8) +
  geom_hline(yintercept = 2,lty=2,col="black",lwd=0.8) +
  
  # Topic
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",legend.title = element_blank())+

  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+

  labs(x="log2 (fold change)",y="-log10 (padj)")

theme(legend.position='none')



