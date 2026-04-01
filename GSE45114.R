library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
# Loading data
GPL10558 <-data.table::fread("C:/GSE45114/GPL5918.txt",skip ="ID")
exprSet<-read.table("C:/GSE45114/GSE45114_series_matrix.txt",comment.char="!",
                    stringsAsFactors=F,
                    sep="\t",
                    header=T)
# Group the samples
sample.type <- read.table("C:/GSE45114/sample_Type.txt",
                          header = F, sep = "\t")
colnames(sample.type) <- c("sample", "type")
sample.type <- data.frame(sample.type,
                          group = str_detect(sample.type$type, pattern = "HCC"))
sample.type$group <- gsub(sample.type$group, pattern = "TRUE", replacement = "Tumor")
sample.type$group <- gsub(sample.type$group, pattern = "FALSE", replacement = "Normal")
group <- sample.type$group
table(group)

# Extract Columns 
probe2symbol_df <- GPL10558 %>% 
  
  dplyr::select(ID,HUGOname)

names(exprSet)[1] <- names(probe2symbol_df)[1]

library(dplyr)
exprSet <- exprSet %>% 
  inner_join(probe2symbol_df,by="ID") %>% 
  
  dplyr::select(-ID) %>% 
  
  dplyr::select(HUGOname, everything()) %>% 
  
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% 
  
  arrange(desc(rowMean)) %>% 
  
  distinct(HUGOname,.keep_all = T) %>% 
  
  dplyr::select(-rowMean) %>% 
  
  tibble::column_to_rownames(colnames(.)[1]) 

#log
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0,arr.ind = T)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)



# Performing differential analysis using the limma method
library(limma)
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design
contrast.matrix <- makeContrasts(Tumor-Normal , levels=design)
contrast.matrix
fit <- lmFit(exprSet, design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',coef=1,number=Inf) 
result<-subset(allDiff,allDiff$adj.P.Val<0.05 & abs(allDiff$logFC) >= 1)
write.csv(result,"C:/GSE45114/GSE45114-Differentially expressed genes.csv", 
          col.names=TRUE, row.names=TRUE, append=FALSE)

# Volcano Chart
library(ggplot2)
library(ggrepel)
library(dplyr)
data <- allDiff
data$gene <- rownames(data)
ggplot(data=data, aes(x=logFC, y =-log10(P.Value)))

data$Significant <- ifelse(data$P.Value < 0.05 & abs(data$logFC) >= 1, 
                           ifelse(data$logFC > 1, "Up", "Down"), "Stable")
ggplot(
  data, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(position="jitter",aes(color = Significant), size=1) +
  scale_color_manual(values = c("#2f5688","#BBBBBB", "#CC0000")) +
  geom_vline(xintercept=c(-1,1),lty=2,col=c("#2f5688","#CC0000"),lwd=0.8) +
  geom_hline(yintercept = 1,lty=2,col="black",lwd=0.8) +
  
  # Topic
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",legend.title = element_blank())+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(legend.position='none')


# Heatmap
library(pheatmap)
heatdata <- exprSet[rownames(result),]
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)
heatdata<-na.omit(heatdata)
pheatmap(heatdata, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col =annotation_col, 
         annotation_legend=TRUE,
         show_rownames = F,
         show_colnames = F,
         scale = "row", 
         color =colorRampPalette(c("#2f5688", "white","#CC0000"))(100))