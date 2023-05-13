#Differential Gene Expression Analysis using DESeq2 

# Load libraries
library("DESeq2")
library(ggplot2)
library(stringr)
library(ggrepel)

#Read in all RNAseq data
CR_12h <- read.csv("../Data inputs/GE_RawData_Night12Day12.csv",row.names = 1, check.names = FALSE)
Night <- read.csv("../Data inputs/GE_RawData_Night12Night2.csv",row.names = 1, check.names = FALSE)
Day <- read.csv("../Data inputs/GE_RawData_Day12Day2.csv",row.names = 1, check.names = FALSE)
CRvsNight12 <- read.csv("../Data inputs/GE_RawData_CRNight12.csv",row.names = 1, check.names = FALSE)


#Read in column data
ColData_CR12 <- read.csv("../Data inputs/DESeq2 ColData_CR12.csv",row.names = 1)
ColData_Night <- read.csv("../Data inputs/DESeq2 ColData_Night.csv",row.names = 1)
ColData_Day <- read.csv("../Data inputs/DESeq2 ColData_Day.csv",row.names = 1)
ColData_CRvsNight <- read.csv("../Data inputs/DESeq2 ColData_CRvsNight12.csv",row.names = 1)


#Generate DESeq dataframe for each comparison

#Night12 vs. Day12
dds_CR12 <- DESeqDataSetFromMatrix(countData = CR_12h,colData = ColData_CR12,design = ~Condition)
dds_CR12 <- DESeq(dds_CR12)

#Night12 vs Night2
dds_Night <- DESeqDataSetFromMatrix(countData = Night,colData = ColData_Night,design = ~Condition)
dds_Night <- DESeq(dds_Night)

#Day12 vs Day2
dds_Day <- DESeqDataSetFromMatrix(countData = Day,colData = ColData_Day,design = ~Condition)
dds_Day <- DESeq(dds_Day)

#CR spread vs Night12
dds_CRvsNight12 <- DESeqDataSetFromMatrix(countData = CRvsNight12,colData = ColData_CRvsNight,design = ~Condition)
dds_CRvsNight12 <- DESeq(dds_CRvsNight12)


#Generate results dataframes for each condition
res_CR12 <- results(dds_CR12)
res_Night <- results(dds_Night)
res_Day <- results(dds_Day)
res_CRvsNight12 <- results(dds_CRvsNight12)


#order by p values
res_CR12 <- res_CR12[order(res_CR12$padj),]
res_Night <- res_Night[order(res_Night$padj),]
res_Day <- res_Day[order(res_Day$padj),]
res_CRvsNight12 <- res_CRvsNight12[order(res_CRvsNight12$padj),]


#Save results to csv (provide descriptive filename)
write.csv(as.data.frame(res_CR12), 
          file="RNA_DESeq_results_Night12Day12.csv")
write.csv(as.data.frame(res_Night), 
          file="RNA_DESeq_results_Night12Night2.csv")
write.csv(as.data.frame(res_Day), 
          file="RNA_DESeq_results_Day12Day2.csv")
write.csv(as.data.frame(res_CRvsNight12), 
          file="RNA_DESeq_results_CRvsNight12.csv")


#Label genes as significant or not and upregulated or downregulated
resultdf_CR12=as.data.frame(res_CR12)
resultdf_Night=as.data.frame(res_Night)
resultdf_Day=as.data.frame(res_Day)
resultdf_CRvsNight12=as.data.frame(res_CRvsNight12)

# add a column of NAs
resultdf_CR12$diffexpressed <- "Not Significant"
resultdf_Night$diffexpressed <- "Not Significant"
resultdf_Day$diffexpressed <- "Not Significant"
resultdf_CRvsNight12$diffexpressed <- "Not Significant"

# if log2Foldchange > 1.5 and pvalue < 0.05, set as "UP" 
resultdf_CR12$diffexpressed[resultdf_CR12$log2FoldChange >
                         1.5 & resultdf_CR12$padj< 0.05] <-"Significant"
resultdf_Night$diffexpressed[resultdf_Night$log2FoldChange >
                       1.5 & resultdf_Night$padj< 0.05] <-"Significant"
resultdf_Day$diffexpressed[resultdf_Day$log2FoldChange >
                         1.5 & resultdf_Day$padj< 0.05] <-"Significant"
resultdf_CRvsNight12$diffexpressed[resultdf_CRvsNight12$log2FoldChange >
                         1.5 & resultdf_CRvsNight12$padj< 0.05] <-"Significant"

# if log2Foldchange < -1.5 and pvalue < 0.05, set as "DOWN"
resultdf_CR12$diffexpressed[resultdf_CR12$log2FoldChange <
                              -1.5 & resultdf_CR12$padj< 0.05] <-"Significant"
resultdf_Night$diffexpressed[resultdf_Night$log2FoldChange <
                               -1.5 & resultdf_Night$padj< 0.05] <-"Significant"
resultdf_Day$diffexpressed[resultdf_Day$log2FoldChange <
                             -1.5 & resultdf_Day$padj< 0.05] <-"Significant"
resultdf_CRvsNight12$diffexpressed[resultdf_CRvsNight12$log2FoldChange <
                                     -1.5 & resultdf_CRvsNight12$padj< 0.05] <-"Significant"

#Label differentially expressed genes with their gene IDs
resultdf_CR12$delabel <- NA
resultdf_CR12$delabel[resultdf_CR12$diffexpressed != "Not Significant"] <- rownames(resultdf_CR12)[resultdf_CR12$diffexpressed != "Not Significant"]
resultdf_CR12$delabel <- str_split(resultdf_CR12$delabel, "_", simplify = TRUE)[,2]

resultdf_Night$delabel <- NA
resultdf_Night$delabel[resultdf_Night$diffexpressed != "Not Significant"] <- rownames(resultdf_Night)[resultdf_Night$diffexpressed != "Not Significant"]
resultdf_Night$delabel <- str_split(resultdf_Night$delabel, "_", simplify = TRUE)[,2]

resultdf_Day$delabel <- NA
resultdf_Day$delabel[resultdf_Day$diffexpressed != "Not Significant"] <- rownames(resultdf_Day)[resultdf_Day$diffexpressed != "Not Significant"]
resultdf_Day$delabel <- str_split(resultdf_Day$delabel, "_", simplify = TRUE)[,2]

resultdf_CRvsNight12$delabel <- NA
resultdf_CRvsNight12$delabel[resultdf_CRvsNight12$diffexpressed != "Not Significant"] <- rownames(resultdf_CRvsNight12)[resultdf_CRvsNight12$diffexpressed != "Not Significant"]
resultdf_CRvsNight12$delabel <- str_split(resultdf_CRvsNight12$delabel, "_", simplify = TRUE)[,2]


#Create Volcano Plot
p_CR12 <- ggplot(data=resultdf_CR12, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = delabel)) +
  geom_point(size=4) +  scale_color_manual(values = 
                                             c("Not Significant" = "black",  "Significant"="red")) +
  geom_text_repel(nudge_y=0.25, size=10) +
  ggtitle("Differential Gene Expression Between 
          \nNighttime Fed (12h) and Daytime Fed (12h) Mice")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("-log10(pvalue)")+
  theme(legend.title = element_blank())+
  theme(text = element_text(size = 35))

p_Night <- ggplot(data=resultdf_Night, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = delabel)) +
  geom_point(size=4) +  scale_color_manual(values = 
                                             c("Not Significant" = "black",  "Significant"="red")) +
  geom_text_repel(nudge_y=0.25, size=10) +
  ggtitle("Differential Gene Expression Between 
          \nNighttime Fed (12h) and Nighttime Fed (2h) Mice")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("-log10(pvalue)")+
  theme(legend.title = element_blank())+
  theme(text = element_text(size = 35))

p_Day <- ggplot(data=resultdf_Day, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = delabel)) +
  geom_point(size=4) +  scale_color_manual(values = 
                                             c("Not Significant" = "black",  "Significant"="red")) +
  geom_text_repel(nudge_y=0.25, size=10) +
  ggtitle("Differential Gene Expression Between 
          \nDaytime Fed (12h) and Daytime Fed (2h) Mice")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("-log10(pvalue)")+
  theme(legend.title = element_blank())+
  theme(text = element_text(size = 35))

p_CRvsNight12 <- ggplot(data=resultdf_CRvsNight12, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = delabel)) +
  geom_point(size=4) +  scale_color_manual(values = 
                                             c("Not Significant" = "black",  "Significant"="red")) +
  geom_text_repel(nudge_y=0.25, size=10) +
  ggtitle("Differential Gene Expression Between 
          \nCR Spread and Nighttime Fed (12h) Mice")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("-log10(pvalue)")+
  theme(legend.title = element_blank())+
  theme(text = element_text(size = 35))


#Save volcano plots as .jpg images
ggsave('volcanoplot_Night12vsDay12_labeled.jpg', plot=p_CR12, width = 20, height = 15,
       dpi = 600, units = "in")
ggsave('volcanoplot_Night12vsNight2_labeled.jpg', plot=p_Night, width = 20, height = 15,
       dpi = 600, units = "in")
ggsave('volcanoplot_Day12vsDay2_labeled.jpg', plot=p_Day, width = 20, height = 15,
       dpi = 600, units = "in")
ggsave('volcanoplot_CRvsNight12_labeled.jpg', plot=p_CRvsNight12, width = 20, height = 15,
       dpi = 600, units = "in")
