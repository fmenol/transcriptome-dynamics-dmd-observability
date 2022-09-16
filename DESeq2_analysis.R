library(htmltools)
library("DESeq2")
library(ggplot2)
library(tximport)
library(tidyverse)
library(edgeR) 

# load the data and metadata
targets <- read_tsv("data/kallisto-results/studydesign.txt")
path <- file.path("data/kallisto-results", targets$sample, "abundance.tsv") # set file paths to your mapped data
all(file.exists(path)) 
sampleLabels <- targets$sample

Tx <- rtracklayer::import("data/kallisto-results/GCA_000009225.1_ASM922v1_genomic.gtf")
Tx <- as_tibble(Tx)
Tx <- dplyr::select (Tx,-c(source,type,score,phase))
Tx <- distinct(Tx, transcript_id,gene_id, .keep_all= TRUE)
# rename column which contains unique transcript identifier to target_id . 
# rename column which contains unique locus tag to gene_name
Tx <- dplyr::rename(Tx, target_id = protein_id, gene_name = locus_tag)
# just select those columns, don't really need the other for this script. 
Tx <- dplyr::select(Tx, "target_id", "gene_name")
# drop rows that contain NA in any column 
Tx <- na.omit(Tx)

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = TRUE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

# saving the TPM for all genes (no low count filtering) 
myTPM <- Txi_gene$abundance
myTPM <- as_tibble(myTPM, rownames = "geneID")
colnames(myTPM) <- c('geneID',targets$sample)
# write.table(myTPM,file="tpm.csv",sep=',') # TPM

# count matrix from kallisto 
cts <- Txi_gene$counts
mode(cts) <- "integer"
colnames(cts) <- c(targets$sample)

# create a DGEList using the edgeR package (using normalized cpm to filter low expression genes)
DGEList.counts <- DGEList(cts)
cpm <- cpm(DGEList.counts) 

# Filter data
cutoff = 5
keepers <- rowSums(cpm>cutoff)>=(length(sampleLabels)-2) # rowSums(cpm>cutoff)>=(length(sampleLabels)/2)
sum(keepers, na.rm = TRUE) 

# subset the cts
cts <- cts[keepers,]

# loading metadata
colData <- data.frame(read_csv("data/kallisto-results/metaData.csv"))
rownames(colData) <- c(colnames(cts))

DESEQ <- function (start,cts,colData) {
  # start by getting the starting datapoint
  # start variable can be "t3" through "t11"
  if ( start == "t3") {
    start <- 1
  } else if ( start == "t4") {
    start <- 5
  } else if ( start == "t5") {
    start <- 9
  } else if ( start == "t6") {
    start <- 13
  } else if ( start == "t7") {
    start <- 17
  } else if ( start == "t8") {
    start <- 21
  } else if ( start == "t9") {
    start <- 25
  } else if ( start == "t10") {
    start <- 29
  } else if ( start == "t11") {
    start <- 33
  }
  
  # for simplicity, let's only keep a single timepoint. 
  # also not using replicate 2
  # subset the data accordingly
  cts <- cts[,c(start,start+1,start+2,start+3)]
  colData <- colData[c(start,start+1,start+2,start+3),]
  print(rownames(colData))
  print(all(colnames(cts) == rownames(colData)))
  
  # construct DESeqDataSet Object
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData=colData,
                                design=~Condition)
  dds <- DESeq(dds)
  res_dds <- results(dds) 
  return(res_dds)
}


# start <- 33
# ctsNew <- cts[,c(start,start+1,start+2,start+3)]
# colDataNew <- colData[c(start,start+1,start+2,start+3),]
# print(rownames(colDataNew))

res_t3 <- DESEQ("t3",cts,colData)
res_t4 <- DESEQ("t4",cts,colData)
res_t5 <- DESEQ("t5",cts,colData)
res_t6 <- DESEQ("t6",cts,colData)
res_t7 <- DESEQ("t7",cts,colData)
res_t8 <- DESEQ("t8",cts,colData)
res_t9 <- DESEQ("t9",cts,colData)
res_t10 <- DESEQ("t10",cts,colData)
res_t11 <- DESEQ("t11",cts,colData)


p_thresh <- 0.05
log2fc_thresh <- 0.26 # fold change of 1.25 or 0.8 as threshold
resSig_t3 <- res_t3[which(res_t3$pvalue < p_thresh,abs(res_t3$log2FoldChange) > log2fc_thresh),]
t3genes <- rownames(resSig_t3)
2^resSig_t3$log2FoldChange
resSig_t4 <- res_t4[which(res_t4$pvalue < p_thresh,abs(res_t4$log2FoldChange) > log2fc_thresh),]
t4genes <- rownames(resSig_t4)
2^resSig_t3$log2FoldChange
resSig_t5 <- res_t5[which(res_t5$pvalue < p_thresh,abs(res_t5$log2FoldChange) > log2fc_thresh),]
t5genes <- rownames(resSig_t5)
2^resSig_t5$log2FoldChange
resSig_t6 <- res_t6[which(res_t6$pvalue < p_thresh,abs(res_t6$log2FoldChange) > log2fc_thresh),]
t6genes <- rownames(resSig_t6)
2^resSig_t6$log2FoldChange
resSig_t7 <- res_t7[which(res_t7$pvalue < p_thresh,abs(res_t7$log2FoldChange) > log2fc_thresh),]
t7genes <- rownames(resSig_t7)
2^resSig_t7$log2FoldChange
resSig_t8 <- res_t8[which(res_t8$pvalue < p_thresh,abs(res_t8$log2FoldChange) > log2fc_thresh),]
t8genes <- rownames(resSig_t8)
2^resSig_t8$log2FoldChange
resSig_t9 <- res_t9[which(res_t9$pvalue < p_thresh,abs(res_t9$log2FoldChange) > log2fc_thresh),]
t9genes <- rownames(resSig_t9)
2^resSig_t9$log2FoldChange
resSig_t10 <- res_t10[which(res_t10$pvalue < p_thresh,abs(res_t10$log2FoldChange) > log2fc_thresh),]
t10genes <- rownames(resSig_t10)
2^resSig_t10$log2FoldChange
resSig_t11 <- res_t11[which(res_t11$pvalue < p_thresh,abs(res_t11$log2FoldChange) > log2fc_thresh),]
t11genes <- rownames(resSig_t11)
2^resSig_t11$log2FoldChange


# write the significant gene names to files
# write.table(t3genes,file="data/deseq_t3genes.csv",sep=',') 
# write.table(t4genes,file="DESeq2_sig_genes/deseq_t4genes.csv",sep=',')
# write.table(t5genes,file="DESeq2_sig_genes/deseq_t5genes.csv",sep=',') 
# write.table(t6genes,file="DESeq2_sig_genes/deseq_t6genes.csv",sep=',') 
# write.table(t7genes,file="DESeq2_sig_genes/deseq_t7genes.csv",sep=',') 
# write.table(t8genes,file="DESeq2_sig_genes/deseq_t8genes.csv",sep=',') 
# write.table(t9genes,file="DESeq2_sig_genes/deseq_t9genes.csv",sep=',') 
# write.table(t10genes,file="DESeq2_sig_genes/deseq_t10genes.csv",sep=',') 
# write.table(t11genes,file="DESeq2_sig_genes/deseq_t11genes.csv",sep=',') 



write.table(resSig_t6,file="data/deseq_t6.csv",sep=',') 
write.table(resSig_t7,file="data/deseq_t7.csv",sep=',') 
write.table(resSig_t8,file="data/deseq_t8.csv",sep=',') 
write.table(resSig_t9,file="data/deseq_t9.csv",sep=',') 
write.table(resSig_t10,file="data/deseq_t10.csv",sep=',') 
write.table(resSig_t11,file="data/deseq_t11.csv",sep=',') 


# informative plots
# plotMA(res_t11, ylim = c(-1, 1) )
# plotDispEsts(dds, ylim = c(1e-6, 1e1) )


# volcano plot
topT <- as.data.frame(resSig.t4)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)


