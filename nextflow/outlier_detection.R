#!/usr/bin/env Rscript
# -P casa

#script created by Benyu Zhou; modified by Jenny Empawi

args <- commandArgs(TRUE)

#--------------------------------gene-----------------------------------#
######Step 1: load data
SE_gene <- readRDS(args[1])
TPM_gene <- SE_gene@assays@data@listData[["TPM"]]

#####Step 2: filter low-expressed genes
library(dplyr)
library(readr)
library(tibble)
library(purrr)

low_expr_TPM <- 0.1
low_expr_TPM_percent <- 0.2

keep_genes_idx <- (rowMeans(TPM_gene>low_expr_TPM)>low_expr_TPM_percent) 
TPM_gene <- TPM_gene[keep_genes_idx,]

keep_genes_idx <- keep_genes_idx[1:nrow(assay(SE_gene))]
SE_gene <- SE_gene[keep_genes_idx, ]

#####Step 3: detect outliers using PCA, RLE, D-statistics
###method 1: PCA 
library('mt')

log_offset <- 0.0001

#prepare gene data for PCA
TPM_gene_t <- as.data.frame(t(TPM_gene))
gene_pca_log <- log10(TPM_gene_t+log_offset)

x <- pca.outlier(gene_pca_log, center = TRUE, scale =FALSE, conf.level = 0.975) 

message(paste0("These samples are marked as candidate outliers in PCA step:") )
print(x$outlier)

pdf(file = "gene_pca.pdf")
x$plot
dev.off()

###method 2: RLE
library(RColorBrewer)
library(ape)
library(reshape2)
library(dplyr)
library(readr)

#set RLE parameters
RLEFilterPercent <- 0.05
DSFilterPercent <- 0.05

#define functions to accomodate rank deficient covariance matrices in https://github.com/cumc/xqtl-pipeline/issues/307
mahalanobis <- function (x, center, cov, inverted = FALSE, ...) {
  x <- if (is.vector(x)) 
    matrix(x, ncol = length(x))
  else as.matrix(x)
  if (!isFALSE(center)) 
    x <- sweep(x, 2L, center)
  if (!inverted) 
    cov <- MASS::ginv(cov, ...)
  setNames(rowSums(x %*% cov * x), rownames(x))
}

RLEFilterLength <- RLEFilterPercent*ncol(TPM_gene)
DSFilter <- DSFilterPercent*ncol(TPM_gene)

logtpm <- log10(TPM_gene + log_offset)
rle <- logtpm-apply(logtpm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.
iqr <- apply(rle,2,IQR)
rle <- melt( rle, variable.name = "Sample",value.name ="TPM", id="ID")
names(rle) <- c("feature","Sample","TPM")
rle_IQR <- rle %>% group_by(Sample) %>% summarise(IQR = IQR(TPM))
rle_IQR_range <- rle_IQR$IQR %>% range %>% abs() %>% max()
rle_IQR_range <- 2*rle_IQR_range %>% ceiling()
bymedian <- with(rle, reorder(Sample, TPM, IQR))  # sort by IQR
par(mar=c(3,3,3,3))

pdf(file = "gene_RLEplot.pdf")
boxplot(TPM ~ bymedian, data=rle, outline=F, ylim = c(-rle_IQR_range, rle_IQR_range), las=2, boxwex=1, col='gray', cex.axis=0.3, main="RLE plot before QC", xlab="", ylab="Residual expression levels", frame=F)
dev.off()

#get the outlier list
RLEfilter_start <- as.integer(ncol(TPM_gene)-round(RLEFilterLength)+1)
RlEfilter_end <- as.integer(ncol(TPM_gene))
RLEFilterList <- levels(bymedian)[RLEfilter_start:RlEfilter_end]
RLEFilterList <- as.character(RLEFilterList)

message(paste0("The right most ", RLEFilterPercent*100, "% samples (N = ", length(RLEFilterList), ") are marked as candidate outliers in RLE (gene):") )
message(RLEFilterList)

###method 3: D-statistics
#get d-statistics outlier list 
sampleDists <- 1 - cor(logtpm, method='spearman')
hc <- hclust(as.dist(sampleDists), method = "complete")
hcphy <- as.phylo(hc)
sampleDists <- 1 - cor(logtpm, method='spearman')
hc <- hclust(as.dist(sampleDists), method = "complete")
hcphy <- as.phylo(hc)
D <- apply(1-sampleDists, 1, median)

pdf(file = "gene_D_stat_hist.pdf")
hist(D, breaks=100, ylab="Number of samples", xlab="D-statistic", main="Histogram of Sample D-statistics before data QC")
dev.off()

DSFilter <- sort(D)[DSFilter]
D<-as.data.frame(D)
D<-data.frame(Sample = rownames(D),D = D$D)
D_filterList = D%>%filter(D <= DSFilter)
D_filterList <- D_filterList$Sample
D_filterList<-as.character(D_filterList)

message(paste0("The right most ", DSFilterPercent*100, "% samples (N=", length(D_filterList), ") are marked as candidate outliers in D-statistics (gene):") )
message(D_filterList)

#####Step 4: save outliers into SE
PCA_outliers_df<-as.data.frame(x$outlier)
PCA_outliers_df$SAMPLE_ID<-rownames(PCA_outliers_df)
PCA_outliers_df$PCAoutliers<-c("Yes")
PCA_outliers_df<-PCA_outliers_df[,-1]
PCA_df<-as.data.frame(SE_gene@colData$SAMPLE_ID)
colnames(PCA_df)<-"SAMPLE_ID" 
PCA_df<-merge(PCA_df,PCA_outliers_df, all=TRUE)
PCA_df[is.na(PCA_df)] <- "No"
SE_gene@colData<-merge(PCA_df,SE_gene@colData,by="SAMPLE_ID")

rle_outliers_df<-as.data.frame(RLEFilterList)
rle_outliers_df$RLEoutliers<-c("Yes")
colnames(rle_outliers_df)<-c("SAMPLE_ID","RLEoutliers")
rle_df<-as.data.frame(SE_gene@colData$SAMPLE_ID)
colnames(rle_df)<-"SAMPLE_ID"
rle_df<-merge(rle_df,rle_outliers_df, all=TRUE)
rle_df[is.na(rle_df)] <- "No"
SE_gene@colData<-merge(rle_df,SE_gene@colData,by="SAMPLE_ID")

ds_outliers_df<-as.data.frame(D_filterList)
ds_outliers_df$Doutliers<-c("Yes")
colnames(ds_outliers_df)<-c("SAMPLE_ID","Doutliers")
ds_df<-as.data.frame(SE_gene@colData$SAMPLE_ID)
colnames(ds_df)<-"SAMPLE_ID"
ds_df<-merge(ds_df,ds_outliers_df, all=TRUE)
ds_df[is.na(ds_df)] <- "No"
SE_gene@colData<-merge(ds_df,SE_gene@colData,by="SAMPLE_ID")

saveRDS(SE_gene,"Gene_Expression.rds")

#--------------------------------isoform---------------------------------------#
######Step 1: load data
SE_isoform <- readRDS(args[2])
TPM_isoform <- SE_isoform@assays@data@listData[["TPM"]]

#####Step 2: filter low-expressed isoforms
keep_isoform_idx <- (rowMeans(TPM_isoform>low_expr_TPM)>low_expr_TPM_percent) 
TPM_isoform <- TPM_isoform[keep_isoform_idx,]

keep_isoform_idx <- keep_isoform_idx[1:nrow(assay(SE_isoform))]
SE_isoform <- SE_isoform[keep_isoform_idx, ]

#####Step 3: detect outliers using PCA, RLE, D-statistics
###method 1: PCA 
#prepare gene data for PCA
TPM_isoform_t<-as.data.frame(t(TPM_isoform))
isoform_pca_log<-log10(TPM_isoform_t+log_offset)

y<-pca.outlier(isoform_pca_log, center = TRUE, scale =FALSE, conf.level = 0.975) 

message(paste0("These samples are marked as candidate outliers in PCA step:") )
message(y$outlier)

pdf(file = "isoform_pca.pdf")
y$plot
dev.off()

###method 2: RLE
RLEFilterLength <- RLEFilterPercent*ncol(TPM_isoform)
DSFilter <- DSFilterPercent*ncol(TPM_isoform)

logtpm <- log10(TPM_isoform + log_offset)
rle <- logtpm-apply(logtpm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.
iqr <- apply(rle,2,IQR)
rle <- melt( rle, variable.name = "Sample",value.name ="TPM", id="ID")
names(rle) <- c("feature","Sample","TPM")
rle_IQR <- rle %>% group_by(Sample) %>% summarise(IQR = IQR(TPM))
rle_IQR_range <- rle_IQR$IQR %>% range %>% abs() %>% max()
rle_IQR_range <- 2*rle_IQR_range %>% ceiling()
bymedian <- with(rle, reorder(Sample, TPM, IQR))  # sort by IQR
par(mar=c(3,3,3,3))

pdf(file = "isoform_RLEplot.pdf")
boxplot(TPM ~ bymedian, data=rle, outline=F, ylim = c(-rle_IQR_range, rle_IQR_range), las=2, boxwex=1, col='gray', cex.axis=0.3, main="RLE plot before QC", xlab="", ylab="Residual expression levels", frame=F)
dev.off()

#get rle outlier list
RLEfilter_start <- as.integer(ncol(TPM_isoform)-round(RLEFilterLength)+1)
RlEfilter_end <- as.integer(ncol(TPM_isoform))
RLEFilterList <- levels(bymedian)[RLEfilter_start:RlEfilter_end]
RLEFilterList <- as.character(RLEFilterList)

message(paste0("The right most ", RLEFilterPercent*100, "% samples (N = ", length(RLEFilterList), ") are marked as candidate outliers in RLE (isoform):") )
message(RLEFilterList)


###method 3: D-stat
#get d-statistics outlier list
sampleDists <- 1 - cor(logtpm, method='spearman')
hc <- hclust(as.dist(sampleDists), method = "complete")
hcphy <- as.phylo(hc)
sampleDists <- 1 - cor(logtpm, method='spearman')
hc <- hclust(as.dist(sampleDists), method = "complete")
hcphy <- as.phylo(hc)
D <- apply(1-sampleDists, 1, median)

pdf(file = "isoform_D_stat_hist.pdf")
hist(D, breaks=100, ylab="Number of samples", xlab="D-statistic", main="Histogram of Sample D-statistics before data QC")
dev.off()

DSFilter <- sort(D)[DSFilter]
D<-as.data.frame(D)
D<-data.frame(Sample = rownames(D),D = D$D)
D_filterList = D%>%filter(D <= DSFilter)
D_filterList <- D_filterList$Sample
D_filterList<-as.character(D_filterList)

message(paste0("The right most ", DSFilterPercent*100, "% samples (N=", length(D_filterList), ") are marked as candidate outliers in D-statistics (isoform):") )
message(D_filterList)

#####Step 4: save result into SE
#isoform
PCA_outliers_df<-as.data.frame(y$outlier)
PCA_outliers_df$SAMPLE_ID<-rownames(PCA_outliers_df)
PCA_outliers_df$PCAoutliers<-c("Yes")
PCA_outliers_df<-PCA_outliers_df[,-1]
PCA_df<-as.data.frame(SE_isoform@colData$SAMPLE_ID)
colnames(PCA_df)<-"SAMPLE_ID" 
PCA_df<-merge(PCA_df,PCA_outliers_df, all=TRUE)
PCA_df[is.na(PCA_df)] <- "No"
SE_isoform@colData<-merge(PCA_df,SE_isoform@colData,by="SAMPLE_ID")

rle_outliers_df<-as.data.frame(RLEFilterList)
rle_outliers_df$RLEoutliers<-c("Yes")
colnames(rle_outliers_df)<-c("SAMPLE_ID","RLEoutliers")
rle_df<-as.data.frame(SE_isoform@colData$SAMPLE_ID)
colnames(rle_df)<-"SAMPLE_ID"
rle_df<-merge(rle_df,rle_outliers_df, all=TRUE)
rle_df[is.na(rle_df)] <- "No"
SE_isoform@colData<-merge(rle_df,SE_isoform@colData,by="SAMPLE_ID")

ds_outliers_df<-as.data.frame(D_filterList)
ds_outliers_df$Doutliers<-c("Yes")
colnames(ds_outliers_df)<-c("SAMPLE_ID","Doutliers")
ds_df<-as.data.frame(SE_isoform@colData$SAMPLE_ID)
colnames(ds_df)<-"SAMPLE_ID"
ds_df<-merge(ds_df,ds_outliers_df, all=TRUE)
ds_df[is.na(ds_df)] <- "No"
SE_isoform@colData<-merge(ds_df,SE_isoform@colData,by="SAMPLE_ID")

saveRDS(SE_isoform,"Isoform_Expression.rds")

