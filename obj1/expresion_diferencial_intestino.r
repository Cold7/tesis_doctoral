library ( magrittr ) # this will allow us to string commands together in a UNIX - pipe - like fashion using % >%

#get the table of read counts
readcounts <- read.table("intestino.count", header = TRUE, sep=",")
#One of the requirements of the assay () slots is that the row.names correspond to the gene IDs and the col.names to the sample names
row.names(readcounts)<- readcounts$Geneid
#str(row.names(readcounts))

##in addition , we need to exclude all columns that do not contain read counts
readcounts <- readcounts [,c(2:7)]
#readcounts
#head(readcounts,n=3)


names(readcounts) <- c("intestino_colitis_1","intestino_colitis_2","intestino_colitis_3","intestino_control_1","intestino_control_2","intestino_control_3")

# make a data.frame with metadata where row.names should match the individual sample names
orig_names <- names(readcounts)
#orig_names
sample_info <- data.frame(condition=gsub("_[0-9]+", "", names(readcounts)), row.names=names(readcounts))
#sample_info

##verifying data
#str(readcounts)
#head(readcounts,n=3)

#importing deseq2 and generating deseqdataset
library ( DESeq2 )
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts,
                                   colData = sample_info,
                                   design = ~ condition)

# remove genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0 , ]

## investigate different library sizes
#colSums(counts(DESeq.ds)) # should be the same as colSums ( readcounts )

#############################################
#DESeq2’s default method to normalize read
#counts to account for differences in 
#sequencing depths is implemented in 
#estimateSizeFactors()
#############################################
#calculate the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors (DESeq.ds)
sizeFactors(DESeq.ds)

#counts() allows you to immediately retrieve the _normalized_ read counts
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE )

#transform size-factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)

# obtain regularized log - transformed values
DESeq.rlog <- rlog(DESeq.ds, blind = FALSE ) #This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. The rlog transformation produces a similar variance stabilizing effect 
rlog.norm.counts <- assay(DESeq.rlog)



#genes DE

###########################
##
## DESeq2
##
###########################

# DESeq2 uses the levels of the condition to determine the order of the comparison
# str(colData(DESeq.ds)$condition)
# set bazo_control as the first-level-factor
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition,"intestino_control")

#Now, running the DGE analysis is very simple:
DESeq.ds <- DESeq(DESeq.ds)

#The results() function lets you extract the base means across samples, moderated log2 fold changes, standard errors, test statistics etc. for every gene.
DGE.results <-results(DESeq.ds, independentFiltering = TRUE , alpha = 0.05)
DGE.volcano <-results(DESeq.ds)
##############################################
#saving data
res <- DGE.volcano[order(DGE.volcano$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(DESeq.ds,normalized =TRUE)), by = 'row.names', sort = FALSE)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj
write.csv(resdata, file="intestino_deseq.csv", sep='\t')

#############################################
#summary (DGE.results)
## the DESeqResult object can basically be handled like a data.frame
#table (DGE.results$padj < 0.05)
#table
#rownames(subset(DGE.results, padj < 0.05))

#library(EnhancedVolcano)
#EnhancedVolcano(DGE.volcano, lab=rownames(DGE.volcano), x="log2FoldChange", y = "padj", xlim = c(-10,10), pCutoff = 0.05, ylab="Padj")

##heatmap
library (NMF)
## aheatmap needs a matrix of values , e.g. , a matrix of DE genes with the transformed read counts for each replicate 
##sort the results according to the adjusted p- value
DGE.results.sorted <- DGE.results[order(DGE.results$padj),]

## identify genes with the desired adjusted p- value cut -off
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))
#this genelist came from the result of deseq, edger and limma, ig the gene was detected by 2 libs, then it form part of this list
geneList <- c('ENSMUSG00000027556', 'ENSMUSG00000048878', 'ENSMUSG00000053977', 'ENSMUSG00000074272', 'ENSMUSG00000022508', 'ENSMUSG00000039220', 'ENSMUSG00000050272', 'ENSMUSG00000013419', 'ENSMUSG00000030364', 'ENSMUSG00000040284', 'ENSMUSG00000035441', 'ENSMUSG00000076471', 'ENSMUSG00000013653', 'ENSMUSG00000076472')
geneList <- c('ENSMUSG00000027556','ENSMUSG00000048878','ENSMUSG00000053977','ENSMUSG00000074272','ENSMUSG00000022508','ENSMUSG00000039220')



##extract the normalized read counts for DE genes into a matrix
hm.mat_DGEgenes <- log.norm.counts[geneList,]
rownames(hm.mat_DGEgenes) <-c("Car1","Hexim1","Cd8a","Ceacam1","Bcl6","Ppp1r10")
rownames(hm.mat_DGEgenes)
## plot the normalized read counts of DE genes sorted by the adjusted p- value
aheatmap (hm.mat_DGEgenes , Rowv = NA , Colv = NA,scale = "row" )

# combine the heatmap with hierarchical clustering and scale ,the read counts per gene to emphasize the sample -type - specific differences

aheatmap (hm.mat_DGEgenes ,Rowv = TRUE , Colv = TRUE , distfun = "euclidean" , hclustfun = "average" ,scale = "row" )
# values are transformed into distances from the center of the row - specific average : ( actual value - mean of the group ) / standard deviation


###########################
##
## EdgeR
##
###########################
#BiocManager::install("edgeR") #
library(edgeR)
sample_info.edger <-factor(c(rep("colitis",3), rep("control",3)))
sample_info.edger <-relevel(sample_info.edger,ref="control")
#Now, DGEList() is the function that converts the count matrix into an edgeR object.
edgeR.DGElist <- DGEList(counts = readcounts , group = sample_info.edger)

#edgeR also recommends removing genes with almost no coverage. In order to determine a sensible cutoff, we
#plot a histogram of counts per million calculated by edgeR’s cpm() function

#hist(log2(rowSums(cpm(edgeR.DGElist))))
#summary(log2(rowSums(cpm(edgeR.DGElist))))

# remove genes that do not have one count per million in the 6 samples

keep <- rowSums(cpm(edgeR.DGElist)>= 1)>=6
edgeR.DGElist <- edgeR.DGElist[keep,]

# recompute library sizes after filtering
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)

#Calculate normalization factors for the library sizes. We use the standard edgeR method here, which is
#the trimmed mean of M-values; if you wanted to use, for example, DESeq’s size factor, you could use
#method = "RLE"). See Table 13 for details of the methods.
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method="TMM")


# specify the design setup - the design matrix looks a bit intimitating , but if
# you just focus on the formula [~ sample_info.edger ] you can see that it 's
# exactly what we used for DESeq2 , too
design <- model.matrix(~sample_info.edger)
# estimate the dispersion for all read counts across all samples
edgeR.DGElist <- estimateDisp(edgeR.DGElist, design)
# fit the negative binomial model
edger_fit <- glmFit(edgeR.DGElist, design)

# perform the testing for every gene using the neg. binomial model
edger_lrt<- glmLRT(edger_fit)

#extract results from edger _lrt$ table plus adjusted p- values
DGE.results_edgeR <- topTags(edger_lrt, n=Inf, # to retrieve all genes
                                   sort.by ="PValue", adjust.method="BH")
write.csv(DGE.results_edgeR, file="intestino_edgeR.csv", sep='\t')


###########################
##
## limma+vomm
##
###########################
library(limma)

#Like DESeq and edgeR, limma starts with a matrix of raw read counts where each gene is represented by a
#row and the columns represent samples. limma assumes that rows with zero or very low counts have been
#removed. In addition, size factors for sequencing depth can be calculated using edgeR’s calcNormFactors()
#function.
sample_info.edger <-factor(c(rep("colitis",3), rep("control",3)))
sample_info.edger <-relevel(sample_info.edger,ref="control")
edgeR.DGElist <- DGEList(counts = readcounts , group = sample_info.edger)
keep <- rowSums(cpm(edgeR.DGElist)>= 1)>=6
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method="TMM")

# limma also needs a design matrix , just like edgeR
design <- model.matrix(~sample_info.edger)

# transform the count data to log2 -counts -per - million and estimate
# the mean - variance relationship , which is used to compute weights
# for each count -- this is supposed to make the read counts
# amenable to be used with linear models

rownames(design) <- colnames(edgeR.DGElist)
voomTransformed <- voom(edgeR.DGElist, design, plot=FALSE)

# fit a linear model for each gene
voomed.fitted <- lmFit(voomTransformed, design=design)
# compute moderated t- statistics , moderated F- statistics ,
# and log - odds of differential expression
voomed.fitted <- eBayes(voomed.fitted)
# extract gene list with logFC and statistical measures
colnames(design) # check how the coefficient is named

DGE.results_limma <- topTable(voomed.fitted, coef = "sample_info.edgercolitis",
                                    number=Inf, adjust.method="BH",
                                    sort.by = "P")
write.csv(DGE.results_limma, file="intestino_limma.csv", sep='\t')

# make a Venn diagram
DE_list <- list ( edger = rownames (subset(DGE.results_edgeR$table , FDR <=0.05)) ,
                   deseq2 = rownames (subset(DGE.results, padj <=0.05)) ,
                   limma = rownames(subset(DGE.results_limma, adj.P.Val <=0.05)))
#install.packages("gplots")
library(gplots)
gplots::venn(DE_list)
