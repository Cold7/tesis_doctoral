library ( magrittr ) # this will allow us to string commands together in a UNIX - pipe - like fashion using % >%
setwd("/home/scontreras/Desktop/rnaseq_raton/ultimo_analisis/conteo/pca_genes_all_all/intestino_bazo_todo")
#get the table of read counts
readcounts <- read.table("all.count", header = TRUE, sep=",")
head(readcounts, n=1)
#One of the requirements of the assay () slots is that the row.names correspond to the gene IDs and the col.names to the sample names
row.names(readcounts)<- readcounts$Geneid
#str(row.names(readcounts))

##in addition , we need to exclude all columns that do not contain read counts
readcounts <- readcounts [,c(2:14)]
#readcounts
head(readcounts,n=3)


names(readcounts) <- c("intestino_colitis_1","intestino_colitis_2","intestino_colitis_3","intestino_control_1","intestino_control_2","intestino_control_3","bazo_colitis_1","bazo_colitis_2","bazo_colitis_3","bazo_colitis_4","bazo_control_1","bazo_control_2","bazo_control_3")

##verifying data
str(readcounts)
head(readcounts,n=3)

# make a data.frame with metadata where row.names should match the individual sample names
orig_names <- names(readcounts)
orig_names
sample_info <- data.frame(condition=gsub("_[0-9]+", "", names(readcounts)), row.names=names(readcounts))
sample_info

##verifying data
str(readcounts)
head(readcounts,n=3)

#importing deseq2 and generating deseqdataset
library ( DESeq2 )
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts,
                                   colData = sample_info,
                                   design = ~ condition)

## you can check the result using the accessors described above :
colData (DESeq.ds ) %>% head
assay(DESeq.ds , "counts" ) %>% head
rowData(DESeq.ds) %>% head

## test what counts() returns
counts(DESeq.ds) %>% str

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

## if you check colData () again , you see that this now contains the sizeFactors
colData(DESeq.ds)

#counts() allows you to immediately retrieve the _normalized_ read counts
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE )
counts.sf_normalized
#transform size-factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)

#You can see how the log2 transformation makes even simple graphs more easily interpretable by generating boxplots of read counts
#par( mfrow =c(2 ,1) ) # to plot the following two images underneath each other

##first, boxplots of non-transformed read counts (one per sample )
#boxplot ( counts.sf_normalized , notch = TRUE ,main = "untransformed read counts" , ylab = "read counts" )
## box plots of log2 - transformed read counts
#boxplot ( log.norm.counts , notch = TRUE ,main = "log2-transformed read counts",ylab= "log2(read counts)")

#Visually exploring normalized read counts
#plot (log.norm.counts[ ,1:2] , cex =.1 , main = "Normalized log2(read counts)" )

# obtain regularized log - transformed values

DESeq.rlog <- rlog(DESeq.ds, blind = FALSE ) 
#This function transforms the count data to the log2 scale in a way which minimizes 
# differences between samples for rows with small counts, and which normalizes with 
#respect to library size. The rlog transformation produces a similar variance 
#stabilizing effect 
rlog.norm.counts <- assay(DESeq.rlog)
head(DESeq.rlog)
## mean -sd plot for rlog - transformed data
##library (vsn)
##library (ggplot2)
##msd_plot <- meanSdPlot ( rlog.norm.counts, ranks = FALSE , # show the data on the original scale
##                          plot = FALSE )
##msd_plot$gg + ggtitle ( " rlog - transformed read counts " ) +ylab ( " standard deviation " )

#doing Hierarchical clustering
# cor () calculates the correlation between columns of a matrix

#metodo pearson
distance.m_rlog <- as.dist (1 - cor(rlog.norm.counts , method = "pearson" ) )
distance.m_rlog
cor(rlog.norm.counts , method = "pearson" )
#plot () can directly interpret the output of hclust ()

#los siguientes algoritmos de clustering jerarquico fueron extraidos de https://rpubs.com/mjimcua/clustering-jerarquico-en-r
plot(hclust(distance.m_rlog, method="single"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
plot(hclust(distance.m_rlog, method="complete"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
plot(hclust(distance.m_rlog, method="average"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
plot(hclust(distance.m_rlog, method="ward"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )

#metodo: kendall
#distance.m_rlog <- as.dist (1 - cor(rlog.norm.counts , method = "kendall" ) )
#distance.m_rlog
##plot () can directly interpret the output of hclust ()

##los siguientes algoritmos de clustering jerarquico fueron extraidos de https://rpubs.com/mjimcua/clustering-jerarquico-en-r
#plot(hclust(distance.m_rlog, method="single"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
#plot(hclust(distance.m_rlog, method="complete"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
#plot(hclust(distance.m_rlog, method="average"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
#plot(hclust(distance.m_rlog, method="ward"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )

library(philentropy)
distance.m_rlog <- distance(cor(rlog.norm.counts , method = "pearson" ), method = "euclidean" )
distance.m_rlog
plot(hclust(distance.m_rlog, method="single"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )

plot(hclust(distance.m_rlog, method="complete"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
plot(hclust(distance.m_rlog, method="average"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )
plot(hclust(distance.m_rlog, method="ward"), labels = colnames(rlog.norm.counts ), main = "rlog transformed read counts\ndistance: Pearson correlation" )


library (NMF)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = cor(rlog.norm.counts , method = "pearson" ), col = col, symm = TRUE )
aheatmap (cor(rlog.norm.counts , method = "pearson" ),Rowv = TRUE , Colv = TRUE , distfun = "correlation" , hclustfun = "average", scale="none", filename="./heatmap_corr_intestino.png" )

#doing Hierarchical clustering
# cor () calculates the correlation between columns of a matrix
distance.m_rlog <- as.dist (1 - cor(out , method = "pearson" ) )
distance.m_rlog
cor(out , method = "pearson" )
#plot () can directly interpret the output of hclust ()
plot(hclust(distance.m_rlog), labels = colnames(out ), main = "rlog transformed read counts\ndistance: Pearson correlation" )


####################################
#Principal Components Analysis (PCA)
####################################
library (vsn)
library (ggplot2)
pc <- prcomp(t(rlog.norm.counts)) #prcomp perform the pca
pc$sdev[1]
pc$sdev[2]
pc$sdev[3]
pc$sdev[4]
total_var <- 0
pc$sdev

for (num in pc$sdev) {
  #num*num due current sdev is the standar deviation, so we are looking for the variance
  total_var <- total_var + (num*num) 
}  

percentage_var <- c()
i <-1
while(i <= 13) {
  print(pc$sdev[i])
  percentage_var[i] <- round(((pc$sdev[i]*pc$sdev[i])/total_var)*100, digits = 6)
  i <- i+1
}
plot(x=percentage_var, main="% de varianza por cada\n componente principal", xlab="PC", ylab="% de varianza", type ="h")
i <- 2
j <- 3
plot(pc$x[,i], pc$x[,j],
     col = colData(DESeq.ds)[,1],
     xlab = paste("PC",i,":",round(((pc$sdev[i]*pc$sdev[i])/total_var)*100, digits = 0),"% of variance"),
     ylab= paste("PC",j,":",round(((pc$sdev[j]*pc$sdev[j])/total_var)*100, digits= 0),"% of variance")      
)
par( mfrow =c(3 ,2) ) # to plot the following two images underneath each other
i <- 1
while (i <= 3) {
  j <- 1
  while(j <= 3){
    if (i!=j){
      plot(pc$x[,i], pc$x[,j],
           col = colData(DESeq.ds)[,1],
           xlab = paste("PC",i,":",round(((pc$sdev[i]*pc$sdev[i])/total_var)*100, digits = 0),"% of variance"),
           ylab= paste("PC",j,":",round(((pc$sdev[j]*pc$sdev[j])/total_var)*100, digits= 0),"% of variance")      
      )
    }
    j <- j+1
    print(i,j)
  }
  i <- i+1
}
library(rgl)
mycolors <- c('royalblue1', 'black', 'green', 'red')
color <- mycolors[as.numeric(colData(DESeq.ds)[,1])]
color
colData(DESeq.ds)[,1]
open3d()
plot3d(pc$x[,1],pc$x[,2],pc$x[,3],
       col = color, 
       type = 's',
       radius = 2,
       xlab = paste("PC1:",round(((pc$sdev[1]*pc$sdev[1])/total_var)*100, digits = 0),"% of variance"),
       ylab= paste("PC2:",round(((pc$sdev[2]*pc$sdev[2])/total_var)*100, digits= 0),"% of variance"),
       zlab = paste("PC3:",round(((pc$sdev[3]*pc$sdev[3])/total_var)*100, digits= 0),"% of variance")
)
legend3d("topright", legend = c('bazo_colitis', 'bazo_control', 'colon_colitis', 'colon_control'), pch = 16, col = mycolors, cex=.8, inset=c(0.02))
snapshot3d(filename = '3dplot.png', fmt = 'png')
writeWebGL( filename="./3dscatter.html" ,  width=600, height=600)

