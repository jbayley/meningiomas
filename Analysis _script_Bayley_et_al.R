###################
#Environment setup#
###################

#Load functions,  packages, and clinical/genomic data used in this analysis
source('Functions_Bayley_et_al.R')
loadPackages()
SampleData <- read.csv('Sample_Data_Bayley_et_al.csv', row.names = 1)
CPG <- read.csv("CPG information.csv", row.names = 1)

#Set colors for figure production
colors_by_group <- c('#6FDE6E', '#235FA4', '#FF4242', '#BDD9BF', '#2E4052', '#E5323B')
names(colors_by_group) <- c('GroupA', 'GroupB', 'GroupC', 'GroupA2', 'GroupB2', 'GroupC2')

#Load methylation signatures we derived (if deriving yourself, these will be replaced once those signatures are defined)
MethSig_cluster1 <- read.csv("Methylation signature MenG A.csv")$CPG_probe
MethSig_cluster2 <- read.csv("Methylation signature MenG B.csv")$CPG_probe
MethSig_cluster3 <- read.csv("Methylation signature MenG C.csv")$CPG_probe
MethSig <- c(MethSig_cluster1, MethSig_cluster2, MethSig_cluster3)

#Directly load the processed data from provided .csv files (see immediately below if you would prefer to process raw data yourself)
#myNorm is the beta matrix obtained from ChAMP pipeline
myNorm <- as.matrix(read.csv("ChAMP beta matrix.csv", row.names = 1))
colnames(myNorm) <- gsub("X", "", gsub("\\.", "-", colnames(myNorm)))

#sesameBetas is the beta matrix obtained from the SeSAMe pipeline 
sesameBetas <- as.matrix(read.csv("SeSAMe beta matrix.csv", row.names = 1))
colnames(sesameBetas) <- gsub("X", "", gsub("\\.", "-", colnames(sesameBetas)))

#expression is the log transformed RNA-seq data from all samples and the three arachnoid controls from GSE139652
expression <- as.matrix(read.csv("Log expression data.csv", row.names = 1))
colnames(expression) <- gsub("X", "", gsub("\\.", "-", colnames(expression)))

####################################
#If desired, use the code provided below to load and process the raw data 
#DNA methylation data can be downloaded from GEO (Series GSE189521) and placed into a folder 'Methylation' in the working directory
#Import methylation data
myLoad <- champ.load(directory = paste0(getwd(), '/Methylation'), arraytype = 'EPIC')

#Can run quality check, but no issues or any bad samples
#champ.QC()

#Normalize methylation data (also sorts probes in decreasing order of overall variance)
myNorm <- normalizeMethylation()

#Can check for batch effect with SVD, but there is none and so no need for COMBAT
#champ.SVD()

#Load data with SeSAMe pipeline if desired
sesameBetas <- loadSesame()

#Raw counts data is provided in a .csv file and can be processed as below, while fastq files are available on GEO if desired (Series GSE189672)
#Load RNA-seq raw counts and generate transformed expression data
counts <- read.csv('Raw counts expression data.csv', row.names = 1)
colnames(counts) <- gsub('X', '', gsub('(\\.)', '-', colnames(counts)))
rownames(counts) <-  gsub('(\\.).*', '', rownames(counts))
colData <- data.frame(row.names = colnames(counts), type = rep(0, dim(counts)[2]))
colData[rownames(SampleData),] <- SampleData$MenG
colData[grep('Arachnoid', rownames(colData)),] <- 'Arachnoid'
colData$type <- factor(colData$type)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~type)
dds <- DESeq(dds[rowSums(counts(dds)) >= 10, ])
expression <- getExpressionData(dds)
####################################

##############################################
#Establishment of meningioma epigenetic types#
##############################################

#Select the ChAMP or Sesame dataset for the following analysis (N.B., we used ChAMP as our baseline, with SeSAMe as a confirmatory processing method)
methData <- myNorm
#or
methData <- sesameBetas

#Minimum condition calculation to determine number of CpGs for the small subset
kappa <- matrix(nrow = 1000, ncol = 2)
for(val in (1:1000)) {
  kappa[val,1] <- val*10
  kappa[val,2] <- kappa(methData[1:(val*10),])
}

#Curve approaches asymptote right around 4k, so use the local min from 3500-5000 as our initial number of DMPs (4,230)
data <- as.data.frame(kappa[50:500,])
theta.0 <- min(data$V2) * 0.5  
model.0 <- lm(log(V2 - theta.0) ~ V1, data=data)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
model <- nls(V2 ~ alpha * exp(beta * V1) + theta , data = data, start = start)
plot(data$V1, data$V2, ylab = 'Condition Number', xlab = 'Number of Differentially Methylated Probes', pch = 1)
lines(data$V1, predict(model, list(x = data$V1)), col = 'skyblue', lwd = 4)
remove(data, theta.0, model.0, alpha.0, beta.0, start, model)
plot(kappa[(350):(500),], type = 'p', ylab = 'Condition Number', xlab = 'Number of Differentially Methylated Probes')
nDMPs <- kappa[kappa[,2]==min(kappa[(350):(500),]),1]


###Identify optimal number of meningioma clusters

#Perform NMF with this set of DMPs to determine optimal rank (note that myNorm is already sorted in order of decreasing CpG site variance)
rank_estimate_initial <- nmf(methData[1:nDMPs,], 2:10, nrun=100,  .opt = 'v')
plot(rank_estimate_initial, main = paste0('Rank estimate for top ', nDMPs, 'most variably methylated probes'))
plot(2:10, rank_estimate_initial$measures$cophenetic, main = paste0('Top ', nDMPs, ' Cophenetic coefficient'), type = 'b', xlab = '', ylab = '', col = 'purple', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
plot(2:10, rank_estimate_initial$measures$residuals, main = paste0('Top ', nDMPs, ' Residuals'), type = 'b', xlab = '', ylab = '', col = 'dark blue', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)

consensusmap(rank_estimate_initial, color = 'blue', main = NULL, annCol = NULL, annRow = NULL, tracks = NULL, info = FALSE, labRow = NA, labCol = NA, legend = NULL)

#Optimal rank appears to be 3, now repeat NMF over the top 10% of DMPs
top10DMPs <- floor(dim(methData)[1]/10)
#Running the following attempts to allocate a 40Gb vector when calculating RSS (causing a fail on my hardware), so run each rank individually
#rank_estimate_top10 <- nmf(methData[1:top10DMPs,], 2:10, nrun=100)
ranks <- 2:10
rank_estimate_top10 <- lapply(ranks, function(r) {
  fit <- nmf(methData[1:top10DMPs,], r, nrun = 100,  .opt = 'v')
  list(fit = fit, consensus = consensus(fit), coph = cophcor(fit), residual = residuals(fit))
})

plot(ranks, sapply(rank_estimate_top10, '[[', 'coph'), main = paste0('Top ', top10DMPs, ' Cophenetic coefficient'), type = 'b', xlab = '', ylab = '', col = 'purple', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
plot(ranks, sapply(rank_estimate_top10, '[[', 'residual'), main = paste0('Top ', top10DMPs, ' Residuals'), type = 'b', xlab = '', ylab = '', col = 'dark blue', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
dispersion <- sapply(sapply(rank_estimate_top10, '[[', 'fit', simplify = FALSE), dispersion, simplify = FALSE)
plot(ranks, dispersion, type = 'b', col = 'dark blue', pch = 19, cex = 2.1, lwd = 3, main = 'Dispersion', xlab = '', ylab = '', xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
silhouettes <- sapply(sapply(rank_estimate_top10, '[[', 'fit', simplify = FALSE), silhouette, simplify = FALSE)
silhouettes[[length(silhouettes)+1]] <- matrix(ncol = 2, nrow = length(silhouettes)) 
for(val in 1:(length(silhouettes)-1)) {
  silhouettes[[length(silhouettes)]][val,] <- c(max(silhouettes[[val]][,'cluster']), mean(silhouettes[[val]][,'sil_width']))
}
plot(silhouettes[[length(silhouettes)]], main = 'Mean Silhouette Width', type = 'b', xlab = '', ylab = '')

consensusmap(sapply(rank_estimate_top10, '[[', 'fit', simplify = FALSE),  color = 'blue', main = NULL, annCol = NULL, annRow = NULL, tracks = NULL, info = FALSE, labRow = NA, labCol = NA, legend = NULL)


#Recheck optimal number of groups using k-means clustering
Kmeans_initial <- ConsensusClusterPlus(methData[1:nDMPs,], maxK=10, reps=500, pItem=0.8, pFeature=1, title='K-means clustering top 4,230 DMPs', clusterAlg='km', distance='euclidean', plot='png')
Kmeans_top10 <- ConsensusClusterPlus(methData[1:top10DMPs,], maxK=10, reps=500, pItem=0.8, pFeature=1, title='K-means clustering top 10 percent of DMPs', clusterAlg='km', distance='euclidean', plot='png')

maps_initial <- list()
maps_top10 <- list()
for(val in 1:9) {
  maps_initial[[val]] <- Kmeans_initial[[val+1]]$consensusMatrix
  maps_top10[[val]] <- Kmeans_top10[[val+1]]$consensusMatrix
}

consensusmap(maps_initial, color = 'blue', main = NULL, legend = FALSE, labRow = NA, labCol = NA, Rowv = FALSE)
consensusmap(maps_top10, color = 'blue', main = NULL, legend = FALSE, labRow = NA, labCol = NA, Rowv = FALSE)

#3 appears to be the optimal number of groups for all analyses
Rank <- 3

###Perform NMF for rank 3 to define initial meningioma clusters
#Run NMF for top DMPs and export results (for Sesame run, switch the 3 and 2 for the cluster number assignment - order 3,2,1)
NMFinitial <- nmf(methData[1:nDMPs,], rank = Rank, nrun = 1000, .opt = 'v')
consensus_initial <- predict(NMFinitial, what = 'consensus')
silh_initial <- silhouette(NMFinitial)
clusters_initial <- cbind(rep(0, length(consensus_initial)), consensus_initial, silh_initial)
colnames(clusters_initial)[1] <- 'Epigenetic_cluster'
#Adjust numbering of epigenetic clusters
clusters_initial[clusters_initial[,'consensus_initial'] == 1,'Epigenetic_cluster'] <- 2
clusters_initial[clusters_initial[,'consensus_initial'] == 2,'Epigenetic_cluster'] <- 3
clusters_initial[clusters_initial[,'consensus_initial'] == 3,'Epigenetic_cluster'] <- 1

#Run NMF for top 10% DMPs and export results
NMFtop10 <- nmf(methData[1:top10DMPs, ], rank = Rank, nrun = 1000, .opt = 'v')
consensus_top10 <- predict(NMFtop10, what = 'consensus')
silh_top10 <- silhouette(NMFtop10)
clusters_top10 <- cbind(rep(0, length(consensus_top10)), consensus_top10, silh_top10)
colnames(clusters_top10)[1] <- 'Epigenetic_cluster'
#Adjust numbering of epigenetic clusters
clusters_top10[clusters_top10[,'consensus_top10'] == 1,'Epigenetic_cluster'] <- 3
clusters_top10[clusters_top10[,'consensus_top10'] == 2,'Epigenetic_cluster'] <- 2
clusters_top10[clusters_top10[,'consensus_top10'] == 3,'Epigenetic_cluster'] <- 1

#Get assignments from these 4 analyses and find samples with concordant clusters
clusters <- data.frame(row.names = rownames(clusters_initial), NMF_initial = clusters_initial[,'Epigenetic_cluster'], NMF_top10 = clusters_top10[,'Epigenetic_cluster'], Kmeans_initial = Kmeans_initial[[3]]$consensusClass, Kmeans_top10 = Kmeans_top10[[3]]$consensusClass)
concordant <- vector()
for(val in 1:dim(clusters)[1]) {
  if(var(as.numeric(clusters[val,])) == 0) {
    concordant <- c(concordant, rownames(clusters)[val])
  }
}

#Idenitfy those samples with concordant classifications among each of the three clusters
concordant_cluster1 <- intersect(concordant, rownames(clusters_initial[clusters_initial[,'Epigenetic_cluster']==1,]))
concordant_cluster2 <- intersect(concordant, rownames(clusters_initial[clusters_initial[,'Epigenetic_cluster']==2,]))
concordant_cluster3 <- intersect(concordant, rownames(clusters_initial[clusters_initial[,'Epigenetic_cluster']==3,]))

#Use champ.DMP to define DMPs, comparing each group against all other samples
concordant_Samples <- data.frame(row.names = concordant, cluster1 = clusters_initial[concordant, 1], cluster2 = clusters_initial[concordant, 1], cluster3 = clusters_initial[concordant, 1])
concordant_Samples$cluster1[concordant_Samples$cluster1!=1] <- 'other'
concordant_Samples$cluster2[concordant_Samples$cluster2!=2] <- 'other'
concordant_Samples$cluster3[concordant_Samples$cluster3!=3] <- 'other'

#Identify DMPs for a given adjusted p-value
pval <- 1e-06
DMP_cluster1 <- champ.DMP(beta = methData[,concordant], pheno=concordant_Samples$cluster1, adjPVal = pval, arraytype = 'EPIC')
DMP_cluster2 <- champ.DMP(beta = methData[,concordant], pheno=concordant_Samples$cluster2, adjPVal = pval, arraytype = 'EPIC')
DMP_cluster3 <- champ.DMP(beta = methData[,concordant], pheno=concordant_Samples$cluster3, adjPVal = pval, arraytype = 'EPIC')

#Chose those DMPs with a difference in mean beta between a group and the other samples > a cutoff (we selected 0.25)
DMP_beta_cutoff <- 0.25
MethSig_cluster1 <- rownames(DMP_cluster1[[1]][abs(DMP_cluster1[[1]]$deltaBeta ) > DMP_beta_cutoff,])
MethSig_cluster2 <- rownames(DMP_cluster2[[1]][abs(DMP_cluster2[[1]]$deltaBeta ) > DMP_beta_cutoff,])
MethSig_cluster3 <- rownames(DMP_cluster3[[1]][abs(DMP_cluster3[[1]]$deltaBeta ) > DMP_beta_cutoff,])

#Remove DMPs found for multiple groups
overlap <- c(intersect(MethSig_cluster1, MethSig_cluster2), intersect(MethSig_cluster1, MethSig_cluster3), intersect(MethSig_cluster2, MethSig_cluster3))
MethSig_cluster1 <- MethSig_cluster1[!MethSig_cluster1 %in% overlap]
MethSig_cluster2 <- MethSig_cluster2[!MethSig_cluster2 %in% overlap]
MethSig_cluster3 <- MethSig_cluster3[!MethSig_cluster3 %in% overlap]
MethSig <- c(MethSig_cluster1, MethSig_cluster2, MethSig_cluster3)

####Now that we have the methylation signatures, perform final NMF using these probes
#Redo rank estimation, although this is largely teleological at this point
rank_estimate_final <- nmf(methData[MethSig,], 2:10, nrun = 100,  .opt = 'v')
plot(rank_estimate_final, main = 'Rank estimate for signature probes')
consensusmap(rank_estimate_final, color = 'blue', main = NULL, annCol = NULL, annRow = NULL, tracks = NULL, info = FALSE, labRow = NA, labCol = NA, legend = NULL)

#Indeed 3 is the optimal rank
NMFfinal <- nmf(methData[MethSig,], rank = Rank, nrun = 1000, .opt = 'v')
consensusmap(NMFfinal)
consensus_final <- predict(NMFfinal, what = 'consensus')
silh_final <- silhouette(NMFfinal)
clusters_final <- cbind(rep(0, length(consensus_final)), consensus_final, silh_final)
colnames(clusters_final)[1] <- 'Epigenetic_cluster'
#Adjust numbering of epigenetic clusters
clusters_final[clusters_final[,'consensus_final'] == 1,'Epigenetic_cluster'] <- 3
clusters_final[clusters_final[,'consensus_final'] == 2,'Epigenetic_cluster'] <- 2
clusters_final[clusters_final[,'consensus_final'] == 3,'Epigenetic_cluster'] <- 1

#Heatmap of signature probes
Group1Dend <- as.dendrogram(hclust(dist(methData[MethSig_cluster1,])))
Group1Dend <- color_branches(Group1Dend, k = 1, col = colors_by_group[4])
Group2Dend <- as.dendrogram(hclust(dist(methData[MethSig_cluster2,])))
Group2Dend <- color_branches(Group2Dend, k = 1, col = colors_by_group[5])
Group3Dend <- as.dendrogram(hclust(dist(methData[MethSig_cluster3,])))
Group3Dend <- color_branches(Group3Dend, k = 1, col = colors_by_group[6])
CPG_dendrogram <- merge(Group1Dend, Group2Dend, Group3Dend)

column_ha <- HeatmapAnnotation(MethType = clusters_final[,'Epigenetic_cluster'],
                               col = list(MethType = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242')))

#Heatmap ordered by epigenetic cluster 
Heatmap(methData[MethSig,], cluster_rows = CPG_dendrogram, row_split = 2, column_split = clusters_final[,'Epigenetic_cluster'], column_order = order.dendrogram(consensushc(NMFfinal, method = 'mcquitty')), show_row_dend = TRUE, top_annotation = column_ha, show_row_names = FALSE, show_column_names = FALSE, column_title = 'Methylation Types', row_title = 'Methylation Signatures', heatmap_legend_param = list(title = 'Beta'), use_raster = FALSE)
#Heatmap clustered using NMF clustering dendrogram
Heatmap(methData[MethSig,], cluster_rows = CPG_dendrogram, row_split = 2, column_split = 3, cluster_columns = consensushc(NMFfinal, method = 'mcquitty'), show_row_dend = TRUE, top_annotation = column_ha, show_row_names = FALSE, show_column_names = FALSE, column_title = 'Methylation Types', row_title = 'Methylation Signatures', heatmap_legend_param = list(title = 'Beta'), use_raster = FALSE)


#r-TSNE of each of our CpG probe subsets (run one of the following three subsets, then the code below to generate a given tSNE map)
tsneSet <- methData[1:nDMPs,]
title <- paste0('t-SNE of ', nDMPs, ' most variable probes')

tsneSet <- methData[1:top10DMPs,]
title <- paste0('t-SNE of ', top10DMPs, ' most variable probes')

tsneSet <- methData[MethSig,]
title <- paste0('t-SNE of methylation signature probes (n = ', dim(tsneSet)[1], ')')

#Run the following after running one of the three preceding subsets
colors <- clusters_final[,'Epigenetic_cluster']
colors[colors == '1'] <- colors_by_group[1]
colors[colors == '2'] <- colors_by_group[2]
colors[colors == '3'] <- colors_by_group[3]
tsne_output <- Rtsne(t(tsneSet), pca = FALSE, perplexity = 8, labels = colnames(methData))
plot(tsne_output$Y, col=colors, pch = 19, asp=0, cex = 2, xlab = '', ylab = '', main = title)
text(tsne_output$Y, labels=colnames(myNorm), col=colors)

###Compare methylation levels of signature probes between epigenetic clusters
signature_meth_levels <- matrix(data = NA, nrow = length(MethSig), ncol = 6, dimnames = list(MethSig, c('Cluster1Mean', 'Cluster2Mean', 'Cluster3Mean', 'Cluster1', 'Cluster2', 'Cluster3')))
signature_meth_levels[,1] <- rowMeans(methData[rownames(signature_meth_levels), rownames(subset(as.data.frame(clusters_final), Epigenetic_cluster == 1))])
signature_meth_levels[,2] <- rowMeans(methData[rownames(signature_meth_levels), rownames(subset(as.data.frame(clusters_final), Epigenetic_cluster == 2))])
signature_meth_levels[,3] <- rowMeans(methData[rownames(signature_meth_levels), rownames(subset(as.data.frame(clusters_final), Epigenetic_cluster == 3))])

#Assign each probe as being either hypermethylated (1) or hypomethylated (-1) in a single cluster.
#We do this by identifying the group for which the mean is most divergent (by first identifying the two groups with the closest means). Then we compare the mean of that most divergent cluster to the other two. 
for(val in 1:dim(signature_meth_levels)[1]) {
  minDelta <- min(abs(signature_meth_levels[val,1]-signature_meth_levels[val,2]), abs(signature_meth_levels[val,1]-signature_meth_levels[val,3]), abs(signature_meth_levels[val,2]-signature_meth_levels[val,3]))
  
  if(minDelta == abs(signature_meth_levels[val,1]-signature_meth_levels[val,2])) {
    signature_meth_levels[val,4] <- 0
    signature_meth_levels[val,5] <- 0
    if(signature_meth_levels[val,3] > signature_meth_levels[val,1] & signature_meth_levels[val,3] > signature_meth_levels[val,2]) {
      signature_meth_levels[val,6] <- 1
    } else if(signature_meth_levels[val,3] < signature_meth_levels[val,1] & signature_meth_levels[val,3] < signature_meth_levels[val,2]) {
      signature_meth_levels[val,6] <- -1
    } else {
      signature_meth_levels[val,6] <- 'Problem'
    }
  } else if(minDelta == abs(signature_meth_levels[val,1]-signature_meth_levels[val,3])) {
    signature_meth_levels[val,4] <- 0
    signature_meth_levels[val,6] <- 0
    if(signature_meth_levels[val,2] > signature_meth_levels[val,1] & signature_meth_levels[val,2] > signature_meth_levels[val,3]) {
      signature_meth_levels[val,5] <- 1
    } else if(signature_meth_levels[val,2] < signature_meth_levels[val,1] & signature_meth_levels[val,2] < signature_meth_levels[val,3]) {
      signature_meth_levels[val,5] <- -1
    } else {
      signature_meth_levels[val,5] <- 'Problem'
    }
  } else if(minDelta == abs(signature_meth_levels[val,2]-signature_meth_levels[val,3])) {
    signature_meth_levels[val,5] <- 0
    signature_meth_levels[val,6] <- 0
    if(signature_meth_levels[val,1] > signature_meth_levels[val,2] & signature_meth_levels[val,1] > signature_meth_levels[val,3]) {
      signature_meth_levels[val,4] <- 1
    } else if(signature_meth_levels[val,1] < signature_meth_levels[val,2] & signature_meth_levels[val,1] < signature_meth_levels[val,3]) {
      signature_meth_levels[val,4] <- -1
    } else {
      signature_meth_levels[val,4] <- 'Problem'
    }
  } else {
    message(val)
  }
}

#Tables identifying the number of probes hypermethylated (1) or hypomethylated (2) in each cluster
#Epigenetic cluser 1
table(signature_meth_levels[,4])
#Epigenetic cluser 2
table(signature_meth_levels[,5])
#Epigenetic cluser 3
table(signature_meth_levels[,6])

#Now look at the subset of CpG promoters (defined as 1st exon, TSS 200, and TSS1500) and those on CpG islands within the signature
CPG_promoter <- rownames(subset(CPG[rownames(methData),], feature %in% c('1stExon', 'TSS200', 'TSS1500')))
CPG_promoter_island <- rownames(subset(CPG[rownames(methData),], feature %in% c('1stExon', 'TSS200', 'TSS1500') & cgi == 'island'))

#Tables for signature probes not on CpG promoters
table(signature_meth_levels[!rownames(signature_meth_levels) %in% CPG_promoter,4])
table(signature_meth_levels[!rownames(signature_meth_levels) %in% CPG_promoter,5])
table(signature_meth_levels[!rownames(signature_meth_levels) %in% CPG_promoter,6])

#Tables for signature probes on CpG promoters
table(signature_meth_levels[intersect(rownames(signature_meth_levels), CPG_promoter),4])
table(signature_meth_levels[intersect(rownames(signature_meth_levels), CPG_promoter),5])
table(signature_meth_levels[intersect(rownames(signature_meth_levels), CPG_promoter),6])

#Tables for signature probes on CpG island promoters
table(signature_meth_levels[intersect(rownames(signature_meth_levels), CPG_promoter_island),4])
table(signature_meth_levels[intersect(rownames(signature_meth_levels), CPG_promoter_island),5])
table(signature_meth_levels[intersect(rownames(signature_meth_levels), CPG_promoter_island),6])

#Now look at all CpG island promoters across the genome
promoterIslands <- matrix(data = NA, nrow = length(CPG_promoter_island), ncol = 7, dimnames = list(CPG_promoter_island, c('Cluster1Mean', 'Cluster2Mean', 'Cluster3Mean', 'Cluster1', 'Cluster2', 'Cluster3', 'MaxDelta')))
promoterIslands[,1] <- rowMeans(methData[rownames(promoterIslands), rownames(subset(as.data.frame(clusters_final), Epigenetic_cluster == 1))])
promoterIslands[,2] <- rowMeans(methData[rownames(promoterIslands), rownames(subset(as.data.frame(clusters_final), Epigenetic_cluster == 2))])
promoterIslands[,3] <- rowMeans(methData[rownames(promoterIslands), rownames(subset(as.data.frame(clusters_final), Epigenetic_cluster == 3))])

for(val in 1:dim(promoterIslands)[1]) {
  minDelta <- min(abs(promoterIslands[val,1]-promoterIslands[val,2]), abs(promoterIslands[val,1]-promoterIslands[val,3]), abs(promoterIslands[val,2]-promoterIslands[val,3]))
  promoterIslands[val,7] <- max(abs(promoterIslands[val,1]-promoterIslands[val,2]), abs(promoterIslands[val,1]-promoterIslands[val,3]), abs(promoterIslands[val,2]-promoterIslands[val,3]))
  
  if(minDelta == abs(promoterIslands[val,1]-promoterIslands[val,2])) {
    promoterIslands[val,4] <- 0
    promoterIslands[val,5] <- 0
    if(promoterIslands[val,3] > promoterIslands[val,1] & promoterIslands[val,3] > promoterIslands[val,2]) {
      promoterIslands[val,6] <- 1
    } else if(promoterIslands[val,3] < promoterIslands[val,1] & promoterIslands[val,3] < promoterIslands[val,2]) {
      promoterIslands[val,6] <- -1
    } else {
      promoterIslands[val,6] <- 'Problem'
    }
  } else if(minDelta == abs(promoterIslands[val,1]-promoterIslands[val,3])) {
    promoterIslands[val,4] <- 0
    promoterIslands[val,6] <- 0
    if(promoterIslands[val,2] > promoterIslands[val,1] & promoterIslands[val,2] > promoterIslands[val,3]) {
      promoterIslands[val,5] <- 1
    } else if(promoterIslands[val,2] < promoterIslands[val,1] & promoterIslands[val,2] < promoterIslands[val,3]) {
      promoterIslands[val,5] <- -1
    } else {
      promoterIslands[val,5] <- 'Problem'
    }
  } else if(minDelta == abs(promoterIslands[val,2]-promoterIslands[val,3])) {
    promoterIslands[val,5] <- 0
    promoterIslands[val,6] <- 0
    if(promoterIslands[val,1] > promoterIslands[val,2] & promoterIslands[val,1] > promoterIslands[val,3]) {
      promoterIslands[val,4] <- 1
    } else if(promoterIslands[val,1] < promoterIslands[val,2] & promoterIslands[val,1] < promoterIslands[val,3]) {
      promoterIslands[val,4] <- -1
    } else {
      promoterIslands[val,4] <- 'Problem'
    }
  } else {
    message(val)
  }
}

#Tables for signature probes on CpG island promoters (limit to those was a variance > 0.01 to exclude trivial probes with nearly identical methylation across clusters)
table(promoterIslands[promoterIslands[,7]>0.01,4])
table(promoterIslands[promoterIslands[,7]>0.01,5])
table(promoterIslands[promoterIslands[,7]>0.01,6])


##############################################################################################
#Analysis of meningiomas from Heidelberg CNS tumor classifier (Capper et al. 2018, GSE109381)#
##############################################################################################
#Load normalized Heidelberg beta matrix, in case you'd like to process the data yourself, we used champ.load, and champ.norm, with default parameters (quality check and SVD did not show any errors warranting correction)
HeidelbergNorm <- as.matrix(read.csv("Heidelberg beta matrix.csv", row.names = 1))

#Heidelberg is 450K while we used EPIC/850K, so create a combined beta matrix of overlapping probes and then sort by decreasing variance
Heidelberg_Patel_Norm <- cbind(HeidelbergNorm[intersect(rownames(HeidelbergNorm), rownames(methData)),], methData[intersect(rownames(HeidelbergNorm), rownames(methData)), ])
Heidelberg_Patel_Norm <- Heidelberg_Patel_Norm[order(rowVars(Heidelberg_Patel_Norm), decreasing = TRUE),]

#######
#First, train a Random Forest model using our epigenetic types to then classify the Heidelberg samples 
trainingData <- methData[intersect(rownames(Heidelberg_Patel_Norm), MethSig),]
MethylationClassifier <- randomForest(as.factor(SampleData[colnames(methData), "MethCluster"]) ~., data=data.frame(t(trainingData)), method='class', ntree = 5000)

predictions_Heidelberg <- cbind(SampleData[colnames(methData), "MethCluster"], predict(MethylationClassifier, newdata=data.frame(t(trainingData)), type="response"), predict(MethylationClassifier, newdata=data.frame(t(trainingData)), type="prob"))
colnames(predictions_Heidelberg) <- c('NMF_type', 'RF_type', 'ProbType1', 'ProbType2', 'ProbType3')
sum(predictions_Heidelberg[,1] != predictions_Heidelberg[,2])
#sum should be 0 if RF model is assigning all training data to the correct types (expected outcome)
predictions_Heidelberg <- cbind(rep('Heidelberg', dim(Heidelberg_Patel_Norm)[2]), 
                                predict(MethylationClassifier, newdata=data.frame(t(Heidelberg_Patel_Norm[intersect(rownames(Heidelberg_Patel_Norm), MethSig),])), type="response"), 
                                predict(MethylationClassifier, newdata=data.frame(t(Heidelberg_Patel_Norm[intersect(rownames(Heidelberg_Patel_Norm), MethSig),])), type="prob"))
predictions_Heidelberg[rownames(SampleData), 1] <- SampleData$MethCluster
colnames(predictions_Heidelberg) <- c('NMF_type', 'RF_type', 'ProbType1', 'ProbType2', 'ProbType3')

########
#Second, NMF to identify clusters from Heidelberg samples using methylation signature probes (hence this analysis is agnostic to our data/samples, aside from the usage of our signature probes)
rank_estimate_Heidelberg <- nmf(HeidelbergNorm[intersect(rownames(HeidelbergNorm), MethSig),], 2:10, nrun=100,  .opt = 'v')
plot(rank_estimate_Heidelberg, main = 'Rank estimate for Heidelberg tumors over signature probes')
plot(2:10, rank_estimate_Heidelberg$measures$cophenetic, main = paste0('Heidelberg NMF Cophenetic coefficient'), type = 'b', xlab = '', ylab = '', col = 'purple', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
plot(2:10, rank_estimate_Heidelberg$measures$residuals, main = paste0('Heidelberg NMF Residuals'), type = 'b', xlab = '', ylab = '', col = 'dark blue', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
consensusmap(rank_estimate_Heidelberg, color = 'blue', main = NULL, annCol = NULL, annRow = NULL, tracks = NULL, info = FALSE, labRow = NA, labCol = NA, legend = NULL)

#3 appears to be optimal rank, now perform final NMF
Rank <- 3
NMF_Heidelberg <- nmf(HeidelbergNorm[intersect(rownames(HeidelbergNorm), MethSig),], rank = Rank, nrun = 1000, .opt = 'v')
consensusmap(NMF_Heidelberg)
consensus_Heidelberg <- predict(NMF_Heidelberg, what = 'consensus')
silh_Heidelberg <- silhouette(NMF_Heidelberg)
clusters_Heidelberg <- cbind(rep(0, length(consensus_Heidelberg)), consensus_Heidelberg, silh_Heidelberg)
colnames(clusters_Heidelberg)[1] <- 'Epigenetic_cluster'
clusters_Heidelberg[clusters_Heidelberg[,'consensus_Heidelberg'] == 1,'Epigenetic_cluster'] <- 3
clusters_Heidelberg[clusters_Heidelberg[,'consensus_Heidelberg'] == 2,'Epigenetic_cluster'] <- 2
clusters_Heidelberg[clusters_Heidelberg[,'consensus_Heidelberg'] == 3,'Epigenetic_cluster'] <- 1

###########
#Now compare assignments made by NMF and the random forest classifier trained on Patel samples 
Heidelberg_combined_classifications <- cbind(predictions_Heidelberg[intersect(rownames(clusters_Heidelberg), rownames(predictions_Heidelberg)), 'RF_type'], 
                                             clusters_Heidelberg[intersect(rownames(clusters_Heidelberg), rownames(predictions_Heidelberg)), 'Epigenetic_cluster'])

#Demonstrate assignments in a dendrogram, NMF is represented by the dendrogram, while leaf color are the random forest classifier assignments
dendrogram_Heidelberg <- consensushc(NMF_Heidelberg, method = 'mcquitty')
colors <- predictions_Heidelberg[labels(dendrogram_Heidelberg), 'RF_type']
colors[colors == 1] <- colors_by_group[1]
colors[colors == 2] <- colors_by_group[2]
colors[colors == 3] <- colors_by_group[3]
labels_colors(dendrogram_Heidelberg) <- colors
dendrogram_Heidelberg <- raise.dendrogram(dendrogram_Heidelberg, 0.02)
dendrogram_Heidelberg <- assign_values_to_leaves_nodePar(dendrogram_Heidelberg, 19, "pch")
dendrogram_Heidelberg <- assign_values_to_leaves_nodePar(dendrogram_Heidelberg, 1, "cex")
dendrogram_Heidelberg <- assign_values_to_leaves_nodePar(dendrogram_Heidelberg, colors, nodePar = "col")
plot(dendrogram_Heidelberg, leaflab = "none", main = "NMF dend w/ Random Forest assignments for Heidelberg", yaxt='n')


############
#tSNE of Patel and Heidelberg samples, colors being the random forest classifier assignments
colors <- cbind(predictions_Heidelberg[, 1:2])
colnames(colors) <- c('HeidelbergBlack', 'HeidelbergAssigned')
colors[colors[,1] != 'Heidelberg',1] <- colors[colors[,1] != 'Heidelberg',2]
colors[colors[,1] == 'Heidelberg',1] <- 'black'
colors[colors[,1] == '1',1] <- colors_by_group[1]
colors[colors[,1] == '2',1] <- colors_by_group[2]
colors[colors[,1] == '3',1] <- colors_by_group[3]
colors[colors[,2] == '1',2] <- colors_by_group[1]
colors[colors[,2] == '2',2] <- colors_by_group[2]
colors[colors[,2] == '3',2] <- colors_by_group[3]

#tSNE of top n most variably methylated probes
n <- 4500
tsne_output <- Rtsne(t(Heidelberg_Patel_Norm[1:n,]), pca = FALSE, perplexity = 8)
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,1], main = paste0('tSNE plot of Heidelberg methylation data: Top ', n, ' most variable probes'))
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,2], main = paste0('tSNE plot of Heidelberg methylation data: Top ', n, ' most variable probes'))

#tSNE of signature probes
tsne_output <- Rtsne(t(Heidelberg_Patel_Norm[intersect(rownames(Heidelberg_Patel_Norm), MethSig),]), pca = FALSE, perplexity = 8)
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,1], main = 'tSNE plot of Heidelberg methylation data, signature probes')
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,2], main = 'tSNE plot of Heidelberg methylation data, signature probes')


##################################
#Compare our data/clusters with probes from Olar et al. 2017
Olar_CPG <- load_Olar_CPG()
Olar_CPG$MMFav <- intersect(Olar_CPG$MMFav, rownames(methData))
Olar_CPG$MMUnFav <- intersect(Olar_CPG$MMUnFav, rownames(methData))
Olar_CPG$MMBoth <- intersect(Olar_CPG$MMBoth, rownames(methData))
Olar_CPG$MMFav64 <- intersect(Olar_CPG$MMFav64, rownames(methData))
Olar_CPG$MMUnFav64 <- intersect(Olar_CPG$MMUnFav64, rownames(methData))
Olar_CPG$MMBoth64 <- intersect(Olar_CPG$MMBoth64, rownames(methData))

#Compare intersection of Olar CpGs to all signature probes and each cluster's signature
length(Olar_CPG$MMBoth)
length(intersect(MethSig, Olar_CPG$MMBoth))
length(intersect(MethSig_cluster1, Olar_CPG$MMBoth))
length(intersect(MethSig_cluster2, Olar_CPG$MMBoth))
length(intersect(MethSig_cluster3, Olar_CPG$MMBoth))

Olar_AvgBeta <- matrix(nrow = length(Olar_CPG$MMBoth), ncol = 6, dimnames = list(Olar_CPG$MMBoth, c('Type1', 'Type2', 'Type3','NotType1', 'NotType2', 'NotType3')))
Olar_AvgBeta[,'Type1'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethCluster == 1))])
Olar_AvgBeta[,'Type2'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethCluster == 2))])
Olar_AvgBeta[,'Type3'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethCluster == 3))])
Olar_AvgBeta[,'NotType1'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethCluster %in% c(2,3)))])
Olar_AvgBeta[,'NotType2'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethCluster %in% c(1,3)))])
Olar_AvgBeta[,'NotType3'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethCluster %in% c(1,2)))])

#Compare average methylation between Meth types for favorable probes (positive favors first type, negative favors second type/group in each comparison)
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMFav,'NotType1']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type2'] - Olar_AvgBeta[Olar_CPG$MMFav,'NotType2']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type3'] - Olar_AvgBeta[Olar_CPG$MMFav,'NotType3']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMFav,'Type2']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMFav,'Type3']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type2'] - Olar_AvgBeta[Olar_CPG$MMFav,'Type3']))

#Compare average methylation between Meth types for unfavorable probes (negative favors first type, positive favors second type/group in each comparison)
table(sign(Olar_AvgBeta[Olar_CPG$MMUnFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMUnFav,'NotType1']))
table(sign(Olar_AvgBeta[Olar_CPG$MMUnFav,'Type2'] - Olar_AvgBeta[Olar_CPG$MMUnFav,'NotType2']))
table(sign(Olar_AvgBeta[Olar_CPG$MMUnFav,'Type3'] - Olar_AvgBeta[Olar_CPG$MMUnFav,'NotType3']))
table(sign(Olar_AvgBeta[Olar_CPG$MMUnFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMUnFav,'Type2']))
table(sign(Olar_AvgBeta[Olar_CPG$MMUnFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMUnFav,'Type3']))
table(sign(Olar_AvgBeta[Olar_CPG$MMUnFav,'Type2'] - Olar_AvgBeta[Olar_CPG$MMUnFav,'Type3']))

#Heatmap of MMFAv and MMUnFav probes
MMFavDend <- as.dendrogram(hclust(dist(methData[Olar_CPG$MMFav,])))
MMFavDend <- color_branches(MMFavDend, k = 1, col = 'orange')
MMUnFavDend <- as.dendrogram(hclust(dist(methData[Olar_CPG$MMUnFav,])))
MMUnFavDend <- color_branches(MMUnFavDend, k = 1, col = 'purple')
Olar_CPG_dendrogram <- merge(MMFavDend, MMUnFavDend)
column_ha <- HeatmapAnnotation(Epigenetic_cluster = SampleData[colnames(methData), 'MethCluster'],
                               col = list(Epigenetic_cluster = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242')))
Heatmap(methData[c(Olar_CPG$MMFav, Olar_CPG$MMUnFav),], cluster_rows = Olar_CPG_dendrogram, row_split = 2, column_split = SampleData[colnames(methData), 'MethCluster'], column_order = order.dendrogram(consensushc(NMFfinal, method = 'mcquitty')), top_annotation = column_ha, column_title = 'Heatmap of Olar probes by methylation types', row_title = 'Olar unfavorable and favorable probes', show_row_dend = TRUE, heatmap_legend_param = list(title = 'Beta'), show_row_names = FALSE, show_column_names = FALSE, use_raster = FALSE)


############
#PLS models#
############
#Load gene-wise CNV data obtained from WES data
CNV <- read.csv("Copy number status by gene.csv", row.names = 1)
colnames(CNV) <- gsub("X", "", gsub("\\.", "-", colnames(CNV)))

#Limit to samples with consistent classifications across each subclassification
PLSMenG_A <- intersect(rownames(subset(SampleData, RNA_numerical == 1 & MethCluster == 1 & Cytogenetic_type_numerical == 1 & NF2_instability_numerical == 1)), colnames(CNV))
PLSMenG_B <- intersect(rownames(subset(SampleData, RNA_numerical == 2 & MethCluster == 2 & Cytogenetic_type_numerical == 2 & NF2_instability_numerical == 2)), colnames(CNV))
PLSMenG_C <- intersect(rownames(subset(SampleData, RNA_numerical == 3 & MethCluster == 3 & Cytogenetic_type_numerical == 3 & NF2_instability_numerical == 3)), colnames(CNV))
PLSGroups <- c(PLSMenG_A, PLSMenG_B, PLSMenG_C)

#Get detailed chromosome location data
ensembl = useDataset('hsapiens_gene_ensembl', mart = useMart('ENSEMBL_MART_ENSEMBL'))
filterType <- 'external_gene_name'
attributeNames <- c('external_gene_name', 'band', 'chromosome_name')
annotations <- getBM(attributes=attributeNames, filters = filterType, values = rownames(expression), mart = ensembl)
annotations <- subset(annotations, band != '')
annotations$ChrLoc <- paste0(annotations$chromosome_name, annotations$band)

#Get all genes with expression, CNV, and promoter CpG data
AllGenes <- intersect(rownames(expression), intersect(rownames(CNV), CPG[CPG_promoter,]$gene))

#Initialze PLS model output matrix
AllGeneR2s <- data.frame(row.names = AllGenes, Chr = rep(0, length(AllGenes)), Components = 0, Rsquared = 0, CNV_importance = 0, CNV_LMcoef = 0, CNV_LMpval = 0, CNV_LMadjR2 = 0,
                         CpGs = 0, TopCPG = 0, TopCPGImp = 0, TopCPG_LMcoef = 0, TopCPG_LMpval = 0, TopCPG_LMadjR2 = 0, 
                         SecondCPG = 0, SecondCPGImp = 0, SecondCPG_LMcoef = 0, SecondCPG_LMpval = 0, SecondCPG_LMadjR2 = 0, 
                         ThirdCPG = 0, ThirdCPGImp = 0, ThirdCPG_LMcoef = 0, ThirdCPG_LMpval = 0, ThirdCPG_LMadjR2 = 0)

#Run PLS modeling of every gene with relevant data, from which subsets can be selected for further analysis
for(val in 1:length(AllGenes)) {
  Model <- getPLS(AllGenes[val], methylation = methData[,PLSGroups], CPGs = CPG[CPG_promoter,], expressionSubset = expression[,PLSGroups])
  if(length(Model)==1) { 
    AllGeneR2s[val, 'Components'] <- Model
    next() 
  }
  
  AllGeneR2s[val, 'Components'] <- Model[[2]]$bestTune$ncomp
  AllGeneR2s[val, 'Rsquared'] <- Model[[2]]$results$Rsquared[Model[[2]]$bestTune$ncomp]
  AllGeneR2s[val, 'CNV_importance'] <- Model[[3]]
  AllGeneR2s[val, 'CpGs'] <- Model[[4]]
  AllGeneR2s[val, 'TopCPG'] <- Model[[5]]
  AllGeneR2s[val, 'TopCPGImp'] <- Model[[6]]
  AllGeneR2s[val, 'SecondCPG'] <- Model[[7]]
  AllGeneR2s[val, 'SecondCPGImp'] <- Model[[8]]
  AllGeneR2s[val, 'ThirdCPG'] <- Model[[9]]
  AllGeneR2s[val, 'ThirdCPGImp'] <- Model[[10]]
  
  data <- as.data.frame(Model[[2]]$trainingData)
  linearModel <- correlateExpression(data[,c('expression', 'CNV')])
  AllGeneR2s[val, 'CNV_LMcoef'] <- linearModel[[2]]$coef
  AllGeneR2s[val, 'CNV_LMpval'] <- linearModel[[2]]$coefPval
  AllGeneR2s[val, 'CNV_LMadjR2'] <- linearModel[[2]]$adjRsq
  
  linearModel <- correlateExpression(data[,c('expression', Model[[5]])])
  AllGeneR2s[val, 'TopCPG_LMcoef'] <- linearModel[[2]]$coef
  AllGeneR2s[val, 'TopCPG_LMpval'] <- linearModel[[2]]$coefPval
  AllGeneR2s[val, 'TopCPG_LMadjR2'] <- linearModel[[2]]$adjRsq
  
  if(Model[[4]] > 1) {
    linearModel <- correlateExpression(data[,c('expression', Model[[7]])])
    AllGeneR2s[val, 'SecondCPG_LMcoef'] <- linearModel[[2]]$coef
    AllGeneR2s[val, 'SecondCPG_LMpval'] <- linearModel[[2]]$coefPval
    AllGeneR2s[val, 'SecondCPG_LMadjR2'] <- linearModel[[2]]$adjRsq
  }
  
  if(Model[[4]] > 2) {
    linearModel <- correlateExpression(data[,c('expression', Model[[9]])])
    AllGeneR2s[val, 'ThirdCPG_LMcoef'] <- linearModel[[2]]$coef
    AllGeneR2s[val, 'ThirdCPG_LMpval'] <- linearModel[[2]]$coefPval
    AllGeneR2s[val, 'ThirdCPG_LMadjR2'] <- linearModel[[2]]$adjRsq
  }
  
  message(val)
}

#Addend mean expression and mean methylation (beta) for the gene and highest ranked promoter site
for(val in 1:dim(AllGeneR2s)[1]) {
  AllGeneR2s$Chr[val] <- annotations[annotations$external_gene_name==rownames(AllGeneR2s)[val],]$ChrLoc
  if(rownames(AllGeneR2s)[val] %in% rownames(expression)) {
    AllGeneR2s$MenG_AAvg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSMenG_A, colnames(CNV))])
    AllGeneR2s$MenG_BAvg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSMenG_B, colnames(CNV))])
    AllGeneR2s$MenG_CAvg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSMenG_C, colnames(CNV))])
    AllGeneR2s$OverallAvg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSGroups, colnames(CNV))])
  } else {
    AllGeneR2s$MenG_AAvg[val] <- NA
    AllGeneR2s$MenG_BAvg[val] <- NA
    AllGeneR2s$MenG_CAvg[val] <- NA
    AllGeneR2s$OverallAvg[val] <- NA
  }
  if(AllGeneR2s$TopCPG[val] %in% rownames(methData)) {
    AllGeneR2s$MenG_ATopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSMenG_A, colnames(CNV))])
    AllGeneR2s$MenG_BTopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSMenG_B, colnames(CNV))])
    AllGeneR2s$MenG_CTopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSMenG_C, colnames(CNV))])
    AllGeneR2s$OverallTopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSGroups, colnames(CNV))])
  } else {
    AllGeneR2s$MenG_ATopCPGBeta[val] <- NA
    AllGeneR2s$MenG_BTopCPGBeta[val] <- NA
    AllGeneR2s$MenG_CTopCPGBeta[val] <- NA
    AllGeneR2s$OverallTopCPGBeta[val] <- NA
  }
  
  message(val)
}

#Write model outputs to a .csv file if desired
write.csv(AllGeneR2s, 'PLS models - All genes.csv')


##################
#If interested in a specific gene, can use the following to generate the individual model and then visualize the relationships between expression, CNV, and the three highest-ranked promoter sites
#Enter gene of interest below
gene <- 'HSPG2'

Model <- getPLS(gene, methylation = methData[,PLSGroups], CPGs = CPG[CPG_promoter,], expressionSubset = expression[,PLSGroups])
data <- as.data.frame(Model[[2]]$trainingData)
colors <- SampleData[rownames(data),]$MenG
colors[colors == '1'] <- colors_by_group[1]
colors[colors == '2'] <- colors_by_group[2]
colors[colors == '3'] <- colors_by_group[3]

#Plot log expression versus CNV
plot(data[,c(2, 1)], col = colors, pch = 19, ylab = 'Log expession', xlab = 'CNV')

#Plot log expression versus top-ranked promoter site
plot(data[,c(Model[[5]], 'expression')], col = colors, pch = 19, ylab = 'Log expession', xlab = 'Top CPG Beta')

#Plot log expression versus second-ranked promoter site
plot(data[,c(Model[[7]], 'expression')], col = colors, pch = 19, ylab = 'Log expession', xlab = 'Second CPG Beta')

#Plot log expression versus third-ranked promoter site
plot(data[,c(Model[[9]], 'expression')], col = colors, pch = 19, ylab = 'Log expession', xlab = 'Third CPG Beta')
###################

####
#Generate matrices and .csv files for analysis of PLS results for gene subset
#Genes on chromosomes 1p and 22q
Ch1pGenes <- unique(subset(annotations[grep('1p', annotations$ChrLoc),], chromosome_name == 1)$external_gene_name)
Ch22qGenes <- unique(subset(annotations, chromosome_name == 22)$external_gene_name)
write.csv(AllGeneR2s[intersect(c(Ch1pGenes, Ch22qGenes), AllGenes),], 'Chromosome 1p and 22q genes.csv')

#Genes used by the transcriptional random forest classifier in Patel et al. 2019
transcriptional_classifier_genes <- read.csv('Patel transcriptional classifier list.csv')
write.csv(AllGeneR2s[intersect(transcriptional_classifier_genes$Gene, AllGenes),], 'Transcriptional classifier genes.csv')


#Identify differentially expressed genes (DEGs) among the PLS cohort 
#Load raw counts data for the PLS subset and prepare for analysis with DESeq
PLScounts <- read.csv('Raw counts expression data.csv', row.names = 1)
colnames(PLScounts) <- gsub('X', '', gsub('(\\.)', '-', colnames(PLScounts)))
rownames(PLScounts) <-  gsub('(\\.).*', '', rownames(PLScounts))
PLScounts <- PLScounts[,c(PLSGroups, 'Arachnoid_1', 'Arachnoid_2', 'Arachnoid_3')]

#Create column data matrix of each MenG and arachnoid samples individually and get expression data
colData <- data.frame(row.names = colnames(PLScounts), type = rep(0, dim(PLScounts)[2]))
colData[rownames(colData) %in% rownames(SampleData),] <- SampleData[rownames(colData)[rownames(colData) %in% rownames(SampleData)], 'MenG']
colData[grep('Arachnoid', rownames(colData)),] <- 'Arachnoid'
colData$type <- factor(colData$type)
PLSdds <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~type)
PLSdds <- DESeq(PLSdds[rowSums(counts(PLSdds)) >= 10, ])
PLSexpression <- getExpressionData(PLSdds)

#Create a column data matrix the will perform differential expression analysis between each MenG/arachnoid and all other sample combined (i.e. MenG A vs all others)
colData <- data.frame(row.names = colnames(PLScounts), groups = factor(c(SampleData[PLSGroups, 'MenG'], rep('Arachnoid',3))))
levels(colData$groups) <- c(levels(colData$groups), 'other')
colData$arachnoid <- colData$groups
colData$arachnoid[colData$arachnoid %in% c(1,2,3)] <- 'other'

ddsGroups <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~groups)
ddsArachnoid <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~arachnoid)

#Process counts data with DEseq, removing anything with low total count numbers the change rownames from ENSG to gene names, removing duplicates
groupDDS <- DESeq(ddsGroups[rowSums(counts(ddsGroups)) >= 10])
arachnoidDDS <- DESeq(ddsArachnoid[rowSums(counts(ddsArachnoid)) >= 10])

ENSGlist <- rownames(groupDDS)
geneList <- getGeneList(ensemblList = ENSGlist)
rownames(groupDDS) <- geneList[match(ENSGlist, geneList[,"ensembl_gene_id"]),"external_gene_name"]
rownames(arachnoidDDS) <- geneList[match(ENSGlist, geneList[,"ensembl_gene_id"]),"external_gene_name"]

for(val in 1:length(ENSGlist)) {
  if(is.na(rownames(groupDDS)[val])) {
    rownames(groupDDS)[val] <- ENSGlist[val]
    rownames(arachnoidDDS)[val] <- ENSGlist[val]
  }
}

groupDDS <- groupDDS[rownames(groupDDS) != "",]
arachnoidDDS <- arachnoidDDS[rownames(arachnoidDDS) != "",]

groupDDS <- groupDDS[!duplicated(rownames(groupDDS)),]
arachnoidDDS <- arachnoidDDS[!duplicated(rownames(arachnoidDDS)),]

#Compare arachnoid to all tumors and each group individually
res.Tumors.Arachnoid <- results(arachnoidDDS, contrast=c('arachnoid','other','Arachnoid'))
res.MenG_A.Arachnoid <- results(groupDDS, contrast=c('groups','1','Arachnoid'))
res.MenG_B.Arachnoid <- results(groupDDS, contrast=c('groups','2','Arachnoid'))
res.MenG_C.Arachnoid <- results(groupDDS, contrast=c('groups','3','Arachnoid'))

#select
pval <- 1E-6
Wald <- 0
DEGs.Tumors.Arachnoid <- subset(res.Tumors.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)
DEGs.MenG_A.Arachnoid <- subset(res.MenG_A.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)
DEGs.MenG_B.Arachnoid <- subset(res.MenG_B.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)
DEGs.MenG_C.Arachnoid <- subset(res.MenG_C.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)

DEGs.Arachnoid <- rownames(DEGs.Tumors.Arachnoid)
DEGs.MenG_A <- rownames(DEGs.MenG_A.Arachnoid)[!rownames(DEGs.MenG_A.Arachnoid) %in% DEGs.Arachnoid]
DEGs.MenG_B <- rownames(DEGs.MenG_B.Arachnoid)[!rownames(DEGs.MenG_B.Arachnoid) %in% DEGs.Arachnoid]
DEGs.MenG_C <- rownames(DEGs.MenG_C.Arachnoid)[!rownames(DEGs.MenG_C.Arachnoid) %in% DEGs.Arachnoid]

AllDEGs <- unique(c(DEGs.MenG_A, DEGs.MenG_B, DEGs.MenG_C, DEGs.Arachnoid))

uniqueMenG_A <- DEGs.MenG_A[!DEGs.MenG_A %in% c(DEGs.MenG_B, DEGs.MenG_C)]
uniqueMenG_B <- DEGs.MenG_B[!DEGs.MenG_B %in% c(DEGs.MenG_A, DEGs.MenG_C)]
uniqueMenG_C <- DEGs.MenG_C[!DEGs.MenG_C %in% c(DEGs.MenG_A, DEGs.MenG_B)]

#Write PLS outputs for each set of DEGs
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), AllDEGs),], 'All DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueMenG_A),], 'MenG A unique DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueMenG_B),], 'Meng B unique DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueMenG_C),], 'MenG C unique DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), DEGs.Arachnoid),], 'Arachnoid DEGs.csv')

#Identify the proportion DEGs on chromosome 1p within the subsets of DEGs unique to each group
length(grep('1p', subset(annotations, external_gene_name %in% uniqueMenG_A & chromosome_name == 1)$ChrLoc))/length(uniqueMenG_A)
length(grep('1p', subset(annotations, external_gene_name %in% uniqueMenG_B & chromosome_name == 1)$ChrLoc))/length(uniqueMenG_B)
length(grep('1p', subset(annotations, external_gene_name %in% uniqueMenG_C & chromosome_name == 1)$ChrLoc))/length(uniqueMenG_C)


#Identify genes with differentially methylatied promoters (DMPs) among the PLS cohort 
#Create a column data matrix the will perform differential methylation analysis between each MenG and all other sample combined (i.e. MenG A vs all others)
DMPpheno <- data.frame(row.names = PLSGroups, MenG = factor(c(SampleData[PLSGroups, 'MenG'])))
levels(DMPpheno$MenG) <- c(levels(DMPpheno$MenG), 'other')
DMPpheno$MenG_A <- DMPpheno$MenG
DMPpheno$MenG_A[DMPpheno$MenG_A %in% c(2,3)] <- 'other'
DMPpheno$MenG_B <- DMPpheno$MenG
DMPpheno$MenG_B[DMPpheno$MenG_B %in% c(1,3)] <- 'other'
DMPpheno$MenG_C <- DMPpheno$MenG
DMPpheno$MenG_C[DMPpheno$MenG_C %in% c(1,2)] <- 'other'

#Use champ.DMP to find statistically significant DMPs (pval < 1e-06, min difference in mean beta 0.1)
OverallMenG_A <- champ.DMP(beta = methData[,PLSGroups], pheno=as.character(DMPpheno$MenG_A), arraytype = 'EPIC')
OverallMenG_B <- champ.DMP(beta = methData[,PLSGroups], pheno=as.character(DMPpheno$MenG_B), arraytype = 'EPIC')
OverallMenG_C <- champ.DMP(beta = methData[,PLSGroups], pheno=as.character(DMPpheno$MenG_C), arraytype = 'EPIC')

pval <- 1E-6
difBeta <- 0.1
OverallMenG_ADMPs <- intersect(CPG_promoter, rownames(subset(OverallMenG_A[[1]], adj.P.Val < pval & abs(deltaBeta) > difBeta)))
OverallMenG_BDMPs <- intersect(CPG_promoter, rownames(subset(OverallMenG_B[[1]], adj.P.Val < pval & abs(deltaBeta) > difBeta)))
OverallMenG_CDMPs <- intersect(CPG_promoter, rownames(subset(OverallMenG_C[[1]], adj.P.Val < pval & abs(deltaBeta) > difBeta)))

MenG_ADMPgenes <- unique(CPG[OverallMenG_ADMPs,]$gene)
MenG_BDMPgenes <- unique(CPG[OverallMenG_BDMPs,]$gene)
MenG_CDMPgenes <- unique(CPG[OverallMenG_CDMPs,]$gene)

DMPgenes <- unique(c(MenG_ADMPgenes, MenG_BDMPgenes, MenG_CDMPgenes))
uniqueDMPgenesMenG_A <- MenG_ADMPgenes[!MenG_ADMPgenes %in% c(MenG_BDMPgenes, MenG_CDMPgenes)]
uniqueDMPgenesMenG_B <- MenG_BDMPgenes[!MenG_BDMPgenes %in% c(MenG_ADMPgenes, MenG_CDMPgenes)]
uniqueDMPgenesMenG_C <- MenG_CDMPgenes[!MenG_CDMPgenes %in% c(MenG_ADMPgenes, MenG_BDMPgenes)]

write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), DMPgenes),], 'All DMP genes.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueDMPgenesMenG_A),], 'Group 1 DMP genes.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueDMPgenesMenG_B),], 'Group 2 DMP genes.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueDMPgenesMenG_C),], 'Group 3 DMP genes.csv')


####################################################################
#Analysis of data from Magill et al. 2020 (GSE151067 and GSE151921)#
####################################################################
#Patel and Magill MenG and subclassifications
Patel_Magill_classifications <- read.csv("Patel and Magill classifications.csv", row.names = 1)

#Load normalized beta matrix for both Patel and Magill samples, in case you'd like to process the data yourself, we used champ.load, and champ.norm, with default parameters (quality check and SVD did not show any errors warranting correction)
Patel_Magill_norm <- read.csv("Patel and Magill combined beta matrix.csv", row.names = 1)
colnames(Patel_Magill_norm) <- gsub("X", "", gsub("(\\.)", "-", colnames(Patel_Magill_norm)))

#Load raw counts data for both Patel and Magill samples, then use DESeq to obtain log transformed expression data
Patel_Magill_counts <- read.csv("Patel and Magill combined raw counts.csv", row.names = 1)
colnames(Patel_Magill_counts) <- gsub("X", "", gsub("(\\.)", "-", colnames(Patel_Magill_counts)))
ddsPatel_Magill <- DESeqDataSetFromMatrix(countData = Patel_Magill_counts, colData = Patel_Magill_classifications[colnames(Patel_Magill_counts),], design = ~RNA)
ddsPatel_Magill <- DESeq(ddsPatel_Magill[rowSums(counts(ddsPatel_Magill)) >= 10, ])
Patel_Magill_expression <- assay(vst(ddsPatel_Magill, blind=TRUE))

####Identify transcriptional and methylation types of Magill samples
#Classify RNA types by training a random forest classifier over Patel samples and then applying to Magill samples using the genes from our previous transcriptional classifier
RNA_classifier_genes <- intersect(read.csv('Patel transcriptional classifier list.csv')$Gene, rownames(Patel_Magill_expression))
trainingData <- Patel_Magill_expression[RNA_classifier_genes, rownames(subset(Patel_Magill_classifications, UCSFnumber == 'Patel' & MenG != 'Unknown'))]
modelRNA <- randomForest(as.factor(Patel_Magill_classifications[colnames(trainingData),]$RNA) ~., data=data.frame(t(trainingData)), method='class', ntree = 5000)
predictionsRNA <- cbind(Patel_Magill_classifications[colnames(Patel_Magill_expression),]$RNA, predict(modelRNA, newdata=data.frame(t(Patel_Magill_expression[RNA_classifier_genes,])), type="response"), predict(modelRNA, newdata=data.frame(t(Patel_Magill_expression[RNA_classifier_genes,])), type="prob"))

#Classify methylation types of Magill samples, similarly to above using a random forest model trained on Patel samples and using our signature probes
Methylation_classifier_probes <- intersect(MethSig, rownames(Patel_Magill_norm))
trainingData <- Patel_Magill_norm[Methylation_classifier_probes, rownames(subset(Patel_Magill_classifications, UCSFnumber == 'Patel' & MenG != 'Unknown'))]
modelMeth <- randomForest(as.factor(Patel_Magill_classifications[colnames(trainingData),]$Meth) ~., data=data.frame(t(trainingData)), method='class', ntree = 5000)
predictionsMeth <- cbind(Patel_Magill_classifications[colnames(Patel_Magill_norm),]$Meth, predict(modelMeth, newdata=data.frame(t(Patel_Magill_norm[Methylation_classifier_probes,])), type="response"), predict(modelMeth, newdata=data.frame(t(Patel_Magill_norm[Methylation_classifier_probes,])), type="prob"))
