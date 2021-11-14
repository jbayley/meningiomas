#Setup and data importing/loading
setwd('C:/Users/james/Desktop/Patel_Methylation/All_Patel_Data/')
source('Functions_Bayley_Final.R')
loadPackages()

#Import methylation data
myLoad <- champ.load(directory = paste0(getwd(), '/Methylation'), arraytype = 'EPIC')

#Can run quality check, but no issues or any bad samples
#champ.QC()

#Normalize methylation data (also sorts probes in decreasing order of overall variance)
myNorm <- normalizeMethylation()

#Can check for batch effect with SVD, but there is none and so no need for COMBAT
#champ.SVD()

#Load information for all CpGs, previously obtained via ChAMP
CPG <- read.csv('CPG information.csv', row.names = 1)

#Load data with SeSAMe pipeline if desired
#sesameBetas <- loadSesame()

#Load clinical, genomic, and data for this cohort
SampleData <- importSampleData()

#Load RNA-seq data and generate transformed expression data
counts <- read.csv('Raw counts expression data.csv', row.names = 1)
colnames(counts) <- gsub('X', '', gsub('(\\.)', '-', colnames(counts)))
rownames(counts) <-  gsub('(\\.).*', '', rownames(counts))
colData <- data.frame(row.names = colnames(counts), type = rep(0, dim(counts)[2]))
colData[rownames(SampleData),] <- SampleData$CompositeGroup
colData[grep('Arachnoid', rownames(colData)),] <- 'Arachnoid'
colData$type <- factor(colData$type)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~type)
dds <- DESeq(dds[rowSums(counts(dds)) >= 10, ])
expression <- getExpressionData(dds)

#Colors for figures
colors_by_group <- c('#6FDE6E', '#235FA4', '#FF4242', '#BDD9BF', '#2E4052', '#E5323B')
names(colors_by_group) <- c('GroupA', 'GroupB', 'GroupC', 'GroupA2', 'GroupB2', 'GroupC2')

###############################################
#Establishment of meningioma epigenetic types#
###############################################

#Select the ChAMP or Sesame dataset
methData <- myNorm
#methData <- sesameBetas

#Minimum condition calculation to determine number of CpGs
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

#Perform NMF with this set of DMPs to determine optimal rank
rank_estimate_initial <- nmf(methData[1:nDMPs,], 2:10, nrun=100,  .opt = 'v')
plot(rank_estimate_initial, main = paste0('Rank estimate for top ', nDMPs, 'most variably methylated probes'))
plot(2:10, rank_estimate_initial$measures$cophenetic, main = paste0('Top ', nDMPs, ' Cophenetic coefficient'), type = 'b', xlab = '', ylab = '', col = 'purple', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
plot(2:10, rank_estimate_initial$measures$residuals, main = paste0('Top ', nDMPs, ' Residuals'), type = 'b', xlab = '', ylab = '', col = 'dark blue', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)

consensusmap(rank_estimate_initial, color = 'blue', main = NULL, annCol = NULL, annRow = NULL, tracks = NULL, info = FALSE, labRow = NA, labCol = NA, legend = NULL)

#Recheck optimal number for groups for top 10% of DMPs
top10DMPs <- floor(dim(methData)[1]/10)
#Running the following attempts to allocate a 40Gb vector when calculating RSS, so run each rank individually
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
Kmeans_top10 <- ConsensusClusterPlus(methData[1:top10DMPs,], maxK=10, reps=500, pItem=0.8, pFeature=1, title='K-means clustering top 10 DMPs', clusterAlg='km', distance='euclidean', plot='png')

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

#Run NMF for top DMPs and export results (for Sesame run, switch the 3 and 2 for the cluster number assignment - order 3,2,1)
NMFinitial <- nmf(methData[1:nDMPs,], rank = Rank, nrun = 1000, .opt = 'v')
consensus_initial <- predict(NMFinitial, what = 'consensus')
silh_initial <- silhouette(NMFinitial)
clusters_initial <- cbind(rep(0, length(consensus_initial)), consensus_initial, silh_initial)
colnames(clusters_initial)[1] <- 'MethGroup'
clusters_initial[clusters_initial[,'consensus_initial'] == 1,'MethGroup'] <- 2
clusters_initial[clusters_initial[,'consensus_initial'] == 2,'MethGroup'] <- 3
clusters_initial[clusters_initial[,'consensus_initial'] == 3,'MethGroup'] <- 1
write.csv(clusters_initial, 'Initial NMF Predictions.csv')

#Run NMF for top 10% DMPs and export results
NMFtop10 <- nmf(methData[1:top10DMPs, ], rank = Rank, nrun = 1000, .opt = 'v')
consensus_top10 <- predict(NMFtop10, what = 'consensus')
silh_top10 <- silhouette(NMFtop10)
clusters_top10 <- cbind(rep(0, length(consensus_top10)), consensus_top10, silh_top10)
colnames(clusters_top10)[1] <- 'MethGroup'
clusters_top10[clusters_top10[,'consensus_top10'] == 1,'MethGroup'] <- 3
clusters_top10[clusters_top10[,'consensus_top10'] == 2,'MethGroup'] <- 2
clusters_top10[clusters_top10[,'consensus_top10'] == 3,'MethGroup'] <- 1
write.csv(clusters_top10, 'Top 10 NMF Predictions.csv')

#Get assignments from these 4 analyses and find samples with discordant clusters
clusters <- data.frame(row.names = rownames(clusters_initial), NMF_initial = clusters_initial[,'MethGroup'], NMF_top10 = clusters_top10[,'MethGroup'], Kmeans_initial = Kmeans_initial[[3]]$consensusClass, Kmeans_top10 = Kmeans_top10[[3]]$consensusClass)
discordant <- vector()
for(val in 1:dim(clusters)[1]) {
  if(var(as.numeric(clusters[val,])) != 0) {
    discordant <- c(discordant, rownames(clusters)[val])
  }
}

#Select those samples concordant between all 3 initial classifications
distilledSamples <- subset(colnames(methData), !colnames(methData) %in% discordant)
distilledSamples1 <- intersect(distilledSamples, rownames(clusters_initial[clusters_initial[,'MethGroup']==1,]))
distilledSamples2 <- intersect(distilledSamples, rownames(clusters_initial[clusters_initial[,'MethGroup']==2,]))
distilledSamples3 <- intersect(distilledSamples, rownames(clusters_initial[clusters_initial[,'MethGroup']==3,]))

#Use champ.DMP to define DMPs, comparing each group against all other samples
distilledGroups <- data.frame(row.names = distilledSamples, pheno1 = clusters[distilledSamples, 1], pheno2 = clusters[distilledSamples, 1], pheno3 = clusters[distilledSamples, 1])
distilledGroups$pheno1[distilledGroups$pheno1!=1] <- 'other'
distilledGroups$pheno2[distilledGroups$pheno2!=2] <- 'other'
distilledGroups$pheno3[distilledGroups$pheno3!=3] <- 'other'

#Identify DMPs for a given adjusted p-value
pval <- 1e-06
DMPgroup1 <- champ.DMP(beta = methData[,distilledSamples], pheno=distilledGroups$pheno1, adjPVal = pval, arraytype = 'EPIC')
DMPgroup2 <- champ.DMP(beta = methData[,distilledSamples], pheno=distilledGroups$pheno2, adjPVal = pval, arraytype = 'EPIC')
DMPgroup3 <- champ.DMP(beta = methData[,distilledSamples], pheno=distilledGroups$pheno3, adjPVal = pval, arraytype = 'EPIC')

#Chose those DMPs with a difference in mean beta between a group and the other samples > a cutoff (we selected 0.25)
cutoff <- 0.25
CPGgroup1 <- rownames(DMPgroup1[[1]][abs(DMPgroup1[[1]]$deltaBeta) > cutoff,])
CPGgroup2 <- rownames(DMPgroup2[[1]][abs(DMPgroup2[[1]]$deltaBeta) > cutoff,])
CPGgroup3 <- rownames(DMPgroup3[[1]][abs(DMPgroup3[[1]]$deltaBeta) > cutoff,])

#Remove DMPs found for multiple groups
overlap1v2 <- intersect(CPGgroup1, CPGgroup2)
overlap1v3 <- intersect(CPGgroup1, CPGgroup3)
overlap2v3 <- intersect(CPGgroup2, CPGgroup3)
CPGgroup1 <- CPGgroup1[!CPGgroup1 %in% c(overlap1v2, overlap1v3)]
CPGgroup2 <- CPGgroup2[!CPGgroup2 %in% c(overlap1v2, overlap2v3)]
CPGgroup3 <- CPGgroup3[!CPGgroup3 %in% c(overlap1v3, overlap2v3)]

#Now that we have the methylation signatures, perform final NMF using these probes
#Redo rank estimation, although this is largely teleological at this point
rank_estimate_final <- nmf(methData[c(CPGgroup1, CPGgroup2, CPGgroup3),], 2:10, nrun=100,  .opt = 'v')
plot(rank_estimate_final, main = 'Rank estimate for signature probes')
consensusmap(rank_estimate_final, color = 'blue', main = NULL, annCol = NULL, annRow = NULL, tracks = NULL, info = FALSE, labRow = NA, labCol = NA, legend = NULL)

#Indeed 3 is the optimal rank
NMFfinal <- nmf(methData[c(CPGgroup1, CPGgroup2, CPGgroup3),], rank = Rank, nrun = 1000, .opt = 'v')
consensusmap(NMFfinal)
consensus_final <- predict(NMFfinal, what = 'consensus')
silh_final <- silhouette(NMFfinal)
clusters_final <- cbind(rep(0, length(consensus_final)), consensus_final, silh_final)
colnames(clusters_final)[1] <- 'MethGroup'
clusters_final[clusters_final[,'consensus_final'] == 1,'MethGroup'] <- 3
clusters_final[clusters_final[,'consensus_final'] == 2,'MethGroup'] <- 2
clusters_final[clusters_final[,'consensus_final'] == 3,'MethGroup'] <- 1
write.csv(clusters_final, 'Final NMF Predictions.csv')

#Heatmap of signature probes
Group1Dend <- as.dendrogram(hclust(dist(methData[CPGgroup1,])))
Group1Dend <- color_branches(Group1Dend, k = 1, col = colors_by_group[4])
Group2Dend <- as.dendrogram(hclust(dist(methData[CPGgroup2,])))
Group2Dend <- color_branches(Group2Dend, k = 1, col = colors_by_group[5])
Group3Dend <- as.dendrogram(hclust(dist(methData[CPGgroup3,])))
Group3Dend <- color_branches(Group3Dend, k = 1, col = colors_by_group[6])
CPG_dendrogram <- merge(Group1Dend, Group2Dend, Group3Dend)
plot(CPG_dendrogram)

column_ha <- HeatmapAnnotation(MethType = clusters_final[,'MethGroup'],
                               col = list(MethType = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242')))

Heatmap(methData[c(CPGgroup1, CPGgroup2, CPGgroup3),], cluster_rows = CPG_dendrogram, row_split = 2, column_split = clusters_final[,'MethGroup'], column_order = order.dendrogram(consensushc(NMFfinal, method = 'mcquitty')), show_row_dend = TRUE, top_annotation = column_ha, show_row_names = FALSE, show_column_names = FALSE, column_title = 'Methylation Types', row_title = 'Methylation Signatures', heatmap_legend_param = list(title = 'Beta'), use_raster = FALSE)
Heatmap(methData[c(CPGgroup1, CPGgroup2, CPGgroup3),], cluster_rows = CPG_dendrogram, row_split = 2, column_split = 3, cluster_columns = consensushc(NMFfinal, method = 'mcquitty'), show_row_dend = TRUE, top_annotation = column_ha, show_row_names = FALSE, show_column_names = FALSE, column_title = 'Methylation Types', row_title = 'Methylation Signatures', heatmap_legend_param = list(title = 'Beta'), use_raster = FALSE)


#r-TSNE of different sets of probes
tsneSet <- methData[1:nDMPs,]
title <- paste0('t-SNE of ', nDMPs, ' most variable probes')
tsneSet <- methData[1:top10DMPs,]
title <- paste0('t-SNE of ', top10DMPs, ' most variable probes')
tsneSet <- methData[c(CPGgroup1, CPGgroup2, CPGgroup3),]
title <- paste0('t-SNE of methylation signature probes (n = ', dim(tsneSet)[1], ')')

colors <- clusters_final[,'MethGroup']
colors[colors == '1'] <- colors_by_group[1]
colors[colors == '2'] <- colors_by_group[2]
colors[colors == '3'] <- colors_by_group[3]
tsne_output <- Rtsne(t(tsneSet), pca = FALSE, perplexity = 8, labels = colnames(methData))
plot(tsne_output$Y, col=colors, pch = 19, asp=0, cex = 2, xlab = '', ylab = '', main = title)
text(tsne_output$Y, labels=colnames(myNorm), col=colors)

#Signature and CpG promoter island methylation levels
signatureMeth <- matrix(data = NA, nrow = length(c(CPGgroup1, CPGgroup2, CPGgroup3)), ncol = 6, dimnames = list(c(CPGgroup1, CPGgroup2, CPGgroup3), c('Type1Mean', 'Type2Mean', 'Type3Mean', 'Type1', 'Type2', 'Type3')))
signatureMeth[,1] <- rowMeans(methData[rownames(signatureMeth), rownames(subset(SampleData, MethGroup == 1))])
signatureMeth[,2] <- rowMeans(methData[rownames(signatureMeth), rownames(subset(SampleData, MethGroup == 2))])
signatureMeth[,3] <- rowMeans(methData[rownames(signatureMeth), rownames(subset(SampleData, MethGroup == 3))])

for(val in 1:dim(signatureMeth)[1]) {
  minDelta <- min(abs(signatureMeth[val,1]-signatureMeth[val,2]), abs(signatureMeth[val,1]-signatureMeth[val,3]), abs(signatureMeth[val,2]-signatureMeth[val,3]))
  
  if(minDelta == abs(signatureMeth[val,1]-signatureMeth[val,2])) {
    signatureMeth[val,4] <- 0
    signatureMeth[val,5] <- 0
    if(signatureMeth[val,3] > signatureMeth[val,1] & signatureMeth[val,3] > signatureMeth[val,2]) {
      signatureMeth[val,6] <- 1
    } else if(signatureMeth[val,3] < signatureMeth[val,1] & signatureMeth[val,3] < signatureMeth[val,2]) {
      signatureMeth[val,6] <- -1
    } else {
      signatureMeth[val,6] <- 'Problem'
    }
  } else if(minDelta == abs(signatureMeth[val,1]-signatureMeth[val,3])) {
    signatureMeth[val,4] <- 0
    signatureMeth[val,6] <- 0
    if(signatureMeth[val,2] > signatureMeth[val,1] & signatureMeth[val,2] > signatureMeth[val,3]) {
      signatureMeth[val,5] <- 1
    } else if(signatureMeth[val,2] < signatureMeth[val,1] & signatureMeth[val,2] < signatureMeth[val,3]) {
      signatureMeth[val,5] <- -1
    } else {
      signatureMeth[val,5] <- 'Problem'
    }
  } else if(minDelta == abs(signatureMeth[val,2]-signatureMeth[val,3])) {
    signatureMeth[val,5] <- 0
    signatureMeth[val,6] <- 0
    if(signatureMeth[val,1] > signatureMeth[val,2] & signatureMeth[val,1] > signatureMeth[val,3]) {
      signatureMeth[val,4] <- 1
    } else if(signatureMeth[val,1] < signatureMeth[val,2] & signatureMeth[val,1] < signatureMeth[val,3]) {
      signatureMeth[val,4] <- -1
    } else {
      signatureMeth[val,4] <- 'Problem'
    }
  } else {
    message(val)
  }
}

table(signatureMeth[,4])
table(signatureMeth[,5])
table(signatureMeth[,6])

CPG_promoter <- rownames(subset(CPG[rownames(methData),], feature %in% c('1stExon', 'TSS200', 'TSS1500')))
CPG_promoter_island <- rownames(subset(CPG[rownames(methData),], feature %in% c('1stExon', 'TSS200', 'TSS1500') & cgi == 'island'))

table(signatureMeth[!rownames(signatureMeth) %in% CPG_promoter,4])
table(signatureMeth[!rownames(signatureMeth) %in% CPG_promoter,5])
table(signatureMeth[!rownames(signatureMeth) %in% CPG_promoter,6])

table(signatureMeth[intersect(rownames(signatureMeth), CPG_promoter),4])
table(signatureMeth[intersect(rownames(signatureMeth), CPG_promoter),5])
table(signatureMeth[intersect(rownames(signatureMeth), CPG_promoter),6])

table(signatureMeth[intersect(rownames(signatureMeth), CPG_promoter_island),4])
table(signatureMeth[intersect(rownames(signatureMeth), CPG_promoter_island),5])
table(signatureMeth[intersect(rownames(signatureMeth), CPG_promoter_island),6])

promoterIslands <- matrix(data = NA, nrow = length(CPG_promoter_island), ncol = 7, dimnames = list(CPG_promoter_island, c('Type1Mean', 'Type2Mean', 'Type3Mean', 'Type1', 'Type2', 'Type3', 'MaxDelta')))
promoterIslands[,1] <- rowMeans(methData[rownames(promoterIslands), rownames(subset(SampleData, MethGroup == 1))])
promoterIslands[,2] <- rowMeans(methData[rownames(promoterIslands), rownames(subset(SampleData, MethGroup == 2))])
promoterIslands[,3] <- rowMeans(methData[rownames(promoterIslands), rownames(subset(SampleData, MethGroup == 3))])

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

table(promoterIslands[,4])
table(promoterIslands[,5])
table(promoterIslands[,6])

table(promoterIslands[promoterIslands[,7]>0.01,4])
table(promoterIslands[promoterIslands[,7]>0.01,5])
table(promoterIslands[promoterIslands[,7]>0.01,6])



#Analysis of Heidelberg classifier meningiomas
setwd('C:/Users/james/Desktop/Patel_Methylation/Heidelberg Data')
HeidelbergLoad <- champ.load(directory = paste0(getwd(), '/Data'), arraytype = "450K")
#All samples are of good quality on QC
champ.QC(beta = HeidelbergLoad$beta, pheno = HeidelbergLoad$pd$Sample_Group)
HeidelbergNorm <- champ.norm(beta = HeidelbergLoad$beta, rgSet = HeidelbergLoad$rgSet, mset = HeidelbergLoad$mset, resultsDir="./CHAMP_Normalization/", cores = 1, arraytype = "450K") 
#SVD did not show any batch effects warranting correction
champ.SVD(beta = HeidelbergNorm, pd = Heidelbergload$pd)

#Sort normalized Heidelberg beta matrix by decreasing overall variance
HeidelbergNorm <- HeidelbergNorm[order(rowVars(HeidelbergNorm), decreasing = TRUE),]

#Heidelberg is 450K while we used EPIC/850K, so create a combined beta matrix of overlapping probes
Heidelberg_Patel_Norm <- cbind(HeidelbergNorm[intersect(rownames(HeidelbergNorm), rownames(methData)),], methData[intersect(rownames(HeidelbergNorm), rownames(methData)), ])
Heidelberg_Patel_Norm <- Heidelberg_Patel_Norm[order(rowVars(Heidelberg_Patel_Norm), decreasing = TRUE),]

#Perform tSNE for the n most variable probes and across our epigentic signature probes (those in common with the 450K array, as our was derived on EPIC 850K)
n <- 4500
tsne_output <- Rtsne(t(HeidelbergNorm[1:n,]), pca = FALSE, perplexity = 8)
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", main = paste0('tSNE plot of Heidelberg methylation data: Top ', n, ' most variable probes'))

tsne_output <- Rtsne(t(HeidelbergNorm[intersect(rownames(HeidelbergNorm), c(CPGgroup1, CPGgroup2, CPGgroup3)),]), pca = FALSE, perplexity = 8)
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", main = 'tSNE plot of Heidelberg methylation data')

#Train a Random Forest model using our epigenetic types to classify the Heidelberg samples 
trainingData <- methData[intersect(rownames(HeidelbergNorm), c(CPGgroup1, CPGgroup2, CPGgroup3)),]
MethylationClassifier <- randomForest(as.factor(SampleData[colnames(methData), "MethGroup"]) ~., data=data.frame(t(trainingData)), method='class', ntree = 5000)

predictions_Heidelberg <- cbind(SampleData[colnames(methData), "MethGroup"], predict(MethylationClassifier, newdata=data.frame(t(trainingData)), type="response"), predict(MethylationClassifier, newdata=data.frame(t(trainingData)), type="prob"))
colnames(predictions_Heidelberg) <- c('NMF_type', 'RF_type', 'ProbType1', 'ProbType2', 'ProbType3')
sum(predictions_Heidelberg[,1] != predictions_Heidelberg[,2])
#sum should be 0 if RF model is assigning all training data to the correct types (expected outcome)
predictions_Heidelberg <- cbind(rep('Heidelberg', dim(Heidelberg_Patel_Norm)[2]), predict(MethylationClassifier, newdata=data.frame(t(Heidelberg_Patel_Norm[intersect(rownames(HeidelbergNorm), c(CPGgroup1, CPGgroup2, CPGgroup3)),])), type="response"), predict(MethylationClassifier, newdata=data.frame(t(Heidelberg_Patel_Norm[intersect(rownames(HeidelbergNorm), c(CPGgroup1, CPGgroup2, CPGgroup3)),])), type="prob"))
predictions_Heidelberg[rownames(SampleData), 1] <- SampleData$MethGroup
colnames(predictions_Heidelberg) <- c('NMF_type', 'RF_type', 'ProbType1', 'ProbType2', 'ProbType3')

#write.csv(predictions_Heidelberg, 'Epigenetic Classifier - Heidelberg assignments.csv')

#Use NMF to identify clusters from Heidelberg samples using methylation signature probes
rank_estimate_Heidelberg <- nmf(HeidelbergNorm[intersect(rownames(HeidelbergNorm), c(CPGgroup1, CPGgroup2, CPGgroup3)),], 2:10, nrun=100,  .opt = 'v')
plot(rank_estimate_Heidelberg, main = 'Rank estimate for Heidelberg tumors over signature probes')
plot(2:10, rank_estimate_Heidelberg$measures$cophenetic, main = paste0('Heidelberg NMF Cophenetic coefficient'), type = 'b', xlab = '', ylab = '', col = 'purple', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
plot(2:10, rank_estimate_Heidelberg$measures$residuals, main = paste0('Heidelberg NMF Residuals'), type = 'b', xlab = '', ylab = '', col = 'dark blue', pch = 19, cex = 2.1, lwd = 3, xaxt='n')
axis(1, xaxp = c(2,10,8), las = 1)
consensusmap(rank_estimate_Heidelberg, color = 'blue', main = NULL, annCol = NULL, annRow = NULL, tracks = NULL, info = FALSE, labRow = NA, labCol = NA, legend = NULL)

#3 appears to be optimal rank, now perform final NMF
Rank <- 3
NMF_Heidelberg <- nmf(HeidelbergNorm[intersect(rownames(HeidelbergNorm), c(CPGgroup1, CPGgroup2, CPGgroup3)),], rank = Rank, nrun = 1000, .opt = 'v')
consensusmap(NMF_Heidelberg)
consensus_Heidelberg <- predict(NMF_Heidelberg, what = 'consensus')
silh_Heidelberg <- silhouette(NMF_Heidelberg)
clusters_Heidelberg <- cbind(rep(0, length(consensus_Heidelberg)), consensus_Heidelberg, silh_Heidelberg)
colnames(clusters_Heidelberg)[1] <- 'MethGroup'
clusters_Heidelberg[clusters_Heidelberg[,'consensus_Heidelberg'] == 1,'MethGroup'] <- 3
clusters_Heidelberg[clusters_Heidelberg[,'consensus_Heidelberg'] == 2,'MethGroup'] <- 2
clusters_Heidelberg[clusters_Heidelberg[,'consensus_Heidelberg'] == 3,'MethGroup'] <- 1
#write.csv(clusters_Heidelberg, 'Heidelberg NMF Predictions.csv')


Kmeans_Heidelberg <- ConsensusClusterPlus(HeidelbergNorm[intersect(rownames(HeidelbergNorm), c(CPGgroup1, CPGgroup2, CPGgroup3)),], maxK=10, reps=500, pItem=0.8, pFeature=1, title='K-means clustering of Heidelberg samples', distance = 'euclidean', clusterAlg='km', plot='png')
maps_Heidleberg <- list()
for(val in 1:9) {
  maps_Heidleberg[[val]] <- Kmeans_Heidelberg[[val+1]]$consensusMatrix
}
consensusmap(maps_Heidleberg, color = 'blue', main = NULL, legend = FALSE, labRow = NA, labCol = NA, Rowv = FALSE)



dendrogram_Heidelberg <- as.dendrogram(Kmeans_Heidelberg[[3]]$consensusTree)
colors <- predictions_Heidelberg[colnames(HeidelbergNorm)[Kmeans_Heidelberg[[3]]$consensusTree$order], 'RF_type']

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



#tSNE of Patel and Heidelberg samples
colors <- cbind(predictions[, 1:2])
colors[colors[,1] == 'Heidelberg',1] <- 'black'
colors[colors[,1] == '1',1] <- colors_by_group[1]
colors[colors[,1] == '2',1] <- colors_by_group[2]
colors[colors[,1] == '3',1] <- colors_by_group[3]
colors[colors[,2] == '1',2] <- colors_by_group[1]
colors[colors[,2] == '2',2] <- colors_by_group[2]
colors[colors[,2] == '3',2] <- colors_by_group[3]

n <- 4500
tsne_output <- Rtsne(t(Heidelberg_Patel_Norm[1:n,]), pca = FALSE, perplexity = 8)
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,1], main = paste0('tSNE plot of Heidelberg methylation data: Top ', n, ' most variable probes'))
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,2], main = paste0('tSNE plot of Heidelberg methylation data: Top ', n, ' most variable probes'))

tsne_output <- Rtsne(t(Heidelberg_Patel_Norm[intersect(rownames(Heidelberg_Patel_Norm), c(CPGgroup1, CPGgroup2, CPGgroup3)),]), pca = FALSE, perplexity = 8)
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,1], main = 'tSNE plot of Heidelberg methylation data, signature probes')
plot(tsne_output$Y, pch = 20, asp=0, cex = 2, xlab = "", ylab = "", col = colors[,2], main = 'tSNE plot of Heidelberg methylation data, signature probes')



#Compare with Olar et al probes
setwd('C:/Users/james/Desktop/Patel_Methylation/All_Patel_Data/')
Olar_CPG <- load_Olar_CPG()
Olar_CPG$MMFav <- intersect(Olar_CPG$MMFav, rownames(methData))
Olar_CPG$MMUnFav <- intersect(Olar_CPG$MMUnFav, rownames(methData))
Olar_CPG$MMBoth <- intersect(Olar_CPG$MMBoth, rownames(methData))
Olar_CPG$MMFav64 <- intersect(Olar_CPG$MMFav64, rownames(methData))
Olar_CPG$MMUnFav64 <- intersect(Olar_CPG$MMUnFav64, rownames(methData))
Olar_CPG$MMBoth64 <- intersect(Olar_CPG$MMBoth64, rownames(methData))

#Compare intersection of Olar CpGs to methylation signatures and also all methylation group DMPs
length(Olar_CPG$MMBoth)
length(intersect(c(CPGgroup1, CPGgroup2, CPGgroup3), Olar_CPG$MMBoth))
length(intersect(CPGgroup1, Olar_CPG$MMBoth))
length(intersect(CPGgroup2, Olar_CPG$MMBoth))
length(intersect(CPGgroup3, Olar_CPG$MMBoth))

Olar_AvgBeta <- matrix(nrow = length(Olar_CPG$MMBoth), ncol = 6, dimnames = list(Olar_CPG$MMBoth, c('Type1', 'Type2', 'Type3','NotType1', 'NotType2', 'NotType3')))
Olar_AvgBeta[,'Type1'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethGroup == 1))])
Olar_AvgBeta[,'Type2'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethGroup == 2))])
Olar_AvgBeta[,'Type3'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethGroup == 3))])
Olar_AvgBeta[,'NotType1'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethGroup %in% c(2,3)))])
Olar_AvgBeta[,'NotType2'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethGroup %in% c(1,3)))])
Olar_AvgBeta[,'NotType3'] <- rowMeans(methData[Olar_CPG$MMBoth, rownames(subset(SampleData, MethGroup %in% c(1,2)))])

#Compare average methylation between Meth types for favorable probes (positive favors first type, negative favors second type)
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMFav,'NotType1']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type2'] - Olar_AvgBeta[Olar_CPG$MMFav,'NotType2']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type3'] - Olar_AvgBeta[Olar_CPG$MMFav,'NotType3']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMFav,'Type2']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type1'] - Olar_AvgBeta[Olar_CPG$MMFav,'Type3']))
table(sign(Olar_AvgBeta[Olar_CPG$MMFav,'Type2'] - Olar_AvgBeta[Olar_CPG$MMFav,'Type3']))

#Compare average methylation between Meth types for unfavorable probes (negative favors first type, positive favors second type)
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
column_ha <- HeatmapAnnotation(MethType = SampleData[colnames(methData), 'MethGroup'],
                               col = list(MethType = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242')))
Heatmap(methData[c(Olar_CPG$MMFav, Olar_CPG$MMUnFav),], cluster_rows = Olar_CPG_dendrogram, row_split = 2, column_split = SampleData[colnames(methData), 'MethGroup'], column_order = order.dendrogram(consensushc(NMFfinal, method = 'mcquitty')), top_annotation = column_ha, column_title = 'Heatmap of Olar probes by methylation types', row_title = 'Olar unfavorable and favorable probes', show_row_dend = TRUE, heatmap_legend_param = list(title = 'Beta'), show_row_names = FALSE, show_column_names = FALSE, use_raster = FALSE)


#Kmeans of our samples over Olar probes
Kmeans_Olar <- ConsensusClusterPlus(methData[c(Olar_CPG$MMFav, Olar_CPG$MMUnFav),], maxK=10, reps=500, pItem=0.8, pFeature=1, title='K-means clustering of Olar probes', clusterAlg='km', distance='euclidean', plot='png')

maps_Olar <- list()
for(val in 1:9) {
  maps_Olar[[val]] <- Kmeans_Olar[[val+1]]$consensusMatrix
}
consensusmap(maps_Olar, color = 'blue', main = NULL, legend = FALSE, labRow = NA, labCol = NA, Rowv = FALSE)

n <- 5
Heatmap(methData[c(Olar_CPG$MMFav, Olar_CPG$MMUnFav),], cluster_rows = Olar_CPG_dendrogram, row_split = 2, cluster_columns = Kmeans_Olar[[n]]$consensusTree, column_split = n, top_annotation = column_ha, column_title = 'Heatmap of Olar probes by methylation types', row_title = 'Olar unfavorable and favorable probes', show_row_dend = TRUE, heatmap_legend_param = list(title = 'Beta'), show_row_names = FALSE, show_column_names = FALSE, use_raster = FALSE)


####################################################################
#Comparison of methylation types to transcriptional and other types#
####################################################################

##############################################################
#Establish 4th criteria to decide overall groups, using a NF2 cutoff and chromosomal losses
merlin <- cbind(expression['NF2', colnames(methData)], SampleData[colnames(methData),]$MethGroup)
colnames(merlin) <- c('Merlin', 'Group')
merlinConcordant <- merlin[rownames(subset(SampleData, AllConcordant == 'Yes')),]

boxplot(Merlin~Group,data=merlin, main='NF2 expression for all samples', xlab='Methylation Type', ylab='Log NF2 expression')
plot(Merlin~Group,data=merlin, main='NF2 expression by methylation type', xlab='Methylation Type', ylab='Log NF2 expression')

boxplot(Merlin~Group,data=merlinConcordant, main='NF2 expression for concordant samples (Meth, RNA, CNV)', xlab='Methylation Type', ylab='Log NF2 expression')
plot(Merlin~Group,data=merlinConcordant, main='NF2 expression for concordant samples (Meth, RNA, CNV)', xlab='Methylation Type', ylab='Log NF2 expression')
#Cutoff of 11 selected on basis of using only concordant samples

#Create Sankey plot for classifications
library(networkD3)
nodes = data.frame('name' = 
                     c('Epigenetic type 1', # Node 0
                       'Epigenetic type 2', # Node 1
                       'Epigenetic type 3', # Node 2
                       'Transcriptional type A', # Node 3
                       'Transcriptional type B', # Node 4
                       'Transcriptional type C', # Node 5
                       'No losses', # Node 6
                       '22 loss', # Node 7
                       '1p/22 loss', # Node 8
                       'NF2-intact', # Node 9
                       'NF2-deficient, low-CNV', # Node 10
                       'NF2-deficient, high-CNV', # Node 11
                       'MenG A', # Node 12
                       'MenG B', # Node 13
                       'MenG C', # Node 14
                       'MenG Unknown')) # Node 15

# Each row represents a link. The first number = node being conntected from,  second number = node connected to, third number = value of the node
links = as.data.frame(matrix(c(
  0, 3, dim(subset(SampleData, MethGroup == '1' & RNA_numerical == '1'))[1],
  0, 4, dim(subset(SampleData, MethGroup == '1' & RNA_numerical == '2'))[1],
  0, 5, dim(subset(SampleData, MethGroup == '1' & RNA_numerical == '3'))[1],
  1, 3, dim(subset(SampleData, MethGroup == '2' & RNA_numerical == '1'))[1],
  1, 4, dim(subset(SampleData, MethGroup == '2' & RNA_numerical == '2'))[1],
  1, 5, dim(subset(SampleData, MethGroup == '2' & RNA_numerical == '3'))[1],
  2, 3, dim(subset(SampleData, MethGroup == '3' & RNA_numerical == '1'))[1],
  2, 4, dim(subset(SampleData, MethGroup == '3' & RNA_numerical == '2'))[1],
  2, 5, dim(subset(SampleData, MethGroup == '3' & RNA_numerical == '3'))[1],
  3, 6, dim(subset(SampleData, RNA_numerical == '1' & CNVGroup == '1'))[1],
  3, 7, dim(subset(SampleData, RNA_numerical == '1' & CNVGroup == '2'))[1],
  3, 8, dim(subset(SampleData, RNA_numerical == '1' & CNVGroup == '3'))[1],
  4, 6, dim(subset(SampleData, RNA_numerical == '2' & CNVGroup == '1'))[1],
  4, 7, dim(subset(SampleData, RNA_numerical == '2' & CNVGroup == '2'))[1],
  4, 8, dim(subset(SampleData, RNA_numerical == '2' & CNVGroup == '3'))[1],
  5, 6, dim(subset(SampleData, RNA_numerical == '3' & CNVGroup == '1'))[1],
  5, 7, dim(subset(SampleData, RNA_numerical == '3' & CNVGroup == '2'))[1],
  5, 8, dim(subset(SampleData, RNA_numerical == '3' & CNVGroup == '3'))[1],
  6, 9, dim(subset(SampleData, CNVGroup == '1' & FourthGroup == '1'))[1],
  6, 10, dim(subset(SampleData, CNVGroup == '1' & FourthGroup == '2'))[1],
  6, 11, dim(subset(SampleData, CNVGroup == '1' & FourthGroup == '3'))[1],
  7, 9, dim(subset(SampleData, CNVGroup == '2' & FourthGroup == '1'))[1],
  7, 10, dim(subset(SampleData, CNVGroup == '2' & FourthGroup == '2'))[1],
  7, 11, dim(subset(SampleData, CNVGroup == '2' & FourthGroup == '3'))[1],
  8, 9, dim(subset(SampleData, CNVGroup == '3' & FourthGroup == '1'))[1],
  8, 10, dim(subset(SampleData, CNVGroup == '3' & FourthGroup == '2'))[1],
  8, 11, dim(subset(SampleData, CNVGroup == '3' & FourthGroup == '3'))[1],
  9, 12, dim(subset(SampleData, FourthGroup == '1' & CompositeGroup == '1'))[1],
  9, 13, dim(subset(SampleData, FourthGroup == '1' & CompositeGroup == '2'))[1],
  9, 14, dim(subset(SampleData, FourthGroup == '1' & CompositeGroup == '3'))[1],
  9, 15, dim(subset(SampleData, FourthGroup == '1' & CompositeGroup == 'Unknown'))[1],
  10, 12, dim(subset(SampleData, FourthGroup == '2' & CompositeGroup == '1'))[1],
  10, 13, dim(subset(SampleData, FourthGroup == '2' & CompositeGroup == '2'))[1],
  10, 14, dim(subset(SampleData, FourthGroup == '2' & CompositeGroup == '3'))[1],
  10, 15, dim(subset(SampleData, FourthGroup == '2' & CompositeGroup == 'Unknown'))[1],
  11, 12, dim(subset(SampleData, FourthGroup == '3' & CompositeGroup == '1'))[1],
  11, 13, dim(subset(SampleData, FourthGroup == '3' & CompositeGroup == '2'))[1],
  11, 14, dim(subset(SampleData, FourthGroup == '3' & CompositeGroup == '3'))[1],
  11, 15, dim(subset(SampleData, FourthGroup == '3' & CompositeGroup == 'Unknown'))[1]),
  byrow = TRUE, ncol = 3))
names(links) = c('source', 'target', 'value')
sankeyNetwork(Links = links, Nodes = nodes,
              Source = 'source', Target = 'target',
              Value = 'value', NodeID = 'name',
              fontSize= 12, nodeWidth = 30)




#Oncoprint
oncoprint_order <- read.csv('Oncoprint_order.csv')
oncoprint_cnv <- read.csv('Oncoprint CNVs.csv', row.names = 1)
colnames(oncoprint_cnv) <- sub('.', '', colnames(oncoprint_cnv))
oncoprint_cnv[is.na(oncoprint_cnv)]<-''
oncoprint_cnv <- t(oncoprint_cnv[oncoprint_order[,1],])
oncoprint_data <- SampleData[oncoprint_order[,1],]

#Set heatmap column settings for A/B/C etc
column_ha <- HeatmapAnnotation(MenG = oncoprint_data$CompositeGroup, MethType = oncoprint_data$MethGroup, RNAType = oncoprint_data$RNA_class, CNVType = oncoprint_data$CNVGroup, NF2Instability = oncoprint_data$FourthGroup, 
                               Gender = oncoprint_data$Gender, WHO.grade = oncoprint_data$WHO.grade,
                               TRAF7 = oncoprint_data$TRAF7, KLF4 = oncoprint_data$KLF4, AKT1 = oncoprint_data$AKT1, NF2 = oncoprint_data$NF2, SMARCB1 = oncoprint_data$SMARCB1, PIK3CA = oncoprint_data$PIK3CA, SMO = oncoprint_data$SMO, POLR2A = oncoprint_data$POLR2A, 
                               Ch1pLoss = oncoprint_data$Ch1pLoss, Ch22Loss = oncoprint_data$Ch22Loss,
                               cbar = anno_oncoprint_barplot(height = unit(1, 'cm')),
                               col = list(MenG = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242', 'Unknown' = 'black'),
                                          MethType = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242'),
                                          RNAType = c('A' = '#6FDE6E', 'B' = '#235FA4', 'C' = '#FF4242'),
                                          CNVType = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242'),
                                          NF2Instability = c('1' = '#6FDE6E', '2' = '#235FA4', '3' = '#FF4242'),
                                          Gender = c('M' = 'blue', 'F' = 'pink'),
                                          WHO.grade = c('WHO I' = 'yellow', 'WHO II' = 'purple'),
                                          Ch1pLoss = c('Y' = 'red', 'N' = 'gray95'),
                                          Ch22Loss = c('Y' = 'red', 'N' = 'gray95'),
                                          TRAF7 = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white'),
                                          KLF4 = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white'),
                                          AKT1 = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white'),
                                          NF2 = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white'),
                                          SMARCB1 = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white'),
                                          PIK3CA = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white'),
                                          SMO = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white'),
                                          POLR2A = c('Y' = 'black', 'N' = 'gray95', 'Unknown' = 'white')
                               ),
                               annotation_legend_param = list(
                                 MenG = list(
                                   title = 'Meningioma Group',
                                   at = c(1, 2, 3, 'Unknown'),
                                   labels = c('MenG A', 'MenG B', 'MenG C', 'Unknown')
                                 ),
                                 MethType = list(
                                   title = 'Epigenetic Type',
                                   at = c(1, 2, 3),
                                   labels = c('Type 1', 'Type 2', 'Type 3')
                                 ),
                                 RNAType = list(
                                   title = 'Transcriptional Type',
                                   at = c('A', 'B', 'C'),
                                   labels = c('Type A', 'Type B', 'Type C')
                                 ),
                                 CNVType = list(
                                   title = 'Chromosomal Type',
                                   at = c(1, 2, 3),
                                   labels = c('No loss', 'Chr22 loss', 'Chr1p loss')
                                 ),
                                 NF2Instability = list(
                                   title = 'NF2/Instability Type',
                                   at = c(1, 2, 3),
                                   labels = c('NF2-intact', 'NF2-deficient, low-CNV', 'NF2-deficient, high-CNV')
                                 ),
                                 Gender = list(
                                   title = 'Sex',
                                   at = c('F', 'M'),
                                   labels = c('Female', 'Male')
                                 ),
                                 WHO.grade = list(
                                   title = 'Histologic grade',
                                   at = c('WHO I', 'WHO II'),
                                   labels = c('WHO grade I', 'WHO grade II')
                                 ),
                                 TRAF7 = list(
                                   title = 'Loss of function mutation',
                                   at = c('Y', 'N'),
                                   labels = c('Yes', 'No')
                                 ),
                                 Ch1pLoss = list(
                                   title = 'Large-scale deletion',
                                   at = c('Y', 'N'),
                                   labels = c('Yes', 'No')
                                 )
                               )
)

col = c('G' = 'green', 'L' = 'red')

alter_fun = list(
  'G' = alter_graphic('rect', fill = 'green'),
  'L' = alter_graphic('rect', fill = 'red')
)

oncoPrint(oncoprint_cnv[c(1,44),], alter_fun = alter_fun, col = col, alter_fun_is_vectorized = FALSE, top_annotation = column_ha, show_pct = FALSE, row_order = 1:2, column_order = 1:ncol(oncoprint_cnv), column_split = oncoprint_data$CompositeGroup, show_column_names = FALSE, column_names_rot = 45)


#Original version
column_ha <- HeatmapAnnotation(MenG = oncoprint_data$CompositeGroup, MethType = oncoprint_data$MethGroup, SesameGroup = oncoprint_data$SesameGroup, RNAGroup = oncoprint_data$RNA_class, CNVGroup = oncoprint_data$CNVGroup, FourthCriteria = oncoprint_data$FourthGroup, Ch1pLoss = oncoprint_data$Ch1pLoss, Ch22Loss = oncoprint_data$Ch22Loss, Gender = oncoprint_data$Gender, WHO.grade = oncoprint_data$WHO.grade, Simpson_Grade = oncoprint_data$Simpson_grade, Recurrence = oncoprint_data$Recurrence,
                               NF2 = oncoprint_data$NF2, SMARCB1 = oncoprint_data$SMARCB1, TRAF7 = oncoprint_data$TRAF7, KLF4 = oncoprint_data$KLF4, AKT1 = oncoprint_data$AKT1, PIK3CA = oncoprint_data$PIK3CA, SMO = oncoprint_data$SMO, POLR2A = oncoprint_data$POLR2A, ARID1A = oncoprint_data$ARID1A, TERTpromoter = oncoprint_data$TERT.Promoter.Mutation,
                               cbar = anno_oncoprint_barplot(height = unit(1, 'cm')),
                               col = list(MenG = c('1' = 'dark green', '2' = 'dark blue', '3' = 'maroon', 'Unknown' = 'dark grey'),
                                          MethType = c('1' = 'green', '2' = 'blue', '3' = 'red'),
                                          SesameGroup = c('1' = 'light green', '2' = 'light blue', '3' = 'pink'),
                                          RNAGroup = c('A' = 'dark green', 'B' = 'dark blue', 'C' = 'maroon', 'Unknown' = 'grey'),
                                          CNVGroup = c('1' = 'green', '2' = 'blue', '3' = 'red'),
                                          FourthCriteria = c('1' = 'dark green', '2' = 'dark blue', '3' = 'maroon'),
                                          Ch1pLoss = c('Y' = 'red', 'N' = 'grey'),
                                          Ch22Loss = c('Y' = 'red', 'N' = 'grey'),
                                          Gender = c('M' = 'blue', 'F' = 'pink'),
                                          WHO.grade = c('WHO I' = 'yellow', 'WHO II' = 'purple'),
                                          Simpson_Grade = c('1' = 'green', '2' = 'orange', '4' =  'blue', 'Unknown' = 'grey'),
                                          Recurrence = c('Y' = 'red', 'N' = 'grey'),
                                          NF2 = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          SMARCB1 = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          TRAF7 = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          KLF4 = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          AKT1 = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          PIK3CA = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          SMO = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          POLR2A = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          ARID1A = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white'),
                                          TERTpromoter = c('Y' = 'black', 'N' = 'grey', 'Unknown' = 'white')
                               )
)

#Merlin for MenG groups and subtypes
merlin <- data.frame(row.names = c(colnames(methData), 'Arachnoid_1', 'Arachnoid_2', 'Arachnoid_3'), 
                                    NF2 = expression['NF2', c(colnames(methData), 'Arachnoid_1', 'Arachnoid_2', 'Arachnoid_3')], 
                                    MethGroup = c(SampleData[colnames(methData),]$MethGroup, rep("Arachnoid",3)),
                                    RNAGroup = c(SampleData[colnames(methData),]$RNA_numerical, rep("Arachnoid",3)),
                                    CNVGroup = c(SampleData[colnames(methData),]$CNVGroup, rep("Arachnoid",3)),
                                    NF2InstabilityGroup = c(SampleData[colnames(methData),]$FourthGroup, rep("Arachnoid",3)),
                                    MenG = c(SampleData[colnames(methData),]$CompositeGroup, rep("Arachnoid",3)),
                                    WHO.grade = c(SampleData[colnames(methData),]$WHO.grade, rep("Arachnoid",3)))
merlin <- subset(merlin, MenG != 'Unknown')
merlin$MethGroup <- factor(merlin$MethGroup , levels=c("Arachnoid", "1", "2", "3"))
merlin$RNAGroup <- factor(merlin$RNAGroup , levels=c("Arachnoid", "1", "2", "3"))
merlin$CNVGroup <- factor(merlin$CNVGroup , levels=c("Arachnoid", "1", "2", "3"))
merlin$NF2InstabilityGroup <- factor(merlin$NF2InstabilityGroup , levels=c("Arachnoid", "1", "2", "3"))
merlin$MenG <- factor(merlin$MenG , levels=c("Arachnoid", "1", "2", "3"))
merlin$WHO.grade <- factor(merlin$WHO.grade , levels=c("Arachnoid", "WHO I", "WHO II"))

boxplot(NF2~MenG,data=merlin, main='NF2 expression in Meningioma Groups', xlab='MenG', ylab='Log NF2 expression', col = c('grey', colors_by_group[1:3]))
boxplot(NF2~drop.levels(MenG), data=subset(merlin, MenG != 'Arachnoid'), main='NF2 expression in Meningioma Groups', xlab='MenG', ylab='Log NF2 expression', col = colors_by_group[1:3])

boxplot(NF2~MethGroup,data=merlin, main='NF2 expression in Methylation Types', xlab='Methylation Type', ylab='Log NF2 expgetwd()ression', col = c('grey', colors_by_group[1:3]))
boxplot(NF2~drop.levels(MethGroup), data=subset(merlin, MenG != "Arachnoid"), main='NF2 expression in Methylation Types', xlab='Methylation Type', ylab='Log NF2 expression', col = colors_by_group[1:3])

boxplot(NF2~RNAGroup,data=merlin, main='NF2 expression in Transcriptional Types', xlab='Transcriptional Type', ylab='Log NF2 expression', col = c('grey', colors_by_group[1:3]))
boxplot(NF2~drop.levels(RNAGroup), data=subset(merlin, MenG != "Arachnoid"), main='NF2 expression in Transcriptional Types', xlab='Transcriptional Type', ylab='Log NF2 expression', col = colors_by_group[1:3])

boxplot(NF2~CNVGroup,data=merlin, main='NF2 expression in Chromosomal Types', xlab='Chromosomal Type', ylab='Log NF2 expression', col = c('grey', colors_by_group[1:3]))
boxplot(NF2~drop.levels(CNVGroup), data=subset(merlin, MenG != "Arachnoid"), main='NF2 expression in Chromosomal Types', xlab='Chromosomal Type', ylab='Log NF2 expression', col = colors_by_group[1:3])

boxplot(NF2~NF2InstabilityGroup,data=merlin, main='NF2 expression in NF2/Instability Types', xlab='NF2/Instability Type', ylab='Log NF2 expression', col = c('grey', colors_by_group[1:3]))
boxplot(NF2~drop.levels(NF2InstabilityGroup), data=subset(merlin, MenG != "Arachnoid"), main='NF2 expression in NF2/Instability Types', xlab='NF2/Instability Type', ylab='Log NF2 expression', col = colors_by_group[1:3])


#MIB for MenG and subtypes
MIB <- data.frame(row.names = colnames(methData), 
                    MIB = as.numeric(SampleData[colnames(methData), 'MIB']), 
                    MethGroup = SampleData[colnames(methData),]$MethGroup,
                    RNAGroup = SampleData[colnames(methData),]$RNA_numerical,
                    CNVGroup = SampleData[colnames(methData),]$CNVGroup,
                    NF2InstabilityGroup = SampleData[colnames(methData),]$FourthGroup,
                    MenG = SampleData[colnames(methData),]$CompositeGroup,
                    WHO.grade = SampleData[colnames(methData),]$WHO.grade)
MIB$MethGroup <- factor(MIB$MethGroup , levels=c("1", "2", "3"))
MIB$RNAGroup <- factor(MIB$RNAGroup , levels=c("1", "2", "3"))
MIB$CNVGroup <- factor(MIB$CNVGroup , levels=c("1", "2", "3"))
MIB$NF2InstabilityGroup <- factor(MIB$NF2InstabilityGroup , levels=c("1", "2", "3"))
MIB$MenG <- factor(MIB$MenG , levels=c("1", "2", "3", "Unknown"))

boxplot(MIB~MenG,data=subset(MIB, MIB != 'NA'), main='MIB-1 index in Meningioma Groups', xlab='MenG', ylab='MIB-1 Index (%)', col = c(colors_by_group[1:3], 'grey'))
boxplot(MIB~drop.levels(MenG),data=subset(MIB, MIB != 'NA' & MenG != 'Unknown'), main='MIB-1 index in Meningioma Groups', xlab='MenG', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])

boxplot(MIB~drop.levels(MenG),data=subset(MIB, MIB != 'NA' & MenG != 'Unknown' & WHO.grade == "WHO I"), main='MIB-1 index in WHO I tumors', xlab='MenG', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])
plot(MIB~drop.levels(MenG),data=subset(MIB, MIB != 'NA' & MenG != 'Unknown' & WHO.grade == "WHO II"), main='MIB-1 index in WHO II tumors', xlab='MenG', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])

plot(MIB~drop.levels(MenG),data=subset(MIB, MIB != 'NA' & MenG != 'Unknown' & WHO.grade == "WHO I"), type = 'p', main='MIB-1 index in WHO I tumors', xlab='MenG', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])
plot(MIB~drop.levels(MenG),data=subset(MIB, MIB != 'NA' & MenG != 'Unknown' & WHO.grade == "WHO II"), main='MIB-1 index in WHO II tumors', xlab='MenG', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])

boxplot(MIB~MethGroup,data=subset(MIB, MIB != 'NA'), main='MIB-1 index in Methylation Types', xlab='Methylation Type', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])
boxplot(MIB~RNAGroup,data=subset(MIB, MIB != 'NA'), main='MIB-1 index in Transcriptional Types', xlab='Transcriptional Type', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])
boxplot(MIB~CNVGroup,data=subset(MIB, MIB != 'NA'), main='MIB-1 index in Chromosomal Types', xlab='Chromosomal Type', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])
boxplot(MIB~NF2InstabilityGroup,data=subset(MIB, MIB != 'NA'), main='MIB-1 index in NF2/Instability Types', xlab='NF2/Instability Type', ylab='MIB-1 Index (%)', col = colors_by_group[1:3])


########################
#Kaplan Meier curves

#All tumors by histology
ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ WHO.grade, data = SampleData), 
  pval=TRUE,
  legend.labs=c('WHO grade I', 'WHO grade II'), legend.title='', legend = 'bottom',
  palette=c('orange', 'purple'),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'All tumors by histologic grade',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

#All tumors by Mengioma Group
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3)), 
  pval=TRUE,
  legend.labs=c('MenG A', 'MenG B', 'MenG C'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'All tumors by Meningioma Group',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3)))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)

#All tumors by histology and MenG
tempData <- rbind(SampleData, SampleData)
tempData[rownames(SampleData), 'WHO.grade'] <- SampleData$CompositeGroup
tempData <- subset(tempData, CompositeGroup != 'Unknown')
ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ WHO.grade, data = tempData), 
  pval=TRUE,
  legend.labs=c('MenG A', 'MenG B', 'MenG C', 'WHO grade I', 'WHO grade II'), legend.title='', legend = 'bottom',
  palette=c('#6FDE6E', '#235FA4', '#FF4242', 'light gray', 'dark gray'),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'All tumors by histologic grade',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

#All tumors by Methylation Type
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ MethGroup, data = SampleData), 
  pval=TRUE,
  legend.labs=c('Methylation Type A', 'Methylation Type B', 'Methylation Type C'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'All tumors by Methylation Type',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ MethGroup, data = SampleData))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)

#All tumors by Transcriptional Type
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ RNA_numerical, data = SampleData), 
  pval=TRUE,
  legend.labs=c('Transcriptional Type A', 'Transcriptional Type B', 'Transcriptional Type C'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'All tumors by Transcriptional Type',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ RNA_numerical, data = SampleData))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)

#All tumors by Chromosomal Type
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ CNVGroup, data = SampleData), 
  pval=TRUE,
  legend.labs=c('No losses', 'Ch22 loss', 'Ch1p/22 loss'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'All tumors by Chromosomal Type',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ CNVGroup, data = SampleData))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)

#All tumors by NF2/Instability Type
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ FourthGroup, data = SampleData), 
  pval=TRUE,
  legend.labs=c('NF2-intact', 'NF2-deficient, low-CNV', 'NF2-deficient, high-CNV'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'All tumors by NF2/Instability Type',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ FourthGroup, data = SampleData))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)


#Gross total resection by histology
ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ WHO.grade, data = subset(SampleData, Simpson_grade %in% 1:3)), 
  pval=TRUE,
  legend.labs=c('WHO grade I', 'WHO grade II'), legend.title='', legend = 'bottom',
  palette=c('orange', 'purple'),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'Gross total resection by histologic grade',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)


#Gross total resection by Meningioma Group
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3 & Simpson_grade %in% 1:3)), 
  pval=TRUE,
  legend.labs=c('MenG A', 'MenG B', 'MenG C'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'Gross total resection by Meningioma Group',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3 & Simpson_grade %in% 1:3)))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)


#WHO grade I by Meningioma Group
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3 & WHO.grade == "WHO I")), 
  pval=TRUE,
  legend.labs=c('MenG A', 'MenG B', 'MenG C'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'WHO grade I by Meningioma Group',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3 & WHO.grade == "WHO I")))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)

#WHO grade II by Meningioma Group
ggsurv <- ggsurvplot(
  fit = survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3 & WHO.grade == "WHO II")), 
  pval=TRUE,
  legend.labs=c('MenG A', 'MenG B', 'MenG C'), legend.title='', legend = 'bottom',
  palette=as.vector(c(colors_by_group[1], colors_by_group[2], colors_by_group[3])),
  xlab = 'Years',
  ylab = 'Recurrence free survival',
  break.time.by = 2,
  xlim = c(0,7),
  surv.scale = 'percent',
  title = 'WHO grade II by Meningioma Group',
  censor = FALSE,
  font.title = 24,
  font.x = 24,
  font.y = 24,
  font.legend = 18,
  pval.size = 12,
  size = 2,
  font.tickslab = 20)

surv_median <- as.vector(summary(survfit(Surv(KM_survival/365.26, KM_outcome) ~ CompositeGroup, data = subset(SampleData, CompositeGroup %in% 1:3 & WHO.grade == "WHO II")))$table[, 'median'])[3]
ggsurv$plot <- ggsurv$plot + geom_segment(aes(x = 0, y = 0.5, xend = surv_median, yend = 0.5), linetype = 'dashed', size = 1.5) +
  geom_segment(aes(x = surv_median, y = 0.5 , xend = surv_median, yend = 0), linetype = 'dashed', size = 1.5)
print(ggsurv)


############
#PLS models#
############

CNV_meth_probes <- read.csv('CNV from methylation probes.csv', row.names = 1)
CNV_meth_segments <- read.csv('CNV from methylation segments.csv', row.names = 1)

setwd('C:/Users/james/Desktop/Patel_Methylation/All_Patel_Data/Meth based CNV - probes/')
setwd('C:/Users/james/Desktop/Patel_Methylation/All_Patel_Data/Meth based CNV - segments/')

#Load and CNV data
CNV <- read.csv('Copy number status by gene.csv', row.names = 1)
colnames(CNV) <- gsub('X', '', gsub('\\.', '-', colnames(CNV)))

##################################################
#TROUBLE SHOOTING
CNV <- CNV_meth_probes
CNV <- CNV_meth_segments
rbind(CNV['POLR2J',], CNV['POLR2J2',], CNV['POLR2J3',])
######################################################

#Limit CpGs to only those desired to be used in PLS modeling
CPG_promoter <- rownames(subset(CPG[rownames(methData),], feature %in% c('1stExon', 'TSS200', 'TSS1500')))

PLSGroup1 <- intersect(rownames(subset(SampleData, RNA_numerical == 1 & MethGroup == 1 & CNVGroup == 1 & FourthGroup == 1)), colnames(CNV))
PLSGroup2 <- intersect(rownames(subset(SampleData, RNA_numerical == 2 & MethGroup == 2 & CNVGroup == 2 & FourthGroup == 2)), colnames(CNV))
PLSGroup3 <- intersect(rownames(subset(SampleData, RNA_numerical == 3 & MethGroup == 3 & CNVGroup == 3 & FourthGroup == 3)), colnames(CNV))
PLSGroups <- c(PLSGroup1, PLSGroup2, PLSGroup3)

#Get detailed chromosome location data
ensembl = useDataset('hsapiens_gene_ensembl', mart = useMart('ENSEMBL_MART_ENSEMBL'))
filterType <- 'external_gene_name'
attributeNames <- c('external_gene_name', 'band', 'chromosome_name')
annotations <- getBM(attributes=attributeNames, filters = filterType, values = rownames(expression), mart = ensembl)
annotations <- subset(annotations, band != '')
annotations$ChrLoc <- paste0(annotations$chromosome_name, annotations$band)

#Get all genes with expression, CNV, and promoter CpG data
AllGenes <- intersect(rownames(expression), intersect(rownames(CNV), CPG[CPG_promoter,]$gene))

AllGeneR2s <- data.frame(row.names = AllGenes, Chr = rep(0, length(AllGenes)), Components = 0, Rsquared = 0, CNV_importance = 0, CNV_LMcoef = 0, CNV_LMpval = 0, CNV_LMadjR2 = 0,
                         CpGs = 0, TopCPG = 0, TopCPGImp = 0, TopCPG_LMcoef = 0, TopCPG_LMpval = 0, TopCPG_LMadjR2 = 0, 
                         SecondCPG = 0, SecondCPGImp = 0, SecondCPG_LMcoef = 0, SecondCPG_LMpval = 0, SecondCPG_LMadjR2 = 0, 
                         ThirdCPG = 0, ThirdCPGImp = 0, ThirdCPG_LMcoef = 0, ThirdCPG_LMpval = 0, ThirdCPG_LMadjR2 = 0)

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

for(val in 1:dim(AllGeneR2s)[1]) {
  AllGeneR2s$Chr[val] <- annotations[annotations$external_gene_name==rownames(AllGeneR2s)[val],]$ChrLoc
  if(rownames(AllGeneR2s)[val] %in% rownames(expression)) {
    AllGeneR2s$Group1Avg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSGroup1, colnames(CNV))])
    AllGeneR2s$Group2Avg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSGroup2, colnames(CNV))])
    AllGeneR2s$Group3Avg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSGroup3, colnames(CNV))])
    AllGeneR2s$OverallAvg[val] <- mean(expression[rownames(AllGeneR2s)[val], intersect(PLSGroups, colnames(CNV))])
  } else {
    AllGeneR2s$Group1Avg[val] <- NA
    AllGeneR2s$Group2Avg[val] <- NA
    AllGeneR2s$Group3Avg[val] <- NA
    AllGeneR2s$OverallAvg[val] <- NA
  }
  if(AllGeneR2s$TopCPG[val] %in% rownames(methData)) {
    AllGeneR2s$Group1TopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSGroup1, colnames(CNV))])
    AllGeneR2s$Group2TopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSGroup2, colnames(CNV))])
    AllGeneR2s$Group3TopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSGroup3, colnames(CNV))])
    AllGeneR2s$OverallTopCPGBeta[val] <- mean(methData[AllGeneR2s$TopCPG[val], intersect(PLSGroups, colnames(CNV))])
  } else {
    AllGeneR2s$Group1TopCPGBeta[val] <- NA
    AllGeneR2s$Group2TopCPGBeta[val] <- NA
    AllGeneR2s$Group3TopCPGBeta[val] <- NA
    AllGeneR2s$OverallTopCPGBeta[val] <- NA
  }
  
  message(val)
}

write.csv(AllGeneR2s, 'PLS models - All genes 7.31.21.csv')

save.image(file = 'After PLS of all genes.RData')










gene <- 'HSPG2'
Model <- getPLS(gene, methylation = methData[,PLSGroups], CPGs = CPG[CPG_promoter,], expressionSubset = expressionPLS)
Model
Model[[2]]$results$Rsquared[Model[[2]]$bestTune$ncomp]
data <- as.data.frame(Model[[2]]$trainingData)
colors <- SampleData[rownames(data),]$CompositeGroup
colors[colors == '2'] <- 4
colors[colors == '3'] <- 2
colors[colors == '1'] <- 3
plot(data[,c(2, 1)], col = colors, pch = 19, ylab = 'Log expession', xlab = 'CNV')
plot(data[,c(Model[[5]], 'expression')], col = colors, pch = 19, ylab = 'Log expession', xlab = 'Top CPG Beta')
plot(data[,c(Model[[7]], 'expression')], col = colors, pch = 19, ylab = 'Log expession', xlab = 'Second CPG Beta')
plot(data[,c(Model[[9]], 'expression')], col = colors, pch = 19, ylab = 'Log expession', xlab = 'Third CPG Beta')
#




abline(linearModel[[1]])
correlation[[2]]
correlation <- correlateExpression(data[,c('expression', 'CNV')])


Norm_promoter <- methData[CPG_promoter,]
Norm_promoter <- methData[order(rowVars(methData), decreasing = TRUE),]
promoterVarGenes <- unique(CPG[rownames(Norm_promoter),]$gene)
write.csv(AllGeneR2s[intersect(promoterVarGenes, rownames(AllGeneR2s)),], 'All Genes sorted by promoter variance.csv')

Ch1pGenes <- unique(subset(annotations[grep('1p', annotations$ChrLoc),], chromosome_name == 1)$external_gene_name)
Ch22Genes <- unique(subset(annotations, chromosome_name == 22)$external_gene_name)

write.csv(AllGeneR2s[intersect(c(Ch1pGenes, Ch22Genes), rownames(AllGeneR2s)),], 'Chromosome 1p and 22 genes.csv')


rfGenes <- read.csv('RF model gene importance.csv')
plot(rfGenes$average)
write.csv(AllGeneR2s[intersect(rfGenes$Gene, rownames(AllGeneR2s)),], 'Random forest genes.csv')

#P-value calculation for PLS results
totalA <- 15398
subsetA <- 1164

totalA <- 1164
subsetA <- 704
totalB <- 141
subsetB <- 15

subsetA/totalA
subsetB/totalB

table <- matrix(data = c(subsetA, subsetB, totalA-subsetA, totalB-subsetB), nrow = 2, ncol = 2)
chisq.test(table)
chisq.test(table)$p.value


pchisq(1187.1, 1, lower.tail=FALSE)

################################3
#DEGs
PLScounts <- counts[,c(PLSGroups, colnames(counts)[grep('Arachnoid', colnames(counts))])]
geneList <- getGeneList(ensemblList = rownames(PLScounts))
rownames(geneList) <- geneList$ensembl_gene_id
geneList <- geneList[!duplicated(geneList$external_gene_name),]
PLScounts <- PLScounts[rownames(PLScounts) %in% geneList$ensembl_gene_id,]
rownames(PLScounts) <- geneList[rownames(PLScounts), 'external_gene_name']

colData <- data.frame(row.names = colnames(PLScounts), type = rep(0, dim(PLScounts)[2]))
colData[rownames(colData) %in% rownames(SampleData),] <- SampleData[rownames(colData)[rownames(colData) %in% rownames(SampleData)], 'CompositeGroup']
colData[grep('Arachnoid', rownames(colData)),] <- 'Arachnoid'
colData$type <- factor(colData$type)
PLSdds <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~type)
PLSdds <- DESeq(PLSdds[rowSums(counts(PLSdds)) >= 10, ])
PLSexpression <- getExpressionData(PLSdds)

colData <- data.frame(row.names = colnames(PLScounts), groups = factor(c(SampleData[PLSGroups, 'CompositeGroup'], rep('Arachnoid',3))))
levels(colData$groups) <- c(levels(colData$groups), 'other')
colData$pheno1 <- colData$groups
colData$pheno1[colData$pheno1 %in% c(2,3)] <- 'other'
colData$pheno2 <- colData$groups
colData$pheno2[colData$pheno2 %in% c(1,3)] <- 'other'
colData$pheno3 <- colData$groups
colData$pheno3[colData$pheno3 %in% c(1,2)] <- 'other'
colData$arachnoid <- colData$groups
colData$arachnoid[colData$arachnoid %in% c(1,2,3)] <- 'other'

ddsGroups <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~groups)
ddsGroup1 <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~pheno1)
ddsGroup2 <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~pheno2)
ddsGroup3 <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~pheno3)
ddsArachnoid <- DESeqDataSetFromMatrix(countData = PLScounts, colData = colData, design = ~arachnoid)

#Process counts data with DEseq, removing anything with low total count numbers
groupDDS <- DESeq(ddsGroups[rowSums(counts(ddsGroups)) >= 10])
group1DDS <- DESeq(ddsGroup1[rowSums(counts(ddsGroup1)) >= 10])
group2DDS <- DESeq(ddsGroup2[rowSums(counts(ddsGroup2)) >= 10])
group3DDS <- DESeq(ddsGroup3[rowSums(counts(ddsGroup3)) >= 10])
arachnoidDDS <- DESeq(ddsArachnoid[rowSums(counts(ddsArachnoid)) >= 10])

#Compare arachnoid to all tumors and each group individually
res.Tumors.Arachnoid <- results(arachnoidDDS, contrast=c('arachnoid','other','Arachnoid'))
res.Group1.Arachnoid <- results(groupDDS, contrast=c('groups','1','Arachnoid'))
res.Group2.Arachnoid <- results(groupDDS, contrast=c('groups','2','Arachnoid'))
res.Group3.Arachnoid <- results(groupDDS, contrast=c('groups','3','Arachnoid'))

#select
pval <- 1E-6
Wald <- 0
DEGs.Tumors.Arachnoid <- subset(res.Tumors.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)
DEGs.Group1.Arachnoid <- subset(res.Group1.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)
DEGs.Group2.Arachnoid <- subset(res.Group2.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)
DEGs.Group3.Arachnoid <- subset(res.Group3.Arachnoid, padj < pval & abs(log2FoldChange) > Wald)

dim(DEGs.Tumors.Arachnoid)
dim(DEGs.Group1.Arachnoid)
dim(DEGs.Group2.Arachnoid)
dim(DEGs.Group3.Arachnoid)

DEGs.Arachnoid <- rownames(DEGs.Tumors.Arachnoid)
DEGs.Group1 <- rownames(DEGs.Group1.Arachnoid)[!rownames(DEGs.Group1.Arachnoid) %in% DEGs.Arachnoid]
DEGs.Group2 <- rownames(DEGs.Group2.Arachnoid)[!rownames(DEGs.Group2.Arachnoid) %in% DEGs.Arachnoid]
DEGs.Group3 <- rownames(DEGs.Group3.Arachnoid)[!rownames(DEGs.Group3.Arachnoid) %in% DEGs.Arachnoid]

AllDEGs <- unique(c(DEGs.Group1, DEGs.Group2, DEGs.Group3, DEGs.Arachnoid))

uniqueGroup1 <- DEGs.Group1[!DEGs.Group1 %in% c(DEGs.Group2, DEGs.Group3)]
uniqueGroup2 <- DEGs.Group2[!DEGs.Group2 %in% c(DEGs.Group1, DEGs.Group3)]
uniqueGroup3 <- DEGs.Group3[!DEGs.Group3 %in% c(DEGs.Group1, DEGs.Group2)]

write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), AllDEGs),], 'All DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueGroup1),], 'Group 1 unique DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueGroup2),], 'Group 2 unique DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueGroup3),], 'Group 3 unique DEGs.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), DEGs.Arachnoid),], 'Arachnoid DEGs.csv')

length(grep('1p', subset(annotations, external_gene_name %in% uniqueGroup1 & chromosome_name == 1)$ChrLoc))/length(uniqueGroup1)
length(grep('1p', subset(annotations, external_gene_name %in% uniqueGroup2 & chromosome_name == 1)$ChrLoc))/length(uniqueGroup2)
length(grep('1p', subset(annotations, external_gene_name %in% uniqueGroup3 & chromosome_name == 1)$ChrLoc))/length(uniqueGroup3)


########################################3
#Compare genes with DMPs unqiue to a given group
DMPpheno <- colData[PLSGroups,c('pheno1', 'pheno2', 'pheno3')]

#Use champ.DMP to find statistically significant DMPs (pval < 1e-06, min difference in mean beta 0.1)
OverallGroup1 <- champ.DMP(beta = methData[,PLSGroups], pheno=as.character(DMPpheno$pheno1), arraytype = 'EPIC')
OverallGroup2 <- champ.DMP(beta = methData[,PLSGroups], pheno=as.character(DMPpheno$pheno2), arraytype = 'EPIC')
OverallGroup3 <- champ.DMP(beta = methData[,PLSGroups], pheno=as.character(DMPpheno$pheno3), arraytype = 'EPIC')

pval <- 1E-6
difBeta <- 0.1
OverallGroup1DMPs <- intersect(CPG_promoter, rownames(subset(OverallGroup1[[1]], adj.P.Val < pval & abs(deltaBeta) > difBeta)))
OverallGroup2DMPs <- intersect(CPG_promoter, rownames(subset(OverallGroup2[[1]], adj.P.Val < pval & abs(deltaBeta) > difBeta)))
OverallGroup3DMPs <- intersect(CPG_promoter, rownames(subset(OverallGroup3[[1]], adj.P.Val < pval & abs(deltaBeta) > difBeta)))

Group1DMPgenes <- unique(CPG[OverallGroup1DMPs,]$gene)
Group2DMPgenes <- unique(CPG[OverallGroup2DMPs,]$gene)
Group3DMPgenes <- unique(CPG[OverallGroup3DMPs,]$gene)

DMPgenes <- unique(c(Group1DMPgenes, Group2DMPgenes, Group3DMPgenes))
uniqueDMPgenesGroup1 <- Group1DMPgenes[!Group1DMPgenes %in% c(Group2DMPgenes, Group3DMPgenes)]
uniqueDMPgenesGroup2 <- Group2DMPgenes[!Group2DMPgenes %in% c(Group1DMPgenes, Group3DMPgenes)]
uniqueDMPgenesGroup3 <- Group3DMPgenes[!Group3DMPgenes %in% c(Group1DMPgenes, Group2DMPgenes)]

write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), DMPgenes),], 'All DMP genes.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueDMPgenesGroup1),], 'Group 1 DMP genes.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueDMPgenesGroup2),], 'Group 2 DMP genes.csv')
write.csv(AllGeneR2s[intersect(rownames(AllGeneR2s), uniqueDMPgenesGroup3),], 'Group 3 DMP genes.csv')


#Boxplot for a given gene
plot_gene <- 'PCLAF'
plot_gene %in% rownames(expression)
plot_data <- data.frame(row.names = c(colnames(methData), 'Arachnoid_1', 'Arachnoid_2', 'Arachnoid_3'), 
                     gene = expression[plot_gene, c(colnames(methData), 'Arachnoid_1', 'Arachnoid_2', 'Arachnoid_3')], 
                     MethGroup = c(SampleData[colnames(methData),]$MethGroup, rep("Arachnoid",3)),
                     RNAGroup = c(SampleData[colnames(methData),]$RNA_numerical, rep("Arachnoid",3)),
                     CNVGroup = c(SampleData[colnames(methData),]$CNVGroup, rep("Arachnoid",3)),
                     NF2InstabilityGroup = c(SampleData[colnames(methData),]$FourthGroup, rep("Arachnoid",3)),
                     MenG = c(SampleData[colnames(methData),]$CompositeGroup, rep("Arachnoid",3)),
                     WHO.grade = c(SampleData[colnames(methData),]$WHO.grade, rep("Arachnoid",3)))

plot_data <- subset(plot_data, MenG != 'Unknown')
plot_data$MethGroup <- factor(plot_data$MethGroup , levels=c("Arachnoid", "1", "2", "3"))
plot_data$RNAGroup <- factor(plot_data$RNAGroup , levels=c("Arachnoid", "1", "2", "3"))
plot_data$CNVGroup <- factor(plot_data$CNVGroup , levels=c("Arachnoid", "1", "2", "3"))
plot_data$NF2InstabilityGroup <- factor(plot_data$NF2InstabilityGroup , levels=c("Arachnoid", "1", "2", "3"))
plot_data$MenG <- factor(plot_data$MenG , levels=c("Arachnoid", "1", "2", "3"))
plot_data$WHO.grade <- factor(plot_data$WHO.grade , levels=c("Arachnoid", "WHO I", "WHO II"))

boxplot(gene~MenG,data=plot_data, main=paste0(plot_gene,  ' expression in Meningioma Groups'), xlab='MenG', ylab=paste0('Log ', plot_gene, ' expression'), col = c('grey', colors_by_group[1:3]))
boxplot(gene~drop.levels(MenG), data=subset(plot_data, MenG != 'Arachnoid'), main=paste0(plot_gene,  ' expression in Meningioma Groups'), xlab='MenG', ylab=paste0('Log ', plot_gene, ' expression'), col = colors_by_group[1:3])

