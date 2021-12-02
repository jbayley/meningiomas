#########################
#Load all used packages #
#########################

loadPackages <- function() {
  
  memory.limit(size=10E10)
  library(ComplexHeatmap)
  library(circlize)
  library(dendsort)
  library(randomForest)
  library(ChAMP)
  library(DESeq2)
  library(biomaRt)
  library(Rtsne)
  library(caret)
  library(plsr)
  library(pls)
  library(sesame)
  library(dendextend)
  library(cluster)
  library(tidyverse)
  library(gdata)
  library(survival)
  library(survminer)
  library(conumee)
  library(minfi)
  library(CopyNeutralIMA)
  library(NMF)
  library(RColorBrewer)
  library(cluster)
  library(ConsensusClusterPlus)
}

###############################################
#QUality check and normalize methylation data #
###############################################

normalizeMethylation <- function(myLoad) {
  
  #Get normalized methylation data
  myNorm <- champ.norm(arraytype="EPIC") 
  
  #Correct sample names
  colnames(myNorm) <- gsub(".*(\\#)", "", colnames(myNorm))
  
  #Sort myNorm by row variance
  myNorm <- myNorm[order(rowVars(myNorm), decreasing = TRUE),]
  
  return(myNorm)
}

##################################
#Import data via SeSAMe pipeline #
##################################

loadSesame <- function() {
  
  #Set directory where methylation raw data located
  idat_dir <- "C:/Users/james/Desktop/Patel_Methylation/All_Patel_Data/Methylation"
  
  #Load data via SeSAMe
  betas <- openSesame(idat_dir)
  
  #Load sample name translations and rename columns
  names <- read.csv("Methylation_to_Sample_Name.csv")
  colnames(betas) <- names[match(colnames(betas), names$Methylation_Name),'Sample_Name']
  
  #Remove any probes with a NA value
  for(val in 1:dim(betas)[2]) {
    betas <- betas[!is.na(betas[,val]),]
  }
  
  #Limit to CpG sites
  betas <- betas[grep('cg', rownames(betas)),]
  
  #Remove probes on sex chromosomes, first by getting data for all probes
  betas <- betas[!rownames(betas) %in% rownames(subset(CPG, CHR %in% c("X", "Y"))),]
  
  #Sort in decreasing order of variance
  betas <- betas[order(rowVars(betas), decreasing = TRUE),]
  
  return(betas)
}

##########################################
#Get variance-stabilized expression data #
##########################################

getExpressionData <- function(dds) {
  
  #Perform variance stabilizing transformation on dds to create expression dataframe
  vsd <- vst(dds, blind=TRUE)
  expression <- assay(vsd)
  ENSGlist <- rownames(expression)
  
  #Get gene names for each ensemble gene ID using getGeneList
  geneList <- getGeneList(ensemblList = ENSGlist)
  
  #Replace ensemble protein IDs with gene name in expression
  rownames(expression) <- geneList[match(ENSGlist, geneList[,"ensembl_gene_id"]),"external_gene_name"]
  
  #Loop through expression and replace missing genes (not found through BioMart) with the gene ID
  for(val in 1:length(ENSGlist)) {
    if(is.na(rownames(expression)[val])) {
      rownames(expression)[val] <- ENSGlist[val]
    }
  }
  
  #Remove any rows with blank rownames
  expression <- expression[rownames(expression) != "",]
  
  #Remove duplicated rownames
  expression <- expression[!duplicated(rownames(expression)),]
  
  return(expression)
}

########################################################
#Get gene and transcript ID information from biomaRt   #
########################################################

getGeneList <- function(ensemblList) {
  
  library(biomaRt)
  
  #Set the dataset from which to pull gene info
  ensembl = useDataset("hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
  
  #If the calls from the database have changed and need to check for updated names, use the following
  #library(tidyverse)
  #listMarts()
  #ensembl = useMart("ENSEMBL_MART_ENSEMBL")
  #list the available datasets (species)
  #listDatasets(ensembl) %>% filter(str_detect(description, "Human"))
  #ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
  #listFilters(ensembl) %>% filter(str_detect(name, "ensembl"))
  #check the available "attributes" - things you can retrieve
  #listAttributes(ensembl) %>%  head(200)
  
  #Ascertain if the ensemblList is made of transcript IDs or gene Ids
  if (substr(ensemblList[1], 1, 4) == "ENSG") {
    filterType <- "ensembl_gene_id"
  } else if (substr(ensemblList[1], 1, 4) == "ENST") {
    filterType <- "ensembl_transcript_id"
  } else {
    message("ensemblList input not valid gene or transcript IDs")
    break()
  }
  
  #Set the list of attributes
  attributeNames <- c(filterType, 'external_gene_name')
  
  #Run the query
  annotations <- getBM(attributes=attributeNames, filters = filterType, values = ensemblList, mart = ensembl)
  
  return(annotations)
}

########################################
#Get Olar et al. MMFav/Unfav CPG lists #
########################################

load_Olar_CPG <- function() {
  
  #Load Olar et al. data
  Olar_data <- read.csv(paste0(getwd(), "/Olar_CPGs.CSV"), fill = FALSE)
  
  #Create vectors for each type of CpG
  MMFav <- Olar_data$MMFav[1:98]
  MMUnFav <- Olar_data$MMUnFav[1:185]
  MMFav64 <- Olar_data$MMFav64[1:24]
  MMUnFav64 <- Olar_data$MMUnFav64[1:40]
  
  #Remove CpGs for which we do not have normalized data (MMFav)
  x <- (1:length(MMFav))
  for (val in x) {
    if(!MMFav[val] %in% rownames(methData)) {
      message(MMFav[val], " is not in our CpGs and was removed from MMFav")
      MMFav[val] <- NA
    }
  }
  MMFav <- MMFav[!is.na(MMFav)]
  
  #Remove CpGs for which we do not have normalized data (MMUnFav)
  x <- (1:length(MMUnFav))
  for (val in x) {
    if(!MMUnFav[val] %in% rownames(methData)) {
      message(MMUnFav[val], " is not in our CpGs and was removed from MMUnFav")
      MMUnFav[val] <- NA
    }
  }
  MMUnFav <- MMUnFav[!is.na(MMUnFav)]
  
  #Remove CpGs for which we do not have normalized data (MMFav64)
  x <- (1:length(MMFav64))
  for (val in x) {
    if(!MMFav64[val] %in% rownames(methData)) {
      message(MMFav64[val], " is not in our CpGs and was removed from MMFav64")
      MMFav64[val] <- NA
    }
  }
  MMFav64 <- MMFav64[!is.na(MMFav64)]
  
  #Remove CpGs for which we do not have normalized data (MMUnFav64)
  x <- (1:length(MMUnFav64))
  for (val in x) {
    if(!MMUnFav64[val] %in% rownames(methData)) {
      message(MMUnFav64[val], " is not in our CpGs and was removed from MMUnFav64")
      MMUnFav64[val] <- NA
    }
  }
  MMUnFav64 <- MMUnFav64[!is.na(MMUnFav64)]
  
  #Create a list of each type of CpG for output
  output <- vector(mode = "list", 6)
  names(output) <- c("MMFav", "MMUnFav",  "MMBoth", "MMFav64", "MMUnFav64", "MMBoth64")
  output[[1]] <- MMFav
  output[[2]] <- MMUnFav
  output[[3]] <- c(MMFav, MMUnFav)
  output[[4]] <- MMFav64
  output[[5]] <- MMUnFav64
  output[[6]] <- c(MMFav64, MMUnFav64)
  
  return(output)
}

###########################################
#Perform partial least squares regression #
###########################################

getPLS <- function(gene, methylation = methData, CPGs = CPG, expressionSubset = expression) {
  if(!gene%in%rownames(expressionSubset)) {
    return("No expression data")
  }
  cpgList <- rownames(CPGs[CPGs$gene == gene,])
  allSamples <- intersect(colnames(methylation), intersect(colnames(expressionSubset), colnames(CNV)))
  
  if(sum(expressionSubset[gene,allSamples]==min(expressionSubset[gene,allSamples])) > length(allSamples)*0.8) {
    return("Vast majority of samples with minimum expression")
  }
  
  if(gene %in% rownames(CNV)) {
    cnv <- CNV[gene, allSamples]
    cnv_status <- "CNV data available"
  } else {
    cnv <- CNV[1, allSamples]
    cnv[1,] <- 0
    cnv_status <-"CNV data NOT available"
  }
  if(length(cpgList)==0) {
    return("No CpGs")
  }
  
  if(length(cpgList) == 1) {
    data <- cbind(cnv = t(cnv), methylation[cpgList, allSamples], expression = expressionSubset[gene, allSamples])
    colnames(data)[1] <- "CNV"
    colnames(data)[2] <- cpgList
  } else {
    #Assemble data of methylation for each CpG and expression for the gene
    data <- cbind(cnv = t(cnv), t(methylation[cpgList, allSamples]),expression = expressionSubset[gene, allSamples])
    colnames(data)[1] <- "CNV"
  } 
  
  #Run partial least squares model and return CpGs by importance and the model
  tryCatch(plsModel <- train(expression ~ ., data = data, method = "pls", scale = TRUE), error = function(e) {})
  if(!exists('plsModel')) {
    message(paste0("Model did not work for ", gene))
    return("PLS model failed due to low variance in data")
  }
  imp <- varImp(plsModel)$importance
  imp <- cbind(imp, CPGs[rownames(imp),])
  imp <- imp[order(-imp$Overall),]
  
  if(cnv_status=="CNV data available") {
    cnv_status <- varImp(plsModel)$importance['CNV','Overall']
  }
  
  justCPGs <- imp[rownames(imp) != "CNV",]
  TopCPG <- rownames(justCPGs)[1]
  TopCPGImp <- justCPGs[1,'Overall']
  if(dim(justCPGs)[1] > 1) {
    SecondCPG <- rownames(justCPGs)[2]
    SecondCPGImp <- justCPGs[2,'Overall']
  } else {
    SecondCPG <- "Less than 2 CpGs"
    SecondCPGImp <- "Less than 2 CpGs"
  }
  if(dim(justCPGs)[1] > 2) {
    ThirdCPG <- rownames(justCPGs)[3]
    ThirdCPGImp <- justCPGs[3,'Overall']
  } else {
    ThirdCPG <- "Less than 3 CpGs"
    ThirdCPGImp <- "Less than 3 CpGs"
  }
  
  return(list(imp, plsModel, cnv_status, dim(justCPGs)[1], TopCPG, TopCPGImp, SecondCPG, SecondCPGImp, ThirdCPG, ThirdCPGImp))
}

###########################################################################
#Function to correlate CNV or a CpG site methylation with gene expression #
###########################################################################

correlateExpression <- function(data) {
  
  linearmodel <- lm(expression ~ ., data = data)
  
  output <- data.frame(row.names = colnames(data)[2])
  output$coef <- summary(linearmodel)$coef[2, 1]
  output$coefPval <- summary(linearmodel)$coef[2, 4]
  output$adjRsq <- summary(linearmodel)$adj.r.squared
  
  return(list(linearmodel, output))
}