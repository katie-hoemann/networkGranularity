
##Load packages 
library("foreign")
library("qgraph")
library("bootnet")
library("dplyr")
library("igraph")
library("R.matlab")


#subj_IDs <- c("PP002", "PP003", "PP006", "PP007")
#subj_IDs <- c("PP2", "PP3", "PP6", "PP7", "PP9", "PP10", "PP11", "PP12", "PP13", "PP14", "PP15", "PP17", "PP18", "PP19", "PP21", "PP22", "PP23", "PP24", "PP25", "PP26", "PP27", "PP28", "PP29", "PP32", "PP33", "PP34", "PP36", "PP37", "PP38", "PP40", "PP43", "PP45", "PP46", "PP47", "PP49", "PP50", "PP51", "PP54", "PP55", "PP56", "PP57", "PP58", "PP59", "PP61", "PP62", "PP63", "PP64", "PP66", "PP67", "PP68")
subj_IDs <- c("PP1", "PP2", "PP3", "PP5", "PP6", "PP7", "PP8", "PP9", "PP10", "PP11", "PP12", "PP13", "PP14", "PP15", "PP16", "PP17", "PP18", "PP19", "PP20", "PP21", "PP22", "PP23", "PP24", "PP25", "PP26", "PP27", "PP28", "PP30", "PP31", "PP32", "PP34", "PP35", "PP36", "PP37", "PP38", "PP39", "PP40", "PP41", "PP42", "PP43", "PP44", "PP45", "PP46", "PP47", "PP48", "PP49", "PP50", "PP51", "PP52", "PP53", "PP54", "PP55", "PP56", "PP58", "PP59", "PP61", "PP62", "PP63", "PP64", "PP66", "PP67", "PP69", "PP70","PP71", "PP72", "PP73", "PP74", "PP75", "PP76", "PP77", "PP79", "PP80", "PP81", "PP82", "PP83", "PP84")

lambda_result_EBIC0 <- vector("numeric")
lambda_result_EBIC.25 <- vector("numeric")
lambda_result_EBIC.5 <- vector("numeric")

for (i in subj_IDs) {
  filename <- paste("/Users/Katie/Documents/R/Input/",i,"_data.csv",sep="")
  PDFname1 <- paste("/Users/Katie/Documents/R/Output/PDFs/",i,"_VariableLambda.pdf",sep="")
  PDFname2 <- paste("/Users/Katie/Documents/R/Output/PDFs/",i,"_BestNetwork_EBIC0.pdf",sep="")
  PDFname3 <- paste("/Users/Katie/Documents/R/Output/PDFs/",i,"_BestNetwork_EBIC.25.pdf",sep="")
  PDFname4 <- paste("/Users/Katie/Documents/R/Output/PDFs/",i,"_BestNetwork_EBIC.5.pdf",sep="")
  Edgelist1 <- paste("/Users/Katie/Documents/R/Output/Matrices/",i,"_Edgelist_EBIC0.mat",sep="")
  Edgelist2 <- paste("/Users/Katie/Documents/R/Output/Matrices/",i,"_Edgelist_EBIC.25.mat",sep="")
  Edgelist3 <- paste("/Users/Katie/Documents/R/Output/Matrices/",i,"_Edgelist_EBIC.5.mat",sep="")
  LambdaResults <- paste("/Users/Katie/Documents/R/Output/LambdaResults.mat",sep="")
  
  # Read in data file (change path to location of the file on your computer)
  data_full <- read.csv(filename)
  
  # Create data frame of variables that we need for analysis: 
  data <- as.data.frame(data_full[,c(1:length(data_full))]) # data for analysis
  
  # Correlation matrix:
  corMat <- cor(data)
  
  #To calculate polychoric correlations:
  #library(polycor)
  #corMat <- polychor(data)
  
  # Labels:
  #Labs <- c("Afraid", "Amused", "Angry", "Bored", "Calm", "Disgusted", "Embarassed", "Excited", "Frustrated", "Grateful", "Happy", "Neutral", "Proud", "Relieved", "Sad", "Serene", "Surprised", "WornOut")
  Labs <- colnames(data)
  
  # EBIC analysis:
  res0 <- EBICglasso(corMat,nrow(data),0,returnAllResults=TRUE,nlambda=100)
  lambda0 <- res0$lambda[which.min(res0$ebic)] #identify lambda that corresponds to the minimum EBIC
  lambda_result_EBIC0[i] <- lambda0
  
  res0.25 <- EBICglasso(corMat,nrow(data),0.25,returnAllResults=TRUE,nlambda=100)
  lambda.25 <- res0.25$lambda[which.min(res0.25$ebic)] #identify lambda that corresponds to the minimum EBIC
  lambda_result_EBIC.25[i] <- lambda.25
  
  res0.5 <- EBICglasso(corMat,nrow(data),0.5,returnAllResults=TRUE,nlambda=100)
  lambda.5 <- res0$lambda[which.min(res0$ebic)] #identify lambda that corresponds to the minimum EBIC
  lambda_result_EBIC.5[i] <- lambda.5
  
  #Estimated networks:
  #pdf(PDFname1,width=5*2,height=2*2)  ##saves a PDF with figure according to specified directory/name 'PDFname'
  #layout(rbind(1:5,6:10))
  #for (i in 1:10){
  #  qgraph(qgraph::wi2net(res0$results$wi[,,i]), labels = Labs, vsize = 15,mar=c(2,2,2,2), theme = "colorblind")
  #  text(0,0.3,paste0("Lambda: ",round(res0$lambda[i],3)),cex=1)
  #  text(0,0.1,paste0("EBIC (gamma = 0): ",round(res0$ebic[i],1)),cex=ifelse(res0$ebic[i] == min(res0$ebic),1,0.6),font = ifelse(res0$ebic[i] == min(res0$ebic),2,1))
  #  text(0,-0.1,paste0("EBIC (gamma = 0.25): ",round(res0.25$ebic[i],1)),cex=ifelse(res0.25$ebic[i] == min(res0.25$ebic),1,0.6),font = ifelse(res0.25$ebic[i] == min(res0.25$ebic),2,1))
  #  text(0,-0.3,paste0("EBIC (gamma = 0.5): ",round(res0.5$ebic[i],1)),cex=ifelse(res0.5$ebic[i] == min(res0.5$ebic),1,0.6),font = ifelse(res0.5$ebic[i] == min(res0.5$ebic),2,1))
  #}
  #dev.off() 

##TO PLOT JUST ONE NETWORK https://cran.r-project.org/web/packages/qgraph/qgraph.pdf
  # Compute graph with tuning = 0 (EBIC)
  EBICgraph1 <- EBICglasso(corMat, nrow(data), 0, nlambda = 100, threshold = TRUE)
  writeMat(Edgelist1, subjectMatrix=EBICgraph1)
  # Plot:
  #pdf(PDFname2,width=2*2,height=2*2)
  #layout(1)
  #EBICgraph1 <- qgraph(EBICgraph1,  title = "EBIC", theme="colorblind", details=TRUE, labels=Labs)
  #dev.off()

  # Compute graph with tuning = 0.25 (EBIC)
  EBICgraph2 <- EBICglasso(corMat, nrow(data), 0.25, nlambda = 100, threshold = TRUE)
  writeMat(Edgelist2, subjectMatrix=EBICgraph2)
  # Plot:
  #pdf(PDFname3,width=2*2,height=2*2)
  #layout(1)
  #EBICgraph2 <- qgraph(EBICgraph2,  title = "EBIC", theme="colorblind", details=TRUE, labels=Labs)
  #dev.off()

  # Compute graph with tuning = 0.5 (EBIC)
  EBICgraph3 <- EBICglasso(corMat, nrow(data), 0.5, nlambda = 100, threshold = TRUE)
  writeMat(Edgelist3, subjectMatrix=EBICgraph3)
  # Plot:
  #pdf(PDFname4,width=2*2,height=2*2)
  #layout(1)
  #EBICgraph3 <- qgraph(EBICgraph3,  title = "EBIC", theme="colorblind", details=TRUE, labels=Labs)
  #dev.off()


}


writeMat(LambdaResults, lambda_EBIC0=cbind(subj_IDs, lambda_result_EBIC0), lambda_EBIC.25=cbind(subj_IDs,lambda_result_EBIC.25), lambda_EBIC.5=cbind(subj_IDs,lambda_result_EBIC.5))

