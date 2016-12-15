biomarkerCategorization <- function(studies,afunction,B=10,DEindex=NULL,seed = 15213){
	if(is.null(DEindex)){
		cat('generate DE index since it is NULL','\n')
		res.obs <- getPvalueAll(studies, afunction)
		pval.obs <- res.obs$p.matrix
		awres.obs <- AWFisher.pvalue(pval.obs)
		fdr.obs <- p.adjust(awres.obs$pvalues, 'BH')
		DEindex <- fdr.obs <= 0.05
	}
		
	es.null <- NULL
	pval.null <- NULL
	weight.null <- NULL
	selfDistDirection <- matrix(0,nrow=sum(DEindex),ncol = sum(DEindex))
	
	cat('permuting analysis:')
	for(b in 1:B){
		cat('.')
		studies.b <- permuteLabels(studies, seed=seed+b)
		res.null <- getPvalueAll(studies.b, afunction)
		pval.null <- res.null$p.matrix
		es.null <- res.null$es.matrix			
		
		awres.null <- AWFisher.pvalue(pval.null)		
		weight.null[[b]] <- awres.null$weights
					
		aDEweightDirection <- (awres.null$weights * sign(es.null))[DEindex,]
		maxDistDirection <- dist(aDEweightDirection,'maximum')
		bselfDistDirection <- ifelse(as.matrix(maxDistDirection) == 0, 1, 0)	
		selfDistDirection <- selfDistDirection + as.matrix(bselfDistDirection)/B			
	} ## b for end of B
	cat('\n', 'calculating variability index','\n')
	
	Tscore <- 0
	aveScore <- Reduce('+',weight.null)/B
	for(b in 1:B){
		Tscore = Tscore + (weight.null[[b]] - aveScore)^2/B
	}

	result <- list(varibility=Tscore, dissimilarity=selfDistDirection, DEindex=DEindex)
	result
	
}

tmp22 <- biomarkerCategorization(studies,function_limma,B=50,DEindex=NULL,seed = 15213)

library(mclust)

distr <- as.dist(tmp22$dissimilarity)
hcr <- hclust(distr, method='ward.D2')
hcr <- hclust(distr, method='average')
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
acut <- cutreeDynamicTree(hcr, maxTreeHeight = 1, deepSplit = TRUE, minModuleSize = 50)
table(acut)
table(distTruth, acut)
adjustedRandIndex(distTruth, acut)


setwd('/home07/xiaoguang/AW/CI/result20160129/code')
source('simu_fun.r')
source('packageAW.R')
source('/home07/xiaoguang/AW/WangLin/20160120_v3/AW.R')
library(limma)
library(dendextend)



workDirectory <- '/home07/xiaoguang/AW/CI/result20160129/biomarkerCategorization'
system(paste('mkdir -p', workDirectory))
setwd(workDirectory)

N0 = 50
K = 4

#### generating the data
set.seed(15213)
generateS = simu.data.partial(K = K, N = N0, G = 10000,typeselection=c(1,2),firstK=2,comb = c(0.04,0.04,0.02,0.2),sigma=1, sigmaRandom = 0.2)
dataExp = generateS$data
truth = generateS$truth
DEtype = generateS$DEtype
sum(truth!=0)  ##	2000

N=N0
ControlLabel = 1:N
caseLabel = (N+1):(2*N)


effectSize <- NULL
p.matrix <- NULL

labels <- NULL

set.seed(15213)
for(i in 1:length(dataExp)){
	print(i)
	aData = dataExp[[i]]
	alabel = c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))
	labels[[i]] <- list(ControlLabel,caseLabel)
	design = cbind(rep(1,length(alabel)),alabel)
	fit <- lmFit(aData,design)
	fit <- eBayes(fit)		

	effectSize <- cbind(effectSize, fit$coefficients[,2])
	p.matrix <- cbind(p.matrix, fit$p.value[,2])
}

effectSize_direction <- sign(effectSize)


AWres = aw.fisher.pvalueC(p.matrix, method="original", weight.matrix=T, core=1)

AWqvalue <- p.adjust(AWres$pvalues,'BH')
DEindex <- AWqvalue < 0.05
sum(DEindex) ## 678


permuteLabels <- function(studies,seed=15213){
	set.seed(seed)
	for(k in 1:length(studies)){
		studies[[k]]$label <- sample(studies[[k]]$label)
	}
	studies
}


Tscore <- 0
res_b <- NULL
selfDistDirection <- matrix(0,nrow=sum(DEindex),ncol = sum(DEindex))

B = 1000


aData = dataExp[[i]]
alabel = c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))

study1 <- list(data=dataExp[[1]], label=alabel)
study2 <- list(data=dataExp[[2]], label=alabel)
study3 <- list(data=dataExp[[3]], label=alabel)

studies <- list(study1,study2,study3)
names(studies) <- c('study1','study2','study3')


function_limma <- function(astudy){
	alabel <- astudy$label
	aData <- astudy$data
	
	adesign = cbind(rep(1,length(alabel)),alabel)
	afit <- lmFit(aData,adesign)
	afit <- eBayes(afit)		

	aeffectsize <- afit$coefficients[,2]
	apvalue <- afit$p.value[,2]
	return(list(pvalue=apvalue, effectSize=aeffectsize))
}

getPvalueSingle <- function(astudy, afunction){
	## input a study and label as a list
	## return p value and effect size
	return(afunction(astudy))
}

tmp <- getPvalueSingle(study1, function_limma)

getPvalueAll <- function(studies, afunction){
	K <- length(studies)
	
	p.matrix <- NULL
	es.matrix <- NULL
	
	for(k in 1:K){
		astudy <- studies[[k]]
		ares <- afunction(astudy)
		p.matrix <- cbind(p.matrix, ares$pvalue)
		es.matrix <- cbind(es.matrix, ares$effectSize)		
	}
	return(list(p.matrix=p.matrix, es.matrix=es.matrix))
}

tmp2 <- getPvalueAll(studies, function_limma)


for(b in 1:B){
	print(b)
	dataExp_b <- permuteX(dataExp,labels,seed=b)
	
	p.matrix_b <- NULL
	effectsize_b <- NULL
	
	for(i in 1:length(dataExp_b)){
		print(i)
		bData = dataExp_b[[i]]
		blabel = c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))
		bdesign = cbind(rep(1,length(blabel)),blabel)
		bfit <- lmFit(bData,bdesign)
		bfit <- eBayes(bfit)		

		effectsize_b <- cbind(effectsize_b, bfit$coefficients[,2])
		p.matrix_b <- cbind(p.matrix_b, bfit$p.value[,2])
	}
	
	
	bAWres = aw.fisher.pvalueC(p.matrix_b, method="original", weight.matrix=T, core=1)
	res_b[[b]] <- bAWres$weights
			
	aDEweightDirection <- (bAWres$weight * sign(effectsize_b))[DEindex,]
	maxDistDirection <- dist(aDEweightDirection,'maximum')
	bselfDistDirection <- ifelse(as.matrix(maxDistDirection) == 0, 1, 0)
	
	selfDistDirection <- selfDistDirection + as.matrix(bselfDistDirection)/B	
}
range(selfDistDirection)

aveScore <- Reduce('+',res_b)/B
for(b in 1:B){
	Tscore = Tscore + (res_b[[b]] - aveScore)^2/B
}

CI <- list(varibility=Tscore, selfDistDirection=selfDistDirection)

save(CI,file='CI.rdata')
if(F){
	load('CI.rdata')
}



## ward distance works perfectly
library(gplots)
png('visualizatinoSelfDistDirectionH.png')
heatmap2Obj <- heatmap.2(CI$selfDistDirection,col=greenred,trace="none",hclustfun = function(x) hclust(x,method = 'ward'),dendrogram='column',
Rowv=TRUE,symm = TRUE, Colv=TRUE,keysize=1.3,labRow='',labCol='')
dev.off()

## try ward distance
k = 7
png('CIDistDirectionHclustTree.png',1000,600)
dend1 <- heatmap2Obj$rowDendrogram
if(F){
	dend1 <- color_branches(dend1, k = 7)
	dend1 <- color_labels(dend1, k = 7)	
}
plot(dend1,main='rowDendrogram from heatmap')
dev.off()

rainbow_fun <- function(n, c=90, l=50, ...) {
if(requireNamespace("colorspace")) {
colorspace::rainbow_hcl(n, c = c, l = l, ...)		
	} else {
		rainbow(n, ...)
	}
}



reorderLabel <- function(alabel,aorder){
  uniLab = unique(alabel)
  resLabel = NULL
  for(i in 1:length(uniLab)){
    resLabel[uniLab[i]==alabel]=aorder[i]
  }
  return(resLabel)
}

clusterMembership0 <- dendextend::cutree(dend1,k=7) # it now works on a dendrogram
# ClusterColors <- rainbow_fun(max(clusterMembership))[colorClusters]
clusterOrder <- unique(clusterMembership0[labels(dend1)])
matchOrder <- match(unique(clusterMembership0),clusterOrder)


reorderLabel <- function(alabel,aorder){
  uniLab = unique(alabel)
  resLabel = NULL
  for(i in 1:length(uniLab)){
    resLabel[uniLab[i]==alabel]=aorder[i]
  }
  return(resLabel)
}

clusterMembership <- reorderLabel(clusterMembership0,matchOrder)

table(clusterMembership0)
table(clusterMembership)

if(F){

	(AWres$weights * effectSize_direction)[DEindex,][clusterMembership0==1,]
	(AWres$weights * effectSize_direction)[DEindex,][clusterMembership0==2,]
	(AWres$weights * effectSize_direction)[DEindex,][clusterMembership0==3,]
	(AWres$weights * effectSize_direction)[DEindex,][clusterMembership0==4,]
	(AWres$weights * effectSize_direction)[DEindex,][clusterMembership0==5,]
	(AWres$weights * effectSize_direction)[DEindex,][clusterMembership0==6,]
	(AWres$weights * effectSize_direction)[DEindex,][clusterMembership0==7,]	
}

colorBar <- NULL
membershipTable <- table(clusterMembership)
allColor <- rainbow_fun(length(membershipTable))
allColor[6] <- '#777777'


library(gplots)
my_palette <- colorRampPalette(c("yellow2", "black", "skyblue"))(n = 300)
png('visualizatinoSelfDistDirection.png',501,500)
heatmap2Obj <- heatmap.2(CI$selfDistDirection,col=my_palette,trace="none",hclustfun = function(x) hclust(x,method = 'ward'),dendrogram='none',ColSideColors=allColor[clusterMembership]
,Rowv=TRUE,symm = TRUE, Colv=TRUE,keysize=1.3,labRow='',labCol='')
dev.off()


## heatmap for each cluster, genes are ordered by variability.
getWsHeatmap <- function(astudy,aCs,ws=NULL, Rowv=NA,geneOrder=NULL,...){
   if(is.null(ws))
     ws=rep(1,nrow(astudy))
   astudy = astudy[ws!=0,]
   coll <- NULL
   resHeatmap <- NULL
   label <- NULL
   uniCs = unique(aCs)
   uniorderCs = sort(uniCs)
   
   for(i in uniorderCs){  
     if(is.null(resHeatmap)){
       resHeatmap = astudy[,i==aCs]
       label = rep(i,sum(i==aCs))
     } 
     else {
       resHeatmap = cbind(resHeatmap,astudy[,i==aCs])
       label = c(label,rep(i,sum(i==aCs)))
     }
   }
   for(alabel in unique(label))
     coll[which(label==alabel)]=palette()[alabel]    
   finalRes = t(scale(t(resHeatmap)))
   if(!is.null(geneOrder)){
     finalRes = finalRes[geneOrder,]
   }
   outIndex = abs(finalRes)>3 
   finalRes[outIndex] = 3*sign(finalRes)[outIndex]
   B=16
   #return(gplots::heatmap.2(finalRes, col="greenred" ,trace="none",Rowv=Rowv,ColSideColors=coll,
   #                         Colv=NA,keysize=1.3,...)  )
   return(heatmap(finalRes,Rowv=Rowv,ColSideColors=coll,
 	  col= rgb(c(rep(0, B), (0:B)/B), c((B:0)/16, rep(0, B)), rep(0, 2*B+1))
   ,scale='none',Colv=NA,...)  )
 }


#

CI$varibility[DEindex,]
rowSums(CI$varibility[DEindex,])/3

sumVaryCut = 0.25
colorFade = 1/2

DEdata <- lapply(dataExp,function(x) x[DEindex,])
labels
memberShips <- unique(clusterMembership)
table(clusterMembership)

allData <- NULL

groups <- c(7,1,4,5,2,3)
table(clusterMembership)[groups]

for(i in groups){
	cat('membership: ')
	cat(i)
	cat('\n')
	
	## gene order based on variability score
	geneOrder <- order(rowSums(CI$varibility[DEindex,])[clusterMembership==i])
	acolor <- rainbow_fun(max(clusterMembership))[i]
	
	rowData <- NULL
	rowLabel <- NULL
	
	## visualize the heatmap
	for(j in 1:length(DEdata)){
		adata <- DEdata[[j]][clusterMembership==i,][geneOrder,]
		bdata <- t(scale(t(adata)))
		alabel <- numeric(length(Reduce(c,labels[[j]])))
		alabel[labels[[j]][[1]]] = 1
		alabel[labels[[j]][[2]]] = 2
		rowData <- cbind(rowData,bdata,NA)
		rowLabel <- c(rowLabel, alabel, NA)		
	}
	allData <- rbind(allData,rowData, matrix(NA,nrow=10,ncol=ncol(rowData)))	
}

colColor <- character(length(rowLabel))
colColor[is.na(rowLabel)] <- 'white'
colColor[rowLabel==1] <- 'black'
colColor[rowLabel==2] <- 'orange'

png('metaPatternVisualizationWithBar.png',width = 800, height = 1200)
heatmap.2(allData,ColSideColors = colColor,col="greenred", trace="none",Rowv=NA,Colv=NA,key=FALSE, labRow='',labCol='',na.color=par("bg"),
lwid=c(0.1,4), lhei=c(0.1,15),margins = c(1, 1))						
dev.off()

png('metaPatternVisualization.png',width = 800, height = 1200)
heatmap.2(allData,col="greenred", trace="none",Rowv=NA,Colv=NA,key=FALSE, labRow='',labCol='',na.color=par("bg"),
lwid=c(0.1,4), lhei=c(0.1,15),lwd=1,	margins = c(1, 1))						
dev.off()

rowDataWeight <- NULL
for(i in groups){
	## visualize the weight with direction
	geneOrder <- order(rowSums(CI$varibility[DEindex,])[clusterMembership==i])
	acolor <- rainbow_fun(max(clusterMembership))[i]
	
	aWeight <- (AWres$weights * effectSize_direction)[DEindex,][clusterMembership==i,][geneOrder,]
	rowDataWeight <- rbind(rowDataWeight,aWeight,matrix(NA,nrow=10,ncol=ncol(aWeight)))	
}


png('AWweightVisualization.png',width = 500, height = 1200)
heatmap.2(rowDataWeight,col=my_palette, trace="none",Rowv=NA,Colv=NA,key=FALSE, labRow='',labCol='',na.color=par("bg"),
lwid=c(0.1,4), lhei=c(0.1,15),margins = c(1, 1))						
dev.off()


rowDataVariability <- NULL
for(i in groups){
	## visualize the weight with direction
	geneOrder <- order(rowSums(CI$varibility[DEindex,])[clusterMembership==i])
	acolor <- rainbow_fun(max(clusterMembership))[i]

	## varibility score
	avaribility <- CI$varibility[DEindex,][clusterMembership==i,][geneOrder,]
	## sum of varibility score
	aRowVaribility <- rowSums(avaribility)/3	

	aRowVariability <- cbind(avaribility,aRowVaribility)
	rowDataVariability <- rbind(rowDataVariability,aRowVariability,matrix(NA,nrow=10,ncol=ncol(aRowVariability)))	
	
}

table(clusterMembership)

B = 100
png('VariabilityVisualization.png',width = 500, height = 1200)
heatmap.2(rowDataVariability,col=gray.colors(B), trace="none",Rowv=NA,Colv=NA,key=FALSE, labRow='',labCol='',na.color=par("bg"),
lwid=c(0.1,4), lhei=c(0.1,15),margins = c(1, 1))						
dev.off()
 
distTruth = DEtype[DEindex]
distTruth[distTruth==''] = 'nonDE'
distTruth[truth[DEindex,1]==1] = 'ssp1+'
distTruth[truth[DEindex,1]==-1] = 'ssp1-'
distTruth[truth[DEindex,2]==1] = 'ssp2+'
distTruth[truth[DEindex,2]==-1] = 'ssp2-'
distTruth[rowSums(truth[DEindex,])==K] = 'homo+'
distTruth[rowSums(truth[DEindex,])==-K] = 'homo-'


consistencyTable <- table(clusterMembership, distTruth)

write.csv(consistencyTable,'consistencyTableOri.csv',quote = FALSE)

## remove features with high variability index




