####
#### Bayes Model for data integration
#### Code by: Chuck Song 10/14/2014

#########################################
#### xiaoguang, 2014/10/18, GSPH 309D
#### ploting trace

#### Part I: MCMC Inference

#### Step 1. Update Y
setwd('/home07/xiaoguang/AW/CI/result20160129/code')
source('simu_fun.r')
source('packageAW.R')
source('/home07/xiaoguang/AW/WangLin/20160120_v3/AW.R')
workingDirectory <- '/home07/xiaoguang/AW/CI/result20160129/varibility'
system(paste('mkdir -p',workingDirectory))
setwd(workingDirectory)
library(MCMCpack)
library(truncnorm)
library(limma)


args = commandArgs(trailingOnly = TRUE)
## args: N, K, sigma1, seed
args_num <- as.numeric(args)

N <- args_num[1]
K <- args_num[2]
sigma1 <- args_num[3]
sigma2 <- args_num[4]
seed <- args_num[5]

G = 10000
B = 1000
pvalueWeight1Cut <- 1e-7

if(F){
	N=20
	K=4
	sigma1=1
	seed=15213
	sigma2=0.2
}


folder <- paste('AW_CI_G_',G,'_N_',N,'_K_',K,'_B_',B,'_sigma1_',sigma1,'_sigma2_',sigma2,sep='')
#folder <- gsub('[.]','',folder)
system(paste('mkdir -p',folder))

fileName0 <- paste('AW_CI_G_',G,'_N_',N,'_K_',K,'_B_',B,'_sigma1_',sigma1,'_sigma2_',sigma2,'_seed_',seed,'.rdata',sep='')
fileName <- paste(folder,fileName0,sep='/')

## generate data
generateS = simu.data.randomeffect(K = K, N = N, G = G, opdirectionRate = 0.01, g=G*0.3,sigma=sigma1, sigmaRandom = sigma2)
dataExp = generateS$data
truth = generateS$truth
sum(truth!=0)  # 7512	
truthRes = truth
sum(truthRes!=0) ##  7512
sum(apply(truthRes,1,function(x) sum(abs(x)))>0) ## 3000 


## get p values
ControlLabel = 1:N
caseLabel = (N+1):(2*N)
labels <- lapply(1:length(dataExp),function(x) list(ControlLabel,caseLabel))
myLimmaP<-function(dataExp,labels,seed=15213){
	set.seed(seed)
	
	
	effectSize <- NULL
	pvalues <- NULL
	
	for(i in 1:length(dataExp)){
		aData = dataExp[[i]]
		alabel = labels[[i]]
		labelData <- numeric(sum(sapply(alabel,length)))
		labelData[alabel[[1]]] = 0
		labelData[alabel[[2]]] = 1
		design = cbind(rep(1,length(labelData)),labelData)
		fit <- lmFit(aData,design)
		fit <- eBayes(fit)		
		pvalues <- cbind(pvalues,fit$p.value[,2])
		effectSize <- cbind(effectSize, fit$coefficients[,2])					
	}	
	result <- list(pvalues=pvalues, effectSize=effectSize)
}


p.matrix <- myLimmaP(dataExp,labels)$pvalues

res <- aw.fisher.pvalueC(p.matrix, method="original", log=F, weight.matrix=T, core=1)


AWqvalue <- p.adjust(res$pvalues,'BH')
DEindex <- AWqvalue < 0.05
sum(DEindex) ## 2736

sum(!apply(truthRes,1,function(x) sum(abs(x)))>0 & DEindex) / sum(DEindex) ## FDR: 0.03508772

permuteX <- function(dataExp,labels,seed=15213){
	set.seed(seed)
    permx <- list()
	for(i in 1:length(dataExp)){
		X <- matrix(,nrow=nrow(dataExp[[i]]),ncol(dataExp[[i]]))
		aData = dataExp[[i]]
		alabel = labels[[i]]

		for(j in 1:length(alabel)){
			aalabel <- alabel[[j]]
			for(k in aalabel){
				X[,k] <- aData[,sample(aalabel,1)]
			}
		}
		permx[[i]] = X
	}	
	permx	
}


Tscore <- 0
res_b <- NULL
selfDistDirection <- matrix(0,nrow=sum(DEindex),ncol = sum(DEindex))

for(b in 1:B){
	print(b)
	dataExp_b <- permuteX(dataExp,labels,seed=b)
	limmaResult <- myLimmaP(dataExp_b,labels)
	p.matrix_b <- limmaResult$pvalues
	effectsize_b <- limmaResult$effectSize
	
	bAWres = aw.fisher.pvalueC(p.matrix_b, method="original", log=F, weight.matrix=T, core=1)
	res_b[[b]] <- bAWres$weights
	res_b[[b]][p.matrix_b<pvalueWeight1Cut] <- 1			
}

aveScore <- Reduce('+',res_b)/B
for(b in 1:B){
	Tscore = Tscore + (res_b[[b]] - aveScore)^2/B
}
print(Tscore[,1])
mean(Tscore[generateS$effectSize!=0])

res <- list(Tscore=abs(Tscore),effectSize=generateS$effectSize,truth=generateS$truth)
save(res,file=fileName)


## 1, why sum of weight are 0?
## 2, why U = 1
