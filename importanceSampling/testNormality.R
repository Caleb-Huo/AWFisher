setwd("~/Desktop/AWFisher")

devtools::install()

library(AWFisher)


for(k in 2:1000){
	G <- 10000
	p.values = matrix(rbeta(K*G, 1,1), ncol=K)
	res = AWFisher.pvalue(p.values)
	ks<-ks.test(res$pvalues, "punif", min=0, max=1, alternative = "two.sided");
	ks
	cat("k:", k, ". ", ks$p.value, "\n")
	
}

Rep <- 1000
L <- 20#bad numbers: 20, 90
for(L in 2:1000){

#p.values = matrix(rbeta(Rep*L, 1,1), ncol=L)
p.values = matrix(runif(Rep*L), ncol=L)
ks0<-ks.test(p.values, "punif", min=0, max=1, alternative = "two.sided");

res = AWFisher.pvalue(p.values)
ks<-ks.test(res$pvalues, "punif", min=0, max=1, alternative = "two.sided");

prange=range(res$pvalues)

if(ks[[2]]<0.05 || prange[2]>1)
{
  print(paste0("Bad L = ", L))
  print(paste0("ks.p= ", ks[[2]]))
  print(prange)
}
}
