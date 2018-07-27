library(AWFisher)


for(k in 2:1000){
	G <- 10000
	p.values = matrix(rbeta(K*G, 1,1), ncol=K)
	res = AWFisher.pvalue(p.values)
	ks<-ks.test(res$pvalues, "punif", min=0, max=1, alternative = "two.sided");
	ks
	cat("k:", k, ". ", ks$p.value, "\n")
	
}

