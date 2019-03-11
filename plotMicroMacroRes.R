tbl = read.table("/Users/danielesilvestro/Desktop/ClownSimulator/abc_windowv_ALL.log",h=T)
selected = which(tbl$Median_nDNA_in_mtDNAspB < 0.2)      # median nDNA in individuals with spB mtDNA: > 80% spA
selected = intersect(selected, which(tbl$mtDNA>0.01))    # at least 1 mtDNA of spB
selected = intersect(selected, which(tbl$nDNA_80spB==0)) # no individuals  with >80% spB nDNA
#selected = intersect(selected, which(tbl$max_fraction_of_larvae_spB<0.2)) # only small fraction of admixture


post = tbl[selected,]
print( c("Acceptance frequency:", length(selected)/dim(tbl)[1]* mean(tbl$accepted)) )

get_prior_range <- function(min_v=0,max_v=1,n_samples=1000,bins=10){
	b = seq(min_v,max_v,length.out=bins)
	m=0
	M=0
	a=0
	for (i in 1:10000){
		x = runif(n_samples,min_v,max_v)
		h = hist(x,breaks=b,plot=F)
		m = c(m,  min(h$counts))
		M = c(M,  max(h$counts))
		a = c(a, mean(h$counts))
	}
	m= sort(m)
	m025 = m[round(0.025*10000)]
	M= sort(M)
	M975 = M[round(0.975*10000)]
	return(list(c(m025,M975),mean(a),b))
}



## ADD PLTS of the prior distrbutions
pr01   = get_prior_range(n_samples=length(selected))
pr1500 = get_prior_range(min_v=0,max_v=1500,n_samples=length(selected),bins=20)

pdf("/Users/danielesilvestro/Desktop/ClownSimulator/output_simulations.pdf",4*3.5, 2*3.5)
par(mfrow=c(2,4))
library(scales)
hist(post$durationHybridization,breaks=pr1500[[3]],main="Duration hybridization phase", xlab="n. generations")
abline(h= pr1500[[2]],col="red")
abline(h= pr1500[[1]],col="red",lty=2)
hist(post$timeSinceHybridization,breaks=pr1500[[3]],main="Time since hybridization", xlab="n. generations")
abline(h= pr1500[[2]],col="red")
abline(h= pr1500[[1]],col="red",lty=2)
hist(post$max_fraction_of_larvae_spB,xlim=c(0,1),main="Fraction of Sp.B immigrants",breaks=pr01[[3]])
abline(h= pr01[[2]],col="red")
abline(h= pr01[[1]],col="red",lty=2)
hist(post$max_death_rate_hybrid, main="Death probability hybrids", xlim = c(0,1))
abline(h= pr01[[2]],col="red")
abline(h= pr01[[1]],col="red",lty=2)


plot(post$timeSinceHybridization, post$max_fraction_of_larvae_spB, xlab="Time since hybridization phase",ylab="Fraction of spB immigrants",pch=19,col=alpha("black",0.1),ylim=c(0,1),xlim=c(0,1500))
plot(post$durationHybridization, post$max_fraction_of_larvae_spB,  xlab="Duration hybridization phase"  ,ylab="Fraction of spB immigrants",pch=19,col=alpha("black",0.1),ylim=c(0,1),xlim=c(0,1500))
hist(post$median_nDNAadmixture, main="Median fraction of spB nDNA",xlim = c(0,1),nclass=3)
#hist(post$Median_nDNA_in_mtDNAspB, main="Median fraction of nDNA: spB",xlim = c(0,1),nclass=3)
hist(post$mtDNA, main="Fraction of individuals with spB mtDNA",xlim = c(0,1),nclass=16)

#dev.new()
par(mfrow=c(2,4))

plot(post$durationHybridization, post$timeSinceHybridization,      xlab="Duration hybridization"  ,ylab="Time since hybridization"  ,pch=19,col=alpha("black",0.1),ylim=c(0,1500),xlim=c(0,1500))
plot(post$durationHybridization, post$mtDNA,                       xlab="Duration hybridization"  ,ylab="Fraction of  mtDNA: spB"  ,pch=19,col=alpha("black",0.1),ylim=c(0,1),xlim=c(0,1500))
plot(post$timeSinceHybridization, post$mtDNA,                      xlab="Time since hybridization",ylab="Fraction of  mtDNA: spB"  ,pch=19,col=alpha("black",0.1),ylim=c(0,1),xlim=c(0,1500))
plot(post$max_fraction_of_larvae_spB, post$mtDNA,                  xlab="Death probability hybrids",ylab="Fraction of  mtDNA: spB"  ,pch=19,col=alpha("black",0.1),ylim=c(0,1),xlim=c(0,1))
plot(post$median_nDNAadmixture, post$mtDNA,                        xlab="Median fraction of nDNA: spB",ylab="Fraction of  mtDNA: spB"  ,pch=19,col=alpha("black",0.1),ylim=c(0,1),xlim=c(0,0.2))

dev.off()