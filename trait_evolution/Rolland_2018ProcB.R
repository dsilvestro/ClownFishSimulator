require(phytools)
#require(diversitree)
library(geiger)
library(TreeSim)
library(mvMORPH)
# RUN UNDER WHITE NOISE

#### GET LOG-UNIF SAMPLES of Ne and generation time
#popsizeBEAST = rgamma(1000000,shape=11.436259,rate=4.919814)
pop_size_range= c(0.6778258, 5.986533) #c(1.0928, 3.7512)  
popsizeBEAST = runif(100000,pop_size_range[1],pop_size_range[2]) #
gen_time     = runif(1000000,3,7)/1000000
popsize      = popsizeBEAST / gen_time
logpopsize = log(popsize)

minPOP = min(logpopsize) #[0.01*length(logpopsize)]
maxPOP = max(logpopsize) #[0.99*length(logpopsize)]

indx = c()
p = minPOP
for ( i in 1:10){ # uniform samples from a discretized 10-bin distribution
	p1 = p+(maxPOP-minPOP)/10
	w_temp = which(logpopsize>p & logpopsize<p1)
	print(c(p,p1))
	indx = c(indx,sample(w_temp,100,replace=F))
	p = p1
	
}
#hist(logpopsize[indx],nclass=100)
#hist(gen_time[indx])
#plot(gen_time[indx]~logpopsize[indx])

logpopsizeVEC = logpopsize[indx]
gen_timeVEC = gen_time[indx]

"
> min(gen_timeVEC)
[1] 3.004135e-06
> max(gen_timeVEC)
[1] 6.999013e-06
> min(logpopsizeVEC)
[1] 11.51334
> max(logpopsizeVEC)
[1] 14.49601


r = c(100041,1974863)


> r = c(100041,1974863)
> r[2]/r[1]
[1] 19.74054
> r[1]/19.74054
[1] 5067.794
> r[2]*19.74054
[1] 38984862
> log(r[2]*19.74054)
[1] 17.47868
> log(r[1]/19.74054)
[1] 8.530661


"





get_macro_sd <- function(micro_rate_WN,empirical_pop_size){
	macro_sd_per_generation = sqrt(micro_rate_WN)/sqrt(empirical_pop_size)
	return(macro_sd_per_generation)
}

run_sim_fast <- function(br_length=1000,gen_time=1,popsize = 100, sd_per_spec_per_gen = 1,WN_type=2 , bounds = c(-Inf,Inf),anc_mean=0){
	n_generations = round(br_length/gen_time)
	mean_by_generation = rnorm(n_generations,0,sd_per_spec_per_gen)	
	M = anc_mean+sum(mean_by_generation)
	return(M)
}

plot_WN_lineage <- function(){
	end_point=c()
	for (i in 1:100){
		M = run_sim()
		print(i)
		if (i==1){
			plot(M,type="l")
		}else{
			lines(M)
		}
		end_point= c(end_point,M[length(M)])
	}
	sd(end_point)/n_generations
}

simBMT_shifts <- function(tree,a=0,mu=0,sig2=1,internal=T,nsim=1,gen_time=1,popsize = 100, WN_sd = 1,EV_model = 2, bounds = c(-Inf,Inf),sd_per_spec_per_gen=0){
	nsp<-length(tree$tip) #nb of species
	
	# first simulate changes along each branch
	nedges<-length(tree$edge.length) #number of edges in the tree
	ndraws<-nedges*nsim #total number of calls to rnorm
	#x<-matrix(data=rnorm(n=ndraws,mean=rep(mu*tree$edge.length,nsim),sd=rep(sqrt(sig2*tree$edge.length),nsim)),nedges,nsim)
	
	# evolve WN
	x<-matrix(data=rep(0,nedges*nsim),nedges,nsim)
	#for (i in 1:nedges){
	#	M = run_sim(tree$edge.length[i],gen_time,popsize,WN_sd,EV_model,bounds)
	#	x[i] = M[length(M)]
	#}
	
	# now add them up
	x_temp<-matrix(data=rep(0,nedges*nsim),nedges,nsim)
	
	y<-array(0,dim=c(nrow(tree$edge),ncol(tree$edge),nsim))
	for(i in 1:nrow(x)){
		#print(c("Edge",i))
		if(tree$edge[i,1]==(nsp+1))
			y[i,1,]<-a
		else
			y[i,1,]<-y[match(tree$edge[i,1],tree$edge[,2]),2,]
		
		#sim_value = run_sim(tree$edge.length[i],gen_time,popsize,WN_sd,EV_model,bounds,anc_mean=y[i,1,])
		#sim_value = sim_value[length(sim_value)]
		sim_value = run_sim_fast(tree$edge.length[i],gen_time,popsize,sd_per_spec_per_gen[i],EV_model,bounds,anc_mean=y[i,1,])
		sim_value[sim_value>max(bounds)] = max(bounds)- (sim_value[sim_value>max(bounds)]-max(bounds))
		sim_value[sim_value<min(bounds)] = min(bounds)+ (min(bounds)-sim_value[sim_value<min(bounds)])
		if (sim_value>max(bounds) | sim_value<min(bounds)){
			sim_value=runif(1,min(bounds),max(bounds))
		}
		
		y[i,2,]<- sim_value
		x_temp[i] = sim_value
	}
	
	x= x_temp
	x<-matrix(data=rbind(y[1,1,],as.matrix(y[,2,])), nedges+1,nsim)
	rownames(x)<-c(nsp+1,tree$edge[,2])
	x<-as.matrix(x[as.character(1:(nsp+tree$Nnode)),])
	rownames(x)[1:nsp]<-tree$tip.label
	
	if(internal==TRUE) {
		return(x[1:nrow(x),]) # include internal nodes
	}
	else {
		return(x[1:length(tree$tip.label),]) # tip nodes only
	}
}

sim_micro_macro_empirical <- function(tree_file, gen_time_range=1,pop_size_range = c(100,10000), micro_rate_range = c(1,1),EV_model = 2,plot_pheno =T,nsim=1,bounds = c(-Inf,Inf),outfile="simulations.txt",use_geiger=F){
	if (use_geiger==T){
		cat(c("pop_size","gen_time","rateBM","rateWN","rateOU","alphaOU","expVarOU","aicBM","aicWN","aicOU","\n"),file=outfile,sep="\t")
	}else{
		cat(c("pop_size","gen_time","rateBM","rateOU","alphaOU","expVarOU","aicBM","aicOU","\n"),file=outfile,sep="\t")
	}
	print("Reading trees...")
	
	deltabounds=max(bounds)-min(bounds)
	soft_b = 0.025
	soft_bounds = c( min(bounds)-soft_b*deltabounds, max(bounds)+soft_b*deltabounds )
	#print(soft_bounds)
	tree_list = read.nexus(tree_file)
	#for (i in 1:nsim){
	i = 0
	while (i < nsim){	
		# popsize = round(runif(1,pop_size_range[1],pop_size_range[2])) #round(mean(pop_size_range))
		tree<- tree_list[[sample(1:length(tree_list),1)]]
		nedges<-length(tree$edge.length)
		#__ # USE VARIABLE POP SIZE		
		#__ #popsize = round(exp(runif(nedges,log(pop_size_range[1]),log(pop_size_range[2])))) # uniform in log space		
		#__ #popsize = round(exp(runif(nedges,log(pop_size_range[1]),log(pop_size_range[2])))) 
		#__ # USE ONLY ONE POP SIZE
		#__ popsize = rep(round(exp(runif(1,log(pop_size_range[1]),log(pop_size_range[2])))) , nedges )
		#__ # USE VARIABLE POP SIZE	1/5 of clarkii - to clarkii's	
		#__ # r_temp = runif(1,pop_size_range[1],pop_size_range[2])
		#__ # popsize = round(runif(nedges,r_temp/5,r_temp))
		#__ #print(log(popsize))
		print("Simulating trait...")
		
		#popsizeBEAST = rep(exp(runif(1,log(pop_size_range[1]),log(pop_size_range[2]))) , nedges )
		#popsizeBEAST = rep(runif(1,pop_size_range[1],pop_size_range[2]) , nedges )
		#gen_time     = runif(1,gen_time_range[1],gen_time_range[2])
		#popsize      = round(popsizeBEAST / gen_time)
		
		j =sample(1:length(logpopsizeVEC),1)
		popsize = exp(logpopsizeVEC[j])
		gen_time  = gen_timeVEC[j]  
		popsize = rep(popsize,nedges)
		
		#__ popsizeBEAST = rep(runif(1,pop_size_range[1],pop_size_range[2]) , nedges )
		#__ gen_time     = rgamma(1,shape=10,rate=2)/1000000
		#__ popsize      = round(popsizeBEAST / gen_time)		
		
		
		micro_rate = runif(1,micro_rate_range[1],micro_rate_range[2])
		macro_rate_per_generation = get_macro_sd(micro_rate,popsize)
		print(log(popsize))
		
		# ADD RANDOM DATING ERROR
		tree$edge.length = tree$edge.length + (runif(1,-0.17,0.17))*tree$edge.length
		
		# ROOT STATE
		if (max(bounds)<Inf){
			anc_state = runif(1,min(bounds),max(bounds))
			print(c("ANC STATE",anc_state))
		}else{
			anc_state=0
		}
		
		
		s =simBMT_shifts(tree,gen_time=gen_time,popsize=popsize, EV_model=EV_model,
			bounds=soft_bounds,sd_per_spec_per_gen=macro_rate_per_generation)
		if (plot_pheno==T){phenogram(tree, s, spread.labels=F)}		
		
		if (use_geiger==T){
			bm = fitContinuous(tree,s[1:length(tree$tip.label)],model = "BM") 
			wn = fitContinuous(tree,s[1:length(tree$tip.label)],model = "white")
			ou = fitContinuous(tree,s[1:length(tree$tip.label)],model = "OU",bounds=list(alpha=c(0,1000000)))
			res=NULL
			res$rateBM   = bm$opt$sigsq
			res$sig2WN   = wn$opt$sigsq
			res$sig2OU   = ou$opt$sigsq
			res$sig2Apha = ou$opt$alpha
			res$sig2alpha= ou$opt$sigsq/(2*ou$opt$alpha)
			res$aicBM    = bm$opt$aicc
			res$aicWN    = wn$opt$aicc
			res$aicOU    = ou$opt$aicc
			print(res) 
			cat(c(mean(popsize),gen_time*1000000,as.numeric(res),"\n"),file=outfile,sep="\t",append=T)
			
		}else{
			run_opt <- function(tree,trait,outfile){
				bm = mvBM(tree,trait) 
				ou = mvOU(tree,trait)
				res=NULL
				res$rateBM   = bm$sigma 
				res$sig2OU   = ou$sigma 
				res$sig2Apha = ou$alpha 
				res$sig2alpha= ou$sigma/(2* ou$alpha)   
				res$aicBM    = bm$AICc  
				res$aicOU    = ou$AICc  
				print(i) 
				cat(c(mean(popsize),gen_time*1000000,as.numeric(res),"\n"),file=outfile,sep="\t",append=T)
			}
			try( run_opt(tree,s[1:length(tree$tip.label)],outfile) )
			i = i+1
			if (i >= nsim){break}
			
		}
		
		# 		
	
		
		
	}
}
 
#### EMPIRICAL DATA
wd = "~/MicroMacro" 
setwd(wd)
tree_file= "~/MicroMacro/macro.trees.11.5803.nex"
# these parameters are not actually used anymore (need to cleanup)
gen_time_range = c(3,7)/1000000
pop_size_range= c(1.0928, 3.7512)  
pop_size_range= c(8, 12)  

plot_pheno = F   # plot phenogram (true anc states)
nsim = 1100     # number of simulations


# get empirical bounds based on OU-macro stationary variance
trait1_optimum  = 0.849
trait2_optimum  = 0.34
trait5_optimum  = -0.54
trait9_optimum  = 0.892
trait10_optimum = 0.283

trait1_stationaryVar  = 3.819 
trait2_stationaryVar  = 0.847 
trait5_stationaryVar  = 0.724 
trait9_stationaryVar  = 0.842 
trait10_stationaryVar = 0.887 

get_bounds <- function(m,v){
	range_size= 1.96*sqrt(v)
	return(c(m-range_size,m+range_size))
	
}

#_ trait 1:  stat_var=0.285, variance=0.291; 
#_ trait 2:  stat_var=0.952, variance=0.969; 
#_ trait 5:  stat_var=0.808, variance=0.823; 
#_ trait 9:  stat_var=0.414, variance=0.422; 
#_ trait 10: stat_var=0.913, variance=0.931.

# Trait 1
wn_micro = c(0.291,0.291)
bounds =c(-0.9964155,3.1070926) 
bounds = get_bounds(trait1_optimum,trait1_stationaryVar)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait1SD.txt",bounds=bounds)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait1noBoundsSD.txt",bounds=c(-1000,1000))


# Trait 2
wn_micro = c(0.969,0.969)
bounds = c(-0.8749808,3.3783053) #c(-2.094159,3.3783053)
bounds = get_bounds(trait2_optimum,trait2_stationaryVar)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait2SD.txt",bounds=bounds)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait2noBoundsSD.txt",bounds=c(-1000,1000))


# Trait 5
wn_micro = c(0.823,0.823)
bounds = c(-2.202651,1.285620)
bounds = get_bounds(trait5_optimum,trait5_stationaryVar)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait5SD.txt",bounds=bounds)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait5noBoundsSD.txt",bounds=c(-1000,1000) )


# Trait 9
wn_micro = c(0.422,0.422)
bounds = c(-1.954580,2.266513)
bounds = get_bounds(trait9_optimum,trait9_stationaryVar)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait9SD.txt",bounds=bounds)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait9noBoundsSD.txt",bounds=c(-1000,1000) )


# Trait 10
wn_micro = c(0.931,0.931)
bounds = c(-1.385118,2.422880)
bounds = get_bounds(trait10_optimum,trait10_stationaryVar)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait10SD.txt",bounds=bounds)

sim_micro_macro_empirical(tree_file,gen_time_range,pop_size_range,micro_rate=wn_micro,EV_model=2,plot_pheno,nsim,outfile="simulations_trait10noBoundsSD.txt",bounds=c(-1000,1000) )



### PLOT RESULTS
bestmodel <-function(vec){
	return(which(vec==min(vec)))
}

plot_sim_output <- function(logfile,trait_n="",BMrate=NA,WNrate=NA,OUrate=NA,OUalpha=NA,Ylim=c(2,2,2,2),plotWN=F,transp=0.2){
	tbl = read.table(logfile,h=T)
	if (plotWN==T){
		AICs = cbind(tbl$aicBM,tbl$aicOU)
		H = hist(apply(AICs,FUN=bestmodel,1),breaks=c(0,1.5,2.5,4),plot=F)$counts
		names(H)=c("BM","OU","WN")		
	}else{
		AICs = cbind(tbl$aicBM,tbl$aicOU)
		H = hist(apply(AICs,FUN=bestmodel,1),breaks=c(0,1.5,2.5),plot=F)$counts
		names(H)=c("BM","OU")	
	}
	barplot(H,main=paste("Trait ",trait_n),col=c("#0868ac" ,"#2ca25f" ))
	Pval = round(length(intersect(which(tbl$rateBM < BMrate[2]),which(tbl$rateBM > BMrate[1])))/length(tbl$rateBM),3)
	plot(tbl$rateBM ~ (tbl$gen_time),pch=19,ylim=c(0,Ylim[1]),main=sprintf("BM sig2 (P = %s)", Pval),col=alpha("#0868ac",transp))
	#boxplot(tbl$rateBM,pch=19,main="BM sig2")
	abline(h=BMrate[1],col="red")
	abline(h=BMrate[2],col="red")
	if (plotWN==T){
		plot(tbl$rateWN ~ (tbl$gen_time),pch=19,ylim=c(0,Ylim[2]),main="WN sig2")
		abline(h=WNrate[1],col="red")
		abline(h=WNrate[2],col="red")
	}
	Pval = round(length(intersect(which(tbl$rateOU < OUrate[2]),which(tbl$rateOU > OUrate[1])))/length(tbl$rateOU),3)
	plot(tbl$rateOU ~ (tbl$gen_time),pch=19,ylim=c(0,Ylim[3]),main=sprintf("OU sig2 (P = %s)",Pval),col=alpha("#2ca25f",transp))
	abline(h=OUrate[1],col="red")
	abline(h=OUrate[2],col="red")
	Pval = round(length(intersect(which(tbl$alphaOU < OUalpha[2]),which(tbl$alphaOU > OUalpha[1])))/length(tbl$alphaOU),3)
	plot(tbl$alphaOU ~ (tbl$gen_time),pch=19,ylim=c(0,Ylim[4]),main=sprintf("OU alpha (P = %s)",Pval),col=alpha("#66c2a4",transp))
	abline(h=OUalpha[1],col="red")
	abline(h=OUalpha[2],col="red")
	#plot(tbl$expVarOU ~ log(tbl$pop_size),pch=19,main="OU exp var")
	#abline(h=OUexp,col="red")
}

plot_sim_output <- function(logfile,trait_n="",BMrate=NA,WNrate=NA,OUrate=NA,OUalpha=NA,Ylim=c(2,2,2,2),plotWN=F,transp=0.2){
	tbl = read.table(logfile,h=T)
	if (plotWN==T){
		AICs = cbind(tbl$aicBM,tbl$aicOU)
		H = hist(apply(AICs,FUN=bestmodel,1),breaks=c(0,1.5,2.5,4),plot=F)$counts
		names(H)=c("BM","OU","WN")		
	}else{
		AICs = cbind(tbl$aicBM,tbl$aicOU)
		H = hist(apply(AICs,FUN=bestmodel,1),breaks=c(0,1.5,2.5),plot=F)$counts
		names(H)=c("BM","OU")	
	}
	barplot(H,main=paste("Trait ",trait_n),col=c("#0868ac" ,"#2ca25f" ))
	Pval = round(length(intersect(which(tbl$rateBM < BMrate[2]),which(tbl$rateBM > BMrate[1])))/length(tbl$rateBM),3)
	plot(tbl$rateBM ~ log(tbl$pop_size),pch=19,ylim=c(0,Ylim[1]),main=sprintf("BM sig2 (P = %s)", Pval),col=alpha("#0868ac",transp))
	#boxplot(tbl$rateBM,pch=19,main="BM sig2")
	abline(h=BMrate[1],col="red")
	abline(h=BMrate[2],col="red")
	if (plotWN==T){
		plot(tbl$rateWN ~ log(tbl$pop_size),pch=19,ylim=c(0,Ylim[2]),main="WN sig2")
		abline(h=WNrate[1],col="red")
		abline(h=WNrate[2],col="red")
	}
	Pval = round(length(intersect(which(tbl$rateOU < OUrate[2]),which(tbl$rateOU > OUrate[1])))/length(tbl$rateOU),3)
	plot(tbl$rateOU ~ log(tbl$pop_size),pch=19,ylim=c(0,Ylim[3]),main=sprintf("OU sig2 (P = %s)",Pval),col=alpha("#2ca25f",transp))
	abline(h=OUrate[1],col="red")
	abline(h=OUrate[2],col="red")
	Pval = round(length(intersect(which(tbl$alphaOU < OUalpha[2]),which(tbl$alphaOU > OUalpha[1])))/length(tbl$alphaOU),3)
	plot(tbl$alphaOU ~ log(tbl$pop_size),pch=19,ylim=c(0,Ylim[4]),main=sprintf("OU alpha (P = %s)",Pval),col=alpha("#66c2a4",transp))
	abline(h=OUalpha[1],col="red")
	abline(h=OUalpha[2],col="red")
	#plot(tbl$expVarOU ~ log(tbl$pop_size),pch=19,main="OU exp var")
	#abline(h=OUexp,col="red")
}

plot_sim_output <- function(logfile,trait_n="",BMrate=NA,WNrate=NA,OUrate=NA,OUalpha=NA,Ylim=c(2,2,2,2),plotWN=F,transp=0.2){
	tbl = read.table(logfile,h=T)
	if (plotWN==T){
		AICs = cbind(tbl$aicBM,tbl$aicOU)
		H = hist(apply(AICs,FUN=bestmodel,1),breaks=c(0,1.5,2.5,4),plot=F)$counts
		names(H)=c("BM","OU","WN")		
	}else{
		AICs = cbind(tbl$aicBM,tbl$aicOU)
		H = hist(apply(AICs,FUN=bestmodel,1),breaks=c(0,1.5,2.5),plot=F)$counts
		names(H)=c("BM","OU")	
	}
	barplot(H,main=paste("Trait ",trait_n),col=c("#0868ac" ,"#2ca25f" ))
	Pval = round(length(intersect(which(tbl$rateBM < BMrate[2]),which(tbl$rateBM > BMrate[1])))/length(tbl$rateBM),3)
	plot(tbl$rateBM,pch=19,ylim=c(0,Ylim[1]),main=sprintf("BM sig2 (P = %s)", Pval),col=alpha("#0868ac",transp))
	#boxplot(tbl$rateBM,pch=19,main="BM sig2")
	abline(h=BMrate[1],col="red")
	abline(h=BMrate[2],col="red")
	if (plotWN==T){
		plot(tbl$rateWN ~ log(tbl$pop_size),pch=19,ylim=c(0,Ylim[2]),main="WN sig2")
		abline(h=WNrate[1],col="red")
		abline(h=WNrate[2],col="red")
	}
	Pval = round(length(intersect(which(tbl$rateOU < OUrate[2]),which(tbl$rateOU > OUrate[1])))/length(tbl$rateOU),3)
	plot(tbl$rateOU,pch=19,ylim=c(0,Ylim[3]),main=sprintf("OU sig2 (P = %s)",Pval),col=alpha("#2ca25f",transp))
	abline(h=OUrate[1],col="red")
	abline(h=OUrate[2],col="red")
	Pval = round(length(intersect(which(tbl$alphaOU < OUalpha[2]),which(tbl$alphaOU > OUalpha[1])))/length(tbl$alphaOU),3)
	plot(tbl$alphaOU,pch=19,ylim=c(0,Ylim[4]),main=sprintf("OU alpha (P = %s)",Pval),col=alpha("#66c2a4",transp))
	abline(h=OUalpha[1],col="red")
	abline(h=OUalpha[2],col="red")
	#plot(tbl$expVarOU ~ log(tbl$pop_size),pch=19,main="OU exp var")
	#abline(h=OUexp,col="red")
}

wd = "~/MicroMacro/" # outout is saved as: "simulations.txt"
setwd(wd)
library(scales)
plotWN=F

if (plotWN==T){
	pdf("results_simulations_CI-SDvarPOPmvMORPHvariance.pdf",width=16, height=8)
}else{
	pdf("results_simulations_CI-SDvarPOPmvMORPHvariance.pdf",width=16, height=8)
}
par(mfrow=c(2,4))
sigBM = c(0.204,0.398)
#sigWN = c(1.891,1.891)
sigOU = c(0.204,0.688)
alpOU = c(0,0.163)
Ylim = c(1,3,1.5,1.5)
plot_sim_output("simulations_trait1SD.txt","1 (BM)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)
plot_sim_output("simulations_trait1noBoundsSD.txt","1 no bounds (BM)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)

par(mfrow=c(2,4))
sigBM = c(0.165,0.407)
#sigWN = c(0.869,0.869)
sigOU = c(0.347,48.922)
alpOU = c(0.206,29.75)
Ylim = c(2,3,10,10)
plot_sim_output("simulations_trait2SD.txt","2 (OU)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)
plot_sim_output("simulations_trait2noBoundsSD.txt","2 no bounds (OU)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)

par(mfrow=c(2,4))
sigBM = c(0.087,0.144)
#sigWN = c(0.712,0.712)
sigOU = c(0.116,0.246)
alpOU = c(0.068,0.187)
Ylim = c(1,5,2.2,2.2)
plot_sim_output("simulations_trait5SD.txt","5 (BM)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)
plot_sim_output("simulations_trait5noBoundsSD.txt","5 no bounds (BM)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)

par(mfrow=c(2,4))
sigBM = c(0.258,0.433)
#sigWN = c(0.783,0.783)
sigOU = c(0.692,3.009) 
alpOU = c(0.42,1.781)
Ylim = c(1,5,10,10)
plot_sim_output("simulations_trait9SD.txt","9 (OU)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)
plot_sim_output("simulations_trait9noBoundsSD.txt","9 no bounds (OU)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)

par(mfrow=c(2,4))
sigBM = c(0.163,0.324)
#sigWN = c(0.906,0.906)
sigOU = c(0.268,1.111)
alpOU = c(0.136,0.631)
Ylim = c(1,5,10,10)
plot_sim_output("simulations_trait10SD.txt","10 (OU)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)
plot_sim_output("simulations_trait10noBoundsSD.txt","10 no bounds (OU,BM)",sigBM,sigWN,sigOU,alpOU,Ylim,plotWN)
dev.off()





## RANDOM STUFF

## CHECK POPSIZE RANGE (Ne)
gen_time_range = c(2,8)/1000000
pop_size_range= c(1.0928, 3.7512)  

popsizeBEAST = runif(100000,pop_size_range[1],pop_size_range[2]) #
popsizeBEAST = exp(runif(100000,log(pop_size_range[1]),log(pop_size_range[2])))
gen_time     = runif(100000,gen_time_range[1],gen_time_range[2])
popsize      = popsizeBEAST / gen_time
hist(log(round(popsize)))

popsizeBEAST = runif(100000,pop_size_range[1],pop_size_range[2])
gen_time     = rgamma(100000,shape=10,rate=2)/1000000
popsize      = popsizeBEAST / gen_time
hist(log(round(popsize)))





