library(MCMCglmm)
library(MASS)
dat<-read.csv("gene_counts.csv")
dat$nfac<-as.factor(paste(dat$divergence, dat$nonsynonymous))
dat$length<-dat$Ln
dat$length[which(dat$synonymous==1)]<-dat$Ls[which(dat$synonymous==1)]
dat$combo<-paste(dat$organism, dat$gene, dat$Duplicate)

dat<-dat[-which(dat$organism=="bombyx" & dat$gene=="vas"),]

#########################################################################################################
# pad with missing observations: not necessary but easier than pissing around getting missing residuals #
#########################################################################################################

missing<-table(dat$combo,dat$nfac)
missing<-missing[which(rowSums(missing)!=4),]

for(i in 1:nrow(missing)){

  combi<-strsplit(rownames(missing)[i], " ")[[1]]
  combj<-which(missing[i,]==0)

  for(j in 1:length(combj)){ 
     dat[nrow(dat)+1,]<-dat[which(dat$organism==combi[1] & dat$gene==combi[2] & dat$Duplicate==combi[3])[1],]
     dat[nrow(dat),"nfac"]<-colnames(missing)[combj[j]]
     dat[nrow(dat),"count"]<-0
     dat[nrow(dat),"length"]<-1
     dat[nrow(dat),"divergence"]<-as.numeric(substr(colnames(missing)[combj[j]],1,1))
     dat[nrow(dat),"nonsynonymous"]<-as.numeric(substr(colnames(missing)[combj[j]],3,3))
  }
}

rownames(dat)<-1:nrow(dat)


m1<-MCMCglmm(count~log(length)+log(length):nonsynonymous+nonsynonymous+divergence+nonsynonymous:divergence+(nonsynonymous+divergence+ nonsynonymous:divergence):(RNAi+Immune+Male), rcov=~us(nfac):gene, family="poisson", data=dat, pl=TRUE, nitt=13000*10, thin=10*5, burnin=3000*10)


# nfac 0 0  synonymous and polymorphic  
# nfac 0 1  synonymous and divergent 
# nfac 1 0  nonsynonymous and polymorphic  
# nfac 1 1  nonsynonymous and divergent  

# Have X = 

#         b   bN  bD   bND 
# 0/0     1   0   0    0
# 0/1     0   0   1    0
# 1/0     0   1   0    0
# 1/1     0   1   1    1

# so that b=solve(X)%*%resid

X<-matrix(0,4,4)
X[1,1]<-1
X[3:4,2]<-1
X[c(2,4),3]<-1
X[4,4]<-1

uncombo<-unique(dat$combo) # unique genes

resid.transform<-model.3A$Liab
Xb<-model.3A$Liab

Xm<-model.matrix(~nonsynonymous:divergence-1+(nonsynonymous:divergence):(RNAi+Immune+Male+Female), data=dat)
# design matrix for nonsynonymous+divergence (ie. selection) effects

sol.hit<-match(colnames(Xm), colnames(model.3A$Sol))
# column positions of relevant terms

for(i in 1:nrow(model.3A$Liab)){
   resid<-model.3A$Liab[i,]-(model.3A$X%*%model.3A$Sol[i,])@x 
   # residuals for each observation at iteration i
   R<-matrix(model.3A$VCV[i,],4,4)
    for(j in 1:length(uncombo)){
       hits<-which(dat$combo==uncombo[j])
       # find positions of residuals gene j
       if(length(hits)==4){
        beta<-solve(X,resid[hits])
        resid.transform[i,hits]<-beta
        # solve for the (random) b effects for each obseravtion 
        Xb[i,hits]<-c(Xm[hits,]%*%model.3A$Sol[i,sol.hit])
        # get nonsynonymous+divergence predictions for each obseravtion 
      }
    }
}

pred<-Xb+resid.transform
# fixed and random effect prediction (only nonsynonymous+divergence predictions make sense)

gen.type<-paste(dat$RNAi, dat$Immune, dat$Male, dat$Female)

gt<-gen.type[which(dat$nonsynonymous==1 & dat$divergence==1)]
# gene type for nonsynonymous+divergence observations

pms<-colMeans(pred)[which(dat$nonsynonymous==1 & dat$divergence==1)][order(gt)]
# posterior means for each nonsynonymous+divergence observation sorted by gene type
ints<-HPDinterval(pred)[which(dat$nonsynonymous==1 & dat$divergence==1)[order(gt)],]
# 95% credible intervals for each nonsynonymous+divergence observation sorted by gene type


par(bty="l")
plot(pms, col=as.factor(gt[order(gt)]), pch=16, cex=0.5,ylim=c(-4,20), xaxt = "n", ylab="Selection Effect", xlab="")
for(i in 1:nrow(ints)){
lines(c(i,i), ints[i,], col=as.factor(gt[order(gt)])[i])
}
# plot each gene's selection effect with 95% credible intervals

dat.red<-dat[which(dat$nonsynonymous==1 & dat$divergence==1),]
dat.red<-dat.red[which(!duplicated(paste(dat.red$RNAi, dat.red$Immune, dat.red$Male, dat.red$Female))),]
# get a data set with each gene type represented once.

X.red<-model.matrix(~nonsynonymous:divergence-1+(nonsynonymous:divergence):(RNAi+Immune+Male+Female), data=dat.red)
red.type<-paste(dat.red$RNAi, dat.red$Immune, dat.red$Male, dat.red$Female)

type.mean<-as.mcmc(t(apply(model.3A$Sol[,sol.hit], 1, function(x){X.red%*%x}))[,order(red.type)])
# get posterior distribution for each gene types *mean effect* 

class.name<-c("All_others","RNAi", "Immune", "Male", "Female")
for(i in 1:ncol(type.mean)){
lines(c(1,cumsum(table(gt)))[i:(i+1)], rep(colMeans(type.mean)[i],2), col=i)
lines(c(1,cumsum(table(gt)))[i:(i+1)], rep(HPDinterval(type.mean)[i,1],2), col=i, lty=2)
lines(c(1,cumsum(table(gt)))[i:(i+1)], rep(HPDinterval(type.mean)[i,2],2), col=i, lty=2)
text(mean(c(1,cumsum(table(gt)))[i:(i+1)]), c(-4.1,-3.5)[1+as.numeric(i==3)], class.name[i], col=i)
}
# plot each gene types mean selection effect with 95% credible intervals




