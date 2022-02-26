#######################################################################
#
#
#   Code to plot Polymorphism data with different sampling strategies and then run asymptotic MK
#
#
#######################################################################




#######################################################################
#plot the age of segregating alleles
#find plot limits
min_Age<-min_density<-1e20
max_Age<-min_density<-0
for(i in 1:length(polymorphismData)){
	dat<-polymorphismData[[i]]
	for(j in 1:length(dat)){
		min_Age<-min(c(min_Age,dat[[j]]$SynAge,dat[[j]]$DelAge))
		max_Age<-max(c(max_Age,dat[[j]]$SynAge,dat[[j]]$DelAge))
		
	}
}
if(min_Age==0)min_Age<-0.5

plotdim<-ceiling(sqrt(length(polymorphismData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)	
)
par(mfrow=ceiling(c(sqrt(length(polymorphismData)),sqrt(length(polymorphismData)))))
for(i in 1:length(polymorphismData)){

	dat<-polymorphismData[[i]]

	com_SynAge<-vector()
	for(j in 1:length(dat)){com_SynAge<-c(com_SynAge,dat[[j]]$SynAge)}
	D_com_SynAge<-density(log10(com_SynAge),cut=TRUE)
	com_DelAge<-vector()
	for(j in 1:length(dat)){com_DelAge<-c(com_DelAge,dat[[j]]$DelAge)}
	D_com_DelAge<-density(log10(com_DelAge),cut=TRUE)	
	
	plot(D_com_SynAge,xlab="",ylab="",type="n",xlim=log10(c(min_Age,max_Age)),ylim=c(0,0.8),main="",axes=FALSE)
	abline(v=log10(signif(10^axTicks(1),1)),col="gray80",lwd=0.5)
	##abline(v=log10(c(300000,600000)))
	polygon(c(min(D_com_DelAge$x),D_com_DelAge$x,max(D_com_DelAge$x)),c(0,D_com_DelAge$y,0),col=adjustcolor("coral4", alpha.f = 0.5),border="coral4")
	polygon(c(min(D_com_SynAge$x),D_com_SynAge$x,max(D_com_SynAge$x)),c(0,D_com_SynAge$y,0),col=adjustcolor("steelblue4", alpha.f = 0.5),border="steelblue4")
	box(col="gray70")
	
	#pretty Plotting
	as.numeric(gsub("_.+","",names(polymorphismData)[i]))->m
	as.numeric(gsub(".+_","",names(polymorphismData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Allele Age\n(x 1000 generations)",side=1,line=5,cex=1.00)
		axis(1,at=log10(signif(10^axTicks(1),1)),labels=signif(10^axTicks(1),1)/1000,cex.axis=0.8)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*Ne*m),2),sep=""),side=2,line=6,cex=1.4)
		mtext("Density",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),las=2,line=1,cex.axis=0.8)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	
	if(i==1)legend("topleft",legend=c("Unconstrained","Weakly deleterious"),text.col=c("steelblue4","coral4"),cex=1.0,bty="n")

}



#######################################################################
#plot the age of segregating alleles WITH POS
#find plot limits
min_Age<-min_density<-1e20
max_Age<-min_density<-0
for(i in 1:length(polymorphismData)){
	dat<-polymorphismData[[i]]
	for(j in 1:length(dat)){
		min_Age<-min(c(min_Age,dat[[j]]$SynAge,dat[[j]]$DelAge,dat[[j]]$PosAge))
		max_Age<-max(c(max_Age,dat[[j]]$SynAge,dat[[j]]$DelAge,dat[[j]]$PosAge))
		
	}
}
if(min_Age==0)min_Age<-0.9

plotdim<-ceiling(sqrt(length(polymorphismData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)
par(mfrow=ceiling(c(sqrt(length(polymorphismData)),sqrt(length(polymorphismData)))))
for(i in 1:length(polymorphismData)){

	dat<-polymorphismData[[i]]

	com_SynAge<-vector()
	for(j in 1:length(dat)){com_SynAge<-c(com_SynAge,dat[[j]]$SynAge)}
	D_com_SynAge<-density(log10(com_SynAge),cut=TRUE)
	#D_com_SynAge$y<-D_com_SynAge$y*length(com_SynAge)
	com_DelAge<-vector()
	for(j in 1:length(dat)){com_DelAge<-c(com_DelAge,dat[[j]]$DelAge)}
	D_com_DelAge<-density(log10(com_DelAge),cut=TRUE)
	#D_com_DelAge$y<-D_com_DelAge$y*length(com_DelAge)
	com_PosAge<-vector()
	for(j in 1:length(dat)){com_PosAge<-c(com_PosAge,dat[[j]]$PosAge)}
	D_com_PosAge<-density(log10(com_PosAge),cut=TRUE)
	#D_com_PosAge$y<-D_com_PosAge$y*length(com_PosAge)	 # can't be seen!
	
	
	plot(D_com_SynAge,xlab="",ylab="",type="n",xlim=log10(c(min_Age,max_Age)),ylim=c(0,1.8),main="",axes=FALSE)
	abline(v=log10(signif(10^axTicks(1),1)),col="gray80",lwd=0.5)
	polygon(c(min(D_com_PosAge$x),D_com_PosAge$x,max(D_com_PosAge$x)),c(0,D_com_PosAge$y,0),col=adjustcolor("olivedrab4", alpha.f = 0.5),border="olivedrab4")
		polygon(c(min(D_com_DelAge$x),D_com_DelAge$x,max(D_com_DelAge$x)),c(0,D_com_DelAge$y,0),col=adjustcolor("coral4", alpha.f = 0.5),border="coral4")
	polygon(c(min(D_com_SynAge$x),D_com_SynAge$x,max(D_com_SynAge$x)),c(0,D_com_SynAge$y,0),col=adjustcolor("steelblue4", alpha.f = 0.5),border="steelblue4")

	box(col="gray70")
	
	#pretty Plotting
	as.numeric(gsub("_.+","",names(polymorphismData)[i]))->m
	as.numeric(gsub(".+_","",names(polymorphismData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Allele Age\n(x 1000 generations)",side=1,line=5,cex=1.00)
		axis(1,at=log10(signif(10^axTicks(1),1)),labels=signif(10^axTicks(1),1)/1000,cex.axis=0.8)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*Ne*m),2),sep=""),side=2,line=6,cex=1.4)
		mtext("Density",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),las=2,line=1,cex.axis=0.8)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	
	if(i==1)legend("topleft",legend=c("Beneficial","Unconstrained","Weakly deleterious"),text.col=c("olivedrab4","steelblue4","coral4"),cex=1.0,bty="n")

}




#######################################################################
#DFE for fixations and polymorphisms

# want to compare to new mutations
# in Dmel BC aims for 2Ns=2000. = 2 x 36 x 1000 x s => 2000/(2*36*500) =0.0556
# initializeMutationType("m2", 0.5, "g", -0.056, 0.3);
# mean/shape=scale

hist(-rgamma(10000,0.056/0.3,0.3)*2*Ne,breaks=c(-1e100,-100,-10,-1,-0.1,0),plot=FALSE)$count->muts ; muts<-muts/sum(muts)
Ne=500
nPops=36

#Plot DFE of Fixations and SNPs
plotdim<-ceiling(sqrt(length(polymorphismData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)
for(i in 1:length(polymorphismData)){

	dat<-polymorphismData[[i]]
	com_snpDFE_Non<-vector()
	for(j in 1:length(dat)){com_snpDFE_Non<-c(com_snpDFE_Non,dat[[j]]$snpDFE_Non)} #work in positive
	hist(com_snpDFE_Non*2*Ne,breaks=c(-1e100,-100,-10,-1,-0.1,0),plot=FALSE)$count->SNPs ; SNPs<-SNPs/sum(SNPs)

	dat<-fixationData[[i]]
	com_DelFixDFE<-vector()
	for(j in 1:length(dat)){com_DelFixDFE<-c(com_DelFixDFE,dat[[j]]$DelFixDFE)} #work in positive
	hist(com_DelFixDFE*2*Ne,breaks=c(-1e100,-100,-10,-1,-0.1,0),plot=FALSE)$count->Diffs; Diffs<-Diffs/sum(Diffs)
		
	rbind(muts,SNPs,Diffs)->toPlot
	##colnames(toPlot)<-c("below -100","-100 to -10","-10 to -1","-1 to 0")
	
	#names(polymorphismData)[i]
	#xlab="2Ns"
	#ylab="Density"

	barplot(toPlot,xlab="",ylab="",main="",col=c("gray70","coral3","coral4"),ylim=c(0,1),beside=TRUE,axes=FALSE,border=NA,space=c(0,1))->mids
	box(col="gray70")
	#axis(1,at=c(1,1e-2,1e-4,1e-6,1e-8,1e-10),labels=c(1,1e-2,1e-4,1e-6,1e-8,1e-10))

	#pretty Plotting
	as.numeric(gsub("_.+","",names(fixationData)[i]))->m
	as.numeric(gsub(".+_","",names(fixationData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		#mtext("Deleterious fitness effect",side=1,line=5,cex=1.02)
		mtext(c("(-inf ,-100]","(-100,-10]","(-10,-1]","(-1,-0.1]","(-0.1,0]"),side=1,line=1,cex=0.7,at=mids[c(2,5,8,11,14)],las=2)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*Ne*m),2),sep=""),side=2,line=5,cex=1.4)
		mtext("Proportion",side=2,line=2,cex=1.0)
		axis(2,at=c(0,1.0),las=2,line=1)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	
	if(i==1){legend("topleft",legend=c("Mutations","Polymorphisms","Fixations"),text.col=c("gray60","coral3","coral4"),cex=1.0,bty="n")}

	
}


#Plot folded SFS

plotdim<-ceiling(sqrt(length(polymorphismData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)	
)
par(mfrow=ceiling(c(sqrt(length(polymorphismData)),sqrt(length(polymorphismData)))))
for(i in 1:length(polymorphismData)){

	dat<-polymorphismData[[i]]
	com_snpSFS_Non<-matrix(0,nrow=length(dat),ncol=2*indivs_sampled)
	com_snpSFS_Syn<-matrix(0,nrow=length(dat),ncol=2*indivs_sampled)
	for(j in 1:length(dat)){
		com_snpSFS_Non[j,as.numeric(names(dat[[j]]$SFS_Non))]<-dat[[j]]$SFS_Non/sum(dat[[j]]$SFS_Non)
		com_snpSFS_Syn[j,as.numeric(names(dat[[j]]$SFS_Syn))]<-dat[[j]]$SFS_Syn/sum(dat[[j]]$SFS_Syn)
	} 
	
	toPlot<-rbind(colMeans(com_snpSFS_Syn),colMeans(com_snpSFS_Non))
	#FOLD it!
	toPlot<-toPlot[,1:indivs_sampled]+toPlot[,(2*indivs_sampled):(indivs_sampled+1)]

	barplot(toPlot,xlab="",ylab="",main="",col=c("steelblue3","coral3"),ylim=c(0,0.6),beside=TRUE,border=NA,space=c(0,0),axes=FALSE)->mids
	box(col="gray70")
	

	
	#pretty Plotting
	as.numeric(gsub("_.+","",names(fixationData)[i]))->m
	as.numeric(gsub(".+_","",names(fixationData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Allele Frequency",side=1,line=5,cex=1.02)
		axis(1,at=colMeans(mids)[seq(from=1,to=indivs_sampled,by=6)],labels=seq(from=1,to=indivs_sampled,by=6),las=2,line=1)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*Ne*m),2),sep=""),side=2,line=6,cex=1.4)
		mtext("Proportion of alleles",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),las=2,line=1)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	
	if(i==1)legend("topleft",legend=c("Unconstrained","Selected"),text.col=c("steelblue3","coral3"),cex=1.0,bty="n")

}




#Plot unfolded, truncated SFS
plotdim<-ceiling(sqrt(length(polymorphismData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)
par(mfrow=ceiling(c(sqrt(length(polymorphismData)),sqrt(length(polymorphismData)))))
for(i in 1:length(polymorphismData)){

	dat<-polymorphismData[[i]]
	com_snpSFS_Non<-matrix(0,nrow=length(dat),ncol=2*indivs_sampled)
	com_snpSFS_Syn<-matrix(0,nrow=length(dat),ncol=2*indivs_sampled)
	for(j in 1:length(dat)){
		com_snpSFS_Non[j,as.numeric(names(dat[[j]]$SFS_Non))]<-dat[[j]]$SFS_Non/sum(dat[[j]]$SFS_Non)
		com_snpSFS_Syn[j,as.numeric(names(dat[[j]]$SFS_Syn))]<-dat[[j]]$SFS_Syn/sum(dat[[j]]$SFS_Syn)
	} 
	
	toPlot<-rbind(colMeans(com_snpSFS_Syn),colMeans(com_snpSFS_Non))[,1:10]
	#Don't FOLD it!
	#toPlot<-toPlot[,1:nPops]+toPlot[,(2*nPops):(nPops+1)]

	barplot(toPlot,xlab="",ylab="",main="",col=c("steelblue3","coral3"),ylim=c(0,0.45),beside=TRUE,border=NA,space=c(0,0.5),axes=FALSE)->mids
	box(col="gray70")
	

	
	#pretty Plotting
	as.numeric(gsub("_.+","",names(fixationData)[i]))->m
	as.numeric(gsub(".+_","",names(fixationData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Allele Count (truncated)",side=1,line=5,cex=1.02)
		axis(1,at=colMeans(mids),labels=1:10,las=1,line=1)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*Ne*m),2),sep=""),side=2,line=6,cex=1.4)
		mtext("Proportion of alleles",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),las=2,line=1)		
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	
	if(i==1)legend("topleft",legend=c("Unconstrained","Selected"),text.col=c("steelblue3","coral3"),cex=1.0,bty="n")

}



#calculate omega_a and alpha, Ka and Ks, based on the actual number of fixations

obs_alpha<-obs_omega_a<-Ka<-Ks<-pis_sample<-pia_het<-pia_sample<-pis_het<-Fst_syn<-Fst_non<-rep(0,length(polymorphismData))
ms<-rep(0,length(polymorphismData))
ss<-rep(0,length(polymorphismData))

for(i in 1:length(polymorphismData)){
	as.numeric(gsub("_.+","",names(polymorphismData)[i]))->ms[i]
	as.numeric(gsub(".+_","",names(polymorphismData)[i]))->ss[i]
	Pdat<-polymorphismData[[i]]
	Fdat<-fixationData[[i]]
	synFix<-0
	nonFix<-0
	posFix<-0
	sP<-vector()
	nP<-vector()
	for(j in 1:length(Pdat)){

		sP<-c(sP,rep(as.numeric(names(Pdat[[j]]$SFS_Syn[-alleles_sampled]))/alleles_sampled,Pdat[[j]]$SFS_Syn[-alleles_sampled]))
		nP<-c(nP,rep(as.numeric(names(Pdat[[j]]$SFS_Non[-alleles_sampled]))/alleles_sampled,Pdat[[j]]$SFS_Non[-alleles_sampled]))
		synFix<-synFix+length(Fdat[[j]]$SynFixTime) #number of fixations,
		nonFix<-nonFix+length(Fdat[[j]]$DelFixTime)+length(Fdat[[j]]$PosFixTime)  #number of fixations,
		#print(c(j,nonFix))
		posFix<-posFix+length(Fdat[[j]]$PosFixTime)
	}
	total_syn_Sites<-round(length(Pdat)*5000 * 50 * 0.24,0) #combine all simulations into one analysis - so multiply length by number of reps
	total_non_Sites<-round(length(Pdat)*5000 * 50 * 0.76,0)	
	posFix/nonFix->obs_alpha[i] 
	round((posFix/total_non_Sites)/(synFix/total_syn_Sites),3)->obs_omega_a[i]

	nonFix/total_non_Sites->Ka[i]
	synFix/total_syn_Sites->Ks[i]

	################################################################################################
	##################### THESE are not sample piS - they are single-fly PiS
	################################################################################################
	
	pis_sample[i]<-sum(1-(sP^2)-((1-sP)^2))/total_syn_Sites
	pis_het[i]<-mean(colMeans(convergenceData[[i]]$Spi_s[,-(1:500)],na.rm=TRUE))
	Fst_syn[i]<-(pis_sample[i]-pis_het[i])/pis_sample[i]
	
	pia_sample[i]<-sum(1-(nP^2)-((1-nP)^2))/total_non_Sites
	pia_het[i]<-mean(colMeans(convergenceData[[i]]$Spi_a[,-(1:500)],na.rm=TRUE))	
	Fst_non[i]<-(pia_sample[i]-pia_het[i])/pia_sample[i]
	
}

#Run assymptoticMK

plotdim<-ceiling(sqrt(length(polymorphismData)))	
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)


#all replicates together (20 slices)
overall_results3<-vector()
Ka_apparant<-Ks_apparant<-rep(0,length(polymorphismData))
for(i in 1:length(polymorphismData)){
	Pdat<-polymorphismData[[i]]
	Fdat<-fixationData[[i]]
	synFix<-0
	nonFix<-0
	posFix<-0
	sP<-vector()
	nP<-vector()
	for(j in 1:length(Pdat)){
		print(paste("synFix is ",synFix,sep=""))
		if(!is.na(Pdat[[j]]$SFS_Syn[as.character(alleles_sampled)])){ #if there is an allele that appears fixed in the sample, 
			print(paste("Adding ",Pdat[[j]]$SFS_Syn[as.character(alleles_sampled)],sep=""))
			synFix<-synFix+Pdat[[j]]$SFS_Syn[as.character(alleles_sampled)] # it needs to be added to the fixed differences
			print(paste("synFix is now ",synFix,sep=""))
			Pdat[[j]]$SFS_Syn<-Pdat[[j]]$SFS_Syn[-length(Pdat[[j]]$SFS_Syn)]; #and removed from the polymorphisms
		}
		if(!is.na(Pdat[[j]]$SFS_Non[as.character(alleles_sampled)])){ #if there is an allele that appears fixed in the sample
			nonFix<-nonFix+Pdat[[j]]$SFS_Non[as.character(alleles_sampled)]
			Pdat[[j]]$SFS_Non<-Pdat[[j]]$SFS_Non[-length(Pdat[[j]]$SFS_Non)];
		}
		sP<-c(sP,rep(as.numeric(names(Pdat[[j]]$SFS_Syn[-alleles_sampled]))/alleles_sampled,Pdat[[j]]$SFS_Syn[-alleles_sampled]))
		nP<-c(nP,rep(as.numeric(names(Pdat[[j]]$SFS_Non[-alleles_sampled]))/alleles_sampled,Pdat[[j]]$SFS_Non[-alleles_sampled]))
		synFix<-synFix+length(Fdat[[j]]$SynFixTime)
		nonFix<-nonFix+length(Fdat[[j]]$DelFixTime)+length(Fdat[[j]]$PosFixTime)
		#print(c(j,nonFix))
		posFix<-posFix+length(Fdat[[j]]$PosFixTime)
	}
	freq_bin<-0.05+(0:9)/10
	synPol<-hist(sP,breaks=c(0:10)/10,plot=FALSE)$count
	nonPol<-hist(nP,breaks=c(0:10)/10,plot=FALSE)$count	
	asymptoticMK(d0=synFix, d=nonFix, df=data.frame(f=freq_bin, p=nonPol, p0=synPol, row.names=NULL), xlow=0.01, xhigh=0.99, output="pdf")->tmp
	overall_results3<-rbind(overall_results3,tmp)
	abline(h=posFix/nonFix,lwd=1,col="red",lty=3)

	nonFix/total_non_Sites->Ka_apparant[i]
	synFix/total_syn_Sites->Ks_apparant[i]	
	
	#pretty Plotting
	as.numeric(gsub("_.+","",names(fixationData)[i]))->m
	as.numeric(gsub(".+_","",names(fixationData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Derived Allele Frequency",side=1,line=5,cex=1.02)
		axis(1,at=axTicks(1),las=1,line=1)
	}
	if(i%%plotdim==1){
		mtext(paste("Fst = ",round(1/(1+4*Ne*m),2),sep=""),side=2,line=6,cex=1.4)
		mtext("Estimated alpha",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),las=2,line=1)		
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	
	
	
}

