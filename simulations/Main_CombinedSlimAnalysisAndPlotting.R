

setwd("G:/Publications/- 2021_Peterfields_MK/SLiM/Simulations_March2021")

###################################################################
#
#   _____                                                   
#  / ____|                                                  
# | |     ___  _ ____   _____ _ __ __ _  ___ _ __   ___ ___ 
# | |    / _ \| '_ \ \ / / _ \ '__/ _` |/ _ \ '_ \ / __/ _ \
# | |___| (_) | | | \ V /  __/ | | (_| |  __/ | | | (_|  __/
#  \_____\___/|_| |_|\_/ \___|_|  \__, |\___|_| |_|\___\___|
#                                  __/ |                    
#                                 |___/                     
#
###################################################################
#
# In which we look at output through the run to check for equilibrium
#
###################################################################


###############################################################################
#
# Read convergence files
#
###############################################################################

#unknown number of files (replicates and conditions) to parse. Simply list them
list.files("Convergence",pattern="convergence_.+txt")->fnames

#create a list to contain the data
convergenceData<-list()

#for each file
for(fname in fnames){
	#simulation parameters are in the file name (migration_sexFrequncy_RunID)
	unlist(strsplit(fname,"_"))[2:3]->sim_pars
	paste(sim_pars,collapse="_")->sp

	#Pull out those lines that start 'Pars'
	#Tab-separated fields are: 
	#	RunID, Generation, migration, SexFrequency, 
	# 	"pi", pi_s_local, pi_a_local,pi_s_global, pi_a_global,Fst,
	#	"D", Synonymous Fixed Differences, Beneficial Fixed Differences, Deleterious or Neutral Fixed differences, Synonymous Polymorphisms, Deleterious or Neutral Fixed Polymorphisms, Beneficial Polymorphisms 

	#read the data as a list of vectors, in case of unexpected row lengths
	data_list<-strsplit(grep("Pars",scan(paste("Convergence/",fname,sep=""),what="character",sep="\n"),value=TRUE),"\t")
	#turn it into a table
	do.call(rbind,data_list)->dat

	if(is.null(convergenceData[[sp]])){
	#if we have not seen this set of parameters before (i.e. first replicate), create an entry in the list
		convergenceData[[sp]]<-list(Spi_a=as.numeric(dat[,8]),Spi_s=as.numeric(dat[,7]),Tpi_a=as.numeric(dat[,10]),Tpi_s=as.numeric(dat[,9]),fst=as.numeric(dat[,11]),Ds=as.numeric(dat[,13])/1000,Dn=as.numeric(dat[,14])/1000,Dna=as.numeric(dat[,15])/1000)
	}else{
	#if we have seen this set of parameters before (i.e. subsequent replicate), add the replicate on to the bottom of the entry for those parameters
		convergenceData[[sp]]$Spi_a<-rbind(convergenceData[[sp]]$Spi_a,as.numeric(dat[,8]))
		convergenceData[[sp]]$Spi_s<-rbind(convergenceData[[sp]]$Spi_s,as.numeric(dat[,7]))
		convergenceData[[sp]]$Tpi_a<-rbind(convergenceData[[sp]]$Tpi_a,as.numeric(dat[,10]))
		convergenceData[[sp]]$Tpi_s<-rbind(convergenceData[[sp]]$Tpi_s,as.numeric(dat[,9]))
		convergenceData[[sp]]$fst<-rbind(convergenceData[[sp]]$fst,as.numeric(dat[,11]))
		convergenceData[[sp]]$Ds<-rbind(convergenceData[[sp]]$Ds,as.numeric(dat[,13])/1000) 
		convergenceData[[sp]]$Dn<-rbind(convergenceData[[sp]]$Dn,as.numeric(dat[,14])/1000)
		convergenceData[[sp]]$Dna<-rbind(convergenceData[[sp]]$Dna,as.numeric(dat[,15])/1000)
	}

}
str(convergenceData)

###############################################################################
#
# set up plotting area to plot local pi_a and pi_s
#
###############################################################################

#pdf output for plotting
pdf(file = "Figures/Slim_Convergence_Figures.pdf", width=8, height=6) 

#plot both pi[s] and pi[a] for sub-pops
plotdim<-ceiling(sqrt(length(convergenceData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)

mu<-1.53e-7		# scaled mutation rate, as used in the simulations. True value 8.96e-9
burnin<-300 	# in 1000's of generations; Statistics seem to converge by this point 

#for each parameter combination
for(i in 1:length(convergenceData)){

	if(is.vector(convergenceData[[i]]$Spi_s)){
	#Catch R's special case of a 1D matrix being a vector	(i.e. 1 replicate)
		upperx<-max(sapply(convergenceData,function(x){length(x[[1]])}))
		uppery<-max(
			sapply(convergenceData,function(x){max(x$Spi_s,na.rm=TRUE)}),
			sapply(convergenceData,function(x){max(x$Spi_a,na.rm=TRUE)})
		)
		convergenceData[[i]]$Spi_s->x;x<-x[!is.na(x)]
		convergenceData[[i]]$Spi_a->z;z<-z[!is.na(z)]
		plot(x,xlab="",ylab="",main="",xlim=c(1,upperx),ylim=c(0,uppery),type="n",axes=FALSE)
		points(smooth(x),type="l",col="steelblue4")		
		points(smooth(z),type="l",col="coral4")
		
	}else{
	#If its a matrix ... (i.e. more than one replicate)
		#find the upper generation number
		upperx<-max(sapply(convergenceData,function(x){ncol(x[[1]])}))
		#find the upper diversity value
		uppery<-max(
			sapply(convergenceData,function(x){max(x$Spi_s,na.rm=TRUE)}),
			sapply(convergenceData,function(x){max(x$Spi_a,na.rm=TRUE)})
		)
		#blank plot
		plot(0,xlab="",ylab="",main="",xlim=c(1,upperx),ylim=c(0,uppery),type="n",axes=FALSE)
		#plot all of the replicates (transparency, thin)
		for(j in 1:(min(10,nrow(convergenceData[[i]]$Spi_s)))){
			points(convergenceData[[i]]$Spi_s[j,],type="l",col=adjustcolor("steelblue4", alpha.f = 0.2))
		}
		for(j in 1:(min(10,nrow(convergenceData[[i]]$Spi_s)))){
			points(convergenceData[[i]]$Spi_a[j,],type="l",col=adjustcolor("coral4", alpha.f = 0.2))
		}
		#calc mean diversity across replicates
		colMeans(convergenceData[[i]]$Spi_s,na.rm=TRUE)->x
		colMeans(convergenceData[[i]]$Spi_a,na.rm=TRUE)->z
		#plot means
		points(smooth(x),type="l",col="steelblue4",lwd=2)
		points(smooth(z),type="l",col="coral4",lwd=2)
		
	}
	box(col="gray70")
	
	# pi[s] is actual; synonymous heterozygosity in a single individual. If mating is at random, this is an estimator of global Ne
	 Ne<-mean(x[-(1:burnin)])/(4*mu) 
	 if(mean(x[-(1:burnin)])<0.008){legend("topright",legend=paste("Ne = ",signif(Ne,3)," ",sep=""), bty="n")}
	 else{legend("bottomright",legend=paste("Ne = ",signif(Ne,3)," ",sep=""), bty="n")}
	
	#pretty Plotting of axes
	as.numeric(gsub("_.+","",names(convergenceData)[i]))->m
	as.numeric(gsub(".+_","",names(convergenceData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Generations (x 1000)",side=1,line=5,cex=1.00)
		axis(1,at=axTicks(1),line=1)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*500*m),2),sep=""),side=2,line=6,cex=1.4) #Assuming inf island model and local Ne of 500
		mtext("Local pi_s (%)",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),labels=round(100*axTicks(2),2),las=2,line=1)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	

	if(i==9)legend("topleft",legend=c("unconstrained","selected"),text.col=c("steelblue4","coral4"),cex=1.0,bty="n")
	
	#parameter estimates, to check
	#print(signif(c(round(1/(1+4*Ne*m),2),s,mean(x[-(1:500)])),2))
}

###############################################################################
#
# set up plotting area to plot empirical Fst
#
###############################################################################

plotdim<-ceiling(sqrt(length(convergenceData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)

mu<-1.53e-7		# scaled mutation rate, as used in the simulations. True value 8.96e-9
burnin<-300 	# in 1000's of generations; choose temporary value for inferred Ne 

#for each parameter combination
for(i in 1:length(convergenceData)){

	if(is.vector(convergenceData[[i]]$Spi_s)){
	#Catch R's special case of a 1D matrix being a vector	(i.e. 1 replicate)
		upperx<-max(sapply(convergenceData,function(x){length(x[[1]])}))
		uppery<-1
		convergenceData[[i]]$fst->x;x<-x[!is.na(x)]
		(convergenceData[[i]]$Tpi_s-convergenceData[[i]]$Spi_s)/convergenceData[[i]]$Tpi_s->y
		(convergenceData[[i]]$Tpi_a-convergenceData[[i]]$Spi_a)/convergenceData[[i]]$Tpi_a->z
		
		plot(x,xlab="",ylab="",main="",xlim=c(1,upperx),ylim=c(0,uppery),type="n",axes=FALSE)
		points(smooth(x),type="l",col="black")	
		points(smooth(y),type="l",col="blue")	
		points(smooth(z),type="l",col="red")	
		
	}else{
	#If its a matrix ... (i.e. more than one replicate)
		convergenceData[[i]]$fst->x

		#find the upper generation number
		upperx<-max(sapply(convergenceData,function(x){ncol(x[[1]])}))
		#find the upper diversity value
		uppery<-1
		#blank plot
		plot(0,xlab="",ylab="",main="",xlim=c(1,upperx),ylim=c(0,uppery),type="n",axes=FALSE)
		#plot all of the replicates (transparency, thin)
		for(j in 1:(min(10,nrow(x)))){
			points(x[j,],type="l",col=adjustcolor("black", alpha.f = 0.1))

		}

		#plot means
		points(smooth(colMeans(x)),type="l",col="black",lwd=2)

		
	}
	box(col="gray70")
	#Write empirical Fst - 
	Fst<-mean(colMeans(x)[-(1:burnin)])
	legend("topright",legend=paste("F_it = ",signif(Fst,2),sep=""), bty="n")# The nominal 'F_st' measured is actually F_it, if there is inbreeding (such as clone selfing)
	
	#pretty Plotting of axes
	as.numeric(gsub("_.+","",names(convergenceData)[i]))->m
	as.numeric(gsub(".+_","",names(convergenceData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Generations (x 1000)",side=1,line=5,cex=1.00)
		axis(1,at=axTicks(1),line=1)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(Fst,2),sep=""),side=2,line=6,cex=1.4) 
		mtext("F_it",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),labels=round(axTicks(2),2),las=2,line=1)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	

	
	#parameter estimates, to check
	#print(signif(c(round(1/(1+4*Ne*m),2),s,mean(x[-(1:500)])),2))
}

###############################################################################
#
# Plotting area to plot fixations over time
# Logic as above 
#
###############################################################################

#plot Ddel and Dpos and Dsyn
par(mfrow=c(sqrt(length(convergenceData)),sqrt(length(convergenceData))))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)	
for(i in 1:length(convergenceData)){

	if(is.vector(convergenceData[[i]]$Spi_s)){
		upperx<-max(sapply(convergenceData,function(x){length(x[[1]])}))
		uppery<-max(
			sapply(convergenceData,function(x){max(x$Ds,na.rm=TRUE)}),
			sapply(convergenceData,function(x){max(x$Dn,na.rm=TRUE)}),
			sapply(convergenceData,function(x){max(x$Dna,na.rm=TRUE)})
		)
		convergenceData[[i]]$Ds->S;S<-S[!is.na(S)]
		convergenceData[[i]]$Dn->N;N<-N[!is.na(N)]
		convergenceData[[i]]$Dna->A;A<-A[!is.na(A)]
		plot(x,xlab="Generation x 1000",ylab="Cumulative Fixations",main=names(convergenceData[i]),xlim=c(1,upperx),ylim=c(0,uppery),type="n")
		
		points(smooth(S),type="l",col="steelblue4")		
		points(smooth(N),type="l",col="coral4")
		points(smooth(A),type="l",col="olivedrab4")
	}else{
		upperx<-max(sapply(convergenceData,function(x){ncol(x[[1]])}))
		uppery<-max(
			sapply(convergenceData,function(x){max(x$Ds,na.rm=TRUE)}),
			sapply(convergenceData,function(x){max(x$Dn,na.rm=TRUE)}),
			sapply(convergenceData,function(x){max(x$Dna,na.rm=TRUE)})
		)
		plot(0,xlab="",ylab="",main="",xlim=c(1,upperx),ylim=c(0,uppery),type="n",axes=FALSE)
		#plot(x,xlab="Generation x 1000",ylab="Cumulative Fixations",main=names(convergenceData[i]),xlim=c(1,upperx),ylim=c(0,uppery),type="n")
		for(j in 1:(min(10,nrow(convergenceData[[i]]$Ds)))){
			points(convergenceData[[i]]$Ds[j,],type="l",col="lightsteelblue3")
		}
		for(j in 1:(min(10,nrow(convergenceData[[i]]$Dn)))){
			points(convergenceData[[i]]$Dn[j,],type="l",col="bisque3")
		}
		for(j in 1:(min(10,nrow(convergenceData[[i]]$Dna)))){
			points(convergenceData[[i]]$Dna[j,],type="l",col="olivedrab3")
		}
		
		colMeans(convergenceData[[i]]$Ds,na.rm=TRUE)->S
		colMeans(convergenceData[[i]]$Dn,na.rm=TRUE)->N
		colMeans(convergenceData[[i]]$Dna,na.rm=TRUE)->A
		
		points(smooth(S),type="l",col="steelblue4",lwd=2)
		points(smooth(N),type="l",col="coral4",lwd=2)
		points(smooth(A),type="l",col="olivedrab4",lwd=2)
	}
	box(col="gray70")
	
	#pretty Plotting
	as.numeric(gsub("_.+","",names(convergenceData)[i]))->m
	as.numeric(gsub(".+_","",names(convergenceData)[i]))->s
	if(i>((plotdim-1)*plotdim)){
		mtext("Generations (x 1000)",side=1,line=5,cex=1.00)
		axis(1,at=axTicks(1),line=1)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*500*m),2),sep=""),side=2,line=6,cex=1.4) #assuming inf island model and simulated local Ne of 500
		mtext("No. fixed (x1000)",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),las=2,line=1)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	

	if(i==1)legend("topleft",legend=c("Unconstrained","Weakly Deleterious","Beneficial"),text.col=c("steelblue4","coral4","olivedrab4"),col=c("steelblue4","coral4","olivedrab4"),cex=1.0,bty="n",lty=1)
}

### Stop plotting convergence!
dev.off()


###################################################################
#   
#  ______ _           _   _                 
# |  ____(_)         | | (_)                
# | |__   ___  ____ _| |_ _  ___  _ __  ___ 
# |  __| | \ \/ / _` | __| |/ _ \| '_ \/ __|
# | |    | |>  < (_| | |_| | (_) | | | \__ \
# |_|    |_/_/\_\__,_|\__|_|\___/|_| |_|___/
#                                           
#                                           
###################################################################
#
# In which we count the number of fixations occurring since equilibrium is reached
#
###################################################################



#Select a burn-in threshold on the basis of the convergence analysis. 
# Pi appears to be stationary after this point
# Fixation rate (fixations per time) appears linear after this point
# SFS may not be at equilibrium, but we will treat this as an acceptable proxy
BurninThreshold<-300 



#get file names
list.files("Fixations", pattern="fixations_*")->fnames

#empty list to hold data
fixationData<-list()
#for each file 
for (i in 1:length(fnames)){
	fname<-fnames[i]

	#extract key parameters from the file name
	unlist(strsplit(fname,"_"))[2:3]->sim_pars
	paste(sim_pars,collapse="_")->sp
	
	#read in the data file, and construct a list of the mutations in it
	#columns are mutation number,  MutationID,Mutation type, location, selection coef, dominance, deme of origin, generation arose, generation fixed
	read.table(paste("Fixations/",fname,sep=""),skip=2,as.is=TRUE)->dat
	
	#only retain those that arose after diversity reached equilibrium 
	dat<-dat[dat[,8]>BurninThreshold,]

	#find selection coefficients and fixation times
	outdata<-list(DelFixDFE=dat[dat[,3]=="m2",5],
				SynFixTime=dat[dat[,3]=="m1",9]-dat[dat[,3]=="m1",8],
				DelFixTime=dat[dat[,3]=="m2",9]-dat[dat[,3]=="m2",8],
				PosFixTime=dat[dat[,3]=="m3",9]-dat[dat[,3]=="m3",8]
				)
						
	if(is.null(fixationData[[sp]])){
	#if we have not yet seen a replicate with this parameter combination, create a new entry
		fixationData[[sp]]<-list()
		fixationData[[sp]][[1]]<-outdata
	}else{
	#if we have seen it, add the replicate	
		fixationData[[sp]][[length(fixationData[[sp]])+1]]<-outdata
	}
	
}


###############################################################################
#
# Plot Fixation time
#
###############################################################################

pdf(file = "Figures/Slim_Fixation_Times.pdf",width=8, height=6) 



#find plot limits
min_FixTime<-min_density<-1e20
max_FixTime<-min_density<-0
for(i in 1:length(fixationData)){
	dat<-fixationData[[i]]
	for(j in 1:length(dat)){
		min_FixTime<-min(c(min_FixTime,dat[[j]]$SynFixTime,dat[[j]]$DelFixTime,dat[[j]]$PosFixTime))
		max_FixTime<-max(c(max_FixTime,dat[[j]]$SynFixTime,dat[[j]]$DelFixTime,dat[[j]]$PosFixTime))
	}
}
min_FixTime;max_FixTime

#set plotting area
plotdim<-ceiling(sqrt(length(fixationData)))
par(
	mfrow=c(plotdim,plotdim),
	mar=c(2,1,0,1),
	oma=c(6,8,4,1)
)

#for each parameter combination
for(i in 1:length(fixationData)){
	#find key parameters
	as.numeric(gsub("_.+","",names(fixationData)[i]))->m
	as.numeric(gsub(".+_","",names(fixationData)[i]))->s
	N<-500
	n<-36
	
	#empirical Ne. This requires information from Convergence data!
	colMeans(convergenceData[[i]]$Spi_s,na.rm=TRUE)->x;
	mean(x[-(1:500)])/(4*1.53e-7)->EmpNe

	#data for this parameter combination
	dat<-fixationData[[i]]

	# Obtain nice KDE density plots for fixation time (on a log10 scale)
	com_SynFixTime<-vector()
	for(j in 1:length(dat)){com_SynFixTime<-c(com_SynFixTime,dat[[j]]$SynFixTime)}
	D_com_SynFixTime<-density(log10(com_SynFixTime),cut=TRUE)
	com_DelFixTime<-vector()
	for(j in 1:length(dat)){com_DelFixTime<-c(com_DelFixTime,dat[[j]]$DelFixTime)}
	D_com_DelFixTime<-density(log10(com_DelFixTime),cut=TRUE)	
	com_PosFixTime<-vector()
	for(j in 1:length(dat)){com_PosFixTime<-c(com_PosFixTime,dat[[j]]$PosFixTime)}
	D_com_PosFixTime<-density(log10(com_PosFixTime),cut=TRUE)
	
	plot(D_com_PosFixTime,xlab="",ylab="",type="n",xlim=log10(c(min_FixTime,max_FixTime)),ylim=c(0,7),main="",axes=FALSE)
	

	#expected fixation time of neutral allele, based on empirical Ne 
	abline(v=log10(4*EmpNe),lty=2,col=adjustcolor("black", alpha.f = 0.5))

	# plot density and means
	abline(v=log10(mean(com_SynFixTime)),col="steelblue4",lty=1)
	polygon(c(min(D_com_SynFixTime$x),D_com_SynFixTime$x,max(D_com_SynFixTime$x)),c(0,D_com_SynFixTime$y,0),col=adjustcolor("steelblue4", alpha.f = 0.5),border="steelblue4")
	abline(v=log10(mean(com_PosFixTime)),col="olivedrab4",lty=1)
	polygon(c(min(D_com_PosFixTime$x),D_com_PosFixTime$x,max(D_com_PosFixTime$x)),c(0,D_com_PosFixTime$y,0),col=adjustcolor("olivedrab4", alpha.f = 0.5),border="olivedrab4")
	abline(v=log10(mean(com_DelFixTime)),col="coral4",lty=1)
	polygon(c(min(D_com_DelFixTime$x),D_com_DelFixTime$x,max(D_com_DelFixTime$x)),c(0,D_com_DelFixTime$y,0),col=adjustcolor("coral4", alpha.f = 0.5),border="coral4")
	box(col="gray70")
	

	#Add legends and labels where appropriate
	if(i>((plotdim-1)*plotdim)){
		mtext("Fixation time\n(x 1000 generations)",side=1,line=5,cex=1.00)
		axis(1,at=log10(signif(10^axTicks(1),1)),labels=signif(10^axTicks(1),1)/1000)
	}
	if(i%%plotdim==1){
		mtext(paste("F_st = ",round(1/(1+4*N*m),2),sep=""),side=2,line=6,cex=1.4)
		mtext("Density",side=2,line=4,cex=1.0)
		axis(2,at=axTicks(2),las=2,line=1)	
	}
	if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.4,line=2)}	
	if(i==9)legend("topright",legend=c("Unconstrained","Weakly deleterious","Beneficial", "4Ne Generations"),text.col=c("steelblue4","coral4","olivedrab4","gray40"),cex=1.0,bty="n")
}
dev.off()

###############################################################################
#
# Tabulation of fixation data for subsequence MK-like analyses
#
###############################################################################

fixations<-vector()
r<-1
for(i in 1:length(fixationData)){
	dat<-fixationData[[i]]
	for(j in 1:length(dat)){
		fixations<-rbind(fixations,c(as.numeric(unlist(strsplit(names(fixationData[i]),"_"))),r,length(dat[[j]]$SynFixTime),length(dat[[j]]$DelFixTime),length(dat[[j]]$PosFixTime),length(dat[[j]]$DelFixTime)+length(dat[[j]]$PosFixTime)))
		r<-r+1
	}
}
colnames(fixations)<-c("migration","sex_freq","Rep","Syn","Nonsyn_del","Nonsyn_pos","Nonsyn")
write.table(fixations,file="SimulatedeFixationData.tsv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)






####################################################################################
#
#  _____      _                                  _     _                   
# |  __ \    | |                                | |   (_)                  
# | |__) |__ | |_   _ _ __ ___   ___  _ __ _ __ | |__  _ ___ _ __ ___  ___ 
# |  ___/ _ \| | | | | '_ ` _ \ / _ \| '__| '_ \| '_ \| / __| '_ ` _ \/ __|
# | |  | (_) | | |_| | | | | | | (_) | |  | |_) | | | | \__ \ | | | | \__ \
# |_|   \___/|_|\__, |_| |_| |_|\___/|_|  | .__/|_| |_|_|___/_| |_| |_|___/
#                __/ |                    | |                              
#               |___/                     |_|                              
#
####################################################################################

								
#Define convenience function to read a mutations file
read_SNPS<-function(fname){
	unlist(strsplit(fname,"_"))[2:3]->sim_pars
	paste(sim_pars,collapse="_")->sp

	#read in the data file, and construct a list of the mutations in it
	f <- file(fname)
	open(f,"r")
	ReadMutations<-FALSE
	outdata<-list()
	#I could not think of a good way to parse this file, so this will have to do!
	while(length(nline <- readLines(f,n=1)) > 0) {

		if(nline=="Mutations:"){
			#set reading flag on and re-enter loop
			ReadMutations<-TRUE
			next
		}
		if(nline=="Genomes:"){
			#set reading flag off
			ReadMutations<-FALSE
			next
		}
		if(ReadMutations){
			#process mutation lines
			temp<-unlist(strsplit(nline," +"))
			# (1) a temporary unique identifier used only in this output section, 
			# (2) a permanent unique identifier kept by SLiM throughout a run, 
			# (3) the mutation type, 
		#1	# (4) the base position, 
		#2	# (5) the selection coefficient (here always 0 since this is a neutral model), 
		#3	# (6) the dominance coefficient (here always 0.5), 
			# (7) the identifier of the subpopulation in which the mutation first arose, 
		#4	# (8) the generation in which it arose, 
			# (9) the prevalence of the mutation (the number of genomes that contain the mutation (genome=haploid genome)
			
			#if we've not seen this mutation before
			if(is.null(outdata[[temp[2]]])){
				#make an entry
				outdata[[temp[2]]]<-list(type=temp[3],count=as.numeric(temp[9]),info=as.numeric(temp[c(4,5,6,8)]))
			}else{ #if we have seen this mutation before
				#update the increment by the number of copies in this individual
				outdata[[temp[2]]]$count<-outdata[[temp[2]]]$count+as.numeric(temp[9])
			}
		}
	}
	close(f)
	return(outdata)
}

###################################################################
#
# Part I: In which we count the number of segregating polymorphisms at the end of the simulation, based on a single individual sampled per deme
#
###################################################################



###############################################################################
#
# Read in the files for a single individual per deme
#
###############################################################################



indivs_sampled<-36 #(one diploid from each of 36 populations)

list.files("Mutations_1from36", pattern="mutations_0.[0-9].+.txt")->fnames

final_generation<-1000000

polymorphismData<-list()
for (i in 1:length(fnames)){

	#final run data 
	outdata<-read_SNPS(paste("Mutations_1from36/",fnames[i],sep=""))
		
	unlist(strsplit(fnames[i],"_"))[2:3]->sim_pars
	paste(sim_pars,collapse="_")->sp
	
	#eventuaally need to  strip out those alleles that are fixed in our sample
	#these are present because they're segregating in the population
	# , but not yet - need to know they appear fixed
	#sapply(outdata,function(x){x$count!=2*indivs_sampled})->not_fixed
	#outdata<-outdata[not_fixed]
	
	#those which are synonymous
	sapply(outdata,function(x){x$type=="m1"})->syn

	#those which are non-synonymous
	sapply(outdata,function(x){x$type!="m1"})->non

	#those which are negatively selected
	sapply(outdata,function(x){x$type=="m2"})->del

	#those which are positively selected
	sapply(outdata,function(x){x$type=="m3"})->pos	

	#SFS count frequency with table, but 1:(2*indivs_sampled) is to force inclusion of frequencies that are not seen
	table(c(sapply(outdata,function(x){x$count})[syn],1:(2*indivs_sampled)))-1->SFS_Syn
	table(c(sapply(outdata,function(x){x$count})[non],1:(2*indivs_sampled)))-1->SFS_Non
	
	#age distributions 
	sapply(outdata,function(x){final_generation-x$info}[4])[syn]->SynAge
	sapply(outdata,function(x){final_generation-x$info}[4])[del]->DelAge
	sapply(outdata,function(x){final_generation-x$info}[4])[pos]->PosAge
	
	#DFE of mutations  
	sapply(outdata,function(x){x$info}[2])[del]->snpDFE_Non

	#add that list to a list of all polymorphism data seen across all datasets
	if(is.null(polymorphismData[[sp]])){
		polymorphismData[[sp]][[1]]<-list(SFS_Syn=SFS_Syn,SFS_Non=SFS_Non,SynAge=SynAge,DelAge=DelAge,PosAge=PosAge,snpDFE_Non=snpDFE_Non)
	}else{
		polymorphismData[[sp]][[length(polymorphismData[[sp]])+1]]<-list(SFS_Syn=SFS_Syn,SFS_Non=SFS_Non,SynAge=SynAge,DelAge=DelAge,PosAge=PosAge,snpDFE_Non=snpDFE_Non)
	}
	print(paste("Finished reading file ",i,sep=""))	
}

polymorphismData->polymorphismData_1from36

#polymorphism counts
polymorphisms<-vector()
r<-1
for(i in 1:length(polymorphismData)){
	dat<-polymorphismData[[i]]
	for(j in 1:length(dat)){
		polymorphisms<-rbind(polymorphisms,c(as.numeric(unlist(strsplit(names(polymorphismData[i]),"_"))),r,length(dat[[j]]$SynAge),length(dat[[j]]$DelAge),length(dat[[j]]$PosAge),length(dat[[j]]$DelAge)+length(dat[[j]]$PosAge)))
		r<-r+1
	}
}
colnames(polymorphisms)<-c("migration","sex_freq","Rep","Syn","Non_del","Non_pos","Non")

polymorphisms->polymorphisms_1from36

####################################################################################
#  
# plotting
#
####################################################################################
indivs_sampled<-36
pdf(file = "Figures/SamplePolymorphism_1Indiv_from_36Demes.pdf",width=8, height=6) 
	alleles_sampled<-36*2
	polymorphismData<-polymorphismData_1from36
	source("AssymptoticMK.R")
	source("PlotPolymorphism_And_PlotAsymptoticMK.R")	
dev.off()

AsymptoticMK_results_1from36<-data.frame(Migration=ms,Sex=ss,obs_alpha=obs_alpha,alpha_original=overall_results3$alpha_original,alpha_assymptotic=overall_results3$alpha_asymptotic,CI_low=overall_results3$CI_low,CI_high=overall_results3$CI_high)
AsymptoticMK_results_1from36<-cbind(AsymptoticMK_results_1from36,obs_omega_a,Ks,Ks_apparant,Ka,Ka_apparant,pis_sample,pis_het,Fst_syn,pia_sample,pia_het,Fst_non)





####################################################################################
#
# Part II: In which we count the number of segregating polymorphisms at the end of the simulation, based on 36 individuals sampled from a single deme
#
###################################################################

#Select a burnin threshold on the basis of the convergence analysis. 
# Pi appears to be stationary after this point
# Fixation rate (fixations per time) appears linear after this point


list.files("Mutations_36from1", pattern="mutations_0.[0-9].+.txt")->fnames

indivs_sampled<-36 #(36 diploids from one populations)

final_generation<-1000000

polymorphismData<-list()
for (i in 1:length(fnames)){

	#final run data 
	outdata<-read_SNPS(paste("Mutations_36from1/",fnames[i],sep=""))
		
	unlist(strsplit(fnames[i],"_"))[2:3]->sim_pars
	paste(sim_pars,collapse="_")->sp
	
	#need to strip out those alleles that are fixed in our sample
	#these are present because they're segregating in the population
	sapply(outdata,function(x){x$count!=2*indivs_sampled})->not_fixed
	outdata<-outdata[not_fixed]
	
	#those which are synonymous
	sapply(outdata,function(x){x$type=="m1"})->syn

	#those which are non-synonymous
	sapply(outdata,function(x){x$type!="m1"})->non

	#those which are negatively selected
	sapply(outdata,function(x){x$type=="m2"})->del

	#those which are positively selected
	sapply(outdata,function(x){x$type=="m3"})->pos	

	#SFS count frequency with table, but 1:(2*indivs_sampled) is to force inclusion of frequencies that are not seen
	table(c(sapply(outdata,function(x){x$count})[syn],1:(2*indivs_sampled)))-1->SFS_Syn
	table(c(sapply(outdata,function(x){x$count})[non],1:(2*indivs_sampled)))-1->SFS_Non
	
	#print(SFS_Non)
	
	#age distributions 
	sapply(outdata,function(x){final_generation-x$info}[4])[syn]->SynAge
	sapply(outdata,function(x){final_generation-x$info}[4])[del]->DelAge
	sapply(outdata,function(x){final_generation-x$info}[4])[pos]->PosAge
	
	#DFE of mutations  
	sapply(outdata,function(x){x$info}[2])[del]->snpDFE_Non

	#add that list to a list of all polymorphism data seen across all datasets
	if(is.null(polymorphismData[[sp]])){
		polymorphismData[[sp]][[1]]<-list(SFS_Syn=SFS_Syn,SFS_Non=SFS_Non,SynAge=SynAge,DelAge=DelAge,PosAge=PosAge,snpDFE_Non=snpDFE_Non)
	}else{
		polymorphismData[[sp]][[length(polymorphismData[[sp]])+1]]<-list(SFS_Syn=SFS_Syn,SFS_Non=SFS_Non,SynAge=SynAge,DelAge=DelAge,PosAge=PosAge,snpDFE_Non=snpDFE_Non)
	}
	print(paste("Finished reading file ",i,sep=""))	
}

polymorphismData->polymorphismData_36from1


#polymorphism counts
polymorphisms<-vector()
r<-1
for(i in 1:length(polymorphismData)){
	dat<-polymorphismData[[i]]
	for(j in 1:length(dat)){
		polymorphisms<-rbind(polymorphisms,c(as.numeric(unlist(strsplit(names(polymorphismData[i]),"_"))),r,length(dat[[j]]$SynAge),length(dat[[j]]$DelAge),length(dat[[j]]$PosAge),length(dat[[j]]$DelAge)+length(dat[[j]]$PosAge)))
		r<-r+1
	}
}
colnames(polymorphisms)<-c("migration","sex_freq","Rep","Syn","Non_del","Non_pos","Non")
polymorphisms->polymorphisms_36from1



####################################################################################
#  
# plotting
#
####################################################################################

indivs_sampled<-36 #(36 diploids from one populations)

pdf(file = "Figures/SamplePolymorphism_36Indiv_from_1Deme.pdf",width=8, height=6) 
	alleles_sampled<-36*2
	polymorphismData<-polymorphismData_1from36
	source("AssymptoticMK.R")
	source("PlotPolymorphism_And_PlotAsymptoticMK.R")	

dev.off()

AsymptoticMK_results_36from1<-data.frame(Migration=ms,Sex=ss,obs_alpha=obs_alpha,alpha_original=overall_results3$alpha_original,alpha_assymptotic=overall_results3$alpha_asymptotic,CI_low=overall_results3$CI_low,CI_high=overall_results3$CI_high)
AsymptoticMK_results_36from1<-cbind(AsymptoticMK_results_36from1,obs_omega_a,Ks,Ks_apparant,Ka,Ka_apparant,pis_sample,pis_het,Fst_syn,pia_sample,pia_het,Fst_non)
signif(AsymptoticMK_results_36from1,4)


# plot(AsymptoticMK_results[,3],AsymptoticMK_results[,5],xlim=c(0,0.5),ylim=c(0,0.5),xlab="Observed alpha",ylab="Asymptotic Estimate of alpha",pch=19);
# arrows(AsymptoticMK_results[,3],AsymptoticMK_results[,5],AsymptoticMK_results[,3],AsymptoticMK_results[,6],length=0.01,angle=90,xlim=c(0,0.5),ylim=c(0,0.5),xlab="Observed alpha",ylab="Asymptotic Estimate of alpha",pch=19);
# arrows(AsymptoticMK_results[,3],AsymptoticMK_results[,5],AsymptoticMK_results[,3],AsymptoticMK_results[,7],length=0.01,angle=90,xlim=c(0,0.5),ylim=c(0,0.5),xlab="Observed alpha",ylab="Asymptotic Estimate of alpha",pch=19);


####################################################################################
#
# Part III: In which we count the number of segregating polymorphisms at the end of the simulation, in this case based on 10 individuals sampled from each of 36 demes
# Edit this code for other sampling strategies
#
###################################################################

list.files("Mutations_10from36", pattern="mutations_0.[0-9].+.txt")->fnames

indivs_sampled<-10*36 #(10 diploids from each of 36 populations)

final_generation<-1000000

polymorphismData<-list()
for (i in 1:length(fnames)){

	#final run data 
	outdata<-read_SNPS(paste("Mutations_10from36/",fnames[i],sep=""))
		
	unlist(strsplit(fnames[i],"_"))[2:3]->sim_pars
	paste(sim_pars,collapse="_")->sp
	
	#need to strip out those alleles that are fixed in our sample
	#these are present because they're segregating in the population
	sapply(outdata,function(x){x$count!=2*indivs_sampled})->not_fixed
	outdata<-outdata[not_fixed]
	
	#those which are synonymous
	sapply(outdata,function(x){x$type=="m1"})->syn

	#those which are non-synonymous
	sapply(outdata,function(x){x$type!="m1"})->non

	#those which are negatively selected
	sapply(outdata,function(x){x$type=="m2"})->del

	#those which are positively selected
	sapply(outdata,function(x){x$type=="m3"})->pos	

	#SFS count frequency with table, but 1:(2*indivs_sampled) is to force inclusion of frequencies that are not seen
	table(c(sapply(outdata,function(x){x$count})[syn],1:(2*indivs_sampled)))-1->SFS_Syn
	table(c(sapply(outdata,function(x){x$count})[non],1:(2*indivs_sampled)))-1->SFS_Non
	
	#print(SFS_Syn)
	
	#age distributions 
	sapply(outdata,function(x){final_generation-x$info}[4])[syn]->SynAge
	sapply(outdata,function(x){final_generation-x$info}[4])[del]->DelAge
	sapply(outdata,function(x){final_generation-x$info}[4])[pos]->PosAge
	
	#DFE of mutations  
	sapply(outdata,function(x){x$info}[2])[del]->snpDFE_Non

	#add that list to a list of all polymorphism data seen across all datasets
	if(is.null(polymorphismData[[sp]])){
		polymorphismData[[sp]][[1]]<-list(SFS_Syn=SFS_Syn,SFS_Non=SFS_Non,SynAge=SynAge,DelAge=DelAge,PosAge=PosAge,snpDFE_Non=snpDFE_Non)
	}else{
		polymorphismData[[sp]][[length(polymorphismData[[sp]])+1]]<-list(SFS_Syn=SFS_Syn,SFS_Non=SFS_Non,SynAge=SynAge,DelAge=DelAge,PosAge=PosAge,snpDFE_Non=snpDFE_Non)
	}
	print(paste("Finished reading file ",i,sep=""))	
}

polymorphismData->polymorphismData_10from36
polymorphismData_10from36->polymorphismData

#polymorphism counts
polymorphisms<-vector()
r<-1
for(i in 1:length(polymorphismData)){
	dat<-polymorphismData[[i]]
	for(j in 1:length(dat)){
		polymorphisms<-rbind(polymorphisms,c(as.numeric(unlist(strsplit(names(polymorphismData[i]),"_"))),r,length(dat[[j]]$SynAge),length(dat[[j]]$DelAge),length(dat[[j]]$PosAge),length(dat[[j]]$DelAge)+length(dat[[j]]$PosAge)))
		r<-r+1
	}
}
colnames(polymorphisms)<-c("migration","sex_freq","Rep","Syn","Non_del","Non_pos","Non")
polymorphisms->polymorphisms_10from36


####################################################################################
#  
# Plotting Polymorphism data and running assymptoticMK
# This uses code in AssymptoticMK.R and PlotPolymorphism_And_PlotAsymptoticMK.R  
#
####################################################################################

pdf(file = "Figures/SamplePolymorphism_10Indiv_from_36Demes.pdf",width=8, height=6) 
	alleles_sampled<-36*2
	polymorphismData<-polymorphismData_1from36
	source("AssymptoticMK.R")
	source("PlotPolymorphism_And_PlotAsymptoticMK.R")	
dev.off()

AsymptoticMK_results_10from36<-data.frame(Migration=ms,Sex=ss,obs_alpha=obs_alpha,alpha_original=overall_results3$alpha_original,alpha_assymptotic=overall_results3$alpha_asymptotic,CI_low=overall_results3$CI_low,CI_high=overall_results3$CI_high)
AsymptoticMK_results_10from36<-cbind(AsymptoticMK_results_10from36,obs_omega_a,Ks,Ks_apparant,Ka,Ka_apparant,pis_sample,pis_het,Fst_syn,pia_sample,pia_het,Fst_non)
signif(AsymptoticMK_results_10from36,2)


#Having edited the code above for each of the sampling strategies, this saves all the results all in on R file

save(list=c("convergenceData","fixationData","fixations",
	"polymorphismData_1from36","polymorphismData_36from1","polymorphismData_10from36",
	"polymorphisms_1from36","polymorphisms_36from1","polymorphisms_10from36",
	"AsymptoticMK_results_1from36","AsymptoticMK_results_36from1","AsymptoticMK_results_10from36"),
file="Fields2021_SlimulationOutput.Rdata")



#######################################################################
 #  _____  ______ ______                  _       _           
 # |  __ \|  ____|  ____|                | |     | |          
 # | |  | | |__  | |__     ______    __ _| |_ __ | |__   __ _ 
 # | |  | |  __| |  __|   |______|  / _` | | '_ \| '_ \ / _` |
 # | |__| | |    | |____           | (_| | | |_) | | | | (_| |
 # |_____/|_|    |______|           \__,_|_| .__/|_| |_|\__,_|
 #                                         | |                
 #                                         |_|                
#######################################################################
#
#	This assumes that Multi-dfe is available locally in a folder called "MultiDFE"
#		and that dfe alpha is available locally in a folder called "dfe-alpha-release-2.16"
# 	These locations are coded into the function below
# 	It also requires the polymorphismData and fixationData lists from above  
#
#######################################################################


# The files involved are listed below
#
# Default file name 				Comment 																																				Corresponding parameter in config file
# directory_config.dat 				Contains a single line of text specifying the directory containing the data files for 1 and 2 epoch models.												data_path_1
# directory_config_three_epoch.dat	Contains a single line of text specifying the directory containing the data files for 3 epoch model. (Must be different from data_path_1).				data_path_2
# sfs.txt 							Input file containing the site frequency spectra for analysis by est_dfe.																				sfs_input_file
# est_dfe.out 						Main results file produced by est_dfe. 
# neut_egf.out 						Neutral gene frequency vector file produced by est_dfe.																									-
# sel_egf.out 						Selected gene frequency vector file produced by est_dfe.																								-
# est_alpha_omega.out 				Results file produced by est_alpha_omega. 																												est_alpha_omega_results_file 
# divergence.txt I					Input file for est_alpha_omega containing counts of nucleotide sites and differences. 																	divergence_file
# Example of input file with a single pair of SFSs with 11 alleles:
# 1
# 11
# 317781 1410 167 30 11 12 14 17 29 25 33 41
# 132554 3031 1411 1001 731 618 504 503 419 402 421 563


#Convenience function to run dfe-alpha and MultiDFE
RunAlpha<-function(polymorphismData=polymorphismData,fixationData=fixationData,alleles_sampled=alleles_sampled){
	combined_DFEs<-vector()
	for(i in 1:length(polymorphismData)){
		MnS<-as.numeric(unlist(strsplit(names(polymorphismData)[i],"_")))
		if(file.exists("SFS.txt")){file.remove("SFS.txt")}
		if(file.exists("SFS.mdfe.txt")){file.remove("SFS.mdfe.txt")}
		if(file.exists("divergence.txt")){file.remove("divergence.txt")}
		Pdat<-polymorphismData[[i]]
		Fdat<-fixationData[[i]]
		synFix<-0
		nonFix<-0
		posFix<-0
		sP<-vector()
		nP<-vector()
		for(j in 1:length(Pdat)){
			if(!is.na(Pdat[[j]]$SFS_Syn[as.character(alleles_sampled)])){ #if there is an allele that appears fixed in the sample
				synFix<-synFix+Pdat[[j]]$SFS_Syn[as.character(alleles_sampled)] #not correcting for apparent fixations in DFEalpha
				Pdat[[j]]$SFS_Syn<-Pdat[[j]]$SFS_Syn[-length(Pdat[[j]]$SFS_Syn)];
			}
			if(!is.na(Pdat[[j]]$SFS_Non[as.character(alleles_sampled)])){ #if there is an allele that appears fixed in the sample
				nonFix<-nonFix+Pdat[[j]]$SFS_Non[as.character(alleles_sampled)] #not correcting for apparent fixations in DFEalpha
				Pdat[[j]]$SFS_Non<-Pdat[[j]]$SFS_Non[-length(Pdat[[j]]$SFS_Non)];
			}
			sP<-c(sP,rep(as.numeric(names(Pdat[[j]]$SFS_Syn[-alleles_sampled])),Pdat[[j]]$SFS_Syn[-alleles_sampled]))
			nP<-c(nP,rep(as.numeric(names(Pdat[[j]]$SFS_Non[-alleles_sampled])),Pdat[[j]]$SFS_Non[-alleles_sampled]))
			synFix<-synFix+length(Fdat[[j]]$SynFixTime)
			nonFix<-nonFix+length(Fdat[[j]]$DelFixTime)+length(Fdat[[j]]$PosFixTime)
			#print(c(j,nonFix))
			posFix<-posFix+length(Fdat[[j]]$PosFixTime)
		}
		sP_t<-table(c(sP,1:(alleles_sampled-1)))-1
		nP_t<-table(c(nP,1:(alleles_sampled-1)))-1

		#fold the SFS's
		sP_t[1:(alleles_sampled/2-1)]<-sP_t[1:(alleles_sampled/2-1)]+sP_t[(alleles_sampled-1):(1+alleles_sampled/2)]
		sP_t[(alleles_sampled-1):(1+alleles_sampled/2)]<-0
		
		nP_t[1:(alleles_sampled/2-1)]<-nP_t[1:(alleles_sampled/2-1)]+nP_t[(alleles_sampled-1):(1+alleles_sampled/2)]
		nP_t[(alleles_sampled-1):(1+alleles_sampled/2)]<-0
		
		#calc the number of sites than are not polymorphic
		total_syn_Sites<-round(length(Pdat)*5000 * 50 * 0.24,0) #combine all simulations into one analysis - so multiply length by number of reps. Assume 24% of mutations are synonymous
		total_non_Sites<-round(length(Pdat)*5000 * 50 * 0.76,0)
		
		#create the sfs input file for est_dfe
		cat("1\n",file="SFS.txt",append=FALSE)
		cat(paste(alleles_sampled,"\n",sep=""),file="SFS.txt",append=TRUE)
		cat(c(total_non_Sites-sum(nP_t),nP_t,0),file="SFS.txt",append=TRUE)
		cat("\n",file="SFS.txt",append=TRUE)
		cat(c(total_syn_Sites-sum(sP_t),sP_t,0),file="SFS.txt",append=TRUE)
		cat("\n",file="SFS.txt",append=TRUE)
		#create the sfs input file for Multi-DFE
		cat(paste(alleles_sampled,"\n",sep=""),file="SFS.mdfe.txt",append=TRUE)
		cat(c(total_non_Sites-sum(nP_t),nP_t,0),file="SFS.mdfe.txt",append=TRUE)
		cat("\n",file="SFS.mdfe.txt",append=TRUE)
		cat(c(total_syn_Sites-sum(sP_t),sP_t,0),file="SFS.mdfe.txt",append=TRUE)
		cat("\n",file="SFS.mdfe.txt",append=TRUE)
		
		#create the divergence input file for est_alpha_omega
		cat(c(1,total_non_Sites,nonFix),file="divergence.txt",append=FALSE)
		cat("\n",file="divergence.txt",append=TRUE)
		cat(c(0,total_syn_Sites,synFix),file="divergence.txt",append=TRUE)
		cat("\n",file="divergence.txt",append=TRUE)

		#run multi-DFE components
		all_pars<-c(2,2,2,5,1,3,5,7,9,1,3,5,7,9) #number of parameters in each model
		mdfe_estimates<-list()
		cat(file="SFS.mdfe.txt.MAXL.out",paste("\nParameterSet ",names(polymorphismData_10from36)[i],"\n",sep=" "),append=TRUE)
	#lognormal	2 parameters
		cat(file="SFS.mdfe.txt.MAXL.out","\nLognormal\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 4 -nspikes 0 -ranrep 5 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["log_norm"]]
		names(mdfe_estimates[["log_norm"]])<-gsub(":.+","",L)
	#gamma		2 parameters
		cat(file="SFS.mdfe.txt.MAXL.out","\nGamma\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 2 -nspikes 0 -ranrep 5 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["gamma"]]
		names(mdfe_estimates[["gamma"]])<-gsub(":.+","",L)
	#beta		2 parameters
		cat(file="SFS.mdfe.txt.MAXL.out","\nBeta\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 3 -nspikes 0 -ranrep 5 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["beta"]]
		names(mdfe_estimates[["beta"]])<-gsub(":.+","",L)
	#6 fixed spikes	5 parameters
		cat(file="SFS.mdfe.txt.MAXL.out","\nFixed 6 Spike\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 5 -nspikes 0 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["F6sp"]]
		names(mdfe_estimates[["F6sp"]])<-gsub(":.+","",L)
			
	#estimated spike distributions	{1,3,5,7,9} parameters
		cat(file="SFS.mdfe.txt.MAXL.out","\nSpike 1\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 0 -nspikes 1 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["sp1"]]
		names(mdfe_estimates[["sp1"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nSpike 2\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 0 -nspikes 2 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["sp2"]]
		names(mdfe_estimates[["sp2"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nSpike 3\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 0 -nspikes 3 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["sp3"]]
		names(mdfe_estimates[["sp3"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nSpike 4\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 0 -nspikes 4 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["sp4"]]
		names(mdfe_estimates[["sp4"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nSpike 5\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 0 -nspikes 5 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["sp5"]]
		names(mdfe_estimates[["sp5"]])<-gsub(":.+","",L)
		
	#estimated step distributions	{1,3,5,7,9} parameters
		cat(file="SFS.mdfe.txt.MAXL.out","\nStep 1\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 1 -nspikes 1 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["st1"]]
		names(mdfe_estimates[["st1"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nStep 2\n",append=TRUE)
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 1 -nspikes 2 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["st2"]]
		names(mdfe_estimates[["st2"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nStep 3\n",append=TRUE)	
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 1 -nspikes 3 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["st3"]]
		names(mdfe_estimates[["st3"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nStep 4\n",append=TRUE)	
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 1 -nspikes 4 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["st4"]]
		names(mdfe_estimates[["st4"]])<-gsub(":.+","",L)
		cat(file="SFS.mdfe.txt.MAXL.out","\nStep 5\n",append=TRUE)	
		system("MultiDFE/MultiDFE -conpop 0 -sfsfold 1 -selmode 1 -nspikes 5 -ranrep 10 -file SFS.mdfe.txt")
		unlist(strsplit(system("tail -n 1 SFS.mdfe.txt.MAXL.out",intern=TRUE),"\t"))->L
		as.numeric(gsub("^.+:| ","",L))->mdfe_estimates[["st5"]]
		names(mdfe_estimates[["st5"]])<-gsub(":.+","",L)

		#remove any entries that failed to run for numerical reasons
		pars<-all_pars[which(lapply(mdfe_estimates,length)>2)]
		mdfe_estimates<-mdfe_estimates[lapply(mdfe_estimates,length)>2]
		
		MultiDFEs<-t(sapply(mdfe_estimates,function(x){x[c('L','fix_prob','Nes_0.0_0.1','Nes_0.1_1.0','Nes_1.0_10.0','Nes_10.0_100.0')]}))
		
		#calc dN & dS for alpha and omega_a (equations 10 & 11 in "A Comparison of Models to Infer the Distribution of Fitness Effects of New Mutations" Kousathanas & Keightley (2013) https://doi.org/10.1534/genetics.112.148023
		# eqtn 10 (dN-dS*ubar)/dN
		# eqtn 11(dN-dS*ubar)/dS
		# where u_bar is the probability of fixation (fix_prob) as provided by multiDFE	
		total_syn_Sites<-round(length(Pdat)*5000 * 50 * 0.24,0) #combine all simulations into one analysis - so multiply length by number of reps
		total_non_Sites<-round(length(Pdat)*5000 * 50 * 0.76,0)	
		dS<-synFix/total_syn_Sites
		dN<-nonFix/total_non_Sites

		overall_DFE<-data.frame(
			migration=MnS[1],sex=MnS[2],
			MultiDFEs,
			Nes_Over_100.0=1-MultiDFEs[,"Nes_0.0_0.1"]-MultiDFEs[,"Nes_0.1_1.0"]-MultiDFEs[,"Nes_1.0_10.0"]-MultiDFEs[,"Nes_10.0_100.0"],
			alpha=(dN-(dS*(MultiDFEs[,"fix_prob"])))/dN,
			omega_a=(dN-(dS*(MultiDFEs[,"fix_prob"])))/dS,
			#calc AIC and Akaike weights for the models
			AIC=(2*pars-2*MultiDFEs[,"L"]),
			delta_AIC=(2*pars-2*MultiDFEs[,"L"])-min((2*pars-2*MultiDFEs[,"L"]),na.rm=TRUE),
			AkaikeWeight=exp(-0.5*((2*pars-2*MultiDFEs[,"L"])-min((2*pars-2*MultiDFEs[,"L"]),na.rm=TRUE)))/sum(exp(-0.5*((2*pars-2*MultiDFEs[,"L"])-min((2*pars-2*MultiDFEs[,"L"]),na.rm=TRUE))),na.rm=TRUE)
		)
		
		#run dfe-apha components
		system("dfe-alpha-release-2.16/est_dfe -c Neutral_DFE.config")
		system("dfe-alpha-release-2.16/est_dfe -c Selected_DFE.config")
		system("dfe-alpha-release-2.16/est_alpha_omega -c Est_AlphaOmega.config") #No JC correction and don't remove polys
		system("dfe-alpha-release-2.16/prop_muts_in_s_ranges -c TempDFE/est_dfe.out -o TempDFE/estimates_DFE.txt")
		#read int the results
		#read.table("TempDFE/neut_exp_obs_allele_freqs.csv",head=TRUE)->syn_SFS
		#read.table("TempDFE/sel_exp_obs_allele_freqs.csv",head=TRUE)->nonsyn_SFS
		alpha_estimates<-read.table("TempDFE/est_alpha_omega.out")
		alpha_DFE<-matrix(scan("TempDFE/estimates_DFE.txt"),nrow=3)[3,1:4]
		file.remove(file.path("TempDFE", list.files("TempDFE")))
			
		overall_DFE<-round(rbind(as.vector(unlist(c(MnS,NA,alpha_estimates[2],NA,alpha_DFE,alpha_estimates[6],alpha_estimates[8],NA,NA,NA))),overall_DFE),4)
		rownames(overall_DFE)<-c("DFE_alpha",rownames(overall_DFE)[-1])
		
		combined_DFEs<-rbind(combined_DFEs,overall_DFE)
	}
	return(combined_DFEs)
} #end of RunAlpha function





#Convenience function to plot dfe-alpha and MultiDFE results

PlotAlpha<-function(polymorphismData=polymorphismData,AsymptoticMK_results=AsymptoticMK_results,combined_DFEs=combined_DFEs){
	

	#Plot DFE of Fixations and SNPs
	plotdim<-ceiling(sqrt(length(polymorphismData)))
	par(
		mfrow=c(plotdim,plotdim),
		mar=c(2,1,0,1),
		oma=c(6,8,4,1)
	)
	
	for(i in 1:length(polymorphismData)){

		EstNe<-AsymptoticMK_results$pis_het[i]/(4*1.53e-7) #Ne~pi_s_local/(4*mu)
		#Simulated DFE (from a large random draw, rather than doing some maths)
		hist(-rgamma(100000,0.056/0.3,0.3)*2*EstNe,breaks=c(-1e500,-100,-10,-1,0),plot=FALSE)$count->muts ; muts<-muts/sum(muts)
		Simulated<-c(rev(muts),AsymptoticMK_results$obs_alpha[i])
		#DFE and alpha from DFE_alpha
		DFEalpha<-combined_DFEs[grep("DFE_alpha",rownames(combined_DFEs))[i],c(6:9,10)]	
		
		#DFE and alpha from MultiDFE	
		simPars<-as.numeric(combined_DFEs[grep("DFE_alpha",rownames(combined_DFEs))[i],1:2])
		temp<-combined_DFEs[combined_DFEs$migration==simPars[1]&combined_DFEs$sex==simPars[2],]
		MultiDFE_best<-temp[which.max(temp$AkaikeWeight),c(5:9,10)]
		model<-substr(rownames(MultiDFE_best),1,3)
		MultiDFE_best[1]<-MultiDFE_best[1]+MultiDFE_best[2]
		MultiDFE_best<-MultiDFE_best[-1]
		
		vals<-c(as.vector(as.matrix(rbind(Simulated[1:4],DFEalpha[1:4],MultiDFE_best[1:4]))),Simulated[5],as.numeric(DFEalpha[5]),as.numeric(MultiDFE_best[5]),AsymptoticMK_results$alpha_assymptotic[i])
		mycols<-c("gray10","slateblue4","slateblue","red4")
		plot(0,xlab="",ylab="",main="",ylim=c(-0.35,1),xlim=c(0,23),axes=FALSE,type="n") #or ylim=c(min(0,min(vals)),1)
		segments(0,0,23,0)
		segments(17,-2,17,2,col="gray70")
		#segments(18,Simulated[5],22,Simulated[5],col="gray10")
		
		xleft<-c(1,2,3,5,6,7,9,10,11,13,14,15,18,19,20,21)
		rect(xleft,0,xleft+1,vals,col=c(rep(mycols[1:3],4),mycols))
		legend("bottomleft",legend=c("DFE"),cex=0.8,bty="n")
		legend("bottomright",legend=c("alpha"),cex=0.8,bty="n")
		legend("topleft",legend=c("Simulated","dfe-alpha",paste("multi-dfe (",model,")",sep=""),"asymptotic DFE"),fill=mycols,cex=0.6,bty="n")
		box(col="gray70")
		#axis(1,at=c(1,1e-2,1e-4,1e-6,1e-8,1e-10),labels=c(1,1e-2,1e-4,1e-6,1e-8,1e-10))

		#pretty Plotting
		as.numeric(gsub("_.+","",names(fixationData)[i]))->m
		as.numeric(gsub(".+_","",names(fixationData)[i]))->s
		if(i>((plotdim-1)*plotdim)){
			#mtext("Deleterious fitness effect",side=1,line=5,cex=1.02)
			mtext(c("0-1","1-10","10-100","100+",""),side=1,line=1,cex=0.8,at=c(2.5,6.5,10.5,14.5,20),las=2)
		}
		if(i%%plotdim==1){
			mtext(paste("F_st = ",round(1/(1+4*500*m),2),sep=""),side=2,line=5,cex=1.2)
			mtext("Proportion",side=2,line=2,cex=0.5)
			axis(2,at=c(0,1.0),las=2,line=1)	
		}
		if(i<=plotdim){mtext(paste("Sex 1 in ",s,sep=""),side=3,cex=1.2,line=2)}	
		
		#if(i==1){legend("topright",legend=c("Simulated","dfe-Alpha","Multi-DFE","Asymptotic DFE"),fill=mycols,cex=1.0,bty="n")}
	}

} #End of PlotAlpha()


#all replicates together 

load("Fields2021_SlimulationOutput.Rdata") 

########################################################################
# Sample 1 diploid individuals from each of the 36 simulated demes

alleles_sampled<-2*36 
polymorphismData<-polymorphismData_1from36 
cat(file="SFS.mdfe.txt.MAXL.out","\n1 from 36 demes\n",append=FALSE)
RunAlpha(polymorphismData=polymorphismData_1from36,fixationData=fixationData,alleles_sampled=alleles_sampled)->combined_DFEs_1from36
write.table(combined_DFEs_1from36,file="Inferred_DFE_1from36.tsv",sep="\t",row.names=TRUE,col.names=TRUE)
pdf(file = "Figures/DFE_Estimates_1From36.pdf",width=8, height=6) 
	PlotAlpha(polymorphismData=polymorphismData_1from36,AsymptoticMK_results=AsymptoticMK_results_1from36,combined_DFEs=combined_DFEs_1from36)
dev.off()



########################################################################
# Sample 36 diploid individuals from a single deme (of 36 simulated)

alleles_sampled<-2*36 
polymorphismData<-polymorphismData_36from1 
cat(file="SFS.mdfe.txt.MAXL.out","\n36 from 1 demes\n",append=FALSE)
RunAlpha(polymorphismData=polymorphismData_36from1,fixationData=fixationData,alleles_sampled=alleles_sampled)->combined_DFEs_36from1
write.table(combined_DFEs_36from1,file="Inferred_DFE_36from1.tsv",sep="\t",row.names=TRUE,col.names=TRUE)
pdf(file = "Figures/DFE_Estimates_1From36.pdf",width=8, height=6) 
	PlotAlpha(polymorphismData=polymorphismData_36from1,AsymptoticMK_results=AsymptoticMK_results_36from1,combined_DFEs=combined_DFEs_36from1)
dev.off()




########################################################################
# Sample 10 diploid individuals from each of the 36 simulated demes

alleles_sampled<-2*10*36 
polymorphismData<-polymorphismData_10from36 
cat(file="SFS.mdfe.txt.MAXL.out","\n10 from 36 demes\n",append=FALSE)
RunAlpha(polymorphismData=polymorphismData_10from36,fixationData=fixationData,alleles_sampled=alleles_sampled)->combined_DFEs_10from36
write.table(combined_DFEs_10from36,file="Inferred_DFE_10from36.tsv",sep="\t",row.names=TRUE,col.names=TRUE)
pdf(file = "Figures/DFE_Estimates_1From36.pdf",width=8, height=6) 
	PlotAlpha(polymorphismData=polymorphismData_10from36,AsymptoticMK_results=AsymptoticMK_results_10from36,combined_DFEs=combined_DFEs_10from36)
dev.off()


save(list=c("combined_DFEs_1from36","combined_DFEs_36from1","combined_DFEs_10from36"),file="Fields2021_Slimulation_DFEinference.Rdata")


















































