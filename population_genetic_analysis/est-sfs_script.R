


source("est-sfs_wrapper.R")

#########################################################################
# Example script for two outgroups  ### DOES NOT WORK!!!
#########################################################################


setwd("/localdisk/home/dobbard/Daphnia/est_sfs/DJO_sfs/Daphnia_Trials/DoubleOutgroup")

max_depth<-72 #maximum depth you want to consider, e.g. all 72 chromosomes sampled
min_depth<-50 #minimum depth you want to consider, e.g. only as low as 50 chromosomes sampled


#find the files
file_names<-list.files(pattern=".fas?(ta)?$")

#for each file, make a sereies of est-sfs input files, according to depth and position
for (f in file_names){
#reads in a fasta as a data matrix. sequinr is bad at this, and I could make it fasta
	#read the data in as a matrix
	data_matrix<-toupper(as.matrix(read.alignment(f,format="fasta",forceToLower = FALSE)))
	#format teh data for est-sfs
	prepEstSFS(data_matrix,out_list=c("similis","lumholtzi"))->locus
	#split into 1+2, and 3rd positions
	positionSplitEstSFS(locus,c(TRUE,TRUE,FALSE))->FirstSecondPos
	positionSplitEstSFS(locus,c(FALSE,FALSE,TRUE))->ThirdPos
	
	#strip out multiallelic sites
	FirstSecondPos<-stripMultiAlleles(FirstSecondPos)
	ThirdPos<-stripMultiAlleles(ThirdPos)
	
	#split by depth of observation #May not actually be necessary it seems - although will need to do it later for dfe-alpha
	depthSplitEstSFS(FirstSecondPos)->DsplitFirstSecondPos
	depthSplitEstSFS(ThirdPos)->DsplitThirdPos
	#write files out
	
	writeEstSFSInput(DsplitFirstSecondPos,filestem=paste(gsub("[.].+","",f),"_FirstSecondPos_",sep=""),out_list=c("similis","lumholtzi"))
	writeEstSFSInput(DsplitThirdPos,filestem=paste(gsub("[.].+","",f),"_ThirdPos_",sep=""),out_list=c("similis","lumholtzi"))
}


#combine the files to make suitable est-sfs input. Also make a record of how many lines from which file, and in which order they were concatenated
for(i in max_depth:min_depth){
	pack_estSFSfiles(paste("FirstSecondPos_",i,sep=""))
	pack_estSFSfiles(paste("ThirdPos_",i,sep=""))
}


# Runs est-sfs with a two outgroups outgroup

#run estSFS on full-depth files and record the GTR parameters
runEstSFS("FirstSecondPos_72_combined.tsv",GTRParms=NA,outgroups=2,relative_path="../../")->EsimatedPos12Parms
runEstSFS("ThirdPos_72_combined.tsv",GTRParms=NA,outgroups=2,relative_path="../../")->EsimatedPos3Parms

#run it on those with less than full depth.
for(i in (max_depth-1):min_depth){
	fn<-paste("FirstSecondPos_",i,"_combined.tsv",sep="")
	if(file.exists(fn)){runEstSFS(fn,GTRParms=EsimatedPos12Parms,outgroups=2,relative_path="../../")}
	fn<-paste("ThirdPos_",i,"_combined.tsv",sep="")
	if(file.exists(fn)){runEstSFS(fn,GTRParms=EsimatedPos3Parms,outgroups=2,relative_path="../../")}
}

#generate per-locus files for the ancestor inferences
unpack_estSFSfiles("(FirstSecondPos|ThirdPos)")

#Classify SNPs as synonymous or non-synonymous
file.create("NonSynonymous_Calls.tsv")
file.create("Synonymous_Calls.tsv")
for (f in file_names){
	classifySNPs(f,c("similis","lumholtzi"),genetic_code="Standard")->SNP_classification ## note that it does nothing about fixed differences! (ignores internal states int eh tree, and only uses the prob of the major allele being ancestral
	#write out all the synonymous SNPs
	write.table(SNP_classification[SNP_classification[,'Effect']=="Synonymous",], file="Synonymous_Calls.tsv", sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
	#write out all the non-synonymous SNPs
	write.table(SNP_classification[SNP_classification[,'Effect']=="Nonsynonymous",], file="NonSynonymous_Calls.tsv", sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
}

# Retreive the unfolded SFS, using teh same approach as PK, but for Syn /Nonsyn SNPs rather than 0-fold 2-fold and 4-fold sites
read.table("Synonymous_Calls.tsv",sep="\t",as.is=TRUE)->Syn
read.table("NonSynonymous_Calls.tsv",sep="\t",as.is=TRUE)->Non


## If you want to get the uSFS for a specific group of genes, subset Non and Syn based on column 1 (gene names) before this loop
u_sfs<-matrix(0,ncol=2,nrow=(max_depth-1));
colnames(u_sfs)<-c("Synonymous","Nonsynonymous");
rownames(u_sfs)<-c(1:(max_depth-1))
u_sfs_List<-list()


for(d in sort(unique(c(Non$V5,Syn$V5)),decreasing=TRUE)){
	u_sfs_List[[as.character(d)]]<-u_sfs
	for(i in 1:(d-1)){
		u_sfs_List[[as.character(d)]][i,1]<-sum(Syn$V3[(Syn$V6==i)&(Syn$V5==d)],1-Syn$V3[(Syn$V6==(d-i))&(Syn$V5==d)]) 
		u_sfs_List[[as.character(d)]][i,2]<-sum(Non$V3[(Non$V6==i)&(Non$V5==d)],1-Non$V3[(Non$V6==(d-i))&(Non$V5==d)]) 
	}
	u_sfs_List[[as.character(d)]]<-round(u_sfs_List[[as.character(d)]],0) #I suspect that most packages will want integers
}

#final result!

u_sfs_List #is a list where each elemetn is named for the sample coverage (element '71' has data for sites where 71 chromosomes were sequenced) and where the list element is a matrix of: rows=uSFS, cols= Syn|NonSyn 

###########################################################################
#special case of ignoring the depth and the % derived regardless of sampling depth, which I think is what AssyptoticMK wants

## If you want to get the uSFS for a specific group of genes, subset Non and Syn based on column 1 (gene names) before this loop
u_sfs<-matrix(0,ncol=2,nrow=(100-1));
colnames(u_sfs)<-c("Synonymous","Nonsynonymous");
rownames(u_sfs)<-c(1:(100-1))
for(i in 1:99){
	u_sfs[i,1]<-sum(Syn$V3[(Syn$V7==i)&(Syn$V5>=min_depth)],1-Syn$V3[(Syn$V7==(100-i))&(Syn$V5>=min_depth)]) 
	u_sfs[i,2]<-sum(Non$V3[(Non$V7==i)&(Non$V5>=min_depth)],1-Non$V3[(Non$V7==(100-i))&(Non$V5>=min_depth)]) 
}
u_sfs<-round(u_sfs,0) #I suspect that most packages will want integers


#########################################################################
# End of two-outgroup example script
#########################################################################







#########################################################################
# Example script For a single Outgroup
#########################################################################


setwd("/localdisk/home/dobbard/Daphnia/est_sfs/DJO_sfs/Daphnia_Trials/Single_outgroup")

#find the files
file_names<-list.files(pattern=".fas?(ta)?$")

#for each file, make a sereies of est-sfs input files, according to depth and position
for (f in file_names){
	#reads in a fasta as a data matrix. sequinr is bad at this, and I could make it fasta
	#read the data in as a matrix
	data_matrix<-toupper(as.matrix(read.alignment(f,format="fasta",forceToLower = FALSE)))
	#format teh data for est-sfs
	prepEstSFS(data_matrix,"similis")->locus
	#split into 1+2, and 3rd positions
	positionSplitEstSFS(locus,c(TRUE,TRUE,FALSE))->FirstSecondPos
	positionSplitEstSFS(locus,c(FALSE,FALSE,TRUE))->ThirdPos
	
	#strip out multiallelic sites
	FirstSecondPos<-stripMultiAlleles(FirstSecondPos)
	ThirdPos<-stripMultiAlleles(ThirdPos)
	
	#split by depth of observation #May not actually be necessary it seems - although will need to do it later for dfe-alpha
	depthSplitEstSFS(FirstSecondPos)->DsplitFirstSecondPos
	depthSplitEstSFS(ThirdPos)->DsplitThirdPos
	#write files out
	
	writeEstSFSInput(DsplitFirstSecondPos,filestem=paste(gsub("[.].+","",f),"_FirstSecondPos_",sep=""),"similis")
	writeEstSFSInput(DsplitThirdPos,filestem=paste(gsub("[.].+","",f),"_ThirdPos_",sep=""),"similis")
}

max_depth<-72
min_depth<-50
#combine the files to make suitable est-sfs input. Also make a record of how many lines from which file, and in which order they were concatenated
for(i in max_depth:min_depth){
	pack_estSFSfiles(paste("FirstSecondPos_",i,sep=""))
	pack_estSFSfiles(paste("ThirdPos_",i,sep=""))
}


# Runs est-sfs with a single outgroup

#run estSFS on full-depth files and record the GTR parameters #specify the relative path to teh executable and seed file
runEstSFS("FirstSecondPos_72_combined.tsv",GTRParms=NA,outgroups=1,relative_path="../../")->EsimatedPos12Parms
runEstSFS("ThirdPos_72_combined.tsv",GTRParms=NA,outgroups=1,relative_path="../../")->EsimatedPos3Parms

#run it on those with less than full depth.
for(i in (max_depth-1):min_depth){
	fn<-paste("FirstSecondPos_",i,"_combined.tsv",sep="")
	if(file.exists(fn)){runEstSFS(fn,GTRParms=EsimatedPos12Parms,outgroups=1,relative_path="../../")}
	fn<-paste("ThirdPos_",i,"_combined.tsv",sep="")
	if(file.exists(fn)){runEstSFS(fn,GTRParms=EsimatedPos3Parms,outgroups=1,relative_path="../../")}
}

#generate per-locus files for the ancestor inferences
unpack_estSFSfiles("(FirstSecondPos|ThirdPos)")

#Classify SNPs as synonymous or non-synonymous
file.create("NonSynonymous_Calls.tsv")
file.create("Synonymous_Calls.tsv")
for (f in file_names){
	classifySNPs(f,"similis",genetic_code="Standard")->SNP_classification
	#write out all the synonymous SNPs
	write.table(SNP_classification[SNP_classification[,'Effect']=="Synonymous",], file="Synonymous_Calls.tsv", sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
	#write out all the non-synonymous SNPs
	write.table(SNP_classification[SNP_classification[,'Effect']=="Nonsynonymous",], file="NonSynonymous_Calls.tsv", sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
}

# Retreive the unfolded SFS, using teh same approach as PK, but for Syn /Nonsyn SNPs rather than 0-fold 2-fold and 4-fold sites
read.table("Synonymous_Calls.tsv",sep="\t",as.is=TRUE)->Syn
read.table("NonSynonymous_Calls.tsv",sep="\t",as.is=TRUE)->Non

##### If we wanted a list of uSFS's based on all depths
## If you want to get the uSFS for a specific group of genes, subset Non and Syn based on column 1 (gene names) before this loop
u_sfs<-matrix(0,ncol=2,nrow=(max_depth-1));
colnames(u_sfs)<-c("Synonymous","Nonsynonymous");
rownames(u_sfs)<-c(1:(max_depth-1))
u_sfs_List<-list()
for(d in sort(unique(c(Non$V5,Syn$V5)),decreasing=TRUE)){
	u_sfs_List[[as.character(d)]]<-u_sfs
	for(i in 1:(d-1)){
		u_sfs_List[[as.character(d)]][i,1]<-sum(Syn$V3[(Syn$V6==i)&(Syn$V5==d)],1-Syn$V3[(Syn$V6==(d-i))&(Syn$V5==d)]) 
		u_sfs_List[[as.character(d)]][i,2]<-sum(Non$V3[(Non$V6==i)&(Non$V5==d)],1-Non$V3[(Non$V6==(d-i))&(Non$V5==d)]) 
	}
	u_sfs_List[[as.character(d)]]<-round(u_sfs_List[[as.character(d)]],0) #I suspect that most packages will want integers
}




#########################################################################
# End of single -outgroup example script
#########################################################################





