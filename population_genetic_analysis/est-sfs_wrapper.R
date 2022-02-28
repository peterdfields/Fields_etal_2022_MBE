require(seqinr)

#####################################################
## define functions
#####################################################

# function to convert a data matrix and returns it in a format suitable for generating est-sfs output 
prepEstSFS<-function(data_matrix,out_list){
	bases<-c("A","C","G","T")
	data_matrix<-toupper(data_matrix)
	#outgroups, in the order specified in out_list
	if(length(out_list)==1){
		out_matrix<-as.matrix(data_matrix[out_list,])
	}else{
		out_matrix<-t(as.matrix(data_matrix[out_list,]))
	}
	
	#ingroup
	in_matrix<-data_matrix[!rownames(data_matrix) %in% out_list,]
	#counts per site
	counts<-t(apply(in_matrix,2,function(x){table(c(x[x%in%bases],c("A","C","G","T")))-1}))
	#sampling depth
	depths<-rowSums(counts)
	#record the absolute locations - needed later to reconsitute files!
	locations<-c(1:ncol(data_matrix))
	#create list to return
	For_estSFS<-list(locations=locations,depths=depths,ingroup=counts)
	#for each outgroup present, ad an element to the list
	for(o in 1:length(out_list)){
		#create an empty matrix 
		matrix(0,nrow=ncol(data_matrix),ncol=4)->temp_matrix
		#populate it
		temp_matrix[,1]<-as.numeric(out_matrix[,o]=="A")
		temp_matrix[,2]<-as.numeric(out_matrix[,o]=="C")
		temp_matrix[,3]<-as.numeric(out_matrix[,o]=="G")
		temp_matrix[,4]<-as.numeric(out_matrix[,o]=="T")
		For_estSFS[[out_list[o]]]<-temp_matrix
	}
	
	#return the list
	return(For_estSFS)	
}

# function to split est-sfs input by position (given as a mask called 'positions', which repeats) 
positionSplitEstSFS<-function(For_estSFS,positions=c(TRUE,TRUE,FALSE)){
	#Error if the arguments are incompatible
	if(length(For_estSFS$locations)%%length(positions)!=0){print("Error, positions mask not a multiple of sequence length");return(NA)}
	#for each element in the list
	for(i in 1:length(For_estSFS)){
		if(is.vector(For_estSFS[[i]])){
			#subset vectors Why, R, why?!!?
			For_estSFS[[i]]<-For_estSFS[[i]][positions]
		}else{
			#subset matrices.
			For_estSFS[[i]]<-For_estSFS[[i]][positions,]
		}
	}
	return(For_estSFS)
}


# function to split est-sfs input by ingroup coverage depth (number of chr's covered). Returns a list of "For_estSFS"-like lists, each one named for the depth
stripMultiAlleles<-function(For_estSFS){
	
	#these sites have <=2 alleles
	apply(For_estSFS$ingroup,1,function(x){sum(x>0)<=2})->positions
	#for each element in the list
	for(i in 1:length(For_estSFS)){
		if(is.vector(For_estSFS[[i]])){
			#subset vectors Why, R, why?!!?
			For_estSFS[[i]]<-For_estSFS[[i]][positions]
		}else{
			#subset matrices.
			For_estSFS[[i]]<-For_estSFS[[i]][positions,]
		}
	}
	return(For_estSFS)
}


# function to split est-sfs input by ingroup coverage depth (number of chr's covered). Returns a list of "For_estSFS"-like lists, each one named for the depth
depthSplitEstSFS<-function(For_estSFS,min_depth=1){
	#for each read depth in the i
	output<-list()
	depths_seen<-sort(unique(For_estSFS$depths),decreasing=TRUE)
	for(d in depths_seen[depths_seen>=min_depth]){	
		output[[as.character(d)]]<-For_estSFS
		positions<-(For_estSFS$depths==d)
		#for each element in the list
		for(i in 1:length(For_estSFS)){
			if(is.vector(For_estSFS[[i]])){
				#subset vectors Why, R, why?!!?
				output[[as.character(d)]][[i]]<-For_estSFS[[i]][positions]
			}else{
				#subset matrices.
				output[[as.character(d)]][[i]]<-For_estSFS[[i]][positions,]
			}
		}
	}
	return(output)
}



#function to write est-sfs input
writeEstSFSInput<-function(AllDepths,filestem="sfs_data_",out_list){
	
	#for each depth recorded in the file
	for(d in names(AllDepths)){	
		out_file_name=paste(filestem,d,".tsv",sep="")
		loc_file_name=paste(filestem,d,".locations.txt",sep="")
		
		
		if(!is.vector(AllDepths[[d]]$ingroup)){ #check in case this is a vector not a matrix
			apply(AllDepths[[d]]$ingroup,1,paste,collapse=",")->output
			for(out in out_list){
				output<-paste(output,"\t",apply(AllDepths[[d]][[out]],1,paste,collapse=","),sep="")
			}
		}else{
			paste(AllDepths[[d]]$ingroup,collapse=",")->output
			for(out in out_list){
				output<-paste(output,"\t",paste(AllDepths[[d]][[out]],collapse=","),sep="")
			}
		
		}
		
		write(output,file=out_file_name)
		write(AllDepths[[d]]$locations,file=loc_file_name,ncolumns=1)
	}
	
}

# function to cat together the input files
# need to keep an unambigous order to reconstruct things later, done via the locations and recordfiles
pack_estSFSfiles<-function(property_string){
	#find the file names
	tsvs<-list.files(pattern=paste("*_",property_string,".tsv",sep=""))
	#only do it if files with that property exist
	if(length(tsvs)>0){
		locations<-list.files(pattern=paste("*_",property_string,".locations.txt",sep=""))
		#simple sanity check
		if(length(tsvs)!=length(locations)){print("something bad happened");return(0)}
		#output file names
		outfile_name<-paste(property_string,"_combined.tsv",sep="")
		mapfile_name<-paste(property_string,"_mapfile.txt",sep="")
		recordfile_name<-paste(property_string,"_recordfile.txt",sep="")
		file.create(outfile_name)
		file.create(mapfile_name)
		file.create(recordfile_name)
		
		for(i in 1:length(tsvs)){
			#paste(c("paste", tsvs, ">",outfile_name),collapse=" ") if only!
			scan(tsvs[i],what="character",sep="\n")->temp1
			scan(locations[i],what="character",sep="\n")->temp2
			write(temp1,file=outfile_name,append=TRUE)
			write(temp2,file=mapfile_name,append=TRUE)
			cat(paste(tsvs[i],"\t",locations[i],"\t",length(temp1),"\n",sep=""),file=recordfile_name,append=TRUE)
		}
	}
}

# function to rus est-SFS
# Default to one outgroup. Always uses the GTR, but may or may not provide the GTR parameters
# Assumes DJO's edited version of est-SFS in which the GTR parameters are setable
runEstSFS<-function(DataFileName,GTRParms=NA,outgroups=1,relative_path=""){

	DataFileStem<-gsub(".tsv$","",DataFileName)
	SfsFileName<-paste(DataFileStem,".sfs",sep="")
	AncestralFileName<-paste(DataFileStem,".anc",sep="")
	#make the cfg file
	#if GTR parameters have already been estimated
	if(length(GTRParms)==6){
		cfgFileName<-paste(DataFileStem,"_FixParms.cfg",sep="")
		cat(paste("n_outgroup ",outgroups,"\n",sep=""),file=cfgFileName)
		cat("model 2\n",file=cfgFileName,append=TRUE)
		cat("nrandom -1\n",file=cfgFileName,append=TRUE)
		cat(paste("r6_0 ",GTRParms[1],"\n",sep=""),file=cfgFileName,append=TRUE)
		cat(paste("r6_1 ",GTRParms[2],"\n",sep=""),file=cfgFileName,append=TRUE)
		cat(paste("r6_2 ",GTRParms[3],"\n",sep=""),file=cfgFileName,append=TRUE)
		cat(paste("r6_3 ",GTRParms[4],"\n",sep=""),file=cfgFileName,append=TRUE)
		cat(paste("r6_4 ",GTRParms[5],"\n",sep=""),file=cfgFileName,append=TRUE)
		cat(paste("r6_5 ",GTRParms[6],sep=""),file=cfgFileName,append=TRUE)
		#run est-SFS
		##est-sfs configfile.cfg inputfile.tsv Randseedfile.txt Output.sfs Output.anc
		#be careful with teh random seed file. You might not want to run this function in parallel
		system(paste(relative_path, "est-sfs ", cfgFileName, " ",DataFileName," ",relative_path,"seedfile.txt ",SfsFileName," ",AncestralFileName,sep=""),intern=TRUE)->pars

	#if GTR parameters have NOT already been estimated	
	}else{
		cfgFileName<-paste(DataFileStem,"_EstParms.cfg",sep="")
		cat(paste("n_outgroup ",outgroups,"\n",sep=""),file=cfgFileName)
		cat("model 2\n",file=cfgFileName,append=TRUE)
		cat("nrandom 10",file=cfgFileName,append=TRUE)
		#run est-SFS
		system(paste(relative_path, "est-sfs ", cfgFileName, " ",DataFileName," ",relative_path,"seedfile.txt ",SfsFileName," ",AncestralFileName,sep=""),intern=TRUE)->pars
	}
	
	#recover the (estimated or set) GTR Parameters
	as.numeric(unlist(strsplit(pars[length(pars)]," "))[c(4,6,8,10,12,14)])->GTRParms
	
	#retrun the estimated parametres in case we need them for hte next run
	return(GTRParms)
}


unpack_estSFSfiles<-function(property_string){

	#find the pertinent file names
	list.files(pattern=paste(property_string,".+recordfile.txt",sep=""))->record_files

	#set up to record a list of the tempfiles we've output
	list_of_files<-vector()
	#for each record file
	for(recf in record_files){
		#get the name of the file with teh inferred ancestral states
		gsub("_recordfile.txt","_combined.anc",recf)->infered_anc_file
		if(!file.exists(infered_anc_file)){print("ERROR!!!! Missing ancestral file");return(0)}
		#open a connection to the ancetral state file - may be large, and want to avoid reading it into memory
		connection<-file(infered_anc_file, open = "r"); readLines(connection, n = 8, warn = FALSE)->sfs_header_we_ignore
		#read in the record file
		read.table(recf,sep="\t",as.is=TRUE)->record
		#for each record in that file
		for(i in 1:nrow(record)){
			#find the base locations represented
			scan(record$V2[i])->locs
			#get teh ancestral states
			gsub(" +","\t",gsub("^ +| +$","",readLines(connection, n = record$V3[i], warn = FALSE)))->ancs
			#decide on a temp file to write it to
			tempfilename<-paste(gsub(paste("_",property_string,".+$",sep=""),"",record$V2[i]),".tmpfile",sep="")
			#write the combiend lcoation and state to that file
			cat(paste(locs,"\t",ancs,sep=""),append=TRUE,sep="\n",file=tempfilename)
			#record the anme of the temp file we used
			list_of_files<-c(list_of_files,tempfilename)
		}
		close(connection)
	}
	
	#reopen them to sort them by base location - not vital, but may make them easier to interpret
	for(f in unique(list_of_files)){
		read.table(f,sep="\t",as.is=TRUE,header=FALSE)->locus
		write.table(locus[order(locus$V1),],sep="\t",col.names=FALSE, row.names=FALSE, quote=FALSE, file=gsub(".tmpfile",".ancestralcall",f))
		file.remove(f)
	}

}


classifySNPs<-function(fasta_file_name,outgroup,genetic_code="Standard"){
	bases<-c("A","C","G","T")
	transtable<-list("TTT"="Phe","TTC"="Phe","TTA"="Leu","TTG"="Leu",
						"CTT"="Leu","CTC"="Leu","CTA"="Leu","CTG"="Leu",
						"ATT"="Ile","ATC"="Ile","ATA"="Ile","ATG"="Met",
						"GTT"="Val","GTC"="Val","GTA"="Val","GTG"="Val",
						"TCT"="Ser","TCC"="Ser","TCA"="Ser","TCG"="Ser",
						"CCT"="Pro","CCC"="Pro","CCA"="Pro","CCG"="Pro",
						"ACT"="Thr","ACC"="Thr","ACA"="Thr","ACG"="Thr",
						"GCT"="Ala","GCC"="Ala","GCA"="Ala","GCG"="Ala",
						"TAT"="Tyr","TAC"="Tyr","TAA"="Stop","TAG"="Stop",
						"CAT"="His","CAC"="His","CAA"="Gln","CAG"="Gln",
						"AAT"="Asn","AAC"="Asn","AAA"="Lys","AAG"="Lys",
						"GAT"="Asp","GAC"="Asp","GAA"="Glu","GAG"="Glu",
						"TGT"="Cys","TGC"="Cys","TGA"="Stop","TGG"="Trp",
						"CGT"="Arg","CGC"="Arg","CGA"="Arg","CGG"="Arg",
						"AGT"="Ser","AGC"="Ser","AGA"="Arg","AGG"="Arg",
						"GGT"="Gly","GGC"="Gly","GGA"="Gly","GGG"="Gly"
					)
					
	if(genetic_code=="Invertebrate_Mitochondrial"){
		transtable[["AGA"]]<-"Ser"
		transtable[["AGG"]]<-"Ser"
		transtable[["ATA"]]<-"Met"
		transtable[["TGA"]]<-"Trp"
	}
				
	#read the sequence data in as a matrix
	raw_matrix<-toupper(as.matrix(read.alignment(fasta_file_name,format="fasta",forceToLower = FALSE)))
	#strip out the outgroup
	data_matrix<-raw_matrix[!rownames(raw_matrix)%in%outgroup,]
	#read est-sfs's call on the ancestral state
	fstem<-gsub(".fas?(ta)?$","",f)
	anc_table<-read.table(gsub(".fas?(ta)?$",".ancestralcall",f),sep="\t",header=FALSE,as.is=TRUE)[,c(1,4)]
	#make a list of minor and major calls. Major first  	
	alleles<-apply(data_matrix,2,function(x){sort(table(x[x%in%bases]),decreasing=TRUE)})
	#and note which are polymorphic and were also called by est-sfs
	SNPs<-which((sapply(alleles,length)>1)&((1:length(alleles))%in%anc_table$V1))
	
	#swap the ordering of teh SNPs int eh list is the ancestral one is not the major one. Triallelic SNPs already stripped before est-sfs 
	for(s in SNPs){ 
		if(anc_table$V4[anc_table$V1==s]<0.5){alleles[[as.character(s)]]<-rev(alleles[[as.character(s)]])}
	}
	
	#for each SNP find the ancestral codon, the derived codon, and test if they differ in their amino-acid
	SNPsClass<-vector()
	for(s in SNPs){
		depth<-sum(alleles[[s]])
		#the probability that the Major call is the ancestral one
		PKcall<-anc_table$V4[anc_table$V1==s]
		#Needs to be mirrored, becasue we swapped the order from major/minor to Ancestral/Derived
		if(PKcall<0.5){PKcall<-1-PKcall}
		#if there is a SNP at position 1, and the ancestral state was called by est-sfs
		if(s%%3==1){
			paste(c(names(alleles[[s]])[1],names(alleles[[s+1]])[1],names(alleles[[s+2]])[1]),collapse="")->anc_codon
			paste(c(names(alleles[[s]])[2],names(alleles[[s+1]])[1],names(alleles[[s+2]])[1]),collapse="")->der_codon
			if(anc_codon%in%names(transtable)&der_codon%in%names(transtable)){
				if(transtable[[anc_codon]]==transtable[[der_codon]]){
					SNPsClass<-rbind(SNPsClass,c(fstem,s,PKcall,"Synonymous",depth,alleles[[s]][2],round(100*alleles[[s]][2]/depth)))	
				}else{
					SNPsClass<-rbind(SNPsClass,c(fstem,s,PKcall,"Nonsynonymous",depth,alleles[[s]][2],round(100*alleles[[s]][2]/depth)))
				}
			}
		}
		#if there is a SNP at position 2, and the ancestral state was called by est-sfs
		if(s%%3==2){
			paste(c(names(alleles[[s-1]])[1],names(alleles[[s]])[1],names(alleles[[s+1]])[1]),collapse="")->anc_codon
			paste(c(names(alleles[[s-1]])[1],names(alleles[[s]])[2],names(alleles[[s+1]])[1]),collapse="")->der_codon
			if(anc_codon%in%names(transtable)&der_codon%in%names(transtable)){
				if(transtable[[anc_codon]]==transtable[[der_codon]]){
					SNPsClass<-rbind(SNPsClass,c(fstem,s,PKcall,"Synonymous",depth,alleles[[s]][2],round(100*alleles[[s]][2]/depth)))	
				}else{
					SNPsClass<-rbind(SNPsClass,c(fstem,s,PKcall,"Nonsynonymous",depth,alleles[[s]][2],round(100*alleles[[s]][2]/depth)))
				}
			}
		}
		#if there is a SNP at position 3, and the ancestral state was called by est-sfs
		if(s%%3==0){
			paste(c(names(alleles[[s-2]])[1],names(alleles[[s-1]])[1],names(alleles[[s]])[1]),collapse="")->anc_codon
			paste(c(names(alleles[[s-2]])[1],names(alleles[[s-1]])[1],names(alleles[[s]])[2]),collapse="")->der_codon
			if(anc_codon%in%names(transtable)&der_codon%in%names(transtable)){
				if(transtable[[anc_codon]]==transtable[[der_codon]]){
					SNPsClass<-rbind(SNPsClass,c(fstem,s,PKcall,"Synonymous",depth,alleles[[s]][2],round(100*alleles[[s]][2]/depth)))	
				}else{
					SNPsClass<-rbind(SNPsClass,c(fstem,s,PKcall,"Nonsynonymous",depth,alleles[[s]][2],round(100*alleles[[s]][2]/depth)))
				}
			}
		}
	}
	colnames(SNPsClass)<-c("Locus","Position","ProbMajorAlleleAncestral","Effect","Chromosomes","DerivedAlleles","DerivedPercAlleles")
	return(SNPsClass)
}














