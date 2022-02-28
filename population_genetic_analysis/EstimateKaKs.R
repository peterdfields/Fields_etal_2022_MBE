unlist(strsplit(commandArgs(trailingOnly = TRUE)," "))[1]->fname
#long for for testing dsimilis_genes/rna8337.fa
name_stem<-gsub(".+/|.fa.*","",fname)
print(name_stem)
tname<-paste(name_stem,".temp.fas",sep="")
aname<-paste(name_stem,".aligned.fas",sep="")


require(seqinr) #package for read/write fasta

paste("dmagna_genes/",name_stem,".fasta",sep="")->m_nam
paste("dsimilis_genes/",name_stem,".fasta",sep="")->s_nam

#read in the files


if(file.exists(s_nam)){
        read.fasta(s_nam)[[1]]->sm
}else{
        stop("No File1")
}

if(file.exists(m_nam)){
        read.fasta(m_nam)[[1]]->mg
}else{
        stop("No File2")
}



#function to return the DNA in the orientation which gives the longest ORF.

#only problem will occur if several are equally good
BestFrame<-function(DNA){
        getAnnot(DNA)->n
                #translate in all frames
        prots<-list(
                pF0=paste(translate(DNA,frame = 0, sens = "F"),collapse=""),
                pF1=paste(translate(DNA,frame = 1, sens = "F"),collapse=""),
                pF2=paste(translate(DNA,frame = 2, sens = "F"),collapse=""),
                pR0=paste(translate(DNA,frame = 0, sens = "R"),collapse=""),
                pR1=paste(translate(DNA,frame = 1, sens = "R"),collapse=""),
                pR2=paste(translate(DNA,frame = 2, sens = "R"),collapse="")
        )
                #record (padded) DNA in all frames
        dns<-list(
                dF0=paste(c(DNA),collapse=""),
                dF1=paste(c("-","-",DNA),collapse=""),
                dF2=paste(c("-",DNA),collapse=""),
                dR0=paste(c(rev(toupper(comp(DNA)))),collapse=""),
                dR1=paste(c("-","-",rev(toupper(comp(DNA)))),collapse=""),
                dR2=paste(c("-",rev(toupper(comp(DNA)))),collapse="")
        )

                #split on stop codons
        lapply(prots,strsplit,"\\*")->peps


                #look through all translations to find the longest, or in case of ties, the one that starts with an "M"
        pep<-""
        dna<-0
        pepl<-0
        fr<-0
        for (i in 1:6){ #for each frame
                tmp<-unlist(peps[[i]]) #get the peptides
                for (j in 1:length(tmp)){ #for each peptide
                        if((nchar(tmp[j])>=pepl)&(substr(pep,1,1)=="M")){ #if its as long or longer than our longest and starts with an M, keep it
							pep<-tmp[j] #keep it
							pepl<-nchar(pep) #record length
							fr<-i #record frame
                        }else{ #otherwise, if it's longer 9regardles of start codon) keep it
							if(nchar(tmp[j])>pepl){
								pep<-tmp[j]
								pepl<-nchar(pep)
								fr<-i
							}
                        }
                }
        }

		if(nchar(dns[[fr]])%%3==0){
			return(dns[[fr]]) #return the padded DNA in teh best frame
		}
		if(nchar(dns[[fr]])%%3==1){
			return(paste(dns[[fr]],"-","-",sep="") )#return the padded DNA in teh best frame
		}
		if(nchar(dns[[fr]])%%3==2){
			return(paste(dns[[fr]],"-",sep="") )#return the padded DNA in teh best frame
		}		
}

BestFrame(mg)->bmg
BestFrame(sm)->bsm

#generate temporaty file to align
cat(paste(">",name_stem,"_magna",sep=""),"\n",file=tname,append=FALSE)
cat(bmg,"\n",file=tname,append=TRUE)
cat(paste(">",name_stem,"_similis",sep="","\n"),file=tname,append=TRUE)
cat(bsm,file=tname,append=TRUE)

#count how many non-terminal stops there are
sum(rev(translate(strsplit(bmg,"")[[1]],frame = 0, sens = "F"))[-1]=="*")->stop_a
sum(rev(translate(strsplit(bsm,"")[[1]],frame = 0, sens = "F"))[-1]=="*")->stop_b

#run PRANK alignments, outputting in format compatible with PAML
if((stop_a+stop_b)==0){ #if our frame has no non-terminal stops, align as codons
		system(paste("prank -codon -noanchors -once -nomafft -f=paml -d=",tname," -o=",name_stem,sep=""))
}else{ #if our frame has non-terminal stops, align as codons
        cat(name_stem,append=TRUE,file="Aligned_For_Stops.txt")
        #align
        system(paste("prank -DNA -noanchors -once -nomafft -f=paml -d=",tname," -o=",name_stem,sep=""))
}
#and make a fasta version
system(paste("prank -convert -d=",name_stem,".best.phy -o=",name_stem,".aligned",sep=""))


#make a PAML ctrl file
#system(paste("cat codeml.ctl | sed 's/XXX/",name_stem,".best.phy/' | sed 's/YYY/",name_stem,".paml_out/' > codeml.ctl",sep=""))

#run paml
system(paste("codeml",sep=""))

#extract the line of interest from paml
system(paste("grep '^t=' ",name_stem,".paml_out",sep=""), intern=TRUE)->dnds_string
print(dnds_string)
c(name_stem,unlist(strsplit(dnds_string,"[ =]+"))[c(2,4,6,10,12,14)])->dnds
cat(dnds,"\n",sep="\t",file="PAML_DnDs_Estimates.tsv",append=TRUE)

#now we need to identify which sites in d.magna are represented in the dn/ds estimation
read.fasta(paste(name_stem,".aligned.fas",sep=""))->aln
aln[[1]]->mag
aln[[2]]->sim
#replace with N's those positions of mag that are not bases in mag or are not bases in sim
mag[!(sim%in%c("a","g","c","t")&mag%in%c("a","g","c","t"))]<-"N"
#replace partial codons with 3N's as codeml will have ignored these
as.vector(apply(t(matrix(mag,nrow=3)),1,function(x){if("N"%in%x){rep("N",3)}else{x}}))->mag #I'll leave figuring how this is done as an exercise to the reader.

#need to check the direction
if(gsub("-","",bmg)==toupper(paste(mg,collapse=""))){ #if it wasn't turned around when finding the best frame
	strsplit(paste(mag,collapse=""),"[Nn]+")[[1]]->chunks #these are the chunks that are used to calculate dnds
	paste(mg,collapse="")->mag_full #this is the full magna sequence we read in, now we need the locations of the 'chunk' substrings
}else{ #if it /was/ turned around then reverse the sequences. Since we don't care what the code actually reads, no need to translate as well
	strsplit(paste(rev(comp(mag)),collapse=""),"[Nn]+")[[1]]->chunks #these are the chunks that are used to calculate dnds
	paste(mg,collapse="")->mag_full #this is the full magna sequence we read in, now we need the locations of the 'chunk' substrings
}
	
	
indexer<-0
for(i in 1:length(chunks)){

	if(nchar(chunks[i])>0){ #want to ignore empty chunks
		#find first match in string as we new they're in order
		regexpr(chunks[i], mag_full, fixed = TRUE)->local_strt
		local_strt+nchar(chunks[i])-1->local_stp
		
		#remove string to that point to avoid finding a shorter match from later in the list early in the string
		mag_full<-substr(mag_full,local_stp+1,nchar(mag_full))

		#record absolute positions
		local_strt+indexer->strt
		local_stp+indexer->stp

		cat(strt,"\t",stp,"\n",sep="",append=TRUE,file=paste(name_stem,".locations.tsv",sep=""))
		
		
		#update indexer
		indexer<-stp
	}
}
