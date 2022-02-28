setwd(".")
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
system(paste("prank -d=",name_stem,".aligned.fas -o=",name_stem,".aligned"," -translate -F",sep=""))
system(paste("perl pal2nal.v14/pal2nal.pl ",name_stem,".aligned.best.pep.fas ",name_stem,".aligned.best.nuc.fas -output paml -nogap > for_paml/",name_stem,".codon",sep=""))



