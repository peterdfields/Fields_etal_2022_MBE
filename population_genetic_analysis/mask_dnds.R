setwd(".")
unlist(strsplit(commandArgs(trailingOnly = TRUE)," "))[1]->fname
#long for for testing dsimilis_genes/rna8337.fa
name_stem<-gsub(".+/|.fa.*","",fname)
print(name_stem)


require(seqinr) #package for read/write fasta

#read in the files

system(paste("prank -d=",name_stem,".fas -o=",name_stem,".aligned"," -translate -F",sep=""))
system(paste("perl ./pal2nal.pl ",name_stem,".aligned.best.pep.fas ",name_stem,".aligned.best.nuc.fas -output paml -nogap > for_paml/",name_stem,".codon",sep=""))
