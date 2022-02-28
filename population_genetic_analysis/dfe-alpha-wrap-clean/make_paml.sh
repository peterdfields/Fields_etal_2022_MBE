#!/bin/bash
for f in ./FAtmp/*
do
	echo ${f%.*}
	echo $f
	perl fasta2phylip.pl $f > ${f%.*}.phy
	perl phylip2paml.pl ${f%.*}.phy > ${f%.*}.seq
done 

