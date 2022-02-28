#!/bin/sh

sample=$1
describer=$(echo ${sample} | sed 's/.aligned_masked.fas//')

mafft --add ${describer}.fa --reorder ${describer}.aligned_masked.fas > ${describer}.fasta
