#Input:
#Number of groups to test
#List of aligned fasta files from different species <- prefix must be the same between polymorphism and divergence data (same number as groups to test)
#Path to polymorphism data
#Path to divergence data
#Perform permutation test between groups?

from __future__ import division
import re
import os
import csv
import math
import sys
import shutil
import subprocess
from Bio.Phylo.PAML import yn00
from optparse import OptionParser
import random
#subprocess.check_call(['sudo yum install -y gsl-devel'],shell = True)


gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R','AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A','GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C','TGT':'C', 'TGA':'_', 'TGG':'W'}


######################OPTIONS########################################
parser = OptionParser()
parser.add_option("-i", action="store_true", default=False, help="Find individual alphas?")
parser.add_option("-r", action="store_true", default=False, help="Round SFS to nearest 10? DFE-alpha has a limit to how many SFS it can take")
parser.add_option("--pf", action="store",help="If individual alphas, this is polymorphism data filename")
parser.add_option("--df", action="store",help="If individual alphas, this is divergence data filename")
parser.add_option("-g", action="store",help="number of groups to compare, if individual alphas then leave blank",type='int')
parser.add_option("-l", action="store", help="Only for groups - list of aligned divergence/polymorphism fasta files - must have same name")
parser.add_option("-p", action="store",help="Only for groups - where are these aligned polymorphism fasta files?")
parser.add_option("-d", action="store", help="Only for groups - where are these aligned divergence fasta files?")
parser.add_option("-b", action="store_true", default=False, help="Only for groups - bootstrap values")
parser.add_option("--pt", action="store_true", default=False, help="Only for multiple groups - do a permutation test")
parser.add_option("--l2", action="store", help="Only for groups if doing a permutation test, second list. This should be the group you expect to have lower alphas/omegas, this is a one sided test")

(options, args) = parser.parse_args()

######################################################################

def fasta_to_sequences(file='test.fa'): #Input fasta and turn it into a list of 2 sequences
	f=open(file)
	f2=f.read()
	f2 = f2.replace('-','N') #RE doesnt like the dash, turn into N's
	f2 = f2.split('\n')
	f2 = filter(None, f2)
	f3 = []
	ran = range(1,len(f2))
	ran = [i for i in ran if i % 2 != 0]
	for y in ran:
		f3.append(f2[y])
	if len(f2[1]) != len(f2[3]):
		f3 = 'ALIGNMENTS NOT SAME LENGTH' #We expect the geneious alignments to be the same length
	return f3

def sliding_window(dna): #Turn sequences into list of codons
        pos1 = 0
        pos2 = 3
        x = range(0,len(dna))
        sliding_seq = []
        for each in x:
                if len(dna[pos1:pos2]) == 3:
                        sliding_seq.append(dna[pos1:pos2])
                        pos1 = pos1+3
                        pos2 = pos2 +3

        return sliding_seq

def dna_to_protein(dna, x = 0, y=0): #Written for another purpose
        if y ==0:			#x = 0,1,2 to start from other frames
                dna = dna[x:]		#y = 0, 1 - 1 if you want reverse complement
        if y == 1:
                DNA_c=dna.replace('A', 'X')
                DNA_c=DNA_c.replace('T', 'Y')
                DNA_c=DNA_c.replace('X','T')
                DNA_c=DNA_c.replace('Y','A')
                DNA_c=DNA_c.replace('C','X')
                DNA_c=DNA_c.replace('G','Y')
                DNA_c=DNA_c.replace('X','G')
                DNA_c=DNA_c.replace('Y','C')
                dna = DNA_c[::-1]
                dna = dna[x:]
        dna = sliding_window(dna)
        protein = []
        for each in dna:
                for nuc, aa in gencode.items():
                        if each == nuc:
                                protein.append(aa)
                if each.count('N') > 0:
                        protein.append('*')
        return "".join(protein)

def multifasta2doublefasta(file):
	f = open(file)
	fasta = f.read()
	fasta =  fasta.split('\n')
	fasta = filter(None, fasta)

	Ncounts = []
	dashcounts = []
	for line in fasta:
		if re.search(r'^>', line):
			t = 3 #placeholder
		else:
			x = line.count('N')
			y = line.count('-')
			Ncounts.append(x+y)
			
	lowN = min(Ncounts)
	if Ncounts.count(lowN) < 2:
		print('Minumum value of Ns in sequence is not in more than one entry')
		print(f)

	f.close()
	fasta_new = []
	for line in fasta:
		if re.search(r'^>', line):
			t = 3 #placeholder
		else:
			if line.count('N') + line.count('-') == lowN:
				fasta_new.append(fasta[fasta.index(line)-1])
				fasta_new.append(line)
				fasta_new.append(fasta[fasta.index(line)-1])
				fasta_new.append(line)
			elif line.count('N') + line.count('-') < lowN+lowN:
				fasta_new.append(fasta[fasta.index(line)-1])
				fasta_new.append(line)

	fasta_new = fasta_new[0:4]
	return fasta_new





def dN_dS(Seq1, Seq2):
	dS = 0
	dN = 0
	if re.search(r"[KMRYSWBVHDXN]", Seq1): #Exclude codons with ambiguities
		print('')
	else:
		if re.search(r"[KMRYSWBVHDXN]", Seq2): #Exclude codons with anbiguities
			print('')
		else:		
			if Seq1 != Seq2: #If they arent the same codon
				if dna_to_protein(Seq1) == dna_to_protein(Seq2): #But they are the same aa
					dS = dS +1	#Add one to dS
				else: #If not same protein, count number of differences
					u = zip(Seq1,Seq2)
					counter = 0			
					for k,j in u:
						if k!=j:
							counter = counter + 1
					if counter == 1: #If only one diff, then increase dN by 1
						dN = dN + 1
					if counter == 2:
						#Now have to figure out whether 2 nonsynonymous or 1 syn and 1 nonsyn
						x = Seq1
						y = Seq2
						a = y[0] + x[1] + y[2]
						b = x[0] + y[1] + y[2]
						c = x[0] + x[1] + y[2]
						d = y[0] + y[1] + x[2]
						e = y[0] + x[1] + x[2]#All one nucleotide differences between Seq1 and Seq2
						f = x[0] + y[1] + x[2]#If any are synonymous, assume this occured, and nonsynonymous as the other
						all_poss = [a,b,c,d,e,f]
						all_poss = set(all_poss)
						all_aa = []
						for codon in all_poss:
							all_aa.append(dna_to_protein(codon))
						if len(all_aa) == len(set(all_aa)):
							dN = dN + 2
						else:
							dN = dN + 1
							dS = dS + 1
	return dN,dS

def to_precision(x,p):
#From http://randlet.com/blog/python-significant-figures-format/
    x = float(x)
    if x == 0.:
        return "0." + "0"*(p-1)
    out = []
    if x < 0:
        out.append("-")
        x = -x
    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)
    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)
    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1
    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1
    m = "%.*g" % (p, n)
    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)
    return "".join(out)

def round_to(n, precision): # Will round those MAF to nearest 1/maxNum of samples
    correction = 0.5 if n >= 0 else -0.5
    return int( n/precision+correction ) * precision

#Function for nonsyn and syn SFS
def fasta_to_SFS(file):
	alignment = fasta_to_sequences(file)
	align_codons = []
	for each in alignment:
		align_codons.append(sliding_window(each))

	align_codons = map(list, zip(*align_codons)) # results in a listOFlists-each element being a posmatch codon
	#Remove codons with N's
	aligncod_wo_Ns = []
	for each in align_codons:
		temp = []
		for all in each:
			if re.search(r'N', all):
				print(all)
			else:
				temp.append(all)
		if len(temp) > len(each)/2: # If more than half are missing data, don't use the codon
			aligncod_wo_Ns.append(temp)
	invariant = 0
	pn_freq = []
	ps_freq = []
	sampled_alleles = []
	for codon in aligncod_wo_Ns:
		sampled_alleles.append(len(codon))
		diff_codon = list(set(codon))
		if len(diff_codon) == 1:
			invariant = invariant + 3
		if len(diff_codon) == 2: #Per codon position, some have 2 or 3 different codons
			invariant = invariant + 2
			MAF = codon.count(diff_codon[0])
			if MAF > len(codon)/2:
				MAF = len(codon)-MAF #If the codon is the major allele, get freq dor minor allele (MAF)
			pNpS = dN_dS(diff_codon[0],diff_codon[1])
			if pNpS[0] == 1:
				pn_freq.append(MAF/len(codon))
			if pNpS[1] == 1:
				ps_freq.append(MAF/len(codon))
			if pNpS[0] == 2: # for now assume that two hit codons don't exist
				print('2 Changes per nonsynonymous codons')
				print(diff_codon)
			if pNpS[1] == 2: # for now assume that two hit codons dont exist
				print('2 Changes per nonsynonymous codons')
				print(diff_codon)		
		if len(diff_codon) == 3: # for now assume that two hit codons don't exist
			print('3 different codons at same position!!!!!!')
			print(diff_codon)
			#y = map(list, zip(*diff_codon))
			#codon_sites = map(list, zip(*codon))
			#for each in y:
			#	y[y.index(each)] = ''.join(each)
			#for each in y:
			#	if len(set(each)) > 1:
			#		SNPs = list(set(each))
			#			for site in SNPs:
			#				codon_sites[y.index(each)].count(site)

	MAFbins = range(1,max(sampled_alleles)+1)
	temp = []
	for each in MAFbins:
		temp.append(each/ max(sampled_alleles))
	MAFbins = temp # bins f
	del temp
	temp = []
	for each in MAFbins: #Folded SFS
		if each < 0.5:
			temp.append(each)
	MAFbins = temp 
	del temp
	countsbins = [0]*len(MAFbins)
	temp = []
	for each in pn_freq: #Round each to nearest bin and output 5 sig figs
		temp.append(to_precision(round_to(each, 1/max(sampled_alleles)), 5))
	pn_freq = temp
	temp = []
	for each in ps_freq: #Round each to nearest bin and output 5 sig figs
		temp.append(to_precision(round_to(each, 1/max(sampled_alleles)), 5))
	ps_freq = temp
	temp = []
	for each in MAFbins:	
		temp.append(to_precision(each, 5))
	MAFbins = temp
	del temp
	PNcountsbins = [0]*len(MAFbins)
	PScountsbins = [0]*len(MAFbins)
	for bin in MAFbins: # Computes SFS of counts
		bincount_pn = pn_freq.count(bin)
		PNcountsbins[MAFbins.index(bin)] = bincount_pn
		bincount_ps = ps_freq.count(bin)
		PScountsbins[MAFbins.index(bin)] = bincount_ps
	#Now for invariant sites
	if file.count('/') > 0:
		fastafile = file.split('/')[-1]
	else: 
		fastafile = file
	shutil.copyfile(file, 'FAtmp/'+fastafile)
	#for file in os.listdir('FAtmp'):#CAREFUL! Replaces fasta files in FAtmp
	name = fastafile.split('.')[0]
	x = multifasta2doublefasta('FAtmp/'+fastafile)
	with open('FAtmp/'+name+'.fasta', 'w') as f:
		f.write('\n'.join(x))
	subprocess.check_call(['sh make_paml.sh'],shell = True)
	yn = yn00.Yn00(alignment='./FAtmp/'+name+'.seq',out_file = "delete.out",working_dir = ".")
	yn.set_options(weighting=1, commonf3x4= 0, icode = 0, verbose = 0, ndata = 1)
	data = yn.run(verbose = True)
	first = data.keys()[0]
	second =data[first].keys()[0]
	datainterest = data[first][second]['YN00']
	#Dn LnD Ds LsD
	#final_data = [name, str(datainterest['dN']),str(datainterest['N']),str(datainterest['dS']), str(datainterest['S'])]
	invariantN = datainterest['N'] - sum(PNcountsbins)
	invariantS = datainterest['S'] - sum(PScountsbins)
	PNcountsbins.insert(0,str(invariantN))
	PScountsbins.insert(0,str(invariantS))
	for files in os.listdir('FAtmp'):
		os.remove('FAtmp/'+ files)
	#Add zeroes to list until its the right size, since SFS is folded 
	z = max(sampled_alleles) - len(PNcountsbins)
	for i in range(0,z):
		PNcountsbins.append(0)
		PScountsbins.append(0)
	return [PNcountsbins, PScountsbins]
	
def round_down(num, divisor):
    return num - (num%divisor)


#Make sfs.txt
def make_sfs(ind_alpha = options.pf, group_alpha = options.l):
	#If individual alpha, then set number of samples to one
	if options.i:
		number_of_samples = 1
	else:
		f = open(options.l)
		number_of_samples = 0
		for line in f:
			number_of_samples = number_of_samples + 1
	if options.i:
		sfstxt = [[number_of_samples]]
		SFS = fasta_to_SFS(ind_alpha)
		alleles_sampled = len(SFS[0]) - 1
		sfstxt.append([alleles_sampled])
		sfstxt.append(SFS[0])
		sfstxt.append(SFS[1])
		sfstxt_final = sfstxt
	else:
		sfstxt = [[number_of_samples]]
		f = open(group_alpha)
		for each in f:
			each=each.rstrip()
			path_file = options.p+each
			SFS = fasta_to_SFS(path_file)
			alleles_sampled = len(SFS[0]) - 1
			sfstxt.append([alleles_sampled])
			sfstxt.append(SFS[0])
			sfstxt.append(SFS[1])
			sfstxt_final = sfstxt
	#Reduce number of SFS with same number of individuals
		number_of_sfslines = int(sfstxt[0][0])*3
		range_alleles = range(1,number_of_sfslines,3)
		sampled_list = []
		for i in range_alleles:
			sampled_list.append(sfstxt[i][0])
		sampled_list = list(set(sampled_list))
		sfstxt_final = [[len(sampled_list)]]
		for each in sampled_list:
			indices = [i for i, x in enumerate(sfstxt) if x == [each]] #Use as an index for all values
			neutral = []
			selected = []
			for i in indices:
				selected.append(sfstxt[i+1])
				neutral.append(sfstxt[i+2])
			selected = map(list, zip(*selected))
			neutral = map(list, zip(*neutral))
			selected_final = []
			neutral_final = []
			for j in range(0, len(selected)):
				selected_final.append(sum(float(k) for k in selected[j]))
			for j in range(0, len(neutral)):
				neutral_final.append(sum(float(k) for k in neutral[j]))
			sfstxt_final.append([each])
			sfstxt_final.append(selected_final)
			sfstxt_final.append(neutral_final)
#Downsample if too many SFS
	if options.r:
		number_of_sfslines = int(sfstxt_final[0][0])*3
		range_alleles = range(1,number_of_sfslines,3)
		sampled_list = []
		for i in range_alleles:
			sampled_list.append(sfstxt_final[i][0])

		for i in range(0,len(sampled_list)):
			if sampled_list[i] > 10:
				down = round_down(sampled_list[i], 10)
				bins = []
				for j in range(0, down+1):
					bins.append(to_precision(j/down,5))

				nonsyn = sfstxt_final[sfstxt_final[1:].index([sampled_list[i]])+2] # index for nonsyn SFS
				syn = sfstxt_final[sfstxt_final[1:].index([sampled_list[i]])+3] # index for syn SFS
				oldbins = []
				for bn in range(1,len(syn)):
					oldbins.append(to_precision(bn/len(syn[1:]), 5))
				bin_ref = []
				for each in oldbins:
					bin_ref.append(to_precision(round_to(float(each), 1/down),5))
				new_syn = []
				new_nonsyn = []
				for p in range(0,len(bin_ref))[::-1]:
					if bin_ref.count(bin_ref[p]) == 1:
						new_syn.insert(0,syn[p+1])
						new_nonsyn.insert(0,nonsyn[p+1])
						del bin_ref[p]
					elif bin_ref.count(bin_ref[p]) > 1:
						indices = [i for i, x in enumerate(bin_ref) if x == bin_ref[p]]
						if p < max(indices):
							summ = 0
							summ2 = 0
							for indeck in indices:
								summ = summ + int(syn[p+1])
								summ2 = summ2 + int(nonsyn[p+1])
							new_syn.insert(0,float(summ))
							new_nonsyn.insert(0,float(summ))
							del bin_ref[p]
				new_syn.insert(0, syn[0])
				new_nonsyn.insert(0,nonsyn[0])
				sfstxt_final[sfstxt_final.index(syn)] = new_syn
				sfstxt_final[sfstxt_final.index(nonsyn)] = new_nonsyn

		for i in range_alleles:
			if sfstxt_final[i][0] > 10:
				sfstxt_final[i] = [round_down(sfstxt_final[i][0], 10)]
		sfstxt = sfstxt_final

		number_of_sfslines = int(sfstxt[0][0])*3
		range_alleles = range(1,number_of_sfslines,3)
		sampled_list = []
		for i in range_alleles:
			sampled_list.append(sfstxt[i][0])
		sampled_list = list(set(sampled_list))
		sfstxt_final = [[len(sampled_list)]]
		for each in sampled_list:
			indices = [i for i, x in enumerate(sfstxt) if x == [each]] #Use as an index for all values
			neutral = []
			selected = []
			for i in indices:
				selected.append(sfstxt[i+1])
				neutral.append(sfstxt[i+2])
			selected = map(list, zip(*selected))
			neutral = map(list, zip(*neutral))
			selected_final = []
			neutral_final = []
			for j in range(0, len(selected)):
				selected_final.append(sum(float(k) for k in selected[j]))
			for j in range(0, len(neutral)):
				neutral_final.append(sum(float(k) for k in neutral[j]))
			sfstxt_final.append([each])
			sfstxt_final.append(selected_final)
			sfstxt_final.append(neutral_final)

	return sfstxt_final


#Make config file for neutral site class
data_path_1='/home/s1348847/Programs/dfe-alpha-release-2.14/data'
sfs_input_file='./sfs.txt'
est_dfe_results_dir='.'
site_class='0'
fold='1'
epochs='2'
search_n2='1'
t2_variable='1'
t2='50'
config0 = [['data_path_1'+'\t'+data_path_1],['sfs_input_file'+'\t'+sfs_input_file],['est_dfe_results_dir'+'\t'+est_dfe_results_dir],['site_class'+'\t'+site_class],['fold'+'\t'+fold],['epochs'+'\t'+epochs],['search_n2'+'\t'+search_n2],['t2_variable'+'\t'+t2_variable],['t2'+'\t'+t2]]
with open('config0.txt','wb') as f:
	writer = csv.writer(f)
	writer.writerows(config0)

#Make config file for selected site class
data_path_1='/home/s1348847/Programs/dfe-alpha-release-2.14/data'
sfs_input_file='./sfs.txt'
est_dfe_results_dir='.'
est_dfe_demography_results_file='est_dfe.out'
site_class='1'
fold='1'
epochs='2'
mean_s_variable='1'
mean_s='-0.1'
beta_variable='1'
beta='0.5'
config1 = [['data_path_1'+'\t'+data_path_1],['sfs_input_file'+'\t'+sfs_input_file],['est_dfe_results_dir'+'\t'+est_dfe_results_dir],['est_dfe_demography_results_file'+'\t'+est_dfe_demography_results_file],['site_class'+'\t'+site_class],['fold'+'\t'+fold],['epochs'+'\t'+epochs],['mean_s_variable'+'\t'+mean_s_variable],['mean_s'+'\t'+mean_s],['beta_variable'+'\t'+beta_variable],['beta'+'\t'+beta]]
with open('config1.txt','wb') as f:
	writer = csv.writer(f)
	writer.writerows(config1)


#Make alpha-omega config file
data_path_1='/home/s1348847/Programs/dfe-alpha-release-2.14/data/'
divergence_file='./div.txt'
est_alpha_omega_results_file='./est_alpha_omega.out'
est_dfe_results_file='./est_dfe.out'
neut_egf_file='./neut_egf.out'
sel_egf_file='./sel_egf.out'
do_jukes_cantor='0'
remove_poly='0'
config_alpha = [['data_path_1'+'\t'+data_path_1],['divergence_file'+'\t'+divergence_file],['est_alpha_omega_results_file'+'\t'+est_alpha_omega_results_file],['est_dfe_results_file'+'\t'+est_dfe_results_file],['neut_egf_file'+'\t'+neut_egf_file],['sel_egf_file'+'\t'+sel_egf_file],['do_jukes_cantor'+'\t'+do_jukes_cantor],['remove_poly'+'\t'+remove_poly]]
with open('config_alpha.txt','wb') as f:
	writer = csv.writer(f)
	writer.writerows(config_alpha)

sfstxt_final=make_sfs()		
with open('sfs.txt','wb') as f:
	writer = csv.writer(f, delimiter=' ')
	writer.writerows(sfstxt_final)

#Run est-dfe
subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config0.txt'],shell = True)
subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config1.txt'],shell = True)

#Make divergence file
def fasta_to_paml(file):
	if file.count('/') > 0:
		fastafile = file.split('/')[-1]
	else: 
		fastafile = file
	shutil.copyfile(file, 'FAtmp/'+fastafile)
	name = fastafile.split('.')[0]
	subprocess.check_call(['sh make_paml.sh'],shell = True)
	yn = yn00.Yn00(alignment='./FAtmp/'+name+'.seq',out_file = "delete.out",working_dir = ".")
	yn.set_options(weighting=1, commonf3x4= 0, icode = 0, verbose = 0, ndata = 1)
	data = yn.run(verbose = True)
	first = data.keys()[0]
	second =data[first].keys()[0]
	datainterest = data[first][second]['YN00']
	#Dn LnD Ds LsD
	final_data = [name, str(datainterest['dN']),str(datainterest['N']),str(datainterest['dS']), str(datainterest['S'])]
	for files in os.listdir('FAtmp'):
		os.remove('FAtmp/'+ files)
	return final_data

def make_divtxt(genelist=options.l):
	f = open(genelist)
	all_div = []
	for line in f:
		#if os.path.isfile('FAtmp/'+line):
		line=line.rstrip()
		all_div.append(fasta_to_paml(options.d+line))
	Ln = []
	Ls = []
	Dn = []
	Ds = []
	for each in all_div:
		Ln.append(float(each[2]))
		Ls.append(float(each[4]))
		Dn.append(float(each[2])*float(each[1]))
		Ds.append(float(each[4])*float(each[3]))
	Ln = sum(Ln)
	Ls = sum(Ls)
	Dn = sum(Dn)
	Ds = sum(Ds)
	divtxt = [['1'+'\t'+str(Ln)+'\t'+str(Dn)],['0'+'\t'+str(Ls)+'\t'+str(Ds)]]
	return divtxt


#If individual then find divergence of only that one gene, if a group given in a list, then add them all up
if options.i:
	div = fasta_to_paml(options.df)
	Ln = div[2]
	Ls = div[4]
	Dn = float(div[2])*float(div[1])
	Ds = float(div[4])*float(div[3])
	divtxt = [['1'+'\t'+str(Ln)+'\t'+str(Dn)],['0'+'\t'+str(Ls)+'\t'+str(Ds)]]
else:
	divtxt = make_divtxt()

with open('div.txt','wb') as f:
	writer = csv.writer(f)
	writer.writerows(divtxt)
subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_alpha_omega -c config_alpha.txt'],shell = True)

alpha_omega_final = []
file = open('est_alpha_omega.out')
for line in file:
	alpha_omega_final.append(line)
file.close()
	
if options.b:
	bootstrap = []
	for z in range(0,1000):
		file = open(options.l)
		file = file.read()
		file = file.split('\n')
		file = filter(None, file)
		random_list = []
		for i in range(0,len(file)):
			random_list.append(random.choice(file))
		with open('random_list.txt', 'wb') as f:
			f.write("\n".join(map(str,random_list)))
		sfstxt_final=make_sfs('random_list.txt')		
		with open('sfs.txt','wb') as f:
			writer = csv.writer(f, delimiter=' ')
			writer.writerows(sfstxt_final)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config0.txt'],shell = True)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config1.txt'],shell = True)
		divtxt = make_divtxt('random_list.txt')
		with open('div.txt','wb') as f:
			writer = csv.writer(f)
			writer.writerows(divtxt)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_alpha_omega -c config_alpha.txt'],shell = True)
		file2 = open('est_alpha_omega.out')
		for line in file2:
			bootstrap.append(line)
		file2.close()
	
	bootstrap2 = []
	for each in bootstrap:
		each = each.rstrip()
		bootstrap2.append(each.split(' '))
	bootstrap2 = map(list, zip(*bootstrap2))
	alphas = bootstrap2[5]
	alphas = sorted(alphas, key=float)
	omegaA = bootstrap2[7]
	omegaA = sorted(omegaA, key=float)
	print("The alpha 95% C.I. is from "+ alphas[24]+ " to " + alphas[-25])
	print("The omega A 95% C.I. is from "+ omegaA[24]+ " to " + omegaA[-25])
	x="The alpha 95% C.I. is from "+ alphas[24]+ " to " + alphas[-25]
	y="The omega A 95% C.I. is from "+ omegaA[24]+ " to " + omegaA[-25]
	if len(bootstrap) < 999:
		x = 'Error in bootstrapping'
		y = 'Error in bootstrapping'
	out = [alpha_omega_final,[x],[y]]
	out_name = options.l.split('.')[0]
	out_name=out_name+'_bootstrap_results.txt'
	with open(out_name, 'wb') as f:
		writer = csv.writer(f)
		writer.writerows(out)

if options.pt:
	l_alpha = float(alpha_omega_final[0].rstrip().split(' ')[5])
	l_omega = float(alpha_omega_final[0].rstrip().split(' ')[7])
	sfstxt_final2=make_sfs(options.l2)
	with open('sfs.txt','wb') as f:
			writer = csv.writer(f, delimiter=' ')
			writer.writerows(sfstxt_final2)
	subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config0.txt'],shell = True)
	subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config1.txt'],shell = True)
	divtxt2 = make_divtxt(options.l2)
	with open('div.txt','wb') as f:
			writer = csv.writer(f)
			writer.writerows(divtxt2)
	subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_alpha_omega -c config_alpha.txt'],shell = True)
	alpha_omega_final2 = []
	file2 = open('est_alpha_omega.out')
	for line in file2:
		alpha_omega_final2.append(line)
	l2_alpha = float(alpha_omega_final2[0].rstrip().split(' ')[5])
	l2_omega = float(alpha_omega_final2[0].rstrip().split(' ')[7])
	alpha_diff = l_alpha - l2_alpha
	omega_diff = l_omega - l2_omega
	permute_differences_alpha = []
	permute_differences_omega = []
	for z in range(0,1000):
		firstlist = open(options.l)
		firstlist = firstlist.read()
		firstlist = firstlist.split('\n')
		firstlist = filter(None, firstlist)
		seclist = open(options.l2)
		seclist = seclist.read()
		seclist = seclist.split('\n')
		seclist = filter(None, seclist)
		wholelist = firstlist + seclist
		random_list1 = []
		random_list2 = []
		for i in range(0,len(firstlist)):
			random_list1.append(random.choice(wholelist)) #Mkae two fake lists the same length as each group but with mixed categories
		with open('random_list1.txt', 'wb') as f:
			f.write("\n".join(map(str,random_list1)))
		for i in range(0,len(seclist)):
			random_list2.append(random.choice(wholelist))
		with open('random_list2.txt', 'wb') as f:
			f.write("\n".join(map(str,random_list2)))

		sfstxt_final1=make_sfs('random_list1.txt')
		with open('sfs.txt','wb') as f:
				writer = csv.writer(f, delimiter=' ')
				writer.writerows(sfstxt_final1)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config0.txt'],shell = True)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config1.txt'],shell = True)
		divtxt1 = make_divtxt('random_list1.txt')
		with open('div.txt','wb') as f:
				writer = csv.writer(f)
				writer.writerows(divtxt1)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_alpha_omega -c config_alpha.txt'],shell = True)
		alpha_omega_final_ran1 = []
		file1 = open('est_alpha_omega.out')
		for line in file1:
			alpha_omega_final_ran1.append(line)
		ran1_alpha = float(alpha_omega_final_ran1[0].rstrip().split(' ')[5])
		ran1_omega = float(alpha_omega_final_ran1[0].rstrip().split(' ')[7])

		sfstxt_final2=make_sfs('random_list2.txt')
		with open('sfs.txt','wb') as f:
				writer = csv.writer(f, delimiter=' ')
				writer.writerows(sfstxt_final2)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config0.txt'],shell = True)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_dfe -c config1.txt'],shell = True)
		divtxt2 = make_divtxt('random_list2.txt')
		with open('div.txt','wb') as f:
				writer = csv.writer(f)
				writer.writerows(divtxt2)
		subprocess.check_call(['/home/s1348847/Programs/dfe-alpha-release-2.14/est_alpha_omega -c config_alpha.txt'],shell = True)
		alpha_omega_final_ran2 = []
		file2 = open('est_alpha_omega.out')
		for line in file2:
			alpha_omega_final_ran2.append(line)
		ran2_alpha = float(alpha_omega_final_ran2[0].rstrip().split(' ')[5])
		ran2_omega = float(alpha_omega_final_ran2[0].rstrip().split(' ')[7])
		permute_differences_alpha.append(ran1_alpha - ran2_alpha)
		permute_differences_omega.append(ran1_omega - ran2_omega)

	permute_differences_alpha = sorted(permute_differences_alpha)
	permute_differences_omega = sorted(permute_differences_omega)
#Now find the position closest to the sorted list


	alpha_pvalue = to_precision(permute_differences_alpha.index(min(permute_differences_alpha, key=lambda x:abs(x-alpha_diff))),4)
	alpha_pvalue = (1000 - int(float(alpha_pvalue)))/1000
	omega_pvalue = to_precision(permute_differences_omega.index(min(permute_differences_omega, key=lambda x:abs(x-omega_diff))),4)
	omega_pvalue = (1000 - int(float(omega_pvalue)))/1000

	x="The alpha pvalue for " + options.l +" and " + options.l2 + " is "+ str(alpha_pvalue)
	y="The omega pvalue for " + options.l + " and " + options.l2 + " is "+ str(omega_pvalue)
	if len(permute_differences_alpha) < 999:
		x = 'Error in bootstrapping'
		y = 'Error in bootstrapping'
	out = [alpha_omega_final,alpha_omega_final2,[x],[y]] #Give alpha values and omega values for both lists, and the p values
	out_name = options.l.split('.')[0] + options.l2.split('.')[0].split('/')[-1]
	out_name=out_name+'_permutation_results.txt'
	with open(out_name, 'wb') as f:
		writer = csv.writer(f)
		writer.writerows(out)
