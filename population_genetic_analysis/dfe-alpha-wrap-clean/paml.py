from Bio.Phylo.PAML import yn00
import os
all_data = []
for each in os.listdir('./FAtmp/'):
	print(each)
	if each.split('.')[1] == 'seq':
		print(each)

		yn = yn00.Yn00(alignment = './FAtmp/' + each, out_file = "delete.out", working_dir = ".")

		yn.set_options(weighting=1, commonf3x4= 0, icode = 0, verbose = 0, ndata = 1)

		data = yn.run(verbose = True)

		first = data.keys()[0]
		second =data[first].keys()[0]
		datainterest = data[first][second]['YN00']

		name = each.split('.')[0]
		#Dn LnD Ds LsD
		final_data = [name, str(datainterest['dN']),str(datainterest['N']),str(datainterest['dS']), str(datainterest['S'])]
		all_data.append(final_data)

import csv
with open("old_rnai_div.csv","wb") as f:
	writer = csv.writer(f)
	writer.writerows(all_data)
