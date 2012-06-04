#!/usr/bin/env python
# encoding: utf-8

import re
import glob
import pandas
import numpy as np
from random import choice

def phased2pandas(fin):
	fin = open(fin,'rU')
	data = []
	dt = None
	header = None

	# Process rows in file
	for count, row in enumerate(fin):
		row = row.strip()
		row = re.split("[\s:]", row)  # split on whitespace and colon

		# process header
		if count == 0:
			renamed_cols = []

			for col_count, col in enumerate(row[2:]):
				if col_count % 2 == 0: col += 'a'
				else: col += 'b'
				renamed_cols.append(col)

			header = ['I','chrm','start'] + renamed_cols
			
			# create appopriate datatypes
			dt = [('I','|S10'),('chrm','|S10'),('start','<i4')]
			dt += [(col,'|S10') for col in renamed_cols]
			dt = np.dtype(dt)

		# Update rows
		if count > 0:
			data.append(tuple(row))

	# create dataframe
	data = np.array(data,dtype=dt)
	result = pandas.DataFrame(data)

	# add header, make
	result.columns = header
	result = result.set_index(['chrm','start'])

	fin.close()

	return result


def panda2fasta(p, chrm_start, total_length=10000):

	pseudo_chrm, pseudo_start = chrm_start

	print p.ix[0,:].name[1]
	intial_snp_start =  p.ix[0,:].name[1]


	snp_positions = p.index.get_level_values('start')
	gap_sizes =  np.ediff1d(snp_positions)


	# subtract one to account for snp
	Ns = ['N' * (gap-1) for gap in gap_sizes] 
	
	# create N's for start and end of sequences

	initial_Ns = 'N' * (intial_snp_start - pseudo_start-1)
	fastas = ""
	
	for count, item in enumerate(p.values.T):
		if count == 0: continue # skip 'I'
		
		sequence = "" 
		for snp, N in zip(item,Ns):
			sequence += snp + N 

		final_Ns = 'N' * ((total_length - (len(initial_Ns) + len(sequence))))

		# create fasta info
		fasta_header = ">" + p.columns[count]

		

		final_seq =  initial_Ns + sequence + final_Ns

		# if len(final_seq) > 10000:
		print len(initial_Ns), len(sequence), len(final_Ns), sum([len(initial_Ns), len(sequence), len(final_Ns)])

		fastas += fasta_header + "\n" + final_seq + "\n"

	return fastas


def get_slice(scaffold,start,stop,conversion_panda):
	conversion_slice = conversion_panda[conversion_panda['Scaffold'] == scaffold]
	chrm = conversion_slice.Chrm.values[0]

	pseudo_start = conversion_slice.Start.values[0] + start
	pseudo_stop = conversion_slice.Start.values[0] + stop
	pseudo_chrm = chrm

	return (pseudo_chrm, pseudo_start, pseudo_stop)


def choose_scaffold(chrm):
	index = chrm.index
	
	possible_row = False
	row = None
	while (possible_row is False):
		idx_value = choice(index)
		row = chrm.ix[idx_value,:]
		possible_row = row['size'] > 10000

	return row


# Read in beagle data
cydno = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle/split_on_species/beagle_output/Cydno.non-multiallelic-variants-only.filtered.gatk.beagle-output.phased"
pachi = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle/split_on_species/beagle_output/Pachi.non-multiallelic-variants-only.filtered.gatk.beagle-output.phased"
melpo = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle/split_on_species/beagle_output/Melpo.non-multiallelic-variants-only.filtered.gatk.beagle-output.phased"
cydno = phased2pandas(cydno)
pachi = phased2pandas(pachi)
melpo = phased2pandas(melpo)


# Read in conversion
# Format: Chrm Scaffold Start Stop
psuedoChrm2Scaffold = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/BGI_Genome/Butterfly.merge.list"
psuedoChrm2Scaffold = pandas.read_csv(psuedoChrm2Scaffold, header=None)
psuedoChrm2Scaffold.columns = ['Chrm', 'Scaffold', 'Start', 'Stop']
 	
# Read in RAD Mapped Chromosomes.
RADmapped = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/Hmel1.1/Hmel_hox_RAD_matepair_updated.111215.chromosomes1_10kbSCAFFOLDSthatMAP.csv"
RADmapped = pandas.read_csv(RADmapped)

chrms = RADmapped.chromosome
chrms = list(set(chrms))

for count, RAD_chrm in enumerate(sorted(chrms)):
	chrm = RADmapped[RADmapped.chromosome == RAD_chrm]
	row = choose_scaffold(chrm)
	
	scaffold_start = choice(range(row['size']-10000))
	scaffold_stop = scaffold_start + 10000
	scaffold = row['scaffold']

	if scaffold.startswith("scf719") == True:
		continue

	conversion_slice = psuedoChrm2Scaffold[psuedoChrm2Scaffold['Scaffold'] == scaffold]
	if len(conversion_slice.Chrm.values) == 0:
		continue

	chrm = conversion_slice.Chrm.values[0]

	pseudo_start = conversion_slice.Start.values[0] + scaffold_start
	pseudo_stop = conversion_slice.Start.values[0] + scaffold_stop
	pseudo_chrm = chrm

	print RAD_chrm, pseudo_chrm, pseudo_start, pseudo_stop, pseudo_stop - pseudo_start
	
	cydno_slice = cydno.ix[(pseudo_chrm,pseudo_start):(pseudo_chrm,pseudo_stop)]
	pachi_slice = pachi.ix[(pseudo_chrm,pseudo_start):(pseudo_chrm,pseudo_stop)]
	melpo_slice = melpo.ix[(pseudo_chrm,pseudo_start):(pseudo_chrm,pseudo_stop)]

	beagle_merge = cydno_slice.join(pachi_slice.drop('I',axis=1)).join(melpo_slice.drop('I',axis=1))


	if np.any(pandas.isnull(beagle_merge).values == True) == True:
		continue

	name = '_'.join([str(item) for item in [RAD_chrm, scaffold, scaffold_start, scaffold_stop, pseudo_chrm, pseudo_start, pseudo_stop,]])


	print beagle_merge
	print beagle_merge.values

	fasta = panda2fasta(beagle_merge, chrm_start=[pseudo_chrm,pseudo_start], total_length=10000)

	fasta_name = name + ".fa"
	fout = open(fasta_name,'w')
	fout.write(fasta)

	pass

# Do example slice:
# chrm, start, stop = get_slice("Yb_superscaffold", 400, 10000, conversion_panda=psuedoChrm2Scaffold)
# test_slice =  beagle.ix[(chrm,start):(chrm,stop)]
# total_length = stop - start
# panda2fasta(test_slice, total_length=total_length)



