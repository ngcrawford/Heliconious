#!/usr/bin/env python
# encoding: utf-8

import re
import os
import gzip
import glob
import pandas
import numpy as np
from random import choice

def phased2pandas(fin):
	fin = gzip.open(fin,'rb')
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


def panda2fasta(p, chrm_start, chrm_stop):

	total_length = chrm_stop - chrm_start

	# print p.is
	intial_snp_start =  p.ix[0].name[1]

	snp_positions = p.index.get_level_values('start')
	gap_sizes =  np.ediff1d(snp_positions)

	# subtract one to account for snp
	Ns = ['N' * (gap-1) for gap in gap_sizes] 
	
	# create N's for start and end of sequences

	initial_Ns = 'N' * (intial_snp_start - chrm_start-1)
	fastas = ""
	
	for count, item in enumerate(p.values.T):
		# if count == 0: continue # skip 'I'
		
		sequence = "" 
		for snp, N in zip(item,Ns):
			if snp == 'null': snp = 'N'
			sequence += snp + N 

		final_Ns = 'N' * ((total_length - (len(initial_Ns) + len(sequence))))

		# create fasta info
		fasta_header = ">" + p.columns[count]		

		final_seq =  initial_Ns + sequence + final_Ns

		fastas += fasta_header + "\n" + final_seq + "\n"

	return fastas

# Read in beagle data
cydno = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle_output/Cydno*.phased.gz"
pachi = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle_output/Melpo*.phased.gz"
melpo = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle_output/Pachi*.phased.gz"

cydno = glob.glob(cydno)
pachi = glob.glob(pachi)
melpo = glob.glob(melpo)

cydno.sort()
pachi.sort()
melpo.sort()

cpm_zip = zip(cydno,pachi,melpo)

for count, filenames in enumerate(cpm_zip):
	genus_set = [phased2pandas(g) for g in filenames]
	c,p,m = [g.drop('I',axis=1) for g in genus_set]
	cp = c.join(p)
	cpm = cp.join(m)
	chrm, start_stop = os.path.split(filenames[0])[1].split(".")[0].split("_")[2].split(":")
	start, stop = map(int, start_stop.split("-"))

	fasta = panda2fasta(cpm, start, stop)
	fout = os.path.split(filenames[0])[1].split('.')[0]
	fout = "_".join(fout.split('_')[1:]) + '.fa'
	fout = open(fout, 'w')
	fout.write(fasta)
	fout.close





