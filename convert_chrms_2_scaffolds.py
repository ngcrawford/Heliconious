import os
import sys
import glob
import shlex
import pandas
import subprocess
import numpy as np

def agp_2_panda(fin):
	"""Convert agp file into panda"""

	agp_array = []
	for count, line in enumerate(open(fin)):
		if line.startswith("#"): continue

		line = line.strip().split("#")[0].strip()
		line_parts = line.split('\t')

		# Handle truncated 'Fragment' lines
		if 'fragment' in line_parts: 
			line_parts = line_parts[:5] + ['fragment','1'] + [str(int(line_parts[5]) + 1)] + ['NA']

		agp_array.append(tuple(line_parts))

	# Create special numpy datatype
	dt = [('chrm','|S16'), ('c_start','<i4'), ('c_stop','<i4'),\
		  ('id','<i4'), ('type','|S50'), \
		  ('scaffold','|S50'), ('s_start','<i4'), ('s_stop','<i4'),\
		  ('orientation','|S50'),]
			
	agp_array = np.array(agp_array,dtype=np.dtype(dt))
	return pandas.DataFrame(agp_array)

def slice_chrm(agp, chrm, start, stop):
	"""Incrementally slice AGP panda to get the slice that contains just
	   the region of interest."""

	start = int(start)
	stop = int(stop)

	chrm_slice = agp[agp['chrm'] == chrm]
	start_slice = chrm_slice[chrm_slice['c_stop'] >= start]
	return start_slice[start_slice['c_start'] <= stop]

def slice_2_scaffold_list(s, chrm, start, stop):
	"""Converts slice into 2D list of scaffold IDs, starts, and stops. 
	   Does appropriate math to update positions the first an last 
	   slices."""
	
	start = int(start)
	stop = int(stop)

	scaffolds = []

	# Process single row edge case
	if s.shape[0] == 1:
		first_row = s.values[0]
		scaffold_start =  (start - first_row[1]) + first_row[6]
		scaffold_stop = (stop - first_row[1]) + first_row[6]
		scaffolds.append([first_row[5], scaffold_start, scaffold_stop])

	# Process multi-row slice
	else:	
		first_row = s.values[0]
		scaffold_start =  (start - first_row[1]) + first_row[6]
		scaffolds.append([first_row[5], scaffold_start, first_row[7]])
		
		for row in s.values[1:-1]:
			scaffolds.append(list(row[5:8]))

		last_row = s.values[-1]
		scaffold_stop =  (stop - last_row[1]) + last_row[6]
		scaffolds.append([last_row[5], last_row[6], scaffold_stop])

	return scaffolds

def get_gff3_slices(s_list, chrm_info):

	bgz = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/Hmel1.1/heliconius_melpomene_v1.1_primaryScaffs_wGeneSymDesc.gff3.sorted.bgz"
	mod_gff3 = []

	for item in s_list:
		values = item + [bgz]
		cli = "tabix {3} {0}:{1}-{2}".format(*values)
		cli = shlex.split(cli)
		result = subprocess.Popen(cli, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
		for row in result.strip().split("\n"):
			row_parts = row.split("\t")
			updated_row = chrm_info + row_parts
			mod_gff3.append(updated_row)

	return mod_gff3


def islands_2_gff3(agp, fin, fout):

	fout = open(fout,'w')
	for line in open(fin,'rU'):
		line_parts = line.strip().split()
		s = slice_chrm(agp, *line_parts)
		s_list = slice_2_scaffold_list(s, *line_parts)
		s_gff3 = get_gff3_slices(s_list, line_parts)
		for line in s_gff3:
			fout.write('\t'.join(line)+"\n")

	fout.close()

def islands_2_gene_ids(agp, fin, fout):

	fout = open(fout,'w')
	for line in open(fin,'rU'):
		line_parts = line.strip().split()
		s = slice_chrm(agp, *line_parts)
		s_list = slice_2_scaffold_list(s, *line_parts)
		s_gff3 = get_gff3_slices(s_list, line_parts)
		for count, line in enumerate(s_gff3):
			if count == 0:
				fout.write("probeset\tSystemCode\n")
			if 'gene' in line:
				info = dict([ pair.split("=") for pair in  line[-1].strip(';').split(";")])
				if "ID" in info:
					final = "{0}\t{1}\n".format(info["ID"],'HM')
					fout.write(final)
	pass

def islands_2_gene_fastas(agp, fin,):


	for line in open(fin,'rU'):
		line_parts = line.strip().split()
		s = slice_chrm(agp, *line_parts)
		s_list = slice_2_scaffold_list(s, *line_parts)
		s_gff3 = get_gff3_slices(s_list, line_parts)
		for count, line in enumerate(s_gff3):
			if 'gene' in line:
				info = dict([ pair.split("=") for pair in  line[-1].strip(';').split(";")])
				d = [os.path.split(fin)[-1].strip(".txt")] + [info["ID"]] + line_parts + [line[3]] + line[6:8]
				print line
				print " ".join(d)

	 			cli = "samtools faidx /Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/Hmel1.1/heliconius_melpomene_v1.1_primaryScaffs_CDS.fna {}".format(info["ID"]+'-RA')
				cli = shlex.split(cli)
				result = subprocess.Popen(cli, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]


agp_in = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/Hmel1.1/Hmel1-1_hox_RAD_matepair_chromosomes.agp"
agp = agp_2_panda(agp_in)

files = glob.glob('/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/Islands/C,P,M/FST_C_P.txt')
for fin in files:
	print 'processing', fin
	#fout = os.path.splitext(fin)[0] + ".gff3"
	islands_2_gene_fastas(agp, fin)
	#mod_gff3 = islands_2_gff3(agp, fin, fout)

