#!/usr/bin/env python
# encoding: utf-8
import re
import pandas
import argparse
import tempfile
import subprocess
import numpy as np
from random import choice

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('--input-vcf', help='Path to input VCF file.')
    p.add_argument('--output-vcf', help='Path to output VCF file.')
    p.add_argument('--input-beagle', help='Path to input beagle file.')
    p.add_argument('--output-beagle', help='Path to input beagle file.')

    args = p.parse_args()
    return args


class LiftOver(object):
	"""docstring for LiftOver"""
	def __init__(self, path_to_pseudo_chrms_coords_to_scaffolds, path_old_scaffolds_to_NCBI_accessions):
		super(LiftOver, self).__init__()

		self.path_to_pseudo_chrms_coords_to_scaffolds = path_to_pseudo_chrms_coords_to_scaffolds
		self.path_old_scaffolds_to_NCBI_accessions = path_old_scaffolds_to_NCBI_accessions
		self.agp_scaffolds_to_final_chrms_panda = self.__agp2Panda__()
		self.pseudo_chrms_coords_to_scaffolds_panda = self.__merge2Panda__()
		
	def __agp2Panda__(self):
		agp_file = open(self.path_old_scaffolds_to_NCBI_accessions,'rU')

		column_names = ['new_scaffold', 'new_start', 'new_end',  "new_part_number", 'type', \
				    'old_scaffold','old_start', 'old_end', 'orientation', \
				    'org_scaffold','org_start','org_end','org_order']

		data = []
		for count, line in enumerate(agp_file):

			line = line.strip()
			items = re.split("[^a-zA-Z0-9_+\-]",line)
			items = [item for item in items if item != ""] # remove empty characters.

			if len(items) < 10: continue
			
			if 'yes' in items: 
				items = items[:5] + ["100 fragment",1,100,'-'] + items[-4:]

			data.append(tuple(items))

		dt = [('new_scaffold','|S16'), ('new_start','<i4'), ('new_stop','<i4'), ('new_part_number','<i4'), ('type','|S50'), \
			  ('old_scaffold','|S50'), ('old_start','<i4'), ('old_stop','<i4'), ('orientation','|S50'), \
			  ('org_scaffold','|S50'), ('org_start','<i4'), ('org_stop','<i4'), ('org_order', '<i4')]
		

		dt = np.dtype(dt)
		data = np.array(data,dtype=dt)
		result = pandas.DataFrame(data,)
		agp_file.close()

		return result

	def __merge2Panda__(self):
		fin = open(self.path_to_pseudo_chrms_coords_to_scaffolds,'rU')

		column_names = ['pseudo_chrm', 'old_scaffold', 'pseudo_start', 'pseudo_stop']
		data = []
		for count, line in enumerate(fin):
			# if count > 1000: break
			line = line.strip()
			items = re.split("[^a-zA-Z0-9_+\-]",line)
			items = [item for item in items if item != ""] # remove empty characters.
			data.append(tuple(items))

		dt = [('pseudo_chrm','|S50'), ('old_scaffold','|S50'), ('pseudo_start','<i4'), ('pseudo_stop','<i4')]
		dt = np.dtype(dt)
		
		data = np.array(data,dtype=dt)
		result = pandas.DataFrame(data)
		#result = result.set_index(['pseudo_chrm','pseudo_start'])
		fin.close()

		return result

	def old_scaffolds_to_NCBI_Accessions(self, old_scaffold, pos):
		ostna = self.agp_scaffolds_to_final_chrms_panda

		row = ostna[(ostna['org_scaffold'] == old_scaffold) \
			  & (ostna['org_start'] <= int(pos)) \
			  & (int(pos) <= ostna['org_stop'])]

		try:								  # deal with possible offset
			scaffold_pos = ((int(pos) + 1) + (row["new_start"].values[0] - row["org_start"].values[0]))

		except:
			print row, pos, old_scaffold

		return (row['new_scaffold'].values[0], scaffold_pos)


	def pseudo_pos_to_old_scaffold_pos(self, pseudo_chrm, pos):

		pcctsp = self.pseudo_chrms_coords_to_scaffolds_panda

		row = pcctsp[(pcctsp['pseudo_chrm'] == pseudo_chrm) \
					  & (pcctsp['pseudo_start'] <= int(pos)) \
					  & (int(pos) <= pcctsp['pseudo_stop'])]

		scaffold_pos = (int(pos) - row["pseudo_start"].values[0])

		return (row['old_scaffold'].values[0], scaffold_pos)

	def pseudo_to_NCBI(self, chrm,pos):
		scaffold, scaff_pos = self.pseudo_pos_to_old_scaffold_pos(chrm,pos)
		new_scaffold, new_scaff_pos = self.old_scaffolds_to_NCBI_Accessions(scaffold, scaff_pos)
		return (new_scaffold, new_scaff_pos)

	def updateVCF(self, vcf, vcf_out):
		vcf = open(vcf,'rU')
		out_vcf = open(vcf_out,'w')

		in_snps = False
		count = 0
		for line in vcf:
			if in_snps == False:
				out_vcf.write(line)
			
			if in_snps == True:
				# if count > 1000: break
				chrm, pos, Id, ref, alt, qual, filtr, info, format = line.split()[:9]
				
				try:
					scaffold, scaff_pos = self.pseudo_to_NCBI(chrm,pos)
					genotypes = line.split()[9:]
					line =  '\t'.join([scaffold, str(scaff_pos), Id, ref, \
						               alt, qual, filtr, info, format] + genotypes) + "\n"
					out_vcf.write(line)

				except:
					sterr_out = open('missing_scaffolds.vcf.csv','a')
					out_line = ','.join([scaffold, str(scaff_pos), chrm, str(pos), Id, ref, alt, qual, filtr, info, format])
					sterr_out.write(out_line)
				
				if count == 1000000:
					print count, chrm, pos, scaffold, scaff_pos

				#if count == 10: break

				count += 1

			if line.startswith('#CHROM'): 
				in_snps = True

		vcf.close()
		out_vcf.close()

	def runBeagle(self,input_beagle, output_beagle):

		fin = open(input_beagle,'rU')
		fout = open(output_beagle, 'w')

		scaffold_count = 0
		current_scaffold = None
		header = None
		for count, line in enumerate(fin):
			
			if count == 0: 
				header = line
				fout.write(header)
				continue

			chrm, pos = line.strip().split(" ")[0].split(":")
			pos = int(pos)
			
			# if count > 226 and count < 233:
			try:
				NBCI_chrm, NCBI_pos = self.pseudo_to_NCBI(chrm,pos)
				new_line = [NBCI_chrm, str(NCBI_pos)] + line.strip().split(" ")[1:]
				fout.write("\t".join(new_line)+"\n")
			
			except:
				sterr_out = open('missing_scaffolds.beagle.csv','a')
				outline = ','.join([NBCI_chrm, str(NCBI_pos), chrm, str(pos)]) +"\n"
				sterr_out.write(outline)

			if count == 10000:
				print count, NBCI_chrm, NCBI_pos, chrm, pos
			# if count > 10: break

		pass




args = get_args()

LO = LiftOver(path_to_pseudo_chrms_coords_to_scaffolds='/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/BGI_Genome/Butterfly.merge.list', \
		 	  path_old_scaffolds_to_NCBI_accessions="/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/Hmel1.1/Hmel1-1_primaryScaffolds.agp")

#vcf = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/multiallelic_matt/32-butterflies.multiallelic.all-chrms.gatk.vcf"
#new_vcf = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/multiallelic_matt/32-test.multiallelic.all-chrms.filtered.gatk.vcf"

# beagle = "beagle/split_on_species/Cydno.non-multiallelic-variants-only.filtered.gatk.beagle-input"
# LO.runBeagle(beagle)

if args.input_beagle !=  None:
	LO.runBeagle(args.input_beagle, args.output_beagle)

if args.input_vcf != None:
	LO.updateVCF(args.input_vcf, args.output_vcf)

