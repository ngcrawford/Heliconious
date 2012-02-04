#!/usr/bin/env python
# encoding: utf-8
"""
Some Description
"""

def count_snps_in_BGI(bases_to_process):
	fin = "/Volumes/PROMISE PEGASUS/Heliconius Genome Reseq Data/MZH11020_Butterfly_re-sequencing/result_variation/snp/snp_original_soapsnp/c511/c511.Chr01.cns"
	snps = 0.0
	total_bp = 0.0
	calls = []
	position = []
	for count, line in enumerate(open(fin, 'rU')):
		line_parts = line.strip().split("\t")
		calls.append(line_parts[3])
		position.append(int(line_parts[1]))
		if int(line_parts[1]) >= bases_to_process: break

	return (position, calls)


def count_snps_in_GATK(bases_to_process):

	ab_dict = {'GT':'K',
		   'AC':'M',
		   'AG':'R',
		   'CT':'Y',
		   'CG':'S',
		   'AT':'W',
		   'CGT':'B',
		   'ACG':'V',
		   'ACT':'H',
		   'AGT':'D'}

	fin = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/32butterflies.Chr01.allsites.gatk.vcf"
	in_header = True
	position = []
	calls = []
	for count, line in enumerate(open(fin, 'rU')):
		line = line.strip()
				
		if in_header == False:
			line_parts = line.split()
			pos, ref, alt, snp = int(line_parts[1]), line_parts[3], line_parts[4], line_parts[9]
			snp = snp.split(":")[0]
			final_snp = ''
			if snp == "0/0": final_snp = ref
			if snp == "./.": final_snp = ref
			if snp == '1/1': final_snp = alt
			
			if snp == '0/1': 
				final_snp = [ref,alt]
				final_snp.sort()
				final_snp = ''.join(final_snp)
				final_snp = ab_dict[final_snp]
			
			if snp == '1/0': 
				final_snp = [ref,alt]
				final_snp.sort()
				final_snp = ''.join(final_snp)
				final_snp = ab_dict[final_snp]
			
						 
			calls.append(final_snp)
			position.append(pos)
			if pos >= bases_to_process: break
		
		if line.startswith("#CHROM") == True:
			in_header = False
	
    return (position, calls)

bases_to_process = 100000		
gatk = count_snps_in_GATK(bases_to_process)
soap = count_snps_in_BGI(bases_to_process)

gatk_zip = zip(gatk[0], gatk[1])
soap_zip = zip(soap[0], soap[1])

print gatk_zip[-15:-1]
print soap_zip[-15:-1]


missmatches = 0.0
total_bases = 0.0
for count, (g_call, s_call) in enumerate(zip(gatk[1],soap[1])):
	g_call, s_call = (g_call, s_call)
	if g_call !=  s_call:
		print 
		missmatches += 1.0
	total_bases += 1.0

print 'missmatches', missmatches
print 'total_bases', total_bases
print missmatches/total_bases
		
		
		
		 