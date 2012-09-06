
def count_snps_in_BGI():
	fin = "/Volumes/PROMISE PEGASUS/Heliconius Genome Reseq Data/MZH11020_Butterfly_re-sequencing/result_variation/snp/snp_original_soapsnp/c511/c511.Chr01.cns"
	snps = 0.0
	total_bp = 0.0
	for count, line in enumerate(open(fin, 'rU')):
		line_parts = line.strip().split("\t")
		if line_parts[2] != line_parts[3]: 
			snps += 1
		total_bp += 1



	print 'snps', snps 
	print 'total_bp', total_bp
	print 'proportion', total_bp / snps