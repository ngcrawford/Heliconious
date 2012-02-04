#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import gzip
import argparse

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, nargs='*', 
                        help='The input BGI SNP file.')
    args = parser.parse_args()
    return args


def makeAltField(fields, ab_dict):
	ref = fields[2]
	GT = None
	alt = ''
	
	if fields[3] in ['K','M','R','Y','S','W']:
		GT = ""
		unsorted_alt = ab_dict[fields[3]]
		
		# new het not matching ref
		if ref not in unsorted_alt:
			GT = "0/1"
			alt = fields[5] 

			# picking the best one and ignoring hets is not really the way to
			# do this, but GATK CombineVariants barfs when I give it multiple 
			# alternate alleles
	
		# het matches ref
		if ref in unsorted_alt:
			GT = "0/1"
			for item in unsorted_alt:
				if item != ref:
					alt = item
					break
	else:
		GT = "1/1"
		alt = fields[3]
		
	return (GT,alt) 
    
def processSNPFile(file_path):    
    
    path, filename = os.path.split(file_path)
    print filename
    filename_no_ext, ext = os.path.splitext(filename)
    
    # handle gzip files
    if ext == '.gz':
       fin = gzip.open(file_path, 'rb') 
    else:
        fin = open(file_path, 'rU')
        
    #fout = open(os.path.split(file_path))
    fout = os.path.join(path,filename_no_ext+".vcf")
    fout = open(fout,'w')


    header =  """##fileformat=VCFv4.0
##fileDate=20090805
##source=SOAPsnp
##reference=~/Desktop/Grad_Students/Nick/butterfly_practice/Butterfly.merge.fa
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50 percent of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	 %s\n""" % (filename_no_ext)
	
    fout.write(header)
    ab_dict =  {'A':'A',
    			'C':'C',
    			'G':'G',
    			'T':'T',
    			'K':['G','T'],
    			'M':['A','C'],
    			'R':['A','G'],
    			'Y':['C','T'],
    			'S':['C','G'],
    			'W':['A','T'],
    			'B':['C','G','T'],
    			'V':['A','C','G'],
    			'H':['A','C','T'],
    			'D':['A','G','T']}

    # Format description (from left to right):
    #         1)      Chromosome
    #         2)      Coordinate
    #         3)      Reference nucleotide
    #         4)      Consensus nucleotide (HETs are reported in abbreviations)
    #         5)      Consensus quality (in PHRED scale)
    #         6)      Best nucleotide covering the position
    #         7)      Average quality score of best nucleotide
    #         8)      Number of best nucleotide (in unique hits only)
    #         9)      Number of best nucleotide (in all hits)
    #         10)     Second best nucleotide covering the position
    #         11)     Average quality score of second best nucleotide
    #         12)     Number of second best nucleotide (in unique hits only)
    #         13)     Number of second best nucleotide (in all hits)
    #         14)     Total depth
    #         15)     Rank sum test statistic
    #         16)     Estimate copy number of the site
    #         17)     dbSNP tag (in imputation model the tag will be 1 if a site is found in dbSNP)

    snp_count = 0
    for count, line in enumerate(fin):
    	fields = line.split("\t")
    	if fields[2] != fields[3]:
    		chrm = fields[0]
    		pos = fields[1]
    		snp_id = '.'
    		ref = fields[2]
    		GT, alt = makeAltField(fields, ab_dict)
    		qual = fields[4]
    		filt = '.'
    		info = '.'
    		format = 'GT:AD:DP:GQ'
    		
    		# AD values (one for each of REF and ALT fields) is the 
    		# count of all reads that carried with them the REF and ALT alleles
    		AD = ",".join([fields[8], fields[12]]) 
    		
    		#  DP read depth at this position for this sample (Integer)
    		DP = str(sum([int(fields[7]), int(fields[11])]))
    		
    		# GQ genotype quality, encoded as a phred quality -10log_10p 
    		#   (genotype call is wrong) (Numeric)
    		GQ = fields[6] 						  

    		call_string = ":".join([GT,AD,DP,GQ]) 
    		if ref == alt: continue
    		out = [chrm, pos, snp_id, ref, alt, qual, filt, info, format, call_string]
    		final_call = '\t'.join(out)+"\n"
    		fout.write(final_call)
	
    	snp_count += 1

    	if snp_count % 100000 == 0:
    		print snp_count

    fout.close()
    fin.close()


args = get_args()
if len(args.input) > 1:
    for fin in args.input:
        print fin
        processSNPFile(fin)
else:
    processSNPFile(args.input[0])




