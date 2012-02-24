#!/usr/bin/env python
# encoding: utf-8
"""
vcf2oneliners.py

Created by Nick Crawford on 2011-11-18.
Copyright (c) 2011

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""

import os
import sys
import argparse
import itertools
import numpy as np
from collections import namedtuple
def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-start', type=int,
                        help="Starting position.")
    parser.add_argument('-stop', type=int,
                        help="Stop position.")
    parser.add_argument('-c','--chromosome', type=str,
                        help='Chromosome Id. First column in VCF.')
    parser.add_argument('-w','--window-size', required=True, type=int,
                        help='Window size')              
    parser.add_argument('-i','--input', required=True, 
                        help='Path to VCF file.')
    args = parser.parse_args()
    return args

def makeDataTuple(vcf):
    """Setup a labeled tuple to store the data."""
    
    chrm_data = {}
    
    for count, line in enumerate(vcf):
        
        if "##contig" in line:
            contigline = line.split("<")[1]
            contigline = contigline.strip(">\n").split(",")
            contigline = [item.split("=")[1] for item in contigline]
            chrm_data[contigline[0]] = int(contigline[1])
                
        if line.startswith("#CHROM"): 
            field_labels = line.strip().strip("#").split("\t")
            field_labels = [item.strip().split('.')[0] for item in field_labels]
            break
            
    position_data = namedtuple('base', field_labels)
    return (position_data, chrm_data)

def array2OnelinerAlignment(info, taxa, bases):
    """Convert array of array of taxa and an array of bases to one-liner."""

    oneliner = 'start=' + str(position) + "|"
    for count, seq in enumerate(bases):
        oneliner += taxa[count]+","+''.join(itertools.chain(bases[count])) + ","
    oneliner = oneliner[:-1] + ";"    
    return oneliner

def process_snp_call(snp_call, ref, alt):
    """Process VCF genotype fields.
        The current version is very basic and 
        doesn't directly take into account the
        quality of the call."""
        
    called_base = ""
    snp_call = snp_call.split(":")
        
    # process blanks
    if snp_call[0] == "./.":
        called_base = "-"
    
    else:
        
        allele1, allele2 = snp_call[0].split("/")
        
        # process "0/0"
        if allele1 == '0' and allele2 == '0':
            called_base = ref
        
        # process "0/N"
        if allele1 == '0' and allele2 != '0':
            called_base = ref
            
        # process "N/N"
        # this is a bit hacked. For example
        # '2/3' will be considered '2/2'
        if allele1 != '0' and allele2 != '0':
            pos = int(allele1) -1
            called_base = alt.split(",")[pos]
                    
    return called_base


def callSNPs(current_base, numb_of_seqs):
    """Call the SNPs. Duh!"""
        
    blanks =  np.zeros(numb_of_seqs, np.string0)
    
    if current_base.FILTER == 'LowQual':
        blanks.fill("-")
    
    if current_base.FORMAT == 'GT':
        blanks.fill("-")
    
    for count, snp_call in enumerate(current_base[9:]):
        base = process_snp_call(snp_call, current_base.REF, current_base.ALT)
        blanks[count] = base
        
    return blanks


def main():
    
    # SETUP ARGS
    args = get_args()
    align_range = (args.start, args.stop)
    window_size = args.window_size
    chrm = args.chromosome
    vcf_file_path = args.input
    
    # OPEN VCF
    vcf = open(vcf_file_path, 'rU')

    # SETUP NAMED TUPLE TO STORE INFO FROM A SINGLE BASE
    field_labels = []
    position_data, chrm_data = makeDataTuple(vcf)

    # SETUP MULTIPLE ALIGNMENT ARRAY
    numb_of_seqs = len(position_data._fields[9:])
    alignment = np.zeros((window_size,numb_of_seqs), np.string0)

    # SETUP COUNTERS
    current_base = None
    current_window = 1
    line_count = 0
    windows = range(0, chrm_data[chrm], window_size)
    current_data = []

    # PARSE VCF FIlE
    snp_count = 0    
    for line in vcf:

        # SKIP HEADER
        if line.startswith("#CHROM"): 
            line_count = 0   
             
        # START PROCESSING ALIGNED BASES
        if line.startswith(chrm):
            current_base = position_data._make(line.strip().split("\t"))     
            base_calls = callSNPs(current_base, numb_of_seqs) 
            
            if int(current_base.POS) >= windows[-1]: break
            
            # ADD DATA TO ALIGNMENT FOR CURRENT WINDOW        
            if int(current_base.POS) <= windows[current_window]:
                current_data.append(base_calls.copy())
                test = np.array(current_data).shape           
            
            # PRINTOUT RESULTING ONELINER
            if int(current_base.POS) > windows[current_window]:
                alignment =  np.array(current_data)
                taxa = current_base._fields[9:]
                info = 'chrm=%s,start=%s,stop=%s' % \
                (current_base.CHROM, windows[current_window-1], windows[current_window])                 
                oneliner = array2OnelinerAlignment(info, taxa, alignment.transpose())
                if ":" in oneliner: # this prevents bad alignments from getting printed
                    print oneliner
                current_data = []
                current_window += 1
            
        if line.startswith(chrm) == False: break
        line_count += 1

    # print 'total snps:', snp_count
    # print 'total bases:', line_count
    # print snp_count/line_count
    vcf.close()

if __name__ == '__main__':
    main()

