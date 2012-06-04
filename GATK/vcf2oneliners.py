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
import shlex
import random
import argparse
import itertools
import numpy as np
from subprocess import Popen, PIPE
from collections import namedtuple

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', type=int,
                        help="Starting position.")
    parser.add_argument('--stop', type=int,
                        help="Stop position.")
    parser.add_argument('-c','--chromosome', type=str,
                        help='Chromosome Id. First column in VCF.')
    parser.add_argument('-w','--window-size', required=True, type=int,
                        help='Window size')
    parser.add_argument('-i','--input', required=True,
                        help='Path to VCF file.')
    parser.add_argument('-b','--bootstraps', type=int,
                        help='Calculate bootstraps.')
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
<<<<<<< HEAD
    oneliner = str(info) + ":"
=======

    oneliner = 'start=' + str(position) + "|"
>>>>>>> dev
    for count, seq in enumerate(bases):
        oneliner += taxa[count]+","+''.join(itertools.chain(bases[count])) + ","
    oneliner = oneliner[:-1] + ";"
    return oneliner

<<<<<<< HEAD

=======
>>>>>>> dev
def process_snp_call(snp_call, ref, alt):
    """Process VCF genotype fields.
        The current version is very basic and
        doesn't directly take into account the
        quality of the call."""
<<<<<<< HEAD

=======
        
>>>>>>> dev
    called_base = ""
    snp_call = snp_call.split(":")

    # process blanks, many are reall
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

<<<<<<< HEAD
def count_informative_sites(alignment_array):
    """Informative Sites must have two different SNPs"""
    informative_sites = 0
    for site in alignment_array:
        unique_sites = set(site)
        if len(unique_sites) >= 3:
            informative_sites += 1
    return informative_sites

def get_subset_vcf(chrm, start, stop):
    base_dir = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice"
    cli = """java -Xmx6g -jar /Users/MullenLab/Source/gatk/dist/GenomeAnalysisTK.jar \
      -R {3}/Butterfly.merge.scaffolds.fa  \
      -T SelectVariants \
      --variant {3}/32_butterflies_vcfs/32.butterflies.good_BGI_snps.combined.vcf \
      -L {0}:{1}-{2}""".format(chrm, start, stop, base_dir)

    cli_parts = shlex.split(cli)
    vcf = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()[0]
    return vcf
=======

def main():
    
    # SETUP ARGS
    args = get_args()
    align_range = (args.start, args.stop)
    window_size = args.window_size
    chrm = args.chromosome
    vcf_file_path = args.input
    
    # OPEN VCF
    vcf = open(vcf_file_path, 'rU')
>>>>>>> dev

def generate_bootstraps(chrm, chrm_len, window_size, numb_of_reps):

    start_sites = [random.randint(0, chrm_len-window_size) for item in range(numb_of_reps)]
    replicates = []
    for count, start in enumerate(start_sites):
        vcf = get_subset_vcf(chrm, start, start+window_size)
        yield vcf

def parse_window_vcf(vcf, start, stop, fout):
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
    informative_sites = []

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
            current_data.append(base_calls.copy())

    alignment =  np.array(current_data)
    inform_sites = count_informative_sites(alignment)
    if current_base == None:
        return 'error'
    else:
        taxa = current_base._fields[9:]
        info = 'chrm={0},start={1},stop={2},inform_sites={3}'.format(current_base.CHROM, start, stop, inform_sites)
        oneliner = array2OnelinerAlignment(info, taxa, alignment.transpose())

    if ":" in oneliner and oneliner[-1] == ';': # this prevents bad alignments from getting printed
        return oneliner
    else:
        return 'error'

def parse_vcf(vcf, window_size):
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
    informative_sites = []

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

            # APPLY FILTERS
            if int(current_base.POS) >= windows[-1]: break
            if int(current_base.POS) < args.start and args.start != None: continue
            if int(current_base.POS) >= args.stop and args.stop != None: break

            # ADD DATA TO ALIGNMENT FOR CURRENT WINDOW
            if int(current_base.POS) <= windows[current_window]:
                current_data.append(base_calls.copy())
                test = np.array(current_data).shape

            # PRINTOUT RESULTING ONELINER
            if int(current_base.POS) > windows[current_window]:
                alignment =  np.array(current_data)
                informative_sites.append(pruneUninformativeSites(alignment))
                taxa = current_base._fields[9:]
                info = 'chrm=%s,start=%s,stop=%s' % \
                    (current_base.CHROM, windows[current_window-1], windows[current_window])
                oneliner = array2OnelinerAlignment(info, taxa, alignment.transpose())
                # if ":" in oneliner: # this prevents bad alignments from getting printed
                #      print oneliner
                current_data = []
                current_window += 1

        if line.startswith(chrm) == False: break
        line_count += 1


def main():

    # SETUP ARGS
    args = get_args()
    align_range = (args.start, args.stop)
    window_size = args.window_size
    chrm = args.chromosome
    vcf_file_path = args.input

    # OPEN VCIF
    vcf = open(vcf_file_path, 'rU')
    parse_vcf(vcf)
    vcf.close()
    return informative_sites

if __name__ == '__main__':

    args = get_args()

    # process boot
    if args.bootstraps:
        fout = open('32.butterflies.1000bootstraps.1000bp.oneliners', 'w')
        chrm = args.chromosome
        chrm_len = 1333113
        window_size = args.window_size
        numb_of_reps = args.bootstraps
        for count, vcf in enumerate(generate_bootstraps(chrm, chrm_len, window_size, numb_of_reps)):
            vcf = vcf.split("\n")
            if vcf != 'error':
                print 'made bootrep: ', count
                oneliner = parse_window_vcf(vcf, 0, window_size, fout)
                fout.write(oneliner+"\n")
            else:
                print 'error', vcf
    # informative_sites = main()
    # print informative_sites

