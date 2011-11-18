import os
import sys
import glob
import pysam
import argparse

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input-dir', required=True, 
                        help='The input directory containing the bam files.')
    parser.add_argument('-l','--sample-libray', required=True, 
                    help='The file containg the sample to run info.')
    args = parser.parse_args()
    return args
    

def make_lib_dicts(sample_lib_path):
    fin = open(sample_lib_path, 'rU')
    flow_cell_lane_dict = {}
    sample_name_dict = {}
    for line in fin:
        line = line.strip().split()
        date, time, flow_cell, lane, sample_code =  line[1].split("_")
        sample_name = line[0]
        flow_cell_lane = '%s:%s' % (flow_cell, lane[1:])
        flow_cell_lane_dict[flow_cell_lane] = [date, time, flow_cell, lane, sample_code, sample_name]
        
        
        if sample_name_dict.has_key(sample_name) == False:
            sample_name_dict[sample_name] = [[date, time, flow_cell, lane, sample_code]]
        elif sample_name_dict.has_key(sample_name) == True:
            sample_name_dict[sample_name].append([date, time, flow_cell, lane, sample_code])
        
    return (flow_cell_lane_dict, sample_name_dict)


def addRG2Header(filename, library_dicts):
    
    # CREATE TEMPLATE
    # Read group. Unordered multiple @RG lines are allowed.
    RG_template  ={ 'ID': '',          # Read group identifier. e.g., Illumina flowcell + lane name and number
            'CN': "BGI",       # GATK Not Required. Name of sequencing center producing the read.
            "DS": "",          # GATK Not Required. Description
            "DT": "",          # GATK Not Required. Date the run was produced (ISO8601 date or date/time)
            "PI": 500,         # GATK Not Required. Predicted median insert size.
            "PU": "",          # GATK Not Required. Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD).
            "SM": "",          # Sample. Use pool name where a pool is being sequenced.
            'PL': 'ILLUMINA'}  # Platform/technology used to produce the reads.
    
    samfile = pysam.Samfile(filename, "rb" )
    new_header = samfile.header.copy()
    BAM_path, filename = os.path.split(filename)
    sample_name = filename.split(".")[0]

    # ADD INFO TO TEMPLATE
    RGs =  []
    for sample in library_dicts[1][sample_name]:
        date, time, flow_cell, lane, sample_code = sample
        RG_template = RG_template.copy()
        RG_template['ID'] = flow_cell + ":" + lane[1:]
        RG_template['LB'] = sample_code
        RG_template['SM'] = sample_name
        RG_template['DS'] = sample_name + "-heli-butterfly"
        RG_template['DT'] = date # not sure if time is really time...
        RG_template["PU"] = flow_cell + "." + lane
        RGs.append(RG_template)
    
    new_header["RG"] = RGs
    samfile.close()
    
    return new_header


def makeRGBAM(filename, library_dicts):
    
    # Step 1: Make Modified Header
    new_RG_header = addRG2Header(filename,library_dicts)

    # Massage paths and make outputfiles
    path, filename = os.path.split(filename)
    name, ext = os.path.splitext(filename)
    new_name = name + '.wRG.' + "bam"
    outfile_name =  os.path.join(path,new_name)
    outfile = pysam.Samfile( outfile_name, "wb", header = new_RG_header )

    # Step 2: Process Samfile adding Read Group to Each Read
    samfile = pysam.Samfile(os.path.join(path, filename), "rb" )
    samfile.fetch()
    for count, read in enumerate(samfile.fetch()):
        name = read.qname
        read_group =  ":".join(name.split(":")[:2])
        new_tags = read.tags
        new_tags.append(("RG", read_group))
        read.tags = new_tags
        outfile.write(read)

    outfile.close()
    return 1

args = get_args()
library_dicts = make_lib_dicts('sample_lib.txt')
for filename in glob.glob(args.input_dir+'*'):
    if 'bam' in filename and 'with-RG' not in filename and 'bai' not in filename:
        print 'Processing:', filename
        makeRGBAM(filename, library_dicts)


