import gzip
from pylab import *
from pandas import *

fin = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/gatk_vs_soap/c511.Chr01.combined2.vcf"


def parseSNP(line, column_names, sample_ids):
    data = dict(zip(column_names,line))
    gatk, soap = sample_ids
    gatk = data[gatk].split(":")
    soap = data[soap].split(":")
    format = data['FORMAT'].split(":")
    gatk_dict = dict(zip(format,gatk))
    soap_dict = dict(zip(format,soap))
    if gatk_dict.has_key('GQ') == True and soap_dict.has_key('GQ') == True:
        return (int(data["POS"]), float(gatk_dict['GQ']), float(soap_dict['GQ']))
    else:
        return (int(data['POS']), np.nan, np.nan)

def main():    
    column_names = None
    sample_names = None 
    in_snps = False
    snp_counter = 0
    result = []
    for count, line in enumerate(open(fin)):
    
        if line.startswith("#CHROM"):
            column_names = line.strip().strip("#").split()
            sample_names = column_names[9:]
            in_snps = True
            continue
    
        if in_snps == True:
            snp_counter += 1
            quality = line
            line = line.strip().strip("#").split()
            GQs = parseSNP(line, column_names, sample_names)
            result.append(GQs)   

    finalGQs = array(result)  
    finalGQs = DataFrame(finalGQs[:,1:], columns=['SOAP', 'GATK'], index=finalGQs[:,0])   
    return finalGQs

def process_soap_file():
    soap_result = []
    soap_gzip = "/Volumes/PROMISE PEGASUS/Heliconius_Genome_Reseq_Data/MZH11020_Butterfly_re-sequencing/result_variation/snp/snp_result/c511.snp.gz"
    for count, line in enumerate(gzip.open(soap_gzip, 'rb')):
        line_parts = line.strip().split("\t")
        pvalue = float(line_parts[4])
        chrm, bp = line_parts[:2]
        bp = int(bp)
        soap_result.append((bp, pvalue))
        if bp >= 1333113: break

    soap_result = array(soap_result)
    soap_result = DataFrame(soap_result[:,1:], columns=['SOAP_Phred'], index=soap_result[:,0])
    return soap_result

vcf_GQs = main()
soap_GQs = process_soap_file()

z = DataFrame.join(soap_GQs, vcf_GQs["GATK"])

# GQs_store = HDFStore('chr1.heli.GQs.h5')
# GQs_store['GQs'] = finalGQs
