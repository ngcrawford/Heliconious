
fin = open("/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/VCF/heliconius.vcf",'rU')
fout = open("/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/VCF/heliconius.chr-scf-ids.vcf",'w')
line_template = None

for count, line in enumerate(fin):
        
    if line.startswith('#') != True:
        line_parts = line.strip().split()
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = line_parts[0:9]
        INFO = dict([item.split("=") for item in INFO.split(";")])
        
        if CHROM.startswith("sch") != True:
            new_CHROM = ".".join((CHROM, INFO["SC"]))
        else:
            new_CHROM = CHROM
            
        new_line = [new_CHROM] + line_parts[1:]
        final_line = '\t'.join(new_line) +"\n"
        fout.write(final_line)
    
    else:
        fout.write(line)
    
    if count % 100000 == 0: print count

fout.close()