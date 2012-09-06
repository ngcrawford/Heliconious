

fin = '/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle/final.vcf'
fin = open(fin,'rU')

c_out = '/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle/Cydno.non-multiallelic-variants-only.4beagle.Chr06.gatk.vcf'
p_out = '/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle/Pachi.non-multiallelic-variants-only.4beagle.Chr06.gatk.vcf'
m_out = '/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/beagle/Melpo.non-multiallelic-variants-only.4beagle.Chr06.gatk.vcf'

c_out = open(c_out,'w')
p_out = open(p_out,'w')
m_out = open(m_out,'w')

in_genotypes = False

for count, line in enumerate(fin):

	line_parts = line.strip().split("\t")
	
	if in_genotypes == True:

		c_out.write('\t'.join(line_parts[:9] + line_parts[9:19]) + "\n")
		m_out.write('\t'.join(line_parts[:9] + line_parts[21:31]) + "\n")
		p_out.write('\t'.join(line_parts[:9] + line_parts[31:]) + "\n")


	if line.startswith("#CHROM"):
		c_out.write('\t'.join(line_parts[:9] + line_parts[9:19]) + "\n")
		m_out.write('\t'.join(line_parts[:9] + line_parts[21:31]) + "\n")
		p_out.write('\t'.join(line_parts[:9] + line_parts[31:]) + "\n")
		
		in_genotypes = True


	if in_genotypes == False:
		c_out.write(line)
		m_out.write(line)		
		p_out.write(line)


	if count % 100000 == 0:
		print count

fin.close()
c_out.close()
m_out.close()
p_out.close()