import gzip
import shlex
import pandas
import random
import subprocess


def create_panda_chrms_to_scaffolds(agp):
	agp_chrms = open(agp)
	data = []
	for count, line in enumerate(agp_chrms):
		line = line.strip().split("#")[0].strip()
		if len(line) > 0:
			line = line.split("\t")
			if len(line) == 8:
				new_line = line[:5] + line[6:]
				new_line.insert(6,line[5])
				new_line.insert(6,'1')
				line = new_line
			new_list = [int(a) if a.isdigit() else a for a in line]
			data.append(tuple(new_list))

	column_names = ['chrm','chrm_start','chrm_stop','chrm_id','chrm_type',\
					'scaffold','scaffold_start','scaffold_stop','scaffold_orientation']

	data = pandas.DataFrame(data, columns=column_names)
	return data

def randomly_select_chrm_chunks(data):
	unique_chrms = [item for item in set(data.chrm) if "_unmapped" not in item]
	unique_chrms.sort()

	random_scaffolds = []	
	for chrm in unique_chrms:
		chrm_set = data[data.chrm == chrm]
		chrm_set = chrm_set[chrm_set.chrm_type == 'D']

		choosing = 10
		while (choosing != 0):
			row_id = random.choice(chrm_set.index)
			row = chrm_set.ix[row_id]
			
			choices = []
			if row['scaffold_stop'] - row['scaffold_start'] > 10000 \
				and row['scaffold'] not in choices:
				
				random_scaffolds.append(row)
				choices.append(row['scaffold'])
				choosing -= 1

	return random_scaffolds

def get_regions(choices, region_length=10000):

	regions = []
	for count, choice in enumerate(choices):
		file_name = choice['chrm'], 
		region_start = None
		region_stop = None
		
		picked_region = False
		while (picked_region == False):
			region_start = random.randint(1,choice['scaffold_stop'])
			if region_start + region_length <= choice['scaffold_stop']:
				region_stop = region_start + region_length
				picked_region = True

		offset = choice['chrm_stop'] - choice['chrm_start']
		chrm_region_start = region_start + offset
		chrm_region_stop = region_stop + offset

		regions.append((choice['chrm'], chrm_region_start, chrm_region_stop, \
			choice['scaffold'], region_start, region_stop))

	return regions


def run_Beagle(region, fin, species):

	chrm, chrm_region_start, chrm_region_stop, scaffold, region_start, region_stop = region 
	
	gfin = gzip.open(fin , 'rb')
	header = gfin.next()
	gfin.close()

	filename_data = (chrm, chrm_region_start, chrm_region_stop, \
					 scaffold, region_start, region_stop, species) 
	file_name = '{6}_{0}:{1}-{2}_{3}:{4}-{5}.slice.beagle'.format(*filename_data)

	# Get Slice
	tabix_cmd = 'tabix {0} {1}:{2}-{3}'.format(fin, scaffold, region_start, region_stop)
	tabix_cmd = shlex.split(tabix_cmd)
	p = subprocess.Popen(tabix_cmd, stdout=subprocess.PIPE)
	text = p.communicate()[0]
	beagle_lines =  [item.split('\t') for item in  text.strip().split("\n") if item != '']

	# Get Header
	print 'Creating', file_name
	fout = open(file_name,'w')
	fout.write(header)
	for line in beagle_lines:
		loc = ":".join(line[:2]) + " " + " ".join(line[2:]) + "\n"
		fout.write(loc)
	fout.close()
	print 'running beagle...'
	beagle_out = file_name + '.out'

	beagle_cmd = "java -Xmx8g -jar /Users/MullenLab/Source/beagle-3.3.2/beagle.jar like={0} out={1}".format(file_name, beagle_out)
	beagle_cmd = shlex.split(beagle_cmd)
	p = subprocess.Popen(beagle_cmd, stdout=subprocess.PIPE)


# Cydno = "beagle/split_on_species/Cydno.non-multiallelic-variants-only.filtered.gatk.renamed.beagle-input.gz"
# Melpo = "beagle/split_on_species/Melpo.non-multiallelic-variants-only.filtered.gatk.renamed.beagle-input.gz"
# Pachi = "beagle/split_on_species/Pachi.non-multiallelic-variants-only.filtered.gatk.renamed.beagle-input.gz"

# data = create_panda_chrms_to_scaffolds()
# choices = randomly_select_chrm_chunks(data)
# regions = get_regions(choices)

# [run_Beagle(region, Cydno, 'Cydno') for region in regions]
# [run_Beagle(region, Melpo, 'Melpo') for region in regions]
# [run_Beagle(region, Pachi, 'Pachi') for region in regions]




