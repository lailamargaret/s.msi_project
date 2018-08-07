import os 
import subprocess
import random

directory = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'
UCEC_annotations = '/home/upload/msi_project/annotations/UCEC_annotations.txt'
COAD_READ_annotations = '/home/upload/msi_project/annotations/COAD_READ_annotations.txt'
annotations = {}


with open (UCEC_annotations, 'r') as f:
	for line in f:
		fields = line.split('\t')
		annotations[fields[4]] = fields[22]
with open (COAD_READ_annotations, 'r') as f:
	for line in f:
		fields = line.split('\t')
		annotations[fields[4]] = fields[22]
'''
tumor_bams = '/home/upload/msi_project/tcga_bam/tumor_bams'

for filename in os.listdir(tumor_bams):
	if filename.endswith('.bam'):
		bam_name = filename.split('/')[-1].replace('A.bam', '')
		if bam_name in annotations:
			subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/' + filename, '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])
			subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/' + filename + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])
'''

bam_files = []
for filename in os.listdir(directory):
	if filename.endswith('.bam'):
		bam_name = filename.replace('A.bam', '')
		if annotations[bam_name] == 'MSS': 
			bam_files.append(os.path.join(directory, filename))
'''
for bam in bam_files:
	bam_name = bam.split('/')[-1].replace('A.bam', '')
	if 
	subprocess.call(['cp', bam, '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])
	subprocess.call(['cp', bam + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])
'''

chosen = random.sample(bam_files, 100)

for file in chosen:
	subprocess.call(['mv', file , '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/mss_training_set'])
	subprocess.call(['mv', file + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/mss_training_set'])

'''
directory = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'
subset_bams = []
for filename in os.listdir(directory):
	subset_bams.append(filename)

print subset_bams
print len(subset_bams)
	
annotated = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'

for filename in os.listdir(annotated):
	if filename.endswith('.bam.bai'):
		if filename in subset_bams:
			filepath = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/%s' % filename
			subprocess.call(['rm', filepath])

'''



