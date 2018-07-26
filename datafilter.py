import os 
import subprocess
import random

directory = '/home/upload/msi_project/tcga_bam/tumor_bams'
UCEC_annotations = '/home/upload/msi_project/UCEC_annotations.txt'
COAD_READ_annotations = '/home/upload/msi_project/COAD_READ_annotations.txt'
annotations = []

with open (UCEC_annotations, 'r') as f:
	for line in f:
		annotations.append(line.split('\t')[4])
with open (COAD_READ_annotations, 'r') as f:
	for line in f:
		annotations.append(line.split('\t')[4])

bam_files = []
for filename in os.listdir(directory):
	if filename.endswith('.bam'):
		bam_name = filename.split('/')[-1].replace('A.bam', '')
		if bam_name in annotations:
			bam_files.append(os.path.join(directory, filename))

for bam in bam_files:
	subprocess.call(['cp', bam, '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])
	subprocess.call(['cp', bam + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])

chosen = random.sample(bam_files, 100)

for file in chosen:
	subprocess.call(['cp', file , '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'])
	subprocess.call(['cp', file + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'])
