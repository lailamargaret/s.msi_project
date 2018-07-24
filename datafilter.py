import os 
import subprocess
import random

directory = '/home/upload/msi_project/tcga_bam/tumor_bams'

bam_files = []
for filename in os.listdir(directory):
	if filename.endswith('.bam'):
		bam_files.append(os.path.join(directory, filename))

chosen = random.sample(bam_files, 100)


for file in chosen:
	subprocess.call(['cp', file , '/home/upload/msi_project/tcga_bam/tumor_bams/subset'])
	subprocess.call(['cp', file + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/subset'])


