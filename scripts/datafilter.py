import os 
import subprocess
import random

directory = '/home/upload/msi_project/tcga_bam/tumor_bams'
UCEC_annotations = '/home/upload/msi_project/annotations/UCEC_annotations.txt'
COAD_READ_annotations = '/home/upload/msi_project/annotations/COAD_READ_annotations.txt'
cr_annotations = {}
u_annotations = {}
'''
ucec_count = 0
coadread_count = 0
with open (UCEC_annotations, 'r') as f:
	for line in f:
		fields = line.split('\t')
		u_annotations[fields[4]] = fields[22]
with open (COAD_READ_annotations, 'r') as f:
	for line in f:
		fields = line.split('\t')
		cr_annotations[fields[4]] = fields[20]
annotations = {}
for key in cr_annotations:
	annotations[key] = cr_annotations[key]
for key in u_annotations:
	annotations[key] = u_annotations[key]


coadread_bams = []
ucec_bams = []

#lists of coadread and ucec bams
for filename in os.listdir('/home/upload/msi_project/tcga_bam/COAD-READ'):
        if filename.endswith('.bam'):
                coadread_bams.append(filename)
for filename in os.listdir('/home/upload/msi_project/tcga_bam/UCEC'):
        if filename.endswith('.bam'):
                ucec_bams.append(filename)


#generate a set of all bam files that have MSS MSIL or MSIH labels for machine learning

coadread_count = 0
ucec_count = 0
total_count = 0

for filename in os.listdir(directory):
	if filename.endswith('.bam'):
		bam_name = filename.replace('A.bam', '')
		if bam_name in annotations:
			 if annotations[bam_name] == 'MSI-H' or annotations[bam_name] == 'MSI-L' or annotations[bam_name] == 'MSS':
				fileloc = '/home/upload/msi_project/tcga_bam/tumor_bams/%s' % filename
				subprocess.call(['cp', fileloc, '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset'])
                                subprocess.call(['cp', fileloc + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset'])
				if filename in coadread_bams:
					coadread_count += 1
				if filename in ucec_bams:
					ucec_count += 1
				total_count += 1			
print coadread_count
print ucec_count
print total_count

mss_msi_fullset = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset'
mode_train = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mode_train'
training_set = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/training_set'
validation_set = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/validation_set'
test_set = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set'

setlist = [mss_msi_fullset, mode_train, training_set, validation_set, test_set]

def count_mss_msi(directory):
	print directory
	msi_count = 0
	mss_count = 0

	for filename in os.listdir(directory):
		if filename.endswith('.bam'):
			bam_name = filename.replace('A.bam', '')
			if annotations[bam_name] == 'MSS':
				mss_count += 1
			if annotations[bam_name] == 'MSI-L' or annotations[bam_name] == 'MSI-H':
				msi_count += 1

	print 'msi: %d' % msi_count
	print 'mss: %d' % mss_count



def count_cancer_type(directory):
	print directory
	crc = 0
	uc = 0

	for filename in os.listdir(directory):
		if filename.endswith('.bam'):
			if filename in ucec_bams:
				uc += 1
			if filename in coadread_bams:
				crc += 1
	print 'coad-read: %d' % crc
	print 'ucec: %d' % uc


for bamset in setlist:
	count_mss_msi(bamset)
	count_cancer_type(bamset)

#list of all bam files in full set
fullset = []
for filename in os.listdir(mss_msi_fullset):
	if filename.endswith('.bam'):	
		fullset.append(filename)	


all_mss = []
for filename in fullset:
	bam_name = filename.replace('A.bam', '')
	if annotations[bam_name] == 'MSS':
		all_mss.append(filename)

#remove 100 mss files to train the mode
mode_train = random.sample(all_mss, 100)

for filename in mode_train:
	subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename, '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mode_train'])
	subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mode_train'])

remaining = [x for x in fullset if x not in mode_train]

#remove 300 files for the training set
training_set = random.sample(remaining, 300)

for filename in training_set:
	 subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename, '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/training_set'])
         subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/training_set'])

remaining = [x for x in remaining if x not in training_set]

#remove 100 files for validation set
validation_set = random.sample(remaining, 100)

for filename in validation_set:
         subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename, '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/validation_set'])
         subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/validation_set'])

remaining = [x for x in remaining if x not in validation_set]

#put the remaining 102 files in unused set
test_set = remaining
for filename in test_set:
	 subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename, '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set'])
         subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset/' + filename + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set'])

for filename in os.listdir('/home/upload/msi_project/tcga_bam/tumor_bams'):
	if filename.endswith('.bam'):
		bam_name = filename.replace('A.bam', '')
		if bam_name in annotations:
			if annotations[bam_name] == 'MSI-H' or annotations[bam_name] == 'MSI-L' or annotations[bam_name] == 'MSS':
				subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/' + filename, '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset'])
				subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/' + filename + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset'])

ml_fullset = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset'

crc = 0
ucecc = 0

print len(ucec_bams)
print len(coadread_bams)

for filename in os.listdir(ml_fullset):
	if filename in ucec_bams:
		ucecc += 1
	if filename in coadread_bams:
		crc += 1

print crc
print ucecc


tumor_bams = '/home/upload/msi_project/tcga_bam/tumor_bams'

for filename in os.listdir(tumor_bams):
	if filename.endswith('.bam'):
		bam_name = filename.split('/')[-1].replace('A.bam', '')
		if bam_name in annotations:
			subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/' + filename, '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])
			subprocess.call(['cp', '/home/upload/msi_project/tcga_bam/tumor_bams/' + filename + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])


bam_files = []
for filename in os.listdir(directory):
	if filename.endswith('.bam'):
		bam_name = filename.replace('A.bam', '')
		if annotations[bam_name] == 'MSS': 
			bam_files.append(os.path.join(directory, filename))

for bam in bam_files:
	bam_name = bam.split('/')[-1].replace('A.bam', '')
	if 
	subprocess.call(['cp', bam, '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])
	subprocess.call(['cp', bam + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated'])


chosen = random.sample(bam_files, 100)

for file in chosen:
	subprocess.call(['mv', file , '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/mss_training_set'])
	subprocess.call(['mv', file + '.bai', '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/mss_training_set'])

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

ucec_bams = []
coadread_bams =[]

ucec_count = 0
coadread_count = 0

for filename in os.listdir('/home/upload/msi_project/tcga_bam/COAD-READ'):
	coadread_bams.append(filename)
for filename in os.listdir('/home/upload/msi_project/tcga_bam/UCEC'):
        ucec_bams.append(filename)

num_ucec = len(ucec_bams)
num_coadread = len(coadread_bams)

print num_ucec
print num_coadread

for filename in os.listdir('/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'):
	if filename in coadread_bams:
		num_coadread -= 1
	if filename in ucec_bams:
		num_ucec -= 1		

for filename in os.listdir('/home/upload/msi_project/tcga_bam/tumor_bams/annotated/mss_training_set'):
        if filename in coadread_bams:
                num_coadread -= 1
        if filename in ucec_bams:
                num_ucec -= 1

print num_ucec
print num_coadread

'''

origdir = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set/'
valfile = '/home/upload/msi_project/ML/test_set_full_EDITED.txt'
editeddir = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set/edited'

bams = []
bais = []
with open (valfile, 'r') as f:
	for line in f:
		if line.startswith('bam') or line.startswith('#'):
			continue
		bam_name = line.split('\t')[0] + 'A.bam'
		bai_name = bam_name + '.bai' 

		bams.append(bam_name)
		bais.append(bai_name)

for bam in bams:
	subprocess.call(['cp', origdir + bam, editeddir])

for bam in bais:
	subprocess.call(['cp', origdir + bam, editeddir])














































