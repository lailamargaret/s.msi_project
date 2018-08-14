#related modules
import bamprocess
import count_reads
import methods
from constants import _MSI_LOCI, _QUALITY_THRESHOLDS, _ANNOTATIONS, _MSS_LOCUS_DATA, _ML_MODES

#imported modules
from ortools.linear_solver import pywraplp
import pandas as pd
import numpy as np
# ----------- Main --------------
#store bamfiles in a list
directory = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'
bamfiles = bamprocess.scan_files(directory)

def dataframe(directory, locus):
	setname = directory.split('/')[-1]
	
	bamfiles = bamprocess.scan_files(directory)

	bam_names = []
	lengths = []
	avg_lengths = []
	nums_lengths = []
	stdevs = []
	dists_mode = []
	msi_statuses = []

	for bam in bamfiles:
		bam_name = bam.split('/')[-1].replace('A.bam', '')
		#store msi status as a boolean, skip if not MSI or MSS
		if _ANNOTATIONS[bam_name] == 'MSI':
			msi_status = 1
		elif _ANNOTATIONS[bam_name] == 'MSS':
			msi_status = 0
		else:
			continue
	

		reads = count_reads.count(bam, locus)
		#if there are NO reads at this locus in the file, skip and go to the next bamfile
		if len(reads) == 0:
			continue
		bam_names.append(bam_name)
				
		#compute the lengths of the reads, store in a list, append to the overall lengths list
		bam_lengths = [len(e) for e in reads]
		lengths.append(bam_lengths)
	
		#compute the average length of a bam's reads, append to the overall avg_lengths list
		avg_length = np.mean(bam_lengths)
		avg_lengths.append(avg_length)
		
		#compute number of different lengths of a bam's reads, append to overall nums_lengths list
		num_lengths = len(set(bam_lengths))
		nums_lengths.append(num_lengths)
		
		#compute the standard deviation of a bam's reads, append to overall stdevs list
		stdev = np.std(bam_lengths)
		stdevs.append(stdev)
		
		#compute the average distance from the mode length 
		distance = 0
		count = 0
		for e in bam_lengths:
			distance += abs(e - float(_ML_MODES[locus]))
			count += 1
		dist_mode = distance / count
		dists_mode.append(dist_mode)
	
		#finally, add the MSI status
		msi_statuses.append(msi_status)
		


	data = {'bam_name' : bam_names, 
		'lengths' : lengths, 
		'average_length' : avg_lengths, 
		'num_lengths' : nums_lengths, 
		'stdev': stdevs, 
		'dist_mode' : dists_mode, 
		'msi_status' : msi_statuses}

	df = pd.DataFrame(data)

	print df

	df.to_csv('/home/upload/msi_project/ML/%s/%s_df.txt' % (locus, setname), sep='\t')

def get_set_data(setlist, locus):
	for dataset in setlist:
		dataframe(dataset, locus)

mss_msi_fullset = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mss_msi_fullset'
mode_train = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/mode_train'
training_set = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/training_set'
validation_set = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/validation_set'
test_set = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set'

setlist = [training_set, validation_set]

get_set_data(setlist, 'MSI-06')











# data frame for all bamfiles for one locus
# colunns:
# bam name 	list of all lengths	average length	number of different lengths	standard deviation	distance from mode	MSI status



