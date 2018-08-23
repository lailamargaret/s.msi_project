from constants import _MSI_LOCI, _QUALITY_THRESHOLDS, _ANNOTATIONS, _MSS_LOCUS_DATA, _ML_MODES
import bamprocess
import methods
import count_reads

import pandas as pd
import numpy as np
import subprocess
import os

#bamname	locus1avglen	locus1numlen	locus1stdev	locus1distmode	locus2..	status


def metrics(directory, loci_list):
	setname = directory.split('/')[-1]
	
	bamfiles = bamprocess.scan_files(directory)
	outfile = '/home/upload/msi_project/ML/%s_full.txt' % setname
	with open(outfile, 'w') as f:
		#f.write('#Full %s dataset\nbam_name\t' % setname)
		for locus in loci_list:
			f.write('%s_avg_len\t%s_num_lens\t%s_stdev\t%s_dist_mode\t' % (locus, locus, locus, locus))
		f.write('msi_status\n')
		for bam in bamfiles:
			bam_name = bam.split('/')[-1].replace('A.bam', '')
			if _ANNOTATIONS[bam_name] == 'MSI':
				msi_status = 1
			elif _ANNOTATIONS[bam_name] == 'MSS':
				msi_status = 0
			else:
				continue
			f.write(bam_name + '\t')
			for locus in loci_list:
				reads = count_reads.count(bam, locus)
			
				if len(reads) == 0:
					avg_length = num_lengths = stdev = dist_mode = 'NaN'
				else:				
					lengths = [len(e) for e in reads]
					avg_length = np.mean(lengths)
					num_lengths = len(set(lengths))
					stdev = np.std(lengths)
					distance = count = 0
					for e in lengths:
						distance += abs(e - float(_ML_MODES[locus]))
						count += 1
					dist_mode = distance / count
				f.write(str(avg_length) + '\t' + str(num_lengths) + '\t' + str(stdev) + '\t' + str(dist_mode) + '\t')
			f.write(str(msi_status) + '\n')
				
			
			
def print_series_names():
	for locus in _MSI_LOCI:
		print '\'%s_avg_len\',' % locus
		print '\'%s_num_lens\',' % locus
		print '\'%s_stdev\',' % locus
		print '\'%s_dist_mode\',' % locus
		
def write_lengths_df(directory, loci_list, num_bins):
	setname = directory.split('/')[-1]
	
	bamfiles = bamprocess.scan_files(directory)
	outfile = '/home/upload/msi_project/ML/%s_top_9_lengths_full.txt' % setname
	with open(outfile, 'w') as f:
		#write the header row
		f.write('bam name\t')
		feature_list = []
		for locus in loci_list:
			for i in range(num_bins):
				idx = i + 1
				feature = '%s_%d' % (locus, idx)
				f.write(feature + '\t')
				feature_list.append(feature)	
		f.write('msi_status\n')
		for bam in bamfiles:
			bam_details = []
			bam_name = bam.split('/')[-1].replace('A.bam', '')
			if _ANNOTATIONS[bam_name] == 'MSI':
				msi_status = 1
			elif _ANNOTATIONS[bam_name] == 'MSS':
				msi_status = 0
			else:	
				continue
			features = get_locus_counts(bam, loci_list, num_bins)
			if features == 'N/A':
				continue
			f.write(bam_name + '\t')
			for feature in features:
				f.write(str(feature) + '\t')
			f.write(str(msi_status) + '\n')				
	
						
def get_locus_counts(bamfile, loci_list, num_bins):
	features = []
	for locus in loci_list:
		reads = count_reads.count(bamfile, locus)
		#if any of the loci have too low read count, return N/A
		if len(reads) < 5:
			return 'N/A'
		lengths = [len(e) for e in reads]
		np_reads = np.asarray(lengths)
		counts = np.bincount(np_reads)
		x = [0]
		while counts.size < (num_bins + 1):
			counts = np.append(counts, x)	
		idx = [0]
		counts = np.delete(counts, idx)
		norm_counts = []
		for element in counts:
			norm_counts.append(float(element) / len(reads))
		features.extend(norm_counts)
	return features

def write_features(loci_list, num_bins):
	with open('/home/upload/msi_project/ML/top_9_histogram_features.txt', 'w') as f:
		for locus in loci_list:
			for i in range(num_bins):
				idx = i + 1
				f.write('%s_%d\n' % (locus, idx))


llist = ['MSI-11', 'MSI-14', 'H-10', 'HSPH1-T17', 'BAT-26', 'BAT-25', 'MSI-04', 'MSI-06', 'MSI-07', 'MSI-01', 'MSI-03', 'MSI-09', 'H-09', 'H-08', 'H-01', 'H-03', 'H-02', 'H-04', 'H-07', 'H-06', 'H-05']
top_7 = ['BAT-26', 'MSI-07', 'MSI-09', 'H-06', 'MSI-06', 'MSI-04', 'HSPH1-T17']
top_9 = ['BAT-26', 'MSI-07', 'MSI-09', 'H-06', 'MSI-06', 'MSI-04', 'HSPH1-T17', 'H-02', 'H-03']


#metrics('/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set', llist)
#print_series_names()
num_bins = 50
write_lengths_df('/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/training_set', top_9, num_bins)
write_features(top_9, num_bins)
