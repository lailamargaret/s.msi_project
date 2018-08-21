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
		#f.write('#Full %s dataset\nbam_bame\t' % setname)
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
		

llist = ['MSI-11', 'MSI-14', 'H-10', 'HSPH1-T17', 'BAT-26', 'BAT-25', 'MSI-04', 'MSI-06', 'MSI-07', 'MSI-01', 'MSI-03', 'MSI-09', 'H-09', 'H-08', 'H-01', 'H-03', 'H-02', 'H-04', 'H-07', 'H-06', 'H-05']

metrics('/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set', llist)
print_series_names()


