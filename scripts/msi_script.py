#related modules
import count_reads
import lstproc
import bamprocess
from constants import _MSI_LOCI, _QUALITY_THRESHOLDS, _ANNOTATIONS, _MSS_LOCUS_DATA

#imported modules
import pysam
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

#from sklearn import metrics
#from ortools.linear_solver import pywraplp
# --------- Functions -----------
def get_z_score(bamfile, locus, mismatch = 2, length = 7):
	if float(_MSS_LOCUS_DATA[locus][1]) == 0:
		return 'error'

	accepted_reads = count_reads.count(bamfile, locus, flank_length = length, flank_mismatch = mismatch)
	if len(accepted_reads) == 0:
		return 'error'
	else:
		lengths = [len(e) for e in accepted_reads]
		std_dev = np.std(lengths)
		z_score = ((float(_MSS_LOCUS_DATA[locus][0]) - float(std_dev)) / float(_MSS_LOCUS_DATA[locus][1]))
		z_score = abs(z_score)
		return z_score

def locus_histogram(bamfiles):
	"""	
	Brief: Produces histograms of the z-scores of bamfiles, separated by annotation (MSI or MSS)
	Args: lst
	Return: None
	""" 
	for locus in _MSI_LOCI:
                msi_avgs = []
                mss_avgs = []
                for bam in bamfiles:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        if bam_name in _ANNOTATIONS:
                                status = _ANNOTATIONS[bam_name]
                                if status != 'MSI' and status != 'MSS':
                                        continue
				z_score = get_z_score(bam, locus)
                                if status == 'MSI':
                                        msi_avgs.append(z_score)
                                elif status == 'MSS':
                                        mss_avgs.append(z_score)
                if len(msi_avgs) != 0 or len(mss_avgs) != 0:
                        plt.hist([msi_avgs, mss_avgs], color = ['yellow', 'orange'], label = ['MSI', 'MSS'])
                        plt.title('%s Z-Score Distribution - Subset A' % locus)
                        plt.legend(loc = 'best')
                        plt.xlabel = ('Average MS length (bp)')
                        plt.ylabel('Number of BAM files')
                        saveloc = '/home/upload/msi_project/diag_analysis/method_3/locus_plots/%s_dist.png' % locus
                        plt.savefig(saveloc)
                        plt.clf()



def confusion_matrix(bamfiles, metric_list, reporting_threshold = .9):
	real_pos = 0
	real_neg = 0
	false_pos = 0
	false_neg = 0
	idx = 0

	for bam in bamfiles:
		bam_name = bam.split('/')[-1].replace('A.bam', '')
		if bam_name in _ANNOTATIONS:
			known_status = _ANNOTATIONS[bam_name]
		else:
			known_status = 'Not reported'
		if metric_list[idx] < reporting_threshold:
			msi_status = 'MSS'
		else:
			msi_status = 'MSI'	
		if msi_status == 'MSI' and known_status == 'MSI':
			real_pos += 1
		elif msi_status == 'MSI' and known_status == 'MSS':
			false_pos += 1
		elif msi_status == 'MSS' and known_status == 'MSI':
			false_neg += 1	
		elif msi_status == 'MSS' and known_status == 'MSS':
			real_neg += 1
		idx += 1

	total = real_pos + false_pos + false_neg + real_neg
	perc_correct = (real_pos + false_pos)/float(total)
	#print 'Threshold: %f' % reporting_threshold
	#print 'True positive: %d' % real_pos
	#print 'False positive: %d' % false_pos
	#print 'True negative: %d' % real_neg
	#print 'False negative: %d' % false_neg
	#print 'Percentage correct: %f' % perc_correct

	return real_pos, false_pos 

'''	
def roc_plot(bamfiles):
	z_scores = []
	for bam in bamfiles:
		z_scores.append(get_z_score(bam, 'H-06')
	#for i in 
	#tpr, fpr = confusion_matrix(bamfiles, z_scores)	
'''

# ----------- Main --------------
#store bamfiles in a list
directory = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'
bamfiles = bamprocess.scan_files(directory)

print(bamfiles[1])
print(get_z_score(bamfiles[1], 'H-09'))









