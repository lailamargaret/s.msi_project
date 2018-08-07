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

from sklearn import metrics as skmetrics
from ortools.linear_solver import pywraplp

# --------- Functions -----------
def get_z_score(bamfile, locus, mismatch = 2, length = 7):
	"""
	Brief: Calculates the standard deviation z-score for the given bam file and locus using the mean stdev and stdev from MSS sample
	Args: str, str, int, int
	Return: str (if no locus data, no accepted reads), or float
	"""
	if float(_MSS_LOCUS_DATA[locus][1]) == 0:
		return 'error'

	accepted_reads = count_reads.count(bamfile, locus, flank_length = length, flank_mismatch = mismatch)
	if len(accepted_reads) == 0:
		return 'error'
	else:
		lengths = [len(e) for e in accepted_reads]
		std_dev = np.std(lengths)
		z_score = ((float(_MSS_LOCUS_DATA[locus][0]) - float(std_dev)) / float(_MSS_LOCUS_DATA[locus][1]))
		return z_score

def get_dist_mode(bamfile, locus, mismatch = 2, length = 7):
	"""
	Brief: For a given bam file and locus, calculates the average distance of the accepted polynucleotide runs from the MSS sample mode length
	Args: str, str, int, int
	Return: str (if no accepted reads at that locus), float 
	"""
	accepted_reads = count_reads.count(bamfile, locus, flank_length = length, flank_mismatch = mismatch)
	if len(accepted_reads) == 0:
		return 'error'	
	
	lengths = [len(e) for e in accepted_reads]
	distance = 0
	count = 0
	for read in lengths:
		distance += abs(read - float(_MSS_LOCUS_DATA[locus][2]))
		count += 1
	return distance / count

def get_avg_length(bamfile, locus, mismatch = 2, length = 7):
	accepted_reads = count.count_reads(bamfile, locus, flank_length = length, flank_mismatch = mismatch)
	if len(accepted_reads) == 0:
		return 'error'
	else:	
			
	

def report_dist_mode(bams, mismatch = 2, length = 7):
        """
        Brief: Reports to file the MSI status of a patient to compare with the known status, determining MSI status
                based on absolute distance from the mode
        Args: list, int, int
        Returns: none
        """
        outfile = '/home/upload/msi_project/diag_analysis/method_2/subsetA_dist_mode.txt'
        
	with open (outfile, 'w') as f:
                f.write('Bam name\n')
		for bam in bams:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        f.write(bam_name + '\t')
                        for locus in _MSI_LOCI: #iterate over all loci
                        	dist_mode = get_dist_mode(bam, locus)
				f.write(str(dist_mode) + '\t')
			if bam_name in _ANNOTATIONS:
                                known_status = _ANNOTATIONS[bam_name]
                        else:
                                known_status = 'Not reported'
                        f.write(known_status + '\n')

def report(bams, method, mismatch = 2, length = 7):
        """
        Brief: Reports to file some metric for each bamfile and locus depending on parameter, known MSI status
        Args: list, str, int, int
        Returns: none
        """
	if method == 'z_score':
	        method_no = 3
        elif method == 'dist_mode':
                method_no = 2
        elif method == 'emd':
                method_no = 4
        else:
                print 'Error: unknown method'
                return

        outfile = '/home/upload/msi_project/diag_analysis/method_%d/subsetA_%s.txt' % (method_no, method)

        with open (outfile, 'w') as f:
                f.write('Bam name\n')
                for bam in bams:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        f.write(bam_name + '\t')
                        for locus in _MSI_LOCI: #iterate over all loci
                                if method_no == 2:
					metric = get_dist_mode(bam, locus)
                                elif method_no == 3:
					metric = get_z_score(bam, locus)
				elif method_no == 4:
					metric = get_emd(bam,locus)
				f.write(str(metric) + '\t')
                        if bam_name in _ANNOTATIONS:
                                known_status = _ANNOTATIONS[bam_name]
                        else:
                                known_status = 'Not reported'
                        f.write(known_status + '\n')


def mode_length(in_list):
        '''
        Brief: returns the mode length of elements in the input list. If there is no mode, returns 'error', if there are multiple modes, returns the mode closest to the arithmetic average. If all modes have the same distance to arithmetic mean, returns the smallest mode
        Args: list
        Return: float if no mode (returns mean), int (mode exists)
        '''
        lengths = []
        dict_counts = {}
        list_counts = []
        if len(in_list) == 0:
                return 'error'
        for i in range(len(in_list)):
                lengths.append(len(in_list[i]))

        for i in lengths:
                counti = lengths.count(i)
                list_counts.append(counti)
                dict_counts[i] = counti
        maxcount = max(list_counts)
        if maxcount == 1: #when there is no mode, use the arithmetic mean
                return avg_value(lengths)
        else:
                modelist = []
                for key, item in dict_counts.iteritems():
                        if item == maxcount:
                                modelist.append(str(key))
                if len(modelist) == 1: #there is exactly 1 mode
                        return modelist[0]
                else: #more than 1 mode, return the one closest to the arithmetic mean
                        average = avg_length(in_list)
                        distances = []
                        for i in modelist:
                                distances.append(abs(float(i) - average))
                        min_index = 0
                        current_min = distances[0]
                        for i in range(len(distances)):
                                if distances[i] < current_min:
                                        min_index = i
                                        current_min = distances[i]


def roc_curve(bamfiles, method, plot = False):
	if method == 'z-score':
		method_no = 3
	elif method == 'dist_mode':
		method_no = 2
	elif method == 'emd':
		method_no = 4
	else:
		print 'Error: unknown method'
		return

	outfile = '/home/upload/msi_project/diag_analysis/method_%d/roc_plots/thresholds.txt' % method_no
	with open (outfile, 'w') as f:
		for locus in _MSI_LOCI:
			f.write(locus + '\t')
		f.write('\n')
		for locus in _MSI_LOCI:
			metrics = []
			statuses = []
			for bam in bamfiles:
				bam_name = bam.split('/')[-1].replace('A.bam', '')
				if bam_name in _ANNOTATIONS:
					status = _ANNOTATIONS[bam_name]
					if status != 'MSI' and status != 'MSS':
						continue
					elif status == 'MSI':
						label = 1
					elif status == 'MSS':
						label = 0
					
					if method_no == 2:
						metric = get_dist_mode(bam, locus)
					elif method_no == 3:
						metric = get_z_score(bam, locus)
					elif method_no == 4:
						metric = get_emd(bam, locus)
					
					if metric == 'error':
						continue
					else:
						metrics.append(metric)
						statuses.append(label)
			
			if len(metrics) < 2 or len(statuses) < 2:
				print 'No ROC produced for %s' % locus
				continue
			y = np.array(statuses)
			scores = np.array(metrics)
			fpr, tpr, thresholds = skmetrics.roc_curve(y, scores)
			roc_auc = skmetrics.auc(fpr, tpr)
		
			if plot:
				plt.title('%s ROC Curve' % locus)
				plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
				plt.legend(loc = 'lower right')
				plt.plot([0, 1], [0, 1], 'r--')
				plt.xlim([0, 1])
				plt.ylim([0, 1]) 
				plt.ylabel('True Positive Rate (Sensitivity)')
				plt.xlabel('False Positive Rate (1 - Specificity)')
				saveloc = '/home/upload/msi_project/diag_analysis/method_%d/roc_plots/ROC_%s.png' % (method_no, locus)
				plt.savefig(saveloc)
				plt.clf()
		
			scores = []
			for i in range(len(tpr)):
				scores.append(tpr[i] - fpr[i])
			max_score = scores[0]
			idx = 0
			for i in range(len(scores)):
				if scores[i] > max_score:
					max_score = scores[i]
					idx = i
			f.write(str(thresholds[idx]) + '\t')
					
	
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
				metric = get_dist_mode(bam, locus)
                                if metric == 'error':
					continue
				if status == 'MSI':
                                        msi_avgs.append(float(metric))
                                elif status == 'MSS':
                                        mss_avgs.append(float(metric))
                if len(msi_avgs) != 0 or len(mss_avgs) != 0:
                        plt.hist([msi_avgs, mss_avgs], color = ['yellow', 'orange'], label = ['MSI', 'MSS'])
                        plt.title('%s Distance from Mode Distribution - Subset A' % locus)
                        plt.legend(loc = 'best')
                        plt.xlabel = ('Average Distance from Mode')
                        plt.ylabel('Number of BAM files')
                        saveloc = '/home/upload/msi_project/diag_analysis/method_2/locus_plots/%s_dist.png' % locus
                        plt.savefig(saveloc)
                        plt.clf()
		else:
			print 'No data to plot for %s' % locus

def report_z_score(bams, mismatch = 2, length = 7):
        '''
        Brief: Reports z-score of standard deviation of a sample compared to mean/stdev of MSS known population
        Args: lst, float, int, int
        Return: none
        '''
        outfile = '/home/upload/msi_project/diag_analysis/method_3/subsetA_statuses_zscore.txt'

        with open (outfile, 'w') as f:
                f.write('#mismatch: %s, flank length: %s\n' % (str(mismatch), str(length)))
                f.write('locus\t')
                for locus in _MSI_LOCI:
                        f.write(locus + '\t')
                f.write('Status\n')
                for bam in bams:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        f.write(bam_name + '\t')
                        for locus in _MSI_LOCI:
                                z_score = get_z_score(bam, locus, mismatch = mismatch, length = length)
                                f.write(str(z_score) + '\t')

                        if bam_name in _ANNOTATIONS:
                                known_status = _ANNOTATIONS[bam_name]
                        else:
                                known_status = 'Not reported'

                        f.write(known_status + '\n')


# ----------- Main --------------
#store bamfiles in a list
directory = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'
bamfiles = bamprocess.scan_files(directory)

roc_curve(bamfiles, 'z-score')



