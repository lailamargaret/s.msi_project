#methods for calculating statistical scores per bamfile/locus, reporting out
from constants import _MSI_LOCI, _QUALITY_THRESHOLDS, _ANNOTATIONS, _MSS_LOCUS_DATA
import count_reads
import bamprocess
import lstproc

import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from sklearn import metrics as skmetrics

def get_stdev_score(bamfile, locus, mismatch = 2, length = 7):
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
        accepted_reads = count_reads.count(bamfile, locus, flank_length = length, flank_mismatch = mismatch)
        if len(accepted_reads) == 0:
                return 'error'

def get_num_lengths(bamfile, locus, mismatch = 2, length = 7):
 	"""
	Brief: Calculates the z-score of the number of different lengths present in the file at the given locus as compared to the mss sample
	Args: str, str, int, int
	Return: str (if no accepted reads), float
	"""
	accepted_reads = count_reads.count(bamfile, locus, flank_length = length, flank_mismatch = mismatch)
	if len(accepted_reads) == 0:
		return 'error'
	else:
		lengths = [len(e) for e in accepted_reads]
		num_lengths = len(set(accepted_reads))
		z_score = ((float(_MSS_LOCUS_DATA[locus][3]) - num_lengths) / float(_MSS_LOCUS_DATA[locus][4]))
		return z_score
	
def mss_num_lengths(bamfiles):
	outfile = '/home/upload/msi_project/diag_analysis/method_5/mss_training_data.txt'
	with open (outfile, 'w') as f:
		f.write('BAM\t')
		for locus in _MSI_LOCI:
			f.write(locus + '\t')
		f.write('\n')
		for bam in bamfiles:
			bam_name = bam.split('/')[-1].replace('A.bam', '')
			f.write(bam_name + '\t')
			for locus in _MSI_LOCI:
				accepted_reads = count_reads.count(bam, locus)
				lengths = [len(e) for e in accepted_reads]
				num_lengths = len(set(lengths)) 
				if num_lengths == 0:
					f.write('n/a\t')
				else:
					f.write(str(num_lengths) + '\t')
			f.write('\n')

def report_mode_length(directory):
	bamfiles = bamprocess.scan_files(directory)
	with open ('/home/upload/msi_project/ML/modes.txt', 'w') as f:
		f.write('locus\t')
		for locus in _MSI_LOCI:
			f.write(locus + '\t')
		f.write('\nmode\t')
		for locus in _MSI_LOCI:
			all_reads = []
			for bam in bamfiles:
				all_reads.extend(count_reads.count(bam, locus))
			mode = lstproc.mode_length(all_reads)
			f.write(mode + '\t')
	
	


def report(bams, method, mismatch = 2, length = 7):
        """
        Brief: Reports to file some metric for each bamfile and locus depending on parameter, known MSI status
        Args: list, str, int, int
        Returns: none
        """
        if method == 'stdev':
                method_no = 3
        elif method == 'dist_mode':
                method_no = 2
        elif method == 'emd':
                method_no = 4
	elif method == 'num_lengths':
		method_no = 1
        else:
                print 'Error: unknown method'
                return

        outfile = '/home/upload/msi_project/diag_analysis/method_%d/subsetA_%s.txt' % (method_no, method)

        with open (outfile, 'w') as f:
                f.write('BAM\t')
		for locus in _MSI_LOCI:
			f.write(locus + '\t')
		f.write('Status\n')
                for bam in bams:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        f.write(bam_name + '\t')
                        for locus in _MSI_LOCI: #iterate over all loci
                                if method_no == 2:
                                        metric = get_dist_mode(bam, locus)
                                elif method_no == 3:
                                        metric = get_stdev_score(bam, locus)
                               # elif method_no == 4:
                                 #       metric = get_emd(bam, locus)
				elif method_no == 1:
					metric = get_num_lengths(bam, locus)
                                f.write(str(metric) + '\t')
                        if bam_name in _ANNOTATIONS:
                                known_status = _ANNOTATIONS[bam_name]
                        else:
                                known_status = 'Not reported'
                        f.write(known_status + '\n')


def roc_curve(bamfiles, method, plot = False):
        if method == 'stdev':
                method_no = 3
        elif method == 'dist_mode':
                method_no = 2
        elif method == 'emd':
                method_no = 4
	elif method == 'num_lengths':
		method_no = 1
        else:
                print 'Error: unknown method'
                return

        outfile = '/home/upload/msi_project/diag_analysis/method_%d/roc_plots/thresholds.txt' % method_no
        with open (outfile, 'w') as f:
		locus_thresholds = []
		used_loci = []
		accuracies = []
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
                                                metric = get_stdev_score(bam, locus) * -1
                                       # elif method_no == 4:
                                           #     metric = get_emd(bam, locus)
					elif method_no == 1:
						metric = get_num_lengths(bam, locus)

                                        if metric == 'error' or metric == '':
                                                continue
                                        
					else:
                                                metrics.append(metric)
                                                statuses.append(label)
                        
			if len(metrics) < 2 or len(statuses) < 2:
                                print 'No ROC produced for %s' % locus
                                continue
			else:
				used_loci.append(locus)
                        y = np.array(statuses)
                        scores = np.array(metrics)
                        fpr, tpr, thresholds = skmetrics.roc_curve(y, scores)
                        roc_auc = skmetrics.auc(fpr, tpr)
			accuracies.append(roc_auc)
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
                        locus_thresholds.append(thresholds[idx])
		f.write('locus')
		for locus in used_loci:
			f.write(locus + '\t')
		f.write('\nthreshold')
		for entry in locus_thresholds:
			f.write(str(entry) + '\t')
		f.write('\naccuracy')
		for auc in accuracies:
			f.write(str(auc) + '\t')


def locus_histogram(bamfiles, method):
        """
        Brief: Produces histograms of the score of some method of bamfiles, separated by annotation (MSI or MSS)
        Args: lst
        Return: None
        """
        if method == 'stdev':
                method_no = 3
        elif method == 'dist_mode':
                method_no = 2
        elif method == 'emd':
                method_no = 4
	elif method == 'num_lengths':
		method_no = 1 
        else:
                print 'Error: unknown method'
                return

        for locus in _MSI_LOCI:
                msi_avgs = []
                mss_avgs = []
                for bam in bamfiles:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        if bam_name in _ANNOTATIONS:
                                status = _ANNOTATIONS[bam_name]
                                if status != 'MSI' and status != 'MSS':
                                        continue

                                if method_no == 2:
                                        metric = get_dist_mode(bam, locus)
                                elif method_no == 3:
                                        metric = get_stdev_score(bam, locus)
                                #elif method_no == 4:
                                        #metric = get_emd(bam, locus)
				elif method_no == 1:
					metric = get_num_lengths(bam, locus)

                                if metric == 'error':
                                        continue
                                if status == 'MSI':
                                        msi_avgs.append(float(metric))
                                elif status == 'MSS':
                                        mss_avgs.append(float(metric))

                if len(msi_avgs) != 0 or len(mss_avgs) != 0:
                        plt.hist([msi_avgs, mss_avgs], color = ['yellow', 'orange'], label = ['MSI', 'MSS'])
                        plt.title['%s %s Distribution - Subset A' % (locus, method)]
                        plt.legend(loc = 'best')
                        plt.xlabel = ('%s' % method)
                        plt.ylabel('Number of BAM files')
                        saveloc = '/home/upload/msi_project/diag_analysis/method_%d/locus_plots/%s_dist.png' % (method_no, locus)
                        plt.savefig(saveloc)
                        plt.clf()
                else:
                        print 'No data to plot for %s' % locus



