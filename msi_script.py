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
def scan_files(directory, tumor_only = False):
	bam_files = []
	for filename in os.listdir(directory):
		if tumor_only:
			if filename.endswith('01A.bam'):
				bam_files.append(os.path.join(directory, filename))
		else:
			if filename.endswith('.bam'):
				bam_files.append(os.path.join(directory, filename))
	return bam_files 



def approximate_match(t, p, maxdist, return_mm = False):
	"""
	Brief: find the closest match to substring p within t with a maximum edit distance maxdist
		return the approximate match substring, true/false depending on if a close enough match found
	Args: string t, string p, int maxdist
	Return: bool is_matched, string match
	"""
        is_matched = False
	match = ''
	partial_matches = []
        for i in range(0, len(t) - len(p) + 1):
                nmm = 0
                for j in range(0, len(p)):
                        if t[i+j] != p[j]:
                                nmm += 1
                                if nmm > maxdist:
                                        break
                if nmm <= maxdist:
                        partial_matches.append((i, nmm))
        if len(partial_matches) == 0:
		if return_mm:
			return is_matched, match, 'nm'
		else:
			return is_matched, match 
	#find the index at which nmm is the lowest
        min = partial_matches[0][1]
	index = partial_matches[0][0]
        for i in range(len(partial_matches)):
                if partial_matches[i][1] < min:
                        index = partial_matches[i][0]
                        min = partial_matches[i][1]
        match = t[index:index+len(p)]
	is_matched = True
	if return_mm:
		return is_matched, match, min
	else:
		return is_matched, match

def get_msi_loci(infile):
	"""
	Brief: Creates a dictionary storing information about each msi locus
	Args: string infile
	Return: dict _MSI_LOCI {locus_name | [chromosome, start_pos, end_pos]}
	"""
	_MSI_LOCI = {}
	with open(infile, 'r') as f:
        	for line in f:
                	if (line.startswith('Name')): #skip the header row
                        	continue
                	fields = line.split('\t') #list of each column
                	fields[3] = fields[3].replace('\n', '')
                	_MSI_LOCI[fields[0]] = [fields[1], fields[2], fields[3]]
	return _MSI_LOCI

def count_reads(runfile_loc, locus, run_quality_threshold = 0.1, flank_length = 7,  flank_mismatch = 2, print_full = False, plot = False, show_reads = False, return_mms = False):
	"""
	Brief: Uses hg38 ref to choose locus flanking regions of size flank_length, looks for matching regions up to flank_mismatch edit distance
		Filters polynucleotide runs that are more than run_quality_threshold different from reference run base
		If toggled, prints detailed information about the run
	Args: string runfile_loc, string locus
	Optional args: float run_quality_threshold, int flank_length, int flank_mismatch, bool print_full
	Return: list accepted_runs

	"""
	runfile = pysam.AlignmentFile(runfile_loc, 'rb')
	
	#Locus information
	start_position = _MSI_LOCI[locus][1]
	end_position = _MSI_LOCI[locus][2]
	chromosome = _MSI_LOCI[locus][0]
	run_length = int(end_position) - int(start_position) + 1	
	#increase quality threshold if short run
	if run_length < 10:
		run_quality_threshold = .2
	#if the quality threshold is defined, use that regardless of run length
	if locus in _QUALITY_THRESHOLDS:
		run_quality_threshold = _QUALITY_THRESHOLDS[locus]


	#find the flanking regions of size specified
	ref_seq = pysam.FastaFile('/StampFileShare/refs/hg38.fa')
	front_flank = (ref_seq.fetch('chr' + str(chromosome), int(start_position) - (flank_length + 1), int(start_position) - 1)).upper()
	back_flank = (ref_seq.fetch('chr' + str(chromosome), int(end_position), int(end_position) + flank_length)).upper()
	nucleotide = (ref_seq.fetch('chr' + str(chromosome), int(start_position), int(start_position) + 1)).upper()
	
	#filter runs
	reads = []
	polynucleotide_runs = []
	accepted_runs = []
	rejected_runs = []
	front_mms = []
	back_mms = []

	#gather reads that have coverage in the desired region
	for read in runfile.fetch("chr" + str(chromosome), int(start_position) - 1, int(end_position)):
	        parts = str(read).split('\t')
	        reads.append(parts[9])

	#filter for reads that contain both flanking regions
	for read in reads:
		if show_reads:
			print '\nOriginal read: %s' % read
	        cropped = []
		if return_mms:
			matched1, matched_front, front_mm = approximate_match(read, front_flank, flank_mismatch, return_mm = True)
		else:
			matched1, matched_front = approximate_match(read, front_flank, flank_mismatch)
		if matched1:
	                cropped = read.split(matched_front)
	        	if show_reads:	
				print 'Front reference flank: %s' % front_flank
				print 'Front flank: %s' % matched_front
				print 'Trim at front flank: %s' % cropped[1]
		if cropped != []:
			if return_mms:
                        	matched2, matched_back, back_mm = approximate_match(cropped[1], back_flank, flank_mismatch, return_mm = True)
			else:
                        	matched2, matched_back = approximate_match(cropped[1], back_flank, flank_mismatch)
	                if matched2:
	                        poly_seq = cropped[1].split(matched_back)
				if show_reads:
					print 'Back reference flank: %s' % back_flank
					print 'Back flank: %s' % matched_back
					print 'Trim at back flank: %s ' % poly_seq[0]
	                        polynucleotide_runs.append(poly_seq[0])
		if matched1 and matched2 and return_mms:
			if front_mm != 'nm' and back_mm != 'nm':
				front_mms.append(front_mm)
				back_mms.append(back_mm)

	#filter for polynucleotide runs that have a fraction of incorrect bases lower than the threshold
	for read in polynucleotide_runs:
	        flag = 0
	        for i in range(len(read)):
	                if read[i] != nucleotide:
	                        flag += 1
	        if (len(read) == 0):
	                pass
	        elif (float(flag) / len(read) < run_quality_threshold) :
       		        accepted_runs.append(read)
        	else:
        	        rejected_runs.append(read)
	
	runfile.close
	
	#print relevant information if toggled
	if print_full:
		print 'BAM: %s' % runfile_loc
		print 'locus: %s' % locus
		print 'total reads: %d' % len(reads)
		print 'total reads w flanking sequences: %d' % len(polynucleotide_runs)
		print 'total accepted reads: %d' % len(accepted_runs)
		if len(accepted_runs) == 0:
			print 'average run length: 0'
		else:
			print 'average run length: %f' % avg_length(accepted_runs)
		print 'total rejected reads: %d' % len(rejected_runs)
		if len(rejected_runs) > 0:
			print 'rejected reads:'
			for read in rejected_runs:
		                print read
		print 'accepted reads: '
		for read in accepted_runs:
			print read
	temp = runfile_loc.replace('.bam', '')
        temp = temp.split('/')
        bam_name = temp[-1]

	while plot:	
		saveloc = runfile_loc.replace('.bam', '.%s.graph.png' % locus).replace('/bam', '/length_distribution_graphs')
		title = bam_name + ' - ' + locus		
		
		if len(accepted_runs) < 10:
			print 'Insufficient number of runs to make a graph for %s' % title
			break
		
	        run_lengths = []
        	for each in accepted_runs:
                	run_lengths.append(len(each))
        	maxval = max(run_lengths)
		if run_length > maxval:
			maxval = run_length		

        	arr = np.array(run_lengths)
        	labels, counts = np.unique(arr, return_counts = True)
        	fig, axs = plt.subplots(tight_layout = True)
        	plt.bar(labels, counts, align = 'center', color = '#859ec6')
        	plt.gca().set_xticks(range(maxval+3))
		plt.axvline(x = run_length, ls = 'dashed')
		label = mpatches.Patch(color = '#859ec6', label = 'Avg length: %f' % avg_length(accepted_runs))
		label2 = mpatches.Patch(color = 'blue', label = 'Reference length: %d' % run_length)
		label3 = mpatches.Patch(color = 'white', label ='Runs: %d' % len(accepted_runs))
		plt.legend(loc = 'upper left', fontsize = 'x-small', handles = [label, label2, label3])
        	plt.xlabel('polynucleotide run length')
		plt.ylabel('runs')
        	fig.suptitle(title, y = 0.995)
        	plt.savefig(saveloc)	
		plot = False	
	
	if return_mms:
		return accepted_runs, avg_value(front_mms), avg_value(back_mms), len(reads)
 	else:	
		return accepted_runs

def avg_value(lst):
	'''
	Brief: calculates the average value from a list
	Args: lst
	Return: float
	'''
	if len(lst) == 0:
                return 'n/a'
	sum = 0
	total = 0
	for each in lst:
		if each == 'n/a':
			continue
		sum += each
		total += 1
	return float(sum)/total

def avg_length(lst):
	"""
	Brief: Calculates the average length of elements in a list
	Args: list lst
	Return: float 
	"""
	if len(lst) == 0:
		return 'Insufficient reads'
	sum = 0
	for each in lst:
               	sum += len(each)
	return float(sum)/len(lst)


def print_mm_depth(bams, mismatch = 2, length = 7):
	'''
	Brief: used for optimizing flank length and mismatch parameters, prints edit distance and read depth to a file indicated below
	Args: list, int, int
	Returns: none
	'''
	outfile = '/home/upload/msi_project/mm_depth_analysis/subsetA-mismatch_depth-%d-%d.txt' % (mismatch, length)
	with open (outfile, 'w') as f:
		f.write('#mismatches: %d, flank length: %d\n' % (mismatch, length))
		f.write('\t')
		for locus in _MSI_LOCI:
                	f.write(locus + '\t\t\t')
                f.write('\n\t')
		for i in range(len(_MSI_LOCI)):
			f.write('f1 mm\tf2mm\t% accepted reads\t')
		f.write('\n')
		for bam in bams:
			f.write(bam.split('/')[-1].replace('.bam', '') + '\t')
			for locus in _MSI_LOCI:
				accepted_reads, f1_mm, f2_mm, num_reads = count_reads(bam, locus, flank_length = length, flank_mismatch = mismatch, return_mms = True)
				if num_reads == 0:
					percent_accepted = "no coverage"
				else:
					percent_accepted = float(len(accepted_reads)) / float(num_reads) * 100 
				f.write(str(f1_mm) + '\t' + str(f2_mm) + '\t' + str(percent_accepted) + '\t')
			f.write('\n')



def report_num_lengths(bams, mismatch = 2, length = 7):
	'''
	Brief: Reports to file the MSI status of a patient to compare with the known status, determining MSI status
		based on number of different lengths
	Args: list, int, int
	Returns: none
	'''
	outfile = '/home/upload/msi_project/diag_analysis/method_1/subsetA_statuses_length.txt'
		
	with open (outfile, 'w') as f:
		f.write('#mismatch: %s, flank length: %s\n' % (str(mismatch), str(length)))
		#f.write('BAM\tNUM DIF ELEMS\tSTATUS\tKNOWN STATUS\tAGREE?\n')
		f.write('locus\t')
		for locus in _MSI_LOCI:
			f.write(locus + '\t')
		f.write('\n')
		for bam in bams:
			status_marker = 0
			bam_name = bam.split('/')[-1].replace('A.bam', '')
			f.write(bam_name + '\t')
			for locus in _MSI_LOCI:
				accepted_reads = (count_reads(bam, locus, flank_length = length, flank_mismatch = mismatch))
				if len(accepted_reads) == 0:
					f.write('n/a\t')
				else:
					f.write(str(len(set(accepted_reads))) + '\t')					
				status_marker += 1
			msi_status = 'MSS'
			if status_marker > 0:
				msi_status = 'MSI'
			if bam_name in _ANNOTATIONS:
				known_status = _ANNOTATIONS[bam_name]
			else:
				known_status = 'Not reported'
			agree = False
			if msi_status == known_status:
				agree == True
			f.write(msi_status + '\t' + known_status + '\t' + str(agree) + '\n')


def report_dist_mode(bams, mismatch = 2, length = 7):
        '''
        Brief: Reports to file the MSI status of a patient to compare with the known status, determining MSI status
                based on absolute distance from the mode
        Args: list, int, int
        Returns: none
	'''
        outfile = '/home/upload/msi_project/diag_analysis/method_2/subsetA_statuses_mode.txt'
	
	all_reads = []
	modes = []

	#Create 2D array of accepted read by locus and bam file
	for bam in bams:
		bam_reads = []
		for locus in _MSI_LOCI:
			accepted_reads = count_reads(bam, locus, flank_length = length, flank_mismatch = mismatch)
			bam_reads.append(accepted_reads)
		all_reads.append(bam_reads)	
	
	#Generate a list of the mode length for each locus
	for i in range(len(_MSI_LOCI)):
		for j in range(len(all_reads)):
			locus = []
			locus.extend(all_reads[j][i])
		mode = mode_length(locus)
		modes.append(mode)
	
	#find average distance from the mode for each bam each locus, average for all loci per bam, correlate with annotations 
        with open (outfile, 'w') as f:
                f.write('BAM\n')
        	for i in range(len(all_reads)): #iterate over all bam files
                        bam_name = bams[i].split('/')[-1].replace('A.bam', '')
                        f.write(bam_name + '\t')
                        for j in range(len(modes)): #iterate over all loci
				if modes[j] == 'error':
					avg_distance = 'low loc covg'
				elif len(all_reads[i][j]) == 0:
					avg_distance = 'low bam covg'
				else:
					total_distance = 0
					mode = modes[j]
					for read in all_reads[i][j]:
						total_distance += abs(float(mode) - len(read))
					avg_distance = float(total_distance) / len(all_reads[i][j])
				f.write(str(avg_distance) + '\t')		               
			if bam_name in _ANNOTATIONS:
				known_status = _ANNOTATIONS[bam_name]
                	else:
				known_status = 'Not reported'
			f.write(known_status + '\n')

def report_std_dev(bams, reporting_threshold = .9, mismatch = 2, length = 7):
	'''
	Brief: Reports to file the MSI status of a patient to compare with the known status, determining MSI status
		based on standard deviation
	Args: lst, int, int
	Return: None
	'''	
	outfile = '/home/upload/msi_project/diag_analysis/method_3/mss_training_set_statuses_stdev.txt'

	with open (outfile, 'w') as f:
                f.write('#mismatch: %s, flank length: %s, reporting_threshold: %s\n' % (str(mismatch), str(length), str(reporting_threshold)))
                f.write('locus\t')
                for locus in _MSI_LOCI:
                        f.write(locus + '\t')
                f.write('Average\tCall\tKnown status\n')
		avg_stdevs = [] #average for each bamfile all loci
                for bam in bams:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        f.write(bam_name + '\t')
			locus_stdevs = []
                        for locus in _MSI_LOCI:
				accepted_reads = (count_reads(bam, locus, flank_length = length, flank_mismatch = mismatch))
                                if len(accepted_reads) == 0:
					std_dev = 'n/a'
					f.write('n/a\t')
                                else:
					lengths = [len(e) for e in accepted_reads]
                                        std_dev = np.std(lengths)	
					f.write(str(std_dev) + '\t')
                        	
				locus_stdevs.append(std_dev)	
			
			bam_stdev = avg_value(locus_stdevs)
			avg_stdevs.append(bam_stdev)
			if len(locus_stdevs) == 0:
				msi_status = 'Indeterminate'
			else:
				if bam_stdev < reporting_threshold:
					msi_status = 'MSS'
				else:
					msi_status = 'MSI'

                        if bam_name in _ANNOTATIONS:
                                known_status = _ANNOTATIONS[bam_name]
                        else:
                                known_status = 'Not reported'
                        
			f.write(str(bam_stdev) + '\t' + msi_status + '\t' + known_status + '\n')
		
		return avg_stdevs

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

def get_z_score(bamfile, locus, mismatch = 2, length = 7):
	if float(_MSS_LOCUS_DATA[locus][1]) == 0:
		return 'error'

	accepted_reads = count_reads(bamfile, locus, flank_length = length, flank_mismatch = mismatch)
	if len(accepted_reads) == 0:
		return 'error'
	else:
		lengths = [len(e) for e in accepted_reads]
		std_dev = np.std(lengths)
		z_score = ((float(_MSS_LOCUS_DATA[locus][0]) - float(std_dev)) / float(_MSS_LOCUS_DATA[locus][1]))
		z_score = abs(z_score)
		return z_score

def locus_histogram(bamfiles):
	'''
	Brief: Produces histograms of the z-scores of bamfiles, separated by annotation (MSI or MSS)
	Args: lst
	Return: None
	''' 
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

			return modelist[min_index]
			

def bw_plot(bams):
	'''
	Brief: Print a candlestick plot of the number of accepted reads at each locus
	Args: lst, dict
	Return: none
	'''
	data = []
	label = []
	for locus in _MSI_LOCI:
		temp = []
		label.append(locus)
		for bam in bams:
			runs, favg, bavg, num_reads = count_reads(bam, locus, return_mms = True)
			temp.append(num_reads)
		data.append(temp)
	plt.boxplot(data, labels = label)
	plt.xticks(rotation = 90)
	plt.title('Subset Read Depth')
	plt.savefig('/home/upload/msi_project/subsetA_depth_plot.png')

def get_msi_annotations():
	'''
	Brief: Generate a dict containing all loci from txt file with their msi status
	Args: none
	Return: dict
	'''
	msi_annotations = {}
	annotations_file = '/home/upload/msi_project/annotations/UCEC_annotations.txt'
	with open (annotations_file, 'r') as f:
		for line in f:
			fields = line.split('\t')
			if fields[0] == 'Cancer Type Detailed':
				continue
			if fields[22] == 'MSI-L' or fields[22] == 'MSI-H':
				msi_status = 'MSI'
 			else:
				msi_status = fields[22]
			msi_annotations[fields[4]] = msi_status

	annotations_file = '/home/upload/msi_project/annotations/COAD_READ_annotations.txt'
	with open (annotations_file, 'r') as f:
		for line in f:
			fields = line.split('\t')
			if fields[0] == 'Cancer Type Detailed':
				continue
			if fields[20] == 'MSI-L' or fields[20] == 'MSI-H':
				msi_status = 'MSI'
			else:
				msi_status = fields[20]
			msi_annotations[fields[4]] = msi_status

	return msi_annotations

def get_mss_locus_data():
	
	mss_input = '/home/upload/msi_project/diag_analysis/method_3/MSS_training_data.txt'

	#{'locus' : ['mean', 'stdev']}
        mss_locus_data = {}

        with open (mss_input, 'r') as f:
                lines = f.readlines()

                fields1 = lines[0].split('\t')
                fields1.pop(0)

                fields2 = lines[1].split('\t')
                fields2.pop(0)

                fields3 = lines[2].split('\t')
                fields3.pop(0)

                for i in range(len(fields1)):
                        fields1[i] = fields1[i].replace('\n', '')
                        fields2[i] = fields2[i].replace('\n', '')
                        fields3[i] = fields3[i].replace('\n', '')
                        mss_locus_data[fields1[i]] = [fields2[i], fields3[i]]	
	return mss_locus_data


def status_plot(bams):
	'''
	Brief: Generate histograms showing average length of MS region depending on MSI status
	Args: lst
	Return: none 
	'''
	for locus in _MSI_LOCI:
		msi_avgs = []
		mss_avgs = []
		for bam in bams:
			bam_name = bam.split('/')[-1].replace('A.bam', '')			
		 	if bam_name in _ANNOTATIONS:
				status = _ANNOTATIONS[bam_name]
				if status != 'MSI' and status != 'MSS':
					continue
				average = avg_length(count_reads(bam, locus))
				if average == 'Insufficient reads':
					continue
				if status == 'MSI':
					msi_avgs.append(float(average))
				elif status == 'MSS':
					mss_avgs.append(float(average))
		if len(msi_avgs) != 0 or len(mss_avgs) != 0:
			plt.hist([msi_avgs, mss_avgs], color = ['red', 'blue'], label = ['MSI', 'MSS'])
			plt.title('%s Distribution (subset)' % locus)
			plt.legend(loc = 'best')
			plt.xlabel = ('Average MS length (bp)')
			plt.ylabel('Number of BAM files')
			saveloc = '/home/upload/msi_project/status_corr_dist/subsetA/%s_dist.png' % locus
			plt.savefig(saveloc)			
			plt.clf()	

# ----------- Main --------------
_MSI_LOCI = get_msi_loci('/home/upload/msi_project/loci/msi_loci_edited.txt')
_QUALITY_THRESHOLDS = {'MSI-11': .25, 'MSI-12': .25, 'MSI-01': .5, 'BAT-25': .18}
_ANNOTATIONS = get_msi_annotations()
_MSS_LOCUS_DATA = get_mss_locus_data()


#store bamfiles in a list
directory = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'
bamfiles = scan_files(directory)

locus_histogram(bamfiles)

'''
i = .7
stdevs = report_std_dev(bamfiles)
while i <= 1:
	confusion_matrix(bamfiles, stdevs, reporting_threshold = i)
	i += .05
'''
