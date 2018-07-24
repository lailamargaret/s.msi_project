import pysam
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
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
	Return: dict msi_loci {locus_name | [chromosome, start_pos, end_pos]}
	"""
	msi_loci = {}
	with open(infile, 'r') as f:
        	for line in f:
                	if (line.startswith('Name')): #skip the header row
                        	continue
                	fields = line.split('\t') #list of each column
                	fields[3] = fields[3].replace('\n', '')
                	msi_loci[fields[0]] = [fields[1], fields[2], fields[3]]
	return msi_loci

def count_reads(runfile_loc, locus, run_quality_threshold = 0.1, flank_length = 7,  flank_mismatch = 3, print_full = False, plot = False, show_reads = False, return_mms = False):
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
	start_position = msi_loci[locus][1]
	end_position = msi_loci[locus][2]
	chromosome = msi_loci[locus][0]
	run_length = int(end_position) - int(start_position) + 1	

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
	        elif (float(flag) / len(read) > run_quality_threshold) :
       		        rejected_runs.append(read)
        	else:
        	        accepted_runs.append(read)
	
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
		if show_reads:
			print 'accepted reads: '
			for read in accepted_runs:
				print read
	temp = runfile_loc.replace('.bam', '')
        temp = temp.split('/')
        bam_name = temp[-1]

	while plot:	
		saveloc = runfile_loc.replace('.bam', '.%s.graph.png' % locus).replace('/bam', '/graphs')
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
	for each in lst:
		sum += each
	return float(sum)/len(lst)

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


def print_mm_depth(bams, mismatch = 3, length = 7):
	'''
	Brief: used for optimizing flank length and mismatch parameters, prints edit distance and read depth to a file indicated below
	Args: list, int, int
	Returns: none
	'''
	infile = '/home/upload/msi_project/msi_loci_edited.txt'
	msi_loci = get_msi_loci(infile)	
	
	outfile = '/home/upload/msi_project/mm_depth_analysis/subset-mismatch_depth-%d-%d.txt' % (mismatch, length)
	with open (outfile, 'w') as f:
		f.write('#mismatches: %d, flank length: %d\n' % (mismatch, length))
		f.write('\t')
		for locus in msi_loci:
                	f.write(locus + '\t\t\t')
                f.write('\n\t')
		for i in range(len(msi_loci)):
			f.write('f1 mm\tf2mm\t% accepted reads\t')
		f.write('\n')
		for bam in bams:
			f.write(bam.split('/')[-1].replace('.bam', '') + '\t')
			for locus in msi_loci:
				accepted_reads, f1_mm, f2_mm, num_reads = count_reads(bam, locus, length, mismatch, return_mms = True)
				if num_reads == 0:
					percent_accepted = "no coverage"
				else:
					percent_accepted = float(len(accepted_reads)) / float(num_reads) * 100 
				f.write(str(f1_mm) + '\t' + str(f2_mm) + '\t' + str(percent_accepted) + '\t')
			f.write('\n')



def report(bams, annotations, mismatch = 3, length = 7):
	'''
	Brief: Reports to file the MSI status of a patient
	Args: list, int, int
	Returns: none
	'''
	infile = '/home/upload/msi_project/msi_loci_edited.txt'
	msi_loci = get_msi_loci(infile)

	outfile = '/home/upload/msi_project/subset_statuses.txt'
		
	with open (outfile, 'w') as f:
		f.write('#mismatch: %s, flank length: %s\n' % (str(mismatch), str(length)))
		f.write('BAM\tCOUNT\tSTATUS\tKNOWN STATUS\n')
		for bam in bams:
			status_marker = 0
			bam_name = bam.split('/')[-1].replace('A.bam', '')
			f.write(bam_name + '\t')
			for locus in msi_loci:
				accepted_reads = (count_reads(bam, locus, length, mismatch))
				average = avg_length(accepted_reads)
				length = int(msi_loci[locus][2]) - int(msi_loci[locus][1]) + 1
				if average == 'Insufficient reads' or len(accepted_reads) < 25:
					pass
				elif average > length:
					status_marker += 1
				else:
					status_marker -= 1
				#f.write(str(length) + '\t' + str(average) + '\t')
			msi_status = 'MSS'
			if status_marker > 0:
				msi_status = 'MSI'
			if bam_name in annotations:
				known_status = annotations[bam_name]
			else:
				known_status = 'Not reported'
			f.write(str(status_marker) + '\t' + msi_status + '\t' + known_status + '\n')

def bw_plot(bams, msi_loci):
	'''
	Brief: Print a candlestick plot of the number of accepted reads at each locus
	Args: lst, dict
	Return: none
	'''
	data = []
	label = []
	for locus in msi_loci:
		temp = []
		label.append(locus)
		for bam in bams:
			runs, favg, bavg, num_reads = count_reads(bam, locus, return_mms = True)
			temp.append(num_reads)
		data.append(temp)
	plt.boxplot(data, labels = label)
	plt.xticks(rotation = 90)
	plt.title('Subset Read Depth')
	plt.savefig('/home/upload/msi_project/depth_plot.png')

def get_msi_annotations():
	'''
	Brief: Generate a dict containing all loci from txt file with their msi status
	Args: none
	Return: dict
	'''
	msi_annotations = {}
	annotations_file = '/home/upload/msi_project/UCEC_annotations.txt'
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

	annotations_file = '/home/upload/msi_project/COAD_READ_annotations.txt'
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

def status_plot(bams):
	'''
	Brief: Generate histograms showing average length of MS region depending on MSI status
	Args: lst
	Return: none 
	'''
	for locus in msi_loci:
		msi_avgs = []
		mss_avgs = []
		for bam in bams:
			bam_name = bam.split('/')[-1].replace('A.bam', '')			
		 	if bam_name in annotations:
				status = annotations[bam_name]
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
			saveloc = '/home/upload/msi_project/status_correlations/all_bamfiles/%s_dist.png' % locus
			plt.savefig(saveloc)			
			plt.clf()	

# ----------- Main --------------
msi_loci = get_msi_loci('/home/upload/msi_project/msi_loci_edited.txt')

#store bamfiles in a list
directory = '/home/upload/msi_project/tcga_bam/tumor_bams/subset'
bamfiles = scan_files(directory)

annotations = get_msi_annotations()

bw_plot(bamfiles, msi_loci)
'''

for bam in bamfiles:
	count_reads(bam, 'Mono-27', print_full = True, show_reads = True)
'''
'''
outfile = '/home/upload/msi_project/msi_data.txt'
with open (outfile, 'w') as f:
	f.write('\t')
	for locus in msi_loci:
		f.write(locus + '\t')
		f.write('Accepted Reads\t')
	f.write('\n')
	f.write('Reference length\t')
	for locus in msi_loci:
		f.write(str(int(msi_loci[locus][2]) - int(msi_loci[locus][1])) + '\t\t') 
	f.write('\n')
	for bam in bamfiles:
		f.write(bam + '\t')
		for locus in msi_loci:
			temp = count_reads(bam, locus)
			f.write(str(avg_length(temp)) + '\t')
			f.write(str(len(temp)) + ' reads \t')
		f.write('\n')
'''

