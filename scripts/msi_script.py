#related modules
import bamprocess
import count_reads
import methods
from constants import _MSI_LOCI, _QUALITY_THRESHOLDS, _ANNOTATIONS, _MSS_LOCUS_DATA, _ML_MODES

#imported modules
from ortools.linear_solver import pywraplp
import pandas as pd
import numpy as np
import subprocess
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
# --------------------------- Archive -------------------------------------------#

#store bamfiles in a list
directory = '/home/upload/msi_project/tcga_bam/tumor_bams/annotated/subset'
bamfiles = bamprocess.scan_files(directory)

def dataframe(directory):
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
cr_dir = '/home/upload/msi_project/tcga_bam/COAD-READ'
u_dir = '/home/upload/msi_project/tcga_bam/UCEC'


#--------------------------------------------------- Current ----------------------------------------------------------------#
def get_weights(infile):
  weights = {}
  with open(infile, 'r') as f:
    for line in f:
      fields = line.split('\t')
      feature = fields[0].split('/')[2]
      weight = fields[1].replace(']', '').replace('[', '').replace('\n', '')
      weights[feature] = weight
  return weights


def calling_function(directory, loci):
	tp = tn = fp = fn = 0 
	upper_threshold = 0.6
	lower_threshold = 0.4
	min_loci = 3
	bamfiles = bamprocess.scan_files(directory)
	correct_guesses = 0
	total_files = 0
	scores = {}
	msi_scores = [] 
	mss_scores = []
	for bam in bamfiles:
		features = {}
		prob_sum = 0
		agree = False
		num_loci = 0
		msicall = 0
		bam_name = bam.split('/')[-1].replace('A.bam', '')
		#store msi status as a boolean, skip if not MSI or MSS
		if _ANNOTATIONS[bam_name] == 'MSI':
			msi_status = 1
		elif _ANNOTATIONS[bam_name] == 'MSS':
			msi_status = 0
		else:
			continue
		print 'Bam file: %s' % bam_name
		for locus in loci:
			reads = count_reads.count(bam, locus)

			#if there are NO reads at this locus in the file, skip and go to the next locus
			if len(reads) == 0:
				continue
			num_loci += 1

			#compute the lengths of the reads, store in a list, append to the overall lengths list
			bam_lengths = [len(e) for e in reads]
      
			#compute the average length of a bam's reads, append to the overall avg_lengths list
			avg_length = np.mean(bam_lengths)
			features['%s_avg_len' % locus] = avg_length
			#compute number of different lengths of a bam's reads, append to overall nums_lengths list
			num_lengths = len(set(bam_lengths))
			features['%s_num_lens' % locus] = num_lengths
			#compute the standard deviation of a bam's reads, append to overall stdevs list
			stdev = np.std(bam_lengths)
       			features['%s_stdev' % locus] = stdev
			#compute the average distance from the mode length
			distance = 0
			count = 0
			for e in bam_lengths:
				distance += abs(e - float(_ML_MODES[locus]))
				count += 1
			dist_mode = distance / count
			features['%s_dist_mode' % locus] = dist_mode
			
	
		#prob = calc_prob(locus, avg_length, dist_mode, num_lengths, stdev) 
		prob = calc_prob(weights, features)
		
		scores[bam_name] = prob
		
		if msi_status:
			msi_scores.append(prob)
		else:
			mss_scores.append(prob)
		print 'Prob: %f' % prob
		#	prob_sum += prob
			#if prob >= 0.5:
				#msicall += 1
				#print 'Locus-level call: MSI'
			#else:
				#print 'Locus-level call: MSS'
		
		#no reads at any locus
		#if num_loci < min_loci:
		#	continue
     		
		#count how many files were called
		total_files += 1
		indeterminate_files = 0 
		
		#make a prediction
		guessed_status = 0
		#print '\nNum loci examined: %d' % num_loci
		#print 'Num MSI calls: %d' % msicall
		#perc_msi = float(msicall) / num_loci
		#print '%%MSI: %f' % perc_msi
		#if perc_msi >= threshold:
			#guessed_status = 1
		#avg_prob = prob_sum / num_loci
		if prob > upper_threshold:
			guessed_status = 1
			print 'Predicted status: MSI'
		elif prob < lower_threshold:
			print 'Predicted status: MSS'
		else:
			guessed_status = -1
			total_files -= 1 
			indeterminate_files += 1
			print 'Predicted status: Indeterminate'
			#print 'Predicted status: MSS'		
		
		print 'Known status: %s' % _ANNOTATIONS[bam_name]
		#decide whether prediction is correct
		if guessed_status != -1:
			if guessed_status == msi_status:
				correct_guesses += 1
				print 'Agree: YES'
			else:
				print 'Agree: NO'
		
		print '\n'
    
		#calculate tp, fp, tn, fn
		if guessed_status == 1 and msi_status == 1:
			tp += 1
		elif guessed_status == 1 and msi_status == 0:
			fp += 1
		elif guessed_status == 0 and msi_status == 1:
			fn += 1
		elif guessed_status == 0 and msi_status == 0:
			tn += 1
	
	print 'Summary:'
	print 'Loci examined: ' + ('\n'.join(loci)) 
	print 'Upper threshold: %f' % upper_threshold
	#print 'Threshold: %f' % upper_threshold
	print 'Lower threshold: %f' % lower_threshold
	#print 'Min no. loci: %d' % min_loci
	print 'Number of predictions: %s' % total_files
	print 'Correct predictions: %s' % correct_guesses
	print 'Indeterminate files: %s' % indeterminate_files
	print 'Accuracy: %f' % (float(correct_guesses) / total_files)
	print 'True pos: %d' % tp
	print 'True neg: %d' % tn
	print 'False pos: %d' % fp
	print 'False neg: %d' % fn
	print 'Sensitivity: %f' % (float(tp) / (tp + fn))
	print 'Specificity: %f' % (float(tn) / (tn + fp))

	bins = []
	i = 0.0
	while i < 1.05:
		bins.append(i)
		i += 0.05
	'''
	plt.hist([msi_scores, mss_scores], bins = bins, color = ['red', 'blue'], label = ['MSI', 'MSS'])
        plt.title('Model-Derrived Probabilities: p(MSI)')
	plt.legend(loc = 'best')
        plt.xlabel = ('p(MSI)')
        plt.ylabel('Number of BAM files')
        saveloc = '/home/upload/msi_project/ML/probability_distribution'
        plt.savefig(saveloc)
        plt.clf()
	'''
	return scores 

def calc_prob(weights, inputs):
  keys = list(inputs.keys())
  yhat = 0
  count = 0
  for key in keys:
    weight = float(weights[str(key)]) 
    value = inputs[str(key)]
    yhat += weight * value
    count += 1
  yhat += float(weights['bias_weights'])
  yhat *= -1
  prob = 1 / (1 + pow(2.718281, yhat))
  return prob	

def make_bam_list(lst):
	new_lst = []
	for each in lst:
		edited = each.split('/')[-1].replace('A.bam', '')
		new_lst.append(edited)

	return new_lst
#------------------------------------------------------------ main ---------------------------------------------------------------#
_CR_BAMS = make_bam_list(bamprocess.scan_files('/home/upload/msi_project/tcga_bam/COAD-READ'))
_U_BAMS = make_bam_list(bamprocess.scan_files('/home/upload/msi_project/tcga_bam/UCEC'))

val_directory = '/home/upload/msi_project/tcga_bam/tumor_bams/ml_set/test_set/edited'
loci = ['MSI-11', 'MSI-14', 'H-10', 'HSPH1-T17', 'BAT-26', 'BAT-25', 'MSI-04', 'MSI-06', 'MSI-07', 'MSI-01', 'MSI-03', 'MSI-09', 'H-09', 'H-08', 'H-01', 'H-03', 'H-02', 'H-04', 'H-07', 'H-06', 'H-05']

infile = '/home/upload/msi_project/ML/100000_0.001000_10_1_model_weights.txt'
weights = get_weights(infile)
scores = calling_function(val_directory, loci)

#print scores

cr_stable = []
cr_unstable = []
u_stable = []
u_unstable = []

for key in scores:
	print key
	cr = 0
	u = 0

	if key in _CR_BAMS:
		cr += 1
		if _ANNOTATIONS[key] == 'MSS':
			cr_stable.append(scores[key])
			print key
		elif _ANNOTATIONS[key] == 'MSI':
			cr_unstable.append(scores[key])
			print key
	elif key in _U_BAMS:
		u += 1
		if _ANNOTATIONS[key] == 'MSS':
			u_stable.append(scores[key])
			print key
		elif _ANNOTATIONS[key] == 'MSI':
			u_unstable.append(scores[key])
			print key
	
#print cr 
#print u
print len(u_stable)
print len(u_unstable)
print len(cr_stable)
print len(cr_unstable)

bins = []
i = 0.0
while i < 1.05:
	bins.append(i)
	i += 0.05

plt.hist([u_unstable, u_stable], bins = bins, color = ['darkorange', 'navajowhite'], label = ['MSI', 'MSS'])
plt.title('Model-Derrived Probabilities: p(MSI)')
plt.legend(loc = 'best')
plt.xlabel = ('p(MSI)')
plt.ylabel('Number of BAM files')
saveloc = '/home/upload/msi_project/ML/UCEC_probability_distribution'
plt.savefig(saveloc)
plt.clf()

plt.hist([cr_unstable, cr_stable], bins = bins, color = ['forestgreen', 'yellowgreen'], label = ['MSI', 'MSS'])
plt.title('Model-Derrived Probabilities: p(MSI)')
plt.legend(loc = 'best')
plt.xlabel = ('p(MSI)')
plt.ylabel('Number of BAM files')
saveloc = '/home/upload/msi_project/ML/COAD-READ_probability_distribution'
plt.savefig(saveloc)
plt.clf()














