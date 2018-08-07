#count reads module
import pysam
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


from constants import _MSI_LOCI, _QUALITY_THRESHOLDS
import lstproc
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


def count(runfile_loc, locus, run_quality_threshold = 0.1, flank_length = 7,  flank_mismatch = 2, print_full = False, plot = False, show_reads = False, return_mms = False):
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
                        print 'average run length: %f' % lstproc.avg_length(accepted_runs)
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
                label = mpatches.Patch(color = '#859ec6', label = 'Avg length: %f' % lstproc.avg_length(accepted_runs))
                label2 = mpatches.Patch(color = 'blue', label = 'Reference length: %d' % run_length)
                label3 = mpatches.Patch(color = 'white', label ='Runs: %d' % len(accepted_runs))
                plt.legend(loc = 'upper left', fontsize = 'x-small', handles = [label, label2, label3])
                plt.xlabel('polynucleotide run length')
                plt.ylabel('runs')
                fig.suptitle(title, y = 0.995)
                plt.savefig(saveloc)
                plot = False

        if return_mms:
                return accepted_runs, lstproc.avg_value(front_mms), lstproc.avg_value(back_mms), len(reads)
        else:
                return accepted_runs




















