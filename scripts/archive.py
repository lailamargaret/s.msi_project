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
