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
                f.write('Known status\n')
                avg_stdevs = [] #average for each bamfile all loci
                for bam in bams:
                        bam_name = bam.split('/')[-1].replace('A.bam', '')
                        f.write(bam_name + '\t')
                        locus_stdevs = []
                        for locus in _MSI_LOCI:
                                accepted_reads = count_reads.count(bam, locus, flank_length = length, flank_mismatch = mismatch)
                                if len(accepted_reads) == 0:
                                        std_dev = 'n/a'
                                        f.write('n/a\t')
                                else:
                                        lengths = [len(e) for e in accepted_reads]
                                        std_dev = np.std(lengths)
                                        f.write(str(std_dev) + '\t')

                                locus_stdevs.append(std_dev)

                        bam_stdev = lstproc.avg_value(locus_stdevs)
                        avg_stdevs.append(bam_stdev)

                        if bam_name in _ANNOTATIONS:
                                known_status = _ANNOTATIONS[bam_name]
                        else:
                                known_status = 'Not reported'

                        f.write(known_status + '\n')

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


def bw_plot(bams):
        """
        Brief: Print a candlestick plot of the number of accepted reads at each locus
        Args: lst, dict
        Return: none
        """
        data = []
        label = []
        for locus in _MSI_LOCI:
                temp = []
                label.append(locus)
                for bam in bams:
                        runs, favg, bavg, num_reads = count_reads.count(bam, locus, return_mms = True)
                        temp.append(num_reads)
                data.append(temp)
        plt.boxplot(data, labels = label)
        plt.xticks(rotation = 90)
        plt.title('Subset Read Depth')
        plt.savefig('/home/upload/msi_project/subsetA_depth_plot.png')

def status_plot(bams):
        """
        Brief: Generate histograms showing average length of MS region depending on MSI status
        Args: lst
        Return: none
        """
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


def earthmover_distance(p1, p2):
    dist1 = {x: count / len(p1) for (x, count) in Counter(p1).items()}
    dist2 = {x: count / len(p2) for (x, count) in Counter(p2).items()}
    solver = pywraplp.Solver('earthmover_distance', pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)
 
    variables = dict()
 
    # for each pile in dist1, the constraint that says all the dirt must leave this pile
    dirt_leaving_constraints = defaultdict(lambda: 0)
 
    # for each hole in dist2, the constraint that says this hole must be filled
    dirt_filling_constraints = defaultdict(lambda: 0)
 
    # the objective
    objective = solver.Objective()
    objective.SetMinimization()
 
    for (x, dirt_at_x) in dist1.items():
        for (y, capacity_of_y) in dist2.items():
            amount_to_move_x_y = solver.NumVar(0, solver.infinity(), 'z_{%s, %s}' % (x, y))
            variables[(x, y)] = amount_to_move_x_y
            dirt_leaving_constraints[x] += amount_to_move_x_y
            dirt_filling_constraints[y] += amount_to_move_x_y
            objective.SetCoefficient(amount_to_move_x_y, euclidean_distance(x, y))
 
    for x, linear_combination in dirt_leaving_constraints.items():
        solver.Add(linear_combination == dist1[x])
 
    for y, linear_combination in dirt_filling_constraints.items():
        solver.Add(linear_combination == dist2[y])
 
    status = solver.Solve()
    if status not in [solver.OPTIMAL, solver.FEASIBLE]:
        raise Exception('Unable to find feasible solution')
 
    return objective.Value()


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

