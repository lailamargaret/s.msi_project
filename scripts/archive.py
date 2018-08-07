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

