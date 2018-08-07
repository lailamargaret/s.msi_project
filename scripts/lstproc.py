#list processing module


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



