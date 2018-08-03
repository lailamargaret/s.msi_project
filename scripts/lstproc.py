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





