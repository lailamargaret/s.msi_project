#fuzzy matching algorithm module

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

