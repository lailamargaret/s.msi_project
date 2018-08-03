#directory processing module

import os

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

