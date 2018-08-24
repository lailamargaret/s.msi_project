#directory processing module

import os

def scan_files(directory, tumor_only = False):
	"""
	Brief: Given a directory, produces a list of paths to every bam file within the directory
	Args: directory, str - path to the directory that will be processed
	      tumor_only, bool - if TRUE, will only process tumor bams. If FALSE, will process tumor and matched normal 
	Return: bam_files, lst - a list of paths to the bam files in the directory
	"""
        bam_files = []
        for filename in os.listdir(directory):
                if tumor_only:
                        if filename.endswith('01A.bam'):
                                bam_files.append(os.path.join(directory, filename))
                else:
                        if filename.endswith('.bam'):
                                bam_files.append(os.path.join(directory, filename))
        return bam_files

