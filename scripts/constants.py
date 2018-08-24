#BRIEF: Produces lists and dicts for use in other functions

#------------------------------------------------------------ Current -------------------------------------------------------------------------------#
def get_msi_loci():
        """
        Brief: Creates a dictionary storing information about each msi locus
        Args: string infile
        Return: dict _MSI_LOCI {locus_name | [chromosome, start_pos, end_pos]}
        """
	infile = '/home/upload/msi_project/loci/msi_loci_edited.txt'
        msi_loci = {}
        with open(infile, 'r') as f:
                for line in f:
                        if (line.startswith('Name')): #skip the header row
                                continue
                        fields = line.split('\t') #list of each column
                        fields[3] = fields[3].replace('\n', '')
                        msi_loci[fields[0]] = [fields[1], fields[2], fields[3]]
        return msi_loci

def get_msi_annotations():
        """
        Brief: Generate a dict containing all bam files from txt file with their msi status
	       NOTE: Saves all 'MSI-L' cases as 'MSS', 'MSI-H' as 'MSI'
        Args: none
        Return: dict
        """
        msi_annotations = {}
        annotations_file = '/home/upload/msi_project/annotations/UCEC_annotations.txt'
        with open (annotations_file, 'r') as f:
                for line in f:
                        fields = line.split('\t')
                        if fields[0] == 'Cancer Type Detailed':
                                continue
                        if fields[22] == 'MSI-L' or fields[22] == 'MSS':
                                msi_status = 'MSS'
			if fields[22] == 'MSI-H':
				msi_status = 'MSI'
                        else:
                                msi_status = fields[22]
                        msi_annotations[fields[4]] = msi_status

        annotations_file = '/home/upload/msi_project/annotations/COAD_READ_annotations.txt'
        with open (annotations_file, 'r') as f:
                for line in f:
                        fields = line.split('\t')
                        if fields[0] == 'Cancer Type Detailed':
                                continue
                        if fields[20] == 'MSI-L' or fields[20] == 'MSS':
                                msi_status = 'MSS'
                        if fields[20] == 'MSI-H':
                                msi_status = 'MSI'
			else:
				msi_status = fields[20]
                        msi_annotations[fields[4]] = msi_status

        return msi_annotations

def get_full_annotations():
	"""
	Brief: Generates a dict containing the RAW msi status from annotations files
	Args: none
	Returns: dict
	"""
	full_annotations = {}
	annotations_file = '/home/upload/msi_project/annotations/UCEC_annotations.txt'
	with open(annotations_file, 'r') as f:
		for line in f:
			fields = line.split('\t')
			if fields[0] == 'Cancer Type Detailed':
				continue
			msi_status = fields[22]
			full_annotations[fields[4]] = msi_status
	annotations_file = '/home/upload/msi_project/annotations/COAD_READ_annotations.txt'
	with open(annotations_file, 'r') as f:
		for line in f:
			fields = line.split('\t')
			if fields[0] == 'Caner Type Detailed':
				continue
			msi_status = fields[20]
			full_annotations[fields[4]] = msi_status
	return full_annotations

def get_ml_modes():
	"""
	Brief: Uses the special ML mode training set to determine the mode of a dataset using a defined mode calling function
	Args: None
	Return: dict
	"""
        mlmodes = '/home/upload/msi_project/ML/modes.txt'
        ml_modes = {}
        with open(mlmodes, 'r') as f:
                lines = f.readlines()

                fields1 = lines[0].split('\t')
                fields1.pop(0)

                fields2 = lines[1].split('\t')
                fields2.pop(0)

                for i in range(len(fields1)):
                        fields1[i] = fields1[i].replace('\n', '')
                        fields2[i] = fields2[i].replace('\n', '')
                        ml_modes[fields1[i]] = fields2[i]
        return ml_modes
		
	
#-------------------------------------------------- Archive ----------------------------------------------------------------#
def get_mss_locus_data():
	"""
	Brief: Stores the average standard deviation, standard deviation of the average standard dev, mode, average length, and standard deviation of a MSS training set to create features
	Args: None
	Return: dict
	"""
        mss_input = '/home/upload/msi_project/diag_analysis/MSS_training_data.txt'

        #{'locus' : ['mean', 'stdev', 'mode']}
        mss_locus_data = {}

        with open (mss_input, 'r') as f:
                lines = f.readlines()

                fields1 = lines[0].split('\t')
                fields1.pop(0)

                fields2 = lines[1].split('\t')
                fields2.pop(0)

                fields3 = lines[2].split('\t')
                fields3.pop(0)
	
		fields4 = lines[3].split('\t')
		fields4.pop(0)

		fields5 = lines[4].split('\t')
		fields5.pop(0)
		
		fields6 = lines[5].split('\t')
		fields6.pop(0)

                for i in range(len(fields1)):
                        fields1[i] = fields1[i].replace('\n', '')
                        fields2[i] = fields2[i].replace('\n', '')
                        fields3[i] = fields3[i].replace('\n', '')
                        fields4[i] = fields4[i].replace('\n', '')
			fields5[i] = fields5[i].replace('\n', '')
			fields6[i] = fields6[i].replace('\n', '')
			mss_locus_data[fields1[i]] = [fields2[i], fields3[i], fields4[i], fields5[i], fields6[i]]
        return mss_locus_data

#-------------------------------------- Main ----------------------------------------------------------------------#
#hard-coded quality thresholds for some non-homogenous loci
_QUALITY_THRESHOLDS =  {'MSI-11': .25, 'MSI-12': .25, 'MSI-01': .5, 'BAT-25': .18}

#populate the lsts and dicts of commonly used information
_MSI_LOCI = get_msi_loci()
_ANNOTATIONS = get_msi_annotations()
_FULL_ANNOTATIONS = get_full_annotations()
_MSS_LOCUS_DATA = get_mss_locus_data()
_ML_MODES = get_ml_modes()


