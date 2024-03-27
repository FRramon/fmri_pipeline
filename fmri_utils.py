###########################################################################
###           		 Provide utils functions                            ###
###########################################################################


import os
import shutil
import csv
import pandas as pd
import re


def isEmpty(path): 
    if os.path.exists(path) and not os.path.isfile(path): 
        #print(os.listdir(path))
        # Checking if the directory is empty or not 
        if not os.listdir(path): 
            return True
        else: 
            return False 
    else: 
        return True


def get_list_sessions(base_dir,groups,session):

	# Get subjects ids which participated in session i

	source_data_dir = os.path.join(base_dir,groups)

	subjects_raw = os.listdir(source_data_dir)
	pattern = re.compile(r'^sub-\d')
	subjects = [s for s in subjects_raw if pattern.match(s)]

	haveSes = []

	for s in subjects:
		ses_id = 'ses-' + session
		ses_path = os.path.join(source_data_dir,s,ses_id)
		if not isEmpty(ses_path):
			haveSes.append(s)

	haveSes = [s[4:] for s in haveSes]

	return haveSes

def get_subfunc_sessions(base_dir,groups):

	# Get subjects ids which participated in at least one functional mri

	source_data_dir = os.path.join(base_dir,groups)

	subjects_raw = os.listdir(source_data_dir)
	pattern = re.compile(r'^sub-\d')
	subjects = [s for s in subjects_raw if pattern.match(s)]

	haveSes1 = []
	haveSes2 = []
	haveSes3 = []

	for s in subjects:
		ses_id1 = 'ses-001'
		ses_id2 = 'ses-002'
		ses_id3 = 'ses-003'

		ses_path1 = os.path.join(source_data_dir,s,ses_id1,'func')
		ses_path2 = os.path.join(source_data_dir,s,ses_id2,'func')
		ses_path3 = os.path.join(source_data_dir,s,ses_id3,'func')

		if not isEmpty(ses_path1):
			haveSes1.append(s)
		elif not isEmpty(ses_path2):
			haveSes2.append(s)
		elif not isEmpty(ses_path3):
			haveSes3.append(s)

	#haveSes = [s[4:] for s in haveSes]

	have_at_least_one = [s for s in subjects if s in haveSes1 or s in haveSes2 or s in haveSes3]

	return have_at_least_one

	# print(len(haveSes))

	# transformed_list = ','.join(haveSes)
	# result_list = [transformed_list]

	# print(result_list)

	# return result_list

def get_list_fmriprep(base_dir,groups,session):

# Get subjects ids which participated in session i

	source_data_dir = os.path.join(base_dir,groups)

	subjects_raw = os.listdir(source_data_dir)
	pattern = re.compile(r'^sub-\d')
	subjects = [s for s in subjects_raw if pattern.match(s)]

	haveSes = []
	fmri_on_ses = []

	for s in subjects:
		ses_id = 'ses-' + session
		ses_path = os.path.join(source_data_dir,s,ses_id)
		if not isEmpty(ses_path) and os.path.exists(os.path.join(base_dir,groups,s,"anat")):
			haveSes.append(s)
			file_path = os.path.join(base_dir,groups,s,ses_id,"func",f"{s}_{ses_id}_task-rest_desc-confounds_timeseries.json")

			if os.path.isfile(os.path.join(base_dir,groups,s,ses_id,"func",f"{s}_{ses_id}_task-rest_desc-confounds_timeseries.json")):
				#print(f"{s} fmriprep done")
				fmri_on_ses.append(s)
		
				#print(f"{s} fmriprep not done")

	fmri_on_ses = [s[4:] for s in fmri_on_ses]

	return fmri_on_ses




# L = get_list_fmriprep("/Volumes/LaCie/derivatives","Patients","003")
# print(L)





