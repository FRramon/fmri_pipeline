# Author : François Ramon
# Creation data : 11mars2024


#compute_correlation_matrix


import os
import numpy as np
import matplotlib.pyplot as plt
import nilearn
from nilearn.input_data import NiftiMasker
from nilearn.interfaces import fmriprep
import pandas as pd
from nilearn import plotting
from nilearn.connectome import ConnectivityMeasure
from nilearn import datasets
from nilearn.maskers import NiftiMapsMasker
from nilearn.maskers import NiftiLabelsMasker

# TR = 2


### Time series  ####

def load_dataset(atlas):
	if atlas == "harvard_oxford":
		dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
	elif atlas == "aal":
		dataset = datasets.fetch_atlas_aal(version='SPM12')
	elif atlas == "msdl":
		dataset = datasets.fetch_atlas_msdl()
	elif atlas == "schaefer":
		dataset = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7,resolution_mm = 2)
	elif atlas == "destrieux":
		dataset = datasets.fetch_atlas_destrieux_2009()
	return dataset

def get_time_series(source_dir,group,sub,ses,atlas):

	data_path = os.path.join(source_dir,group,sub,ses,'func')
	funcim = f"{sub}_{ses}_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"

	## Default strategy is motion + high pass + wm_csf
	confounds = fmriprep.load_confounds(os.path.join(data_path,funcim),ica_aroma = 'basic')

	dataset = load_dataset(atlas)
	atlas_filename = dataset["maps"]
	labels = dataset["labels"]

	#masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=False)
	if atlas == "harvard_oxford" or atlas == "aal":
		masker = NiftiMapsMasker(
	    	maps_img=atlas_filename,
	    	standardize="zscore_sample",
	    	standardize_confounds="zscore_sample",
	    	memory="nilearn_cache",
	    	verbose=5,
		)

	elif atlas == "schaefer" or atlas == "destrieux":
		masker = NiftiLabelsMasker(
    	labels_img=atlas_filename,
    	standardize="zscore_sample",
    	standardize_confounds="zscore_sample",
    	memory="nilearn_cache",
    	verbose=5,
	)

	time_series = masker.fit_transform(os.path.join(data_path,funcim),confounds[0])

	return time_series,dataset


def compute_correlation_matrix(source_dir,group,sub,ses,atlas,kind):

	time_series,dataset = get_time_series(source_dir,group,sub,ses,atlas)
	labels = dataset["labels"]



	correlation_measure = ConnectivityMeasure(
	    kind=kind,
	    standardize="zscore_sample",
	)
	correlation_matrix = correlation_measure.fit_transform([time_series])[0]
	np.fill_diagonal(correlation_matrix, 0)

	# if atlas == "harvard_oxford":
	# 	np.fill_diagonal(correlation_matrix, 0)
	# 	plotting.plot_matrix(
	# 	    correlation_matrix, labels=labels[1:], colorbar=True, vmax=0.8, vmin=-0.8,reorder = True
	# 	)
	# 	#plotting.show()
	# elif atlas == "aal":
	# 	np.fill_diagonal(correlation_matrix, 0)
	# 	plotting.plot_matrix(
	# 	    correlation_matrix, labels=labels, colorbar=True, vmax=0.8, vmin=-0.8,reorder = True
	# 	)
		#plotting.show()


	filename = f"{sub}_{ses}_{atlas}_correlation_mat.csv"
	np.savetxt(os.path.join(source_dir,group,sub,ses,"functional_connectivity",filename), correlation_matrix, delimiter=",")

	
# atlas = "msdl"
# dataset = load_dataset(atlas)

# labels = dataset.labels
# print(labels)





