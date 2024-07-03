#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:05:49 2024

@author: francoisramon
"""



import numpy as np
import nibabel as nib
from nilearn.input_data import NiftiLabelsMasker
from nilearn import plotting
import matplotlib.pyplot as plt
import subprocess
from nilearn.connectome import ConnectivityMeasure
import seaborn as sns
from nilearn.interfaces import fmriprep
import matplotlib.pyplot as plt


# Load the NIfTI label file

filepath = '/Volumes/LaCie/pipe_patients_10/pipe_patients_10/main_workflow/connectome/_ses_id_001_subject_id_01/labelconvert/mapflow/_labelconvert2/aparc.nii.gz'
mat = "/Volumes/LaCie/derivatives/Patients/sub-01/anat/sub-01_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5"
outfile = "/Users/francoisramon/Desktop/These/MNI_aparc.nii.gz"
label_img = 'aparc.nii.gz'

# applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
#           --in=/Volumes/LaCie/pipe_patients_10/pipe_patients_10/main_workflow/connectome/_ses_id_001_subject_id_01/labelconvert/mapflow/_labelconvert2/aparc.nii.gz \
#           --warp=/Volumes/LaCie/derivatives/Patients/sub-01/anat/sub-01_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5 \
#           --out=/Users/francoisramon/Desktop/These/MNI_aparc.nii.gz


command =f"antsApplyTransforms -d 3 -i {filepath} -r /Users/francoisramon/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -t {mat} -n NearestNeighbor -o {outfile} "

#subprocess.run(command, shell = True)

label_img = outfile

# # Load a functional image to extract signals from (example functional image)
func_img = '/Volumes/LaCie/derivatives/Patients/sub-01/ses-001/func/sub-01_ses-001_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
confounds = fmriprep.load_confounds(func_img,ica_aroma = 'basic')

# # Initialize the NiftiLabelsMasker
masker = NiftiLabelsMasker(labels_img=label_img, standardize="zscore_sample", standardize_confounds="zscore_sample",
ory="nilearn_cache",
bose=5, memory='nilearn_cache', verbose=5)

# # Extract time series from the functional image using the labels
time_series = masker.fit_transform(func_img,confounds[0])

correlation_measure = ConnectivityMeasure(kind="correlation",standardize="zscore_sample")
correlation_matrix = correlation_measure.fit_transform([time_series])[0]
np.fill_diagonal(correlation_matrix, 0)
print(correlation_matrix.shape[0])

sns.heatmap(correlation_matrix)
plt.show()

transformed_parcellation_img = nib.load(outfile)
plotting.view_img_on_surf(transformed_parcellation_img, surf_mesh='fsaverage', threshold=None)
plotting.show()

# Interactive 3D plot of the transformed parcellation image
view = plotting.view_img(transformed_parcellation_img, threshold=None)
view.open_in_browser()


# Visualize the label image (optional)
plotting.plot_roi(filepath, title='Label image')
plotting.show()

plotting.plot_roi(outfile, title='Label image')
plotting.show()




