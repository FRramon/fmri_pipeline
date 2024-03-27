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


# TR = 2
source_dir = "/Volumes/LaCie/derivatives/Patients"
data_path = os.path.join(source_dir,'sub-01','ses-003','func')


atlas = datasets.fetch_atlas_schaefer_2018(n_rois=1000, yeo_networks=7,resolution_mm = 2)
#atlas = datasets.fetch_atlas_destrieux_2009()
atlas_filename = atlas["maps"]
labels = atlas["labels"]

print(type(atlas_filename))
print(atlas_filename)

### Time series  ####

funcim = "sub-01_ses-003_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"

confounds = fmriprep.load_confounds(os.path.join(data_path,"sub-01_ses-003_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"))
from nilearn.maskers import NiftiMapsMasker, NiftiLabelsMasker

masker = NiftiLabelsMasker(
    labels_img=atlas_filename,
    standardize="zscore_sample",
    standardize_confounds="zscore_sample"
    #memory="nilearn_cache",
    # verbose=5,
)


time_series = masker.fit_transform(os.path.join(data_path,funcim),confounds[0])# os.path.join(data_path,"sub-01_ses-001_task-rest_desc-confounds_timeseries.tsv"))

print(time_series.shape)
print(confounds[0].shape)


####### CORRELATION MATRIX ############


correlation_measure = ConnectivityMeasure(
    kind="correlation",
    standardize="zscore_sample",
)

correlation_matrix = correlation_measure.fit_transform([time_series])[0]

print(correlation_matrix)
print(correlation_matrix.shape)

# Display the correlation matrix

# Mask out the major diagonal
np.fill_diagonal(correlation_matrix, 0)
plotting.plot_matrix(
    correlation_matrix, labels=labels, colorbar=True, vmax=0.8, vmin=-0.8,reorder = True
)

np.savetxt(os.path.join(source_dir,"correlation_mat_sub-01_ses-003.csv"), correlation_matrix, delimiter=",")
with open(os.path.join(source_dir,"labels.csv"), "w") as output:
    output.write(str(labels))

plotting.show()



