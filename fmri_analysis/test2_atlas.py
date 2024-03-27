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
source_dir = "/Volumes/LaCie/derivatives"
data_path = os.path.join(source_dir,'sub-01','ses-003','func')

### Time series  ####

funcim = "sub-01_ses-003_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"

confounds = fmriprep.load_confounds(os.path.join(data_path,"sub-01_ses-003_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"))
## Default strategy is motion + high pass + wm_csf

from nilearn.maskers import NiftiMapsMasker

from nilearn import datasets
dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
#dataset = datasets.fetch_atlas_destrieux_2009()
dataset = datasets.fetch_atlas_aal(version='SPM12')
#dataset = datasets.fetch_coords_seitzman_2018()

atlas_filename = dataset["maps"]
labels = dataset["labels"]

print(f"Atlas ROIs are located in nifti image (4D) at: {atlas_filename}")

from nilearn.maskers import NiftiLabelsMasker
masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=False)



time_series = masker.fit_transform(os.path.join(data_path,funcim),confounds[0])# os.path.join(data_path,"sub-01_ses-001_task-rest_desc-confounds_timeseries.tsv"))




correlation_measure = ConnectivityMeasure(
    kind="correlation",
    standardize="zscore_sample",
)
correlation_matrix = correlation_measure.fit_transform([time_series])[0]

print(correlation_matrix.shape)
#L = [x[1] for x in labels]
#print(L[1:])

#print(dataset["description"])
print(len(labels))
print(labels)


np.fill_diagonal(correlation_matrix, 0)
plotting.plot_matrix(
    correlation_matrix, labels=labels, colorbar=True, vmax=0.8, vmin=-0.8,reorder = True
)

np.savetxt(os.path.join(source_dir,"correlation_mat_aal_sub-01_ses-003.csv"), correlation_matrix, delimiter=",")
with open(os.path.join(source_dir,"labels.csv"), "w") as output:
    output.write(str(labels))

plotting.show()

