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


from nilearn import datasets


ho_atlas = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
schaefer_atlas = datasets.fetch_atlas_schaefer_2018(n_rois = 400,yeo_networks = 7)


# The parcellation is already loaded into memory
#parcellation = schaefer_atlas.maps

# Retrieve fsaverage5 surface dataset for the plotting background. It contains
# the surface template as pial and inflated version and a sulcal depth maps
# which is used for shading
#fsaverage = datasets.fetch_surf_fsaverage()

# Display Destrieux parcellation on fsaverage5 pial surface using nilearn
from nilearn import plotting

plotting.plot_roi(schaefer_atlas.maps,draw_cross = False)


# plotting.plot_surf_roi(fsaverage['infl_left'], roi_map=parcellation,
#                         hemi='left', view='lateral',
#                         bg_map=fsaverage['sulc_left'], bg_on_data=True,
#                         darkness=.5)
plotting.show()

#print(destrieux_atlas.labels)
