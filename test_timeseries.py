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


funcim = "sub-01_ses-003_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
confounds = fmriprep.load_confounds(os.path.join(data_path,"sub-01_ses-003_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"))
from nilearn.maskers import NiftiMapsMasker

hpp_coords = [(-30, -28, -8)]

from nilearn.maskers import NiftiSpheresMasker

seed_masker = NiftiSpheresMasker(
    hpp_coords,
    radius=8,
    detrend=True,
    standardize="zscore_sample",
    standardize_confounds="zscore_sample",
    low_pass=0.1,
    high_pass=0.01,
    t_r=2,
    memory="nilearn_cache",
    memory_level=1,
    verbose=0,
)

seed_time_series = seed_masker.fit_transform(
    os.path.join(data_path,funcim), confounds=confounds[0]
)

from nilearn.maskers import NiftiMasker

brain_masker = NiftiMasker(
    smoothing_fwhm=6,
    detrend=True,
    standardize="zscore_sample",
    standardize_confounds="zscore_sample",
    low_pass=0.1,
    high_pass=0.01,
    t_r=2,
    memory="nilearn_cache",
    memory_level=1,
    verbose=0,
)


brain_time_series = brain_masker.fit_transform(
    os.path.join(data_path,funcim), confounds=confounds[0]
)

print(f"Seed time series shape: ({seed_time_series.shape})")
print(f"Brain time series shape: ({brain_time_series.shape})")

import matplotlib.pyplot as plt

plt.plot(seed_time_series)
plt.title("Seed time series (Left Anterior Hippocampus)")
plt.xlabel("Scan number")
plt.ylabel("Normalized signal")
plt.tight_layout()
plt.show()


import numpy as np

seed_to_voxel_correlations = (
    np.dot(brain_time_series.T, seed_time_series) / seed_time_series.shape[0]
)

print(
    "Seed-to-voxel correlation shape: (%s, %s)"
    % seed_to_voxel_correlations.shape
)
print(
    "Seed-to-voxel correlation: min = %.3f; max = %.3f"
    % (seed_to_voxel_correlations.min(), seed_to_voxel_correlations.max())
)

seed_to_voxel_correlations_img = brain_masker.inverse_transform(
    seed_to_voxel_correlations.T
)
display = plotting.plot_stat_map(
    seed_to_voxel_correlations_img,
    threshold=0.5,
    vmax=1,
    cut_coords=hpp_coords[0],
    title="Seed-to-voxel correlation (HPP seed)",
)
display.add_markers(
    marker_coords=hpp_coords, marker_color="g", marker_size=300
)
# At last, we save the plot as pdf.
from pathlib import Path

output_dir = Path.cwd() / "results" / "plot_seed_to_voxel_correlation"
output_dir.mkdir(exist_ok=True, parents=True)
print(f"Output will be saved to: {output_dir}")

display.savefig(output_dir / "hpp_seed_correlation.png")





