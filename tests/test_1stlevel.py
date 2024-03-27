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



import numpy as np

from nilearn import datasets, plotting
from nilearn.glm.first_level import (
    FirstLevelModel,
    make_first_level_design_matrix,
)
from nilearn.maskers import NiftiSpheresMasker

# %%
# Prepare data and analysis parameters
# ------------------------------------
# Prepare the data.
#adhd_dataset = datasets.fetch_adhd(n_subjects=1)

# Prepare timing
t_r = 2.0
slice_time_ref = 0.0
n_scans = 178

# Prepare seed
pcc_coords = (0, -53, 26)

# %%
# Extract the seed region's time course
# -------------------------------------
# Extract the time course of the seed region.
seed_masker = NiftiSpheresMasker(
    [pcc_coords],
    radius=10,
    detrend=True,
    standardize="zscore_sample",
    low_pass=0.1,
    high_pass=0.01,
    t_r=2.0,
    memory="nilearn_cache",
    memory_level=1,
    verbose=0,
)
seed_time_series = seed_masker.fit_transform(os.path.join(data_path,funcim))
frametimes = np.linspace(0, (n_scans - 1) * t_r, n_scans)

# %%
# Plot the time course of the seed region.
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(9, 3))
ax = fig.add_subplot(111)
ax.plot(frametimes, seed_time_series, linewidth=2, label="seed region")
ax.legend(loc=2)
ax.set_title("Time course of the seed region")
plt.show()

# %%
# Estimate contrasts
# ------------------
# Specify the contrasts.
design_matrix = make_first_level_design_matrix(
    frametimes,
    hrf_model="spm",
    add_regs=seed_time_series,
    add_reg_names=["pcc_seed"],
)
dmn_contrast = np.array([1] + [0] * (design_matrix.shape[1] - 1))
contrasts = {"seed_based_glm": dmn_contrast}

# %%
# Perform first level analysis
# ----------------------------
# Setup and fit GLM.
first_level_model = FirstLevelModel(t_r=t_r, slice_time_ref=slice_time_ref)
first_level_model = first_level_model.fit(
    run_imgs=os.path.join(data_path,funcim), design_matrices=design_matrix
)

# %%
# Estimate the contrast.
print("Contrast seed_based_glm computed.")
z_map = first_level_model.compute_contrast(
    contrasts["seed_based_glm"], output_type="z_score"
)

# Saving snapshots of the contrasts
filename = "dmn_z_map.png"
display = plotting.plot_stat_map(
    z_map, threshold=3.0, title="Seed based GLM", cut_coords=pcc_coords
)
display.add_markers(
    marker_coords=[pcc_coords], marker_color="g", marker_size=300
)
display.savefig(filename)
print(f"Save z-map in '{filename}'.")

# %%
# Generating a report
# -------------------
# It can be useful to quickly generate a
# portable, ready-to-view report with most of the pertinent information.
# This is easy to do if you have a fitted model and the list of contrasts,
# which we do here.
from nilearn.reporting import make_glm_report

report = make_glm_report(
    first_level_model,
    contrasts=contrasts,
    title="ADHD DMN Report",
    cluster_threshold=15,
    min_distance=8.0,
    plot_type="glass",
)

# %%
# We have several ways to access the report:

# report  # This report can be viewed in a notebook
# report.open_in_browser()

# or we can save as an html file
# from pathlib import Path
# output_dir = Path.cwd() / "results" / "plot_adhd_dmn"
# output_dir.mkdir(exist_ok=True, parents=True)
# report.save_as_html(output_dir / 'report.html')
