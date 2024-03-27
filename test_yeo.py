from nilearn import datasets
from nilearn.input_data import NiftiMasker


yeo = datasets.fetch_atlas_yeo_2011()
print(
    "Yeo atlas nifti image (3D) with 17 parcels and liberal mask "
    f" is located at: {yeo['thick_17']}"
)

data = datasets.fetch_development_fmri(n_subjects=10)

print(
    "Functional nifti images (4D, e.g., one subject) "
    f"are located at : {data.func[0]!r}"
)
print(
    "Counfound csv files (of same subject) are located "
    f"at : {data['confounds'][0]!r}"
)


from nilearn.connectome import ConnectivityMeasure
from nilearn.maskers import MultiNiftiLabelsMasker

# ConenctivityMeasure from Nilearn uses simple 'correlation' to compute
# connectivity matrices for all subjects in a list
connectome_measure = ConnectivityMeasure(
    kind="correlation",
    standardize="zscore_sample",
)

# useful for plotting connectivity interactions on glass brain
from nilearn import plotting
print(yeo["thick_17"])
# create masker using MultiNiftiLabelsMasker to extract functional data within
# atlas parcels from multiple subjects using parallelization to speed up the
# computation
masker = MultiNiftiLabelsMasker(
    labels_img=yeo["thick_17"],
    standardize="zscore_sample",
    standardize_confounds="zscore_sample",
    memory="nilearn_cache",
    n_jobs=2,
)



# # extract time series from all subjects
# time_series = masker.fit_transform(data.func, confounds=data.confounds)

# # calculate correlation matrices across subjects and display
# correlation_matrices = connectome_measure.fit_transform(time_series)

# # Mean correlation matrix across 10 subjects can be grabbed like this,
# # using connectome measure object
# mean_correlation_matrix = connectome_measure.mean_

# # grab center coordinates for atlas labels
# coordinates = plotting.find_parcellation_cut_coords(labels_img=yeo["thick_17"])

# # plot connectome with 80% edge strength in the connectivity
# plotting.plot_connectome(
#     mean_correlation_matrix,
#     coordinates,
#     edge_threshold="80%",
#     title="Yeo Atlas 17 thick (func)",
# )
# plotting.show()