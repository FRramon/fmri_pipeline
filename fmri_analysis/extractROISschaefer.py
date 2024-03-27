
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from fmri_utils import *
from nilearn import datasets


atlas_schaefer = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7,resolution_mm = 2)
#atlas = datasets.fetch_atlas_destrieux_2009()
atlas_filename = atlas_schaefer["maps"]
labels_schaefer = atlas_schaefer["labels"]
labels_index_schaefer = range(1,len(labels_schaefer)+1)

new = [str(x)[2:-1] for x in labels_schaefer]
# print(new)

with open(("/Users/francoisramon/Desktop/These/fMRI_pipeline/labels_schaefer.csv"), 'w') as myfile:
    wr = csv.writer(myfile)#, quoting=csv.QUOTE_ALL)
    wr.writerow(new)

# print(type(labels_schaefer[1]))
# print(labels_schaefer[1][])
