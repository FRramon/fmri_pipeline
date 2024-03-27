#################################################
####            Create ROI to ROI file       ####
#################################################

import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from fmri_utils import *
from nilearn import datasets

# subject_raw = sys.argv[1]
# session_raw = sys.argv[2]
# source_dir = sys.argv[3]
# base_dir = sys.argv[4]

# subject_list = subject_raw.split(',')
# ses_list = session_raw.split(',')

# sys.path.append(source_dir + '/code')
# from run_parameters import *

# if group_raw == 'HealthyVolunteers':
# 	result_dir = os.path.join(source_dir , 'results_test2')
# elif group_raw == 'Patients':
# 	result_dir = os.path.join(source_dir , 'pipe_patients')
# elif group_raw == 'Controls':
# 	result_dir = os.path.join(source_dir , 'pipe_controls')

labels_aal = ['Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Olfactory_L', 'Olfactory_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', 'Rectus_L', 'Rectus_R', 'Insula_L', 'Insula_R', 'Cingulum_Ant_L', 'Cingulum_Ant_R', 'Cingulum_Mid_L', 'Cingulum_Mid_R', 'Cingulum_Post_L', 'Cingulum_Post_R', 'Hippocampus_L', 'Hippocampus_R', 'ParaHippocampal_L', 'ParaHippocampal_R', 'Amygdala_L', 'Amygdala_R', 'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', 'Lingual_L', 'Lingual_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R', 'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Precuneus_L', 'Precuneus_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Caudate_L', 'Caudate_R', 'Putamen_L', 'Putamen_R', 'Pallidum_L', 'Pallidum_R', 'Thalamus_L', 'Thalamus_R', 'Heschl_L', 'Heschl_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R', 'Cerebelum_Crus1_L', 'Cerebelum_Crus1_R', 'Cerebelum_Crus2_L', 'Cerebelum_Crus2_R', 'Cerebelum_3_L', 'Cerebelum_3_R', 'Cerebelum_4_5_L', 'Cerebelum_4_5_R', 'Cerebelum_6_L', 'Cerebelum_6_R', 'Cerebelum_7b_L', 'Cerebelum_7b_R', 'Cerebelum_8_L', 'Cerebelum_8_R', 'Cerebelum_9_L', 'Cerebelum_9_R', 'Cerebelum_10_L', 'Cerebelum_10_R', 'Vermis_1_2', 'Vermis_3', 'Vermis_4_5', 'Vermis_6', 'Vermis_7', 'Vermis_8', 'Vermis_9', 'Vermis_10']
labels_index_aal = range(1,len(labels_aal)+1)

labels_ho = ['Frontal Pole', 'Insular Cortex', 'Superior Frontal Gyrus', 'Middle Frontal Gyrus', 'Inferior Frontal Gyrus, pars triangularis', 'Inferior Frontal Gyrus, pars opercularis', 'Precentral Gyrus', 'Temporal Pole', 'Superior Temporal Gyrus, anterior division', 'Superior Temporal Gyrus, posterior division', 'Middle Temporal Gyrus, anterior division', 'Middle Temporal Gyrus, posterior division', 'Middle Temporal Gyrus, temporooccipital part', 'Inferior Temporal Gyrus, anterior division', 'Inferior Temporal Gyrus, posterior division', 'Inferior Temporal Gyrus, temporooccipital part', 'Postcentral Gyrus', 'Superior Parietal Lobule', 'Supramarginal Gyrus, anterior division', 'Supramarginal Gyrus, posterior division', 'Angular Gyrus', 'Lateral Occipital Cortex, superior division', 'Lateral Occipital Cortex, inferior division', 'Intracalcarine Cortex', 'Frontal Medial Cortex', 'Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex)', 'Subcallosal Cortex', 'Paracingulate Gyrus', 'Cingulate Gyrus, anterior division', 'Cingulate Gyrus, posterior division', 'Precuneous Cortex', 'Cuneal Cortex', 'Frontal Orbital Cortex', 'Parahippocampal Gyrus, anterior division', 'Parahippocampal Gyrus, posterior division', 'Lingual Gyrus', 'Temporal Fusiform Cortex, anterior division', 'Temporal Fusiform Cortex, posterior division', 'Temporal Occipital Fusiform Cortex', 'Occipital Fusiform Gyrus', 'Frontal Opercular Cortex', 'Central Opercular Cortex', 'Parietal Opercular Cortex', 'Planum Polare', "Heschl's Gyrus (includes H1 and H2)", 'Planum Temporale', 'Supracalcarine Cortex', 'Occipital Pole']
labels_index_ho = range(1,len(labels_ho)+1)


labels_msdl = ['L Aud', 'R Aud', 'Striate', 'L DMN', 'Med DMN', 'Front DMN', 'R DMN', 'Occ post', 'Motor', 'R DLPFC', 'R Front pol', 'R Par', 'R Post Temp', 'Basal', 'L Par', 'L DLPFC', 'L Front pol', 'L IPS', 'R IPS', 'L LOC', 'Vis', 'R LOC', 'D ACC', 'V ACC', 'R A Ins', 'L STS', 'R STS', 'L TPJ', 'Broca', 'Sup Front S', 'R TPJ', 'R Pars Op', 'Cereb', 'Dors PCC', 'L Ins', 'Cing', 'R Ins', 'L Ant IPS', 'R Ant IPS']
labels_index_msdl = range(1,len(labels_msdl)+1)


atlas_schaefer = datasets.fetch_atlas_schaefer_2018()
#atlas = datasets.fetch_atlas_destrieux_2009()
atlas_filename = atlas_schaefer["maps"]
labels_schaefer = atlas_schaefer["labels"]
clean_labels_schaefer = [str(x)[2:-1] for x in labels_schaefer]
labels_index_schaefer = range(1,len(labels_schaefer)+1)

atlas_destrieux = datasets.fetch_atlas_destrieux_2009()
labels_destrieux = [x[1] for x in atlas_destrieux.labels]
labels_index_destrieux = range(1,len(labels_destrieux)+1)




df_labelconvert_aal = pd.DataFrame(
    {'index': labels_index_aal,
     'labelname': labels_aal
    })

df_labelconvert_ho = pd.DataFrame(
    {'index': labels_index_ho,
     'labelname': labels_ho
    })

df_labelconvert_msdl = pd.DataFrame(
    {'index': labels_index_msdl,
     'labelname': labels_msdl
    })

df_labelconvert_schaefer = pd.DataFrame(
    {'index': labels_index_schaefer,
     'labelname': clean_labels_schaefer
    })

df_labelconvert_destrieux = pd.DataFrame(
    {'index': labels_index_destrieux,
     'labelname': labels_destrieux
    })


# df_labelconvert_seitzman = pd.DataFrame(
#     {'index': labels_index_seitzman,
#      'labelname': clean_labels_seitzman
#     })


## Define a sub list in a ses

## Define a pandas dataframe of the labels given the chosen atlas
#subject_list = ["01","02","03","04","05"]
base_dir = "/Volumes/LaCie/derivatives"
group = "Patients"
ses_list = ["001","002","003"]
atlas = "destrieux"


if atlas == "harvard_oxford":
	df_labelconvert = df_labelconvert_ho
elif atlas == "aal":
	df_labelconvert = df_labelconvert_aal
elif atlas == "msdl":
	df_labelconvert = df_labelconvert_msdl
elif atlas == "schaefer":
	df_labelconvert = df_labelconvert_schaefer
elif atlas == "seitzman":
	df_labelconvert = df_labelconvert_seitzman
elif atlas == "destrieux":
	df_labelconvert = df_labelconvert_destrieux

all_non_zero_entries = []
result_dir = os.path.join(base_dir,group)

for ses in ses_list:

	subject_list = get_list_fmriprep(base_dir,group,ses)
	#subject_list = ["01"]
	print(subject_list)

	for sub in subject_list:

		connectome_dir = os.path.join(base_dir,group,f'sub-{sub}',f'ses-{ses}',"functional_connectivity")
		print(connectome_dir)
		if not os.path.exists(connectome_dir):
			sys.exit(
				f'Error File not Found: Pipeline has not created connectome matrix for {sub} on {ses}')
		else:
	 		 print(f'RUNNING : {sub} - {ses}')

 		
		print(f'--- Creating ROI file')
		if atlas == "harvard_oxford":
			input_file = os.path.join(connectome_dir, f'sub-{sub}_ses-{ses}_harvard_oxford_correlation_mat.csv')
		elif atlas == "aal":
			input_file = os.path.join(connectome_dir, f'sub-{sub}_ses-{ses}_aal_correlation_mat.csv')
		elif atlas == "msdl":
			input_file = os.path.join(connectome_dir, f'sub-{sub}_ses-{ses}_msdl_correlation_mat.csv')
		elif atlas == "schaefer":
			input_file = os.path.join(connectome_dir, f'sub-{sub}_ses-{ses}_schaefer_correlation_mat.csv')
		elif atlas == "destrieux":
			input_file = os.path.join(connectome_dir, f'sub-{sub}_ses-{ses}_destrieux_correlation_mat.csv')


		df = pd.read_csv(input_file)
			# Create a list to store the non-zero entries for each subject and session
		non_zero_entries = []
		for i in range(df.shape[0]):
			for j in range(i+1, df.shape[1]):
				if df.iloc[i, j] != 0 and i != 0 and j != 0:
					non_zero_entries.append({'subject': sub, 'session': ses, 'i': i, 'j': j, 'PearsonCorr': df.iloc[i, j]})

		# Extend the list of all_non_zero_entries with the current subject and session entri
			
		current_df = pd.DataFrame(non_zero_entries)
		print(current_df.head())

 			
		df_withlabels = pd.merge(current_df, df_labelconvert[['index', 'labelname']], left_on='i', right_on='index', how='left')
		df_withlabels.rename(columns={'labelname': 'ROI1'}, inplace=True)
		df_withlabels.drop(columns='index', inplace=True)
		# # Merge again for 'j' to get ROI2
		df_withlabels = pd.merge(df_withlabels, df_labelconvert[['index', 'labelname']], left_on='j', right_on='index', how='left')
		df_withlabels.rename(columns={'labelname': 'ROI2'}, inplace=True)
		df_withlabels.drop(columns='index', inplace=True)

		output_current_file = connectome_dir + f"/{sub}_{ses}_SC_ROIs.csv" 
		df_withlabels.to_csv(output_current_file,mode = 'w', index=False)

		all_non_zero_entries.append(df_withlabels)

result_df = pd.concat(all_non_zero_entries)

group_dir = base_dir + "/grouped_results"
if not os.path.exists(group_dir):  
    os.makedirs(group_dir) 

# # Write the result DataFrame to a new CSV file
output_file = group_dir + f"/FC_{atlas}.csv" 
print(output_file) 
result_df.to_csv(output_file,mode = 'w', index=False)


		# print(result_dir)

# |Group | Sub-id | Ses-id | Metric | ROI1 | ROI2 |


