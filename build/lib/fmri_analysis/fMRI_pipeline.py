import os
import subprocess
import concurrent.futures

from fmri_analysis.fmri_utils import *
from fmri_analysis.fmri_connectivity import *


### Choose the base directory where is the bids database

source_dir = "/Volumes/LaCie/nifti3"
group = "HealthyVolunteers"
output_dir = f"/Volumes/LaCie/derivatives"
### Subject select

subject_list = [x for x in os.listdir(os.path.join(source_dir,group)) if "sub" in x]
#print(subject_list)

#### Task select

clean_hidden_files = False
run_preprocessing = False
createMatrixes = True

### Parameters selection

#atlas = "harvard_oxford"
atlas = "aal"
#atlas = "msdl"
#atlas = "seitzman"
#atlas = "harvard_oxford"
#atlas = "craddock"

kindlist = ["covariance", "correlation", "partial correlation", "tangent", "precision"]
kind = "correlation"

## OPTIONAL Dataset prepare : delete hidden json files

if clean_hidden_files:	
	for sub in subject_list:
		ses001_anat_hidden_files = [f"._{sub}_ses-001_T1w.nii.gz",f"._{sub}_ses-001_T1w.json"]
		ses001_func_hidden_files = [f"._{sub}_ses-001_task-rest_bold.json",f"._{sub}_ses-001_task-rest_bold.nii.gz"]
		ses001_hidden_d_files = [f"._{sub}_participants.tsv",f"._participants.json","._scans.json"]

		ses002_hidden_files = [f"._{sub}_ses-002_task-rest_bold.json",f"._{sub}_ses-002_T1w.nii.gz"]
		ses003_hidden_files = [f"._{sub}_ses-003_task-rest_bold.json",f"._{sub}_ses-002_T1w.nii.gz"]

		for hidd in ses001_func_hidden_files:
			#print(hidd)
			if os.path.isfile(os.path.join(source_dir,group,sub,"ses-001","func",hidd)):
				os.remove(os.path.join(source_dir,group,sub,"ses-001","func",hidd))
		

		for hidd in ses001_anat_hidden_files:
			#print(hidd)
			if os.path.isfile(os.path.join(source_dir,group,sub,"ses-001","anat",hidd)):
				os.remove(os.path.join(source_dir,group,sub,"ses-001","anat",hidd))

		for hidd in ses001_hidden_d_files:
			#print(hidd)
			if os.path.isfile(os.path.join(source_dir,group,hidd)):
				os.remove(os.path.join(source_dir,group,hidd))

		for hidd in ses002_hidden_files:
			if os.path.isfile(os.path.join(source_dir,group,sub,"ses-002","anat",hidd)):
				os.remove(os.path.join(source_dir,group,sub,"ses-002","anat",hidd))
	
		for hidd in ses003_hidden_files:
			if os.path.isfile(os.path.join(source_dir,group,sub,"ses-003","anat",hidd)):
				os.remove(os.path.join(source_dir,group,sub,"ses-003","anat",hidd))
			

########################################################
########                fmriprep              ##########
########################################################
## Je pense que comme on utilise un pipe nipype, il y a le systeme de caching automatique


subject_list_CLI = [x[4:] for x in subject_list]

#subject_list_CLI = " ".join([item for item in subject_list_CLI])
#print(subject_list_CLI)

## Def a function that gets all patients having had at least one ses of func

if run_preprocessing: 
	#subject_list_CLI = ['03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44']
	#subject_list_CLI = " ".join([item for item in subject_list_CLI])
	#subject_list_CLI = ["02"]

	for sub in subject_list_CLI:
		print(sub)
		command = f"docker run --rm -e DOCKER_VERSION_8395080871=24.0.7 -it -v /opt/freesurfer/license.txt:/opt/freesurfer/license.txt:ro -v {source_dir}/{group}:/data:ro -v {output_dir}/{group}:/out nipreps/fmriprep:latest /data /out participant --participant-label {sub} --skip-bids-validation --fs-no-reconall"
		print(command)
		subprocess.run(command, shell = True)



## Create Matrixes with atlas choice

if createMatrixes:
	for ses in ["001"]:
		sub_in_ses = get_list_fmriprep(output_dir,group,ses)
		print(sub_in_ses)
		for sub in sub_in_ses:
			sub_id = "sub-" + sub
			ses_id = "ses-" + ses

			### Create a "functional connectivity repository"
			os.makedirs(os.path.join(output_dir,group,sub_id,ses_id,"functional_connectivity"),exist_ok = True)
			#if not os.path.isfile(os.path.join(output_dir,group,sub_id,ses_id,"functional_connectivity",f"{sub}_{ses}_{atlas}_correlation_mat.csv")):
			print(f"running on subject {sub} session {ses}")
			compute_correlation_matrix(output_dir,group,sub_id,ses_id,atlas,kind)

### Create one folder "connectivity matrix" if the folder derivatives/sub-/ses-/func exists AND is not empty

### Run a nilearn function that take into argument the derivatives folder. sub id ses id, atlas name. confounds ?






# docker run --privileged --rm -e DOCKER_VERSION_8395080871=24.0.7 -it -v /opt/freesurfer/license.txt:/opt/freesurfer/license.txt:ro -v /Volumes/LaCie/nifti3/nifti3/Patients:/data:ro -v /Volumes/LaCie/derivatives/HealthyVolunteers:/out nipreps/fmriprep:latest /data /out participant --participant-label 02 --skip-bids-validation --fs-no-reconall

# docker run -v /opt/freesurfer/license.txt:/opt/freesurfer/license.txt:ro -v /Volumes/LaCie/nifti3/nifti3/Patients:/data:ro -v /Volumes/LaCie/derivatives/HealthyVolunteers:/out nipreps/fmriprep:latest /data /out participant --participant-label 02 --skip-bids-validation --fs-no-reconall


















## Cache avec nilearn???


