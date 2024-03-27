# import os
# from python_on_whales import docker

# freesurfer_bin_path = '/Users/francoisramon/freesurfer'
# FS_license_file = os.path.join(freesurfer_bin_path, 'license.txt')

# input_synb0 = '/Users/francoisramon/Desktop/Thèse/dMRI_pipeline/test_synthb0/inputs'
# output_synb0 = '/Users/francoisramon/Desktop/Thèse/dMRI_pipeline/test_synthb0/outputs'
           
# output_generator = docker.run(
#     "nipreps/fmriprep",
#     ["--user", str(os.getuid()) + ":" + str(os.getgid())],
#     volumes=[(input_synb0, "/INPUTS"), (output_synb0, "/OUTPUTS"), (FS_license_file, "/extra/freesurfer/license.txt")],
#     remove=True, stream=True,
#     )
# for stream_type, stream_content in output_generator:
#     print(f"Stream type: {stream_type}, stream content: {stream_content}")



# fmriprep-docker /Volumes/LaCie/nifti3/nifti3/Patients /Volumes/LaCie/outfmri participant --participant-label 01
# # RUNNING: docker run --rm -it -v /path/to/data/dir:/data:ro \
# #     -v /path/to_output/dir:/out poldracklab/fmriprep:1.0.0 \
# #     /data /out participant


# docker run --rm -it -v /Volumes/LaCie/nifti3/nifti3/Patients:/data:ro -v /Volumes/LaCie/outfmri:/out -u $UID nipreps/fmriprep:latest /data /out participant --participant-label 

fmriprep-docker /Volumes/LaCie/nifti3/nifti3/Patients /Volumes/LaCie/derivatives \
    participant \
    --participant-label 01 \
    --skip-bids-validation \
    --fs-license-file /opt/freesurfer/license.txt \
    --fs-no-reconall \


#docker run --rm -e DOCKER_VERSION_8395080871=24.0.7 -it -u 501 -v /opt/freesurfer/license.txt:/opt/freesurfer/license.txt:ro -v /Volumes/LaCie/nifti3/nifti3/Patients:/data:ro -v /Volumes/LaCie/derivatives:/out nipreps/fmriprep:23.2.0 /data /out participant --participant-label 01 --fs-no-reconall