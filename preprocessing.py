import pandas as pd
import numpy as np
# import pydicom 
import glob
import os
import sys
import SimpleITK as sitk
from pydicom import dcmread
# from pydicom.data import get_testdata_file
# from pathlib import Path 
import nibabel as nib
# import dicom2nifti.settings as settings
# import dicom2nifti  # to convert DICOM files to the NIftI format


def get_infos(path_to_patients_ctas):
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(path_to_patients_ctas)
    info = dcmread(dicom_names[0])

    return(info)

def get_foldernames_of_phases(path_to_patient_series):
    acquisitions = os.listdir(path_to_patient_series)
    if 'S10' in acquisitions:
        acquisitions.remove('S10')
    if 'S40' in  acquisitions:
        acquisitions.remove('S40')
   
    times = []

    for acquisition in acquisitions:
        path_patient_phases = path_patient_series + acquisition + '\\'
        info = get_infos(path_patient_phases)
        # exec(f'info_{acquisition} = info')
        times.append(int(info.AcquisitionTime))

    sorted_times = sorted(times)

    phase1_folder = path_patient_series + acquisitions[times.index(sorted_times[1])] 
    phase2_folder = path_patient_series + acquisitions[times.index(sorted_times[2])] 
    phase3_folder = path_patient_series + acquisitions[times.index(sorted_times[3])] 

    return(phase1_folder, phase2_folder, phase3_folder)


# def load_data (path_to_patients):

#     reader = sitk.ImageSeriesReader()
#     dicom_names = reader.GetGDCMSeriesFileNames(path_to_patients)
#     reader.SetFileNames(dicom_names)
#     img = reader.Execute()
#     info = dcmread(dicom_names[0])

#     return (img, info)

# def get_phases (path_to_patient_series):

#     # times = []
#     acquisitions = os.listdir(path_to_patient_series)
#     if len(acquisitions) < 6:
#         print("Error! Some of the acquisitions are missing!")
#         return(0,0,0,0,0,0)
#     else:
#         # for acquisition in acquisitions:
#         #     path_patient_phases = path_patient_series + acquisition + '\\'
#         #     info = get_infos(path_patient_phases)
#         #      # exec(f'info_{acquisition} = info')
#         #     times.append(int(info.AcquisitionTime))

#         # sorted_times = sorted(times)
#         dir_phase1, dir_phase2, dir_phase3 = get_foldernames_of_phases(path_patient_series)

#         phase1_img, info_of_phase1 = load_data (dir_phase1 + '\\')
#         phase2_img, info_of_phase2 = load_data (dir_phase2 + '\\')
#         phase3_img, info_of_phase3 = load_data (dir_phase3 + '\\')   
 
#     return(phase1_img, info_of_phase1, phase2_img, info_of_phase2, phase3_img, info_of_phase3)

def convert_to_nifti(input_directory, output_directory, output_file_name):
    cmd = f'dcm2niix -z y -f {output_file_name} -o "{output_directory}" "{input_directory}"'
    os.system(cmd)
    return()

def crop_first_phase_and_save_to_nifti(path_to_phase_1, output_path_to_cropped_phase1):
    image_ct_sitk = sitk.ReadImage(path_to_phase_1, sitk.sitkInt16)
    # image_ct = nib.load(path_to_phase_1)
    origin_indx = [0, 0, 335] #crop area start
    size_bounding_box = [(np.shape(image_ct_sitk)[0]), (np.shape(image_ct_sitk)[1]), (np.shape(image_ct_sitk)[2]-335)]
    cropped_volume = sitk.RegionOfInterest(image_ct_sitk, size_bounding_box, origin_indx)
    sitk.WriteImage(cropped_volume,  output_path_to_cropped_phase1 )
    return()

def resave_from_niigz_to_nii(path_to_image):
    file_paths = glob.glob(path_to_image+'*.nii')
    for file_path in file_paths:
        img = sitk.ReadImage(file_path)
        sitk.WriteImage(img, file_path + '.gz')
    return()

def resave_from_nii_to_niigz(path_to_image):
    file_paths = glob.glob(path_to_image+'thrombus_segmentation.nii.gz')
    for file_path in file_paths:
        img = sitk.ReadImage(file_path)
        sitk.WriteImage(img, path_to_image + 'thrombus_segmentation.nii')
    return()

def conversion_of_matlab_nii(path_matlab, path_python, output_dir):
    matlab_scan = sitk.ReadImage(path_matlab)
    python_scan = sitk.ReadImage(path_python)
    matlab_scan.SetDirection(python_scan.GetDirection())
    sitk.WriteImage(matlab_scan, os.path.join(output_dir, '.nii.gz'))

def register_phases_to_native(path_to_data, patient, parametric_map, output_directory_of_registration):
    registeredPhase1 = sitk.Elastix(sitk.ReadImage(path_to_data + 'nativ' + str(patient) + '.nii.gz'),  \
                           sitk.ReadImage(path_to_data + 'cropped_phase1_' + str(patient) + '.nii.gz'), parametric_map)
    registeredPhase2 = sitk.Elastix(sitk.ReadImage(path_to_data + 'nativ' + str(patient) + '.nii.gz'),  \
                           sitk.ReadImage(path_to_data + 'phase2_' + str(patient) + '.nii.gz'), parametric_map)
    registeredPhase3 = sitk.Elastix(sitk.ReadImage(path_to_data + 'nativ' + str(patient) + '.nii.gz'),  \
                           sitk.ReadImage(path_to_data + 'phase3_' + str(patient) + '.nii.gz'), parametric_map)
    
    sitk.WriteImage(registeredPhase1, output_directory_of_registration + str(patient) + '_phase1_registered.nii.gz')
    sitk.WriteImage(registeredPhase2, output_directory_of_registration + str(patient) + '_phase2_registered.nii.gz')
    sitk.WriteImage(registeredPhase3, output_directory_of_registration + str(patient) + '_phase3_registered.nii.gz')

    return()

def Create_tMIP(path_to_phases, patient):
    phase1 = sitk.GetArrayFromImage(sitk.ReadImage(path_to_phases + str(patient) + '_phase1_registered.nii.gz'))
    phase2 = sitk.GetArrayFromImage(sitk.ReadImage(path_to_phases + str(patient) + '_phase2_registered.nii.gz'))
    phase3 = sitk.GetArrayFromImage(sitk.ReadImage(path_to_phases + str(patient) + '_phase3_registered.nii.gz'))

    tMIP = np.maximum(phase1, phase2, phase3)
    tMIP_img = sitk.GetImageFromArray(tMIP)
    tMIP_img.CopyInformation(sitk.ReadImage(path_to_phases + str(patient) + '_phase1_registered.nii.gz'))
    sitk.WriteImage(tMIP_img, path_to_phases + str(patient) + '_fused.nii.gz')
    return()

def Create_patients_folders_and_subfolders(patient):
    output_dir_nifti_files = 'D:\\Projects\\Stroke\\Nifti_data\\' + str(patient)
    output_dir_registration = "D:\\Projects\\Stroke\\Registered_to_radiologic_view\\" + str(patient) + '\\'
    os.makedirs(output_dir_registration)
    os.makedirs(output_dir_nifti_files)
    return()


#%%

# patients = os.listdir('D:\\Projects\\Stroke\\thrombi_data_FNUSA\\')
patients = ['17', '19', '21', '22', '41', '57', '69', '70', '72']
# list_id = patients.index('62')
# patients = patients[list_id:len(patients)]
path_to_patients_native = "D:\\Projects\\Stroke\\Radiologic_view\\result_batch\\dicoms\\"
path_to_patients_phases = "D:\\Projects\\Stroke\\thrombi_data_FNUSA\\"
# pf = sitk.ReadParameterFile(r"D:\Projects\Stroke\Affine_parameter-file_nii.txt")
pf = sitk.ReadParameterFile(r"D:\Projects\Stroke\Rigid_parameter-file_nii.txt")

# patients_folders = os.listdir(path_to_patients_native)
# input_dir = 'D:\\Projects\\Stroke\\SmartBrained\\31'
# output_dir = 'D:\\Projects\\Stroke\\Nifti_data'
# filename = 'pokus'
# cmd = f'dcm2niix -z y -f {filename} -o "{output_dir}" "{input_dir}"'

# os.system(cmd)


# cmd = 'dcm2niix -z y -f pokus2 -o "D:\\Projects\\Stroke\\Nifti_data" "D:\\Projects\\Stroke\\SmartBrained\\31\\"'
# os.system(cmd)
# dicom2nifti.convert_directory(path_to_patients_native, "D:\\Projects\\Stroke\\Nifti_data\\")

# settings.disable_validate_slice_increment()
# dicom2nifti.dicom_series_to_nifti(path_to_patients_native, "D:\\Projects\\Stroke\\Nifti_data\\")


for patient in patients:
    path_patient_cta = path_to_patients_phases + str(patient) + '\\Export\\DICOM\\'
    patient_series = os.listdir(path_patient_cta)
    path_patient_series = path_patient_cta + patient_series[0] + '\\'

    # resave_from_niigz_to_nii('D:\\Projects\\Stroke\\Thrombi\\' + str(patient) + '\\')
    # resave_from_nii_to_niigz('D:\\Projects\\Stroke\\Thrombi\\'+ str(patient)+'\\')
    input_dir_native = path_to_patients_native + str(patient)
    acquisitions = os.listdir(path_patient_series)
    if len(acquisitions) < 4:
        print("Error! Some of the acquisitions are missing!")
        continue
    else:
        input_dir_phase1, input_dir_phase2, input_dir_phase3 = get_foldernames_of_phases(path_patient_series)
        Create_patients_folders_and_subfolders(patient)
        output_dirs = 'D:\\Projects\\Stroke\\Nifti_data\\' + str(patient)
        output_dir_registration = "D:\\Projects\\Stroke\\Registered_to_radiologic_view\\" + str(patient) + '\\'
        # os.makedirs(output_dir_registration + '\\1\\' + str(patient) + 'registered.nii.gz')
        # os.makedirs(output_dir_registration)
        # os.makedirs(output_dirs)

    #     # save to nifti
        convert_to_nifti(input_dir_native, output_dirs, 'nativ' + str(patient))
        convert_to_nifti(input_dir_phase1, output_dirs, 'phase1_' + str(patient))
        convert_to_nifti(input_dir_phase2, output_dirs, 'phase2_' + str(patient))
        convert_to_nifti(input_dir_phase3, output_dirs, 'phase3_' + str(patient))
        crop_first_phase_and_save_to_nifti(output_dirs+ '\\phase1_' + str(patient) + '.nii.gz', output_dirs + '\\cropped_phase1_' + str(patient) + '.nii.gz')

    # registration
    register_phases_to_native(output_dirs + '\\', patient, pf, output_dir_registration)

    # fusion - tMIP
    
    Create_tMIP(output_dir_registration, patient)

    
    

        

    
# %%
