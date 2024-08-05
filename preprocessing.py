import pandas as pd
import numpy as np
import glob
import os
import sys
import SimpleITK as sitk
from pydicom import dcmread
import nibabel as nib

def get_infos(path_to_patients_ctas):
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(path_to_patients_ctas)
    info = dcmread(dicom_names[0])
    return info

def get_foldernames_of_phases(path_to_patient_series):
    acquisitions = os.listdir(path_to_patient_series)
    if 'S10' in acquisitions:
        acquisitions.remove('S10')
    if 'S40' in acquisitions:
        acquisitions.remove('S40')

    times = []

    for acquisition in acquisitions:
        path_patient_phases = path_patient_series + acquisition + '\\'
        info = get_infos(path_patient_phases)
        times.append(int(info.AcquisitionTime))

    sorted_times = sorted(times)

    phase1_folder = path_patient_series + \
        acquisitions[times.index(sorted_times[1])]
    phase2_folder = path_patient_series + \
        acquisitions[times.index(sorted_times[2])]
    phase3_folder = path_patient_series + \
        acquisitions[times.index(sorted_times[3])]

    return phase1_folder, phase2_folder, phase3_folder

def convert_to_nifti(input_directory, output_directory, output_file_name):
    cmd = '"D:\\Projects\\Stroke\\Python_code\\dcm2niix.exe"' + ' -z y -f ' + output_file_name + ' -o ' + output_directory + ' ' + input_directory
    os.system('"' + cmd + '"')
    return output_directory + output_file_name


def crop_first_phase_and_save_to_nifti(path_to_phase_1, output_path_to_cropped_phase1):
    image_ct_sitk = sitk.ReadImage(path_to_phase_1, sitk.sitkInt16)
    origin_indx = [0, 0, 335]  # crop area start
    size_bounding_box = [(np.shape(image_ct_sitk)[0]), (np.shape(
        image_ct_sitk)[1]), (np.shape(image_ct_sitk)[2]-335)]
    cropped_volume = sitk.RegionOfInterest(
        image_ct_sitk, size_bounding_box, origin_indx)
    sitk.WriteImage(cropped_volume,  output_path_to_cropped_phase1)
    return ()


def resave_from_niigz_to_nii(path_to_image):
    file_paths = glob.glob(path_to_image+'*.nii')
    for file_path in file_paths:
        img = sitk.ReadImage(file_path)
        sitk.WriteImage(img, file_path + '.gz')
    return ()


def resave_from_nii_to_niigz(path_to_image):
    file_paths = glob.glob(path_to_image+'thrombus_segmentation.nii.gz')
    for file_path in file_paths:
        img = sitk.ReadImage(file_path)
        sitk.WriteImage(img, path_to_image + 'thrombus_segmentation.nii')
    return ()


def conversion_of_matlab_nii(path_matlab, path_python, output_dir):
    matlab_scan = sitk.ReadImage(path_matlab)
    python_scan = sitk.ReadImage(path_python)
    matlab_scan.SetDirection(python_scan.GetDirection())
    sitk.WriteImage(matlab_scan, os.path.join(output_dir, '.nii.gz'))
    return()

def register_phases_to_native(path_to_fixed, path_to_moving, moving_label, patient, parametric_map, output_directory_of_registration):

    resultImage = sitk.Elastix(sitk.ReadImage(path_to_fixed + '.nii.gz'),
                               sitk.ReadImage(path_to_moving),
                               parametric_map)
    
    sitk.WriteImage(resultImage, output_directory_of_registration +
                    str(patient) + '_' + moving_label + '_registered.nii.gz')

    return ()


def Create_tMIP(path_to_phases, patient):
    phase1 = sitk.GetArrayFromImage(sitk.ReadImage(
        path_to_phases + str(patient) + '_phase1_registered.nii.gz'))
    phase2 = sitk.GetArrayFromImage(sitk.ReadImage(
        path_to_phases + str(patient) + '_phase2_registered.nii.gz'))
    phase3 = sitk.GetArrayFromImage(sitk.ReadImage(
        path_to_phases + str(patient) + '_phase3_registered.nii.gz'))

    tMIP = np.maximum(phase1, phase2, phase3)
    tMIP_img = sitk.GetImageFromArray(tMIP)
    tMIP_img.CopyInformation(sitk.ReadImage(
        path_to_phases + str(patient) + '_phase1_registered.nii.gz'))
    sitk.WriteImage(tMIP_img, path_to_phases + str(patient) + '_fused.nii.gz')
    return ()


def Create_patients_folders_and_subfolders(output_path, patient):

    try:
        output_dir_nifti_files = output_path + str(patient) + '\\'
        output_dir_registration = output_path + \
            'Registered\\' + str(patient) + '\\'
        folders_to_create = [output_dir_nifti_files, output_dir_registration]
        for folder_to_create in folders_to_create:
            os.makedirs(folder_to_create)
    except:
        pass

    return output_dir_nifti_files, output_dir_registration


def Rotate_to_radiologic_view(input_folder_path, output_folder_path, SmartBrain_path):
    cmd = SmartBrain_path + ' ' + input_folder_path + ' ' + output_folder_path
    os.system('"' + cmd + '"')

# %%
if __name__ == "__main__":

    input_folder_path_fixed_images = 'D:\\Projects\\Stroke\\Data\\pokus\\'
    input_folder_path_moving_images = 'D:\\Projects\\Stroke\\Data\\pokus\\'
    output_folder_path = 'D:\\Projects\\Stroke\\Data\\output_folder_pokus\\'
    path_to_SmartBrain = '"C:\\Program Files\\Brno University of Technology\\SmartBrainBatch_cmd\\application\\SmartBrainBatch_cmd.exe"'

    patients = os.listdir(input_folder_path_fixed_images)
    pf_affine = sitk.ReadParameterFile(
        r"D:\Projects\Stroke\Documents\Affine_parameter-file_nii.txt")
    pf_rigid = sitk.ReadParameterFile(
        r"D:\Projects\Stroke\Documents\Rigid_parameter-file_nii_new.txt")
    transformation_types = ['affine', 'rigid']

    Rotate_to_radiologic_view(input_folder_path = input_folder_path_fixed_images, output_folder_path = output_folder_path, SmartBrain_path = path_to_SmartBrain)

    for patient in patients:
        path_patient_cta = input_folder_path_moving_images + \
            str(patient) + '\\Export\\DICOM\\'
        patient_series = os.listdir(path_patient_cta)
        path_patient_series = path_patient_cta + patient_series[0] + '\\'

        # resave_from_niigz_to_nii('D:\\Projects\\Stroke\\Thrombi\\' + str(patient) + '\\')
        # resave_from_nii_to_niigz('D:\\Projects\\Stroke\\Thrombi\\'+ str(patient)+'\\')
        input_dir_native = input_folder_path_fixed_images + str(patient) + '\\'
        acquisitions = os.listdir(path_patient_series)
        if len(acquisitions) < 4:
            print("Error! Some of the acquisitions are missing!")
            continue
        else:
            input_dir_phase1, input_dir_phase2, input_dir_phase3 = get_foldernames_of_phases(
                path_patient_series)
            output_path_nifti_files, output_path_registration = Create_patients_folders_and_subfolders(
                output_folder_path, patient)

        #     # save to nifti
            directions_to_acquisitions = {
                'native': input_dir_native + 'S20\\', 'phase1': input_dir_phase1, 'phase2': input_dir_phase2, 'phase3': input_dir_phase3}

            for i, direction_to_acquisition_key in enumerate(directions_to_acquisitions.keys()):
                nifti_path_acquisition_patient = convert_to_nifti(input_directory=directions_to_acquisitions[direction_to_acquisition_key],
                                                                  output_directory=output_path_nifti_files, output_file_name=direction_to_acquisition_key + '_' + str(patient))
            crop_first_phase_and_save_to_nifti(output_path_nifti_files + 'phase1_' + str(patient) + '.nii.gz', output_path_nifti_files + 'cropped_phase1_' + str(patient) + '.nii.gz')
            
            for i, direction_to_acquisition_key in enumerate(directions_to_acquisitions.keys()):
            # registration
                if i > 0:
                    for transformation_type in transformation_types:
                        if transformation_type == 'affine':
                            register_phases_to_native(path_to_fixed=input_dir_native, path_to_moving=nifti_path_acquisition_patient, moving_label=direction_to_acquisition_key, patient=patient, parametric_map = pf_affine, output_directory_of_registration = output_path_registration)
                        elif transformation_type == 'rigid':
                            register_phases_to_native(path_to_fixed=input_dir_native, path_to_moving=nifti_path_acquisition_patient, moving_label=direction_to_acquisition_key, patient=patient, parametric_map = pf_rigid, output_directory_of_registration = output_path_registration)


        # fusion - tMIP

        Create_tMIP(output_path_registration, patient)


# %%
