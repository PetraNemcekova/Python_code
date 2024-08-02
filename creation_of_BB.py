import SimpleITK as sitk
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# path_to_patient_masks = 'D:\\Projects\\Stroke\\Thrombi\\'
# path_to_patient_native = 'D:\\Projects\\Stroke\\Nifti_data\\'
# path_to_patient_phases = 'D:\\Projects\\Stroke\\Registered_to_radiologic_view\\'
# output_directory_of_bb = 'D:\\Projects\\Stroke\\BB\\'
# patient = 112

# native = sitk.GetArrayFromImage(sitk.ReadImage(path_to_patient_native + str(patient) + '\\' + 'nativ' + str(patient) + '.nii.gz'))
# mask = sitk.GetArrayFromImage(sitk.ReadImage(path_to_patient_masks + str(patient) + '\\' + str(patient) + '_KV.nii.gz'))
# phase1 = sitk.GetArrayFromImage(sitk.ReadImage(path_to_patient_phases + str(patient) + '\\' + str(patient) + '_phae1_registered.nii.gz'))

def create_bounding_box(thr_mask, scan):

    thrombus_positions = np.where(thr_mask == 1)
    begining_of_thr =np.array([np.min(thrombus_positions[0]), np.min(thrombus_positions[1]), np.min(thrombus_positions[2])])
    end_of_thr =np.array([np.max(thrombus_positions[0]), np.max(thrombus_positions[1]), np.max(thrombus_positions[2])])

    cropped_dato = scan[begining_of_thr[0]-5:end_of_thr[0]+5, begining_of_thr[1]-5:end_of_thr[1]+5, begining_of_thr[2]-5:end_of_thr[2]+5]
    positions = [begining_of_thr, end_of_thr]

    return(cropped_dato, positions)

path_to_heterogeneity_maps = 'D:\\Projects\\Stroke\\Radiomics\\Heterogeneity_Features_Radiomics\\Voxel_based\\'
path_to_patient_mask = 'D:\\Projects\\Stroke\\Thrombi\\'
path_to_patient_native = 'D:\\Projects\\Stroke\\Nifti_data\\'
path_to_patient_phases = 'D:\\Projects\\Stroke\\Registered_to_radiologic_view\\'
output_directory_of_bb = 'D:\\Projects\\Stroke\\BB\\'

patient = 112
heterogeneity_maps = os.listdir(path_to_heterogeneity_maps + str(patient))
acquisitions = os.listdir(path_to_heterogeneity_maps + str(patient))

for acquisition in acquisitions:
    if acquisition[0:5]=='phase':
        dato = sitk.GetArrayFromImage(sitk.ReadImage(path_to_patient_phases + str(patient) + '\\' + str(patient) + '_' + acquisition + '_registered.nii.gz'))
    else: 
        dato = sitk.GetArrayFromImage(sitk.ReadImage(path_to_patient_native + str(patient) + '\\' + 'nativ' + str(patient) + '.nii.gz'))

    mask = sitk.GetArrayFromImage(sitk.ReadImage(path_to_patient_mask + str(patient) + '\\' + str(patient) + '_KV.nii.gz'))

    BB_dato, BB_positions = create_bounding_box(mask, dato)
    patients_folder_for_BB = output_directory_of_bb + str(patient) + '\\' + acquisition + '\\'
    os.makedirs(patients_folder_for_BB)
    sitk.WriteImage(sitk.GetImageFromArray(BB_dato), patients_folder_for_BB + str(patient) + '_' + acquisition + '_BB.nii.gz')

    heterogeneity_maps = os.listdir(path_to_heterogeneity_maps + str(patient)+ '\\' + acquisition )
    
    for map_name in heterogeneity_maps[1:]:
        heterogeneity_map = sitk.GetArrayFromImage(sitk.ReadImage(path_to_heterogeneity_maps + str(patient)+ '\\' + acquisition + '\\' + map_name))
        BB_heterogeneity_map, BB_positions = create_bounding_box(mask, heterogeneity_map)

        sitk.WriteImage(sitk.GetImageFromArray(BB_heterogeneity_map), patients_folder_for_BB + map_name + '_BB.nii.gz')

# acquisition_type = 'nativ'
# BB_native, BB_positions = create_bounding_box(mask, native)
# sitk.WriteImage(sitk.GetImageFromArray(BB_native), output_directory_of_bb + str(patient) + '\\' + acquisition_type + '_BB.nii.gz')

# acquisition_type = '1st_phase'
# BB_phase1, BB_positions = create_bounding_box(mask, phase1)
# sitk.WriteImage(sitk.GetImageFromArray(BB_phase1), output_directory_of_bb + str(patient) + '\\' + acquisition_type + '_BB.nii.gz')



# BB_patch = patches.Rectangle((BB_positions[0][1], BB_positions[0][2]), BB_positions[1][1]-BB_positions[0][1],  BB_positions[1][2]-BB_positions[0][2], linewidth=1, edgecolor='r', facecolor='none')

# fig, ax = plt.subplots()

# ax.imshow(dato[82,:,:], cmap = 'gray')
# ax.add_patch(BB_patch)
# plt.show()
