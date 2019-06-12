import numpy as np
import nibabel as nib
import os
from . import file_finder

def create_nifti(path_2dseq, path_visu_pars, path_method, flip_rpe, save_path):
	# Read in relevant header info
	visu_pars = read_bruker_header_file(path_visu_pars)
	method = read_bruker_header_file(path_method)
	offset = get_offset(visu_pars)
	resolution = get_resolution(visu_pars)
	matrix = get_matrix(visu_pars)

	if (flip_rpe and ('RPE' in visu_pars['VisuAcquisitionProtocol'][0])):
		flip_rpe = True
	else:
		flip_rpe = False

	affine = get_affine(visu_pars)

	# Create bvals and bvecs if applicable
	if 'PVM_DwEffBval' in method:
		bval = np.asfarray(method['PVM_DwEffBval'], dtype='float')
		bval[bval < 50] = 0.0
		bvals = np.zeros((1,len(bval)))
		bvals[0,:] = bval

	if 'PVM_DwGradVec' in method:
		bvec = np.asfarray(method['PVM_DwGradVec'], dtype='float')
		for i in range(bvec.shape[0]):
			if np.linalg.norm(bvec[i,:]) == 0:
				bvec[i,:] = 0
			else:
				bvec[i,:] /= np.sqrt(np.sum(bvec[i,:] ** 2))

		bvecs = bvec.T

	# Read in and shape the image filename
	word_size = visu_pars['VisuCoreWordType']

	if word_size == '_16BIT_SGN_INT':
		word_size = np.int16

	elif word_size == '_32BIT_SGN_INT':
		word_size = np.int32

	img = np.fromfile(path_2dseq, dtype=word_size)
	img.byteswap(False)
	img = np.reshape(img, matrix, order='F')

	if flip_rpe:
		# Remove all non-bzero images
		b0_count = 0
		for i in range(bvals.shape[1]):
			if bvals[0,i] < 50:
				b0_count += 1

		img_new = np.zeros((img.shape[0], img.shape[1], img.shape[2], b0_count))
		b0_count = 0
		for i in range(bvals.shape[1]):
			if bvals[0,i] < 50:
				img_new[:,:,:,b0_count] = img[:,:,:,i]
				b0_count += 1

		img = img_new

		img = np.flip(img, 1)

	# Create the NIFTI
	nii = nib.Nifti1Image(img, affine)
	nii.qoffset_x = offset[0]
	nii.qoffset_y = offset[1]
	nii.qoffset_z = offset[2]

	nii.affine[0:3,3] = [offset[0], offset[1], offset[2]]

	# Save the NIFTI
	study_number = path_2dseq.split('/')[-4]
	series = path_2dseq.split('/')[-2]
	if os.path.isdir(save_path + visu_pars['VisuStudyId'][0]):
		save_filepath = (save_path + visu_pars['VisuStudyId'][0] + '/' +
						visu_pars['VisuAcquisitionProtocol'][0] + '_' +
						study_number + '_' + series)

		nib.save(nii, save_filepath)

		if(('bval' in locals()) and (flip_rpe == False)):
			np.savetxt(save_path + visu_pars['VisuStudyId'][0] + '/' +
					   visu_pars['VisuAcquisitionProtocol'][0] + '_' +
					   study_number +  '_' + series + '.bval', bvals, fmt='%f')

		if(('bvec' in locals()) and (flip_rpe == False)):
			np.savetxt(save_path + visu_pars['VisuStudyId'][0] + '/' +
			 	       visu_pars['VisuAcquisitionProtocol'][0] + '_' +
					   study_number +  '_' + series + '.bvec', bvecs, fmt='%f')
	else:
		os.makedirs(save_path + visu_pars['VisuStudyId'][0])
		save_filepath = (save_path + visu_pars['VisuStudyId'][0] + '/' +
						visu_pars['VisuAcquisitionProtocol'][0] + '_' + study_number)
		nib.save(nii, save_filepath)

		# Only save bvals and bvecs if it is not an RPE image
		if(('bval' in locals()) and (flip_rpe == False)):
			np.savetxt(save_path + visu_pars['VisuStudyId'][0] + '/' +
					   visu_pars['VisuAcquisitionProtocol'][0] + '_' +
					   study_number  + '_' + series + '.bval', bvals, fmt='%f')

		if (('bvec' in locals()) and (flip_rpe == False)):
			np.savetxt(save_path + visu_pars['VisuStudyId'][0] + '/' +
					   visu_pars['VisuAcquisitionProtocol'][0] + '_' +
					   study_number + '_' + series + '.bvec', bvecs, fmt='%f')

def read_bruker_header_file(filename):

	f = open(filename, 'r')

	line = f.readline()
	if (line[0:7] == "##TITLE"):
		header = {}

		while line:
			skip_line = 0

			if (line[0:3] == "##$"):
				attribute = line[3:].split('=')

				# If parameter value is just a number
				if (attribute[1][0] != '('):
					header[attribute[0]] = attribute[1][:-1]

				# If parameter value is a string, vector, or array
				if (attribute[1][0] == '('):
					attribute[1] = attribute[1].replace("(", "")
					attribute[1] = attribute[1].replace(")", "")
					attribute[1] = attribute[1].replace(" ", "")

					field_size = attribute[1][:-1].split(",")

					for i in range(len(field_size)):
						try:
							field_size[i] = int(field_size[i])
						except:
							skip_line = 1

					if (skip_line == 0):
						field_size = tuple(field_size)
						header[attribute[0]] = []

						# Save spot in file
						prev_spot = f.tell()
						next_char = f.read(1)
						info = ''

						# Read file until next attribute is reached
						while (next_char != '#') and (next_char != '$'):
							info += next_char

							next_char = f.read(1)
							if (next_char != '#') and (next_char != '$'):
								prev_spot = f.tell()

						# Split Info Field based on type of information
						if len(info) > 0:
							# Handle Strings
							if info[0] == '<':
								info = info.replace('\n','').split('> <')

								for i in range(len(info)):
									info[i] = info[i].replace('<', '').replace('>', '')

							# Handle Cells
							elif info[0] == '(':
								info = info.replace('\n', '').split(') (')

								for i in range(len(info)):
									info[i] = info[i].replace('(', '').replace(')', '').replace('<','').replace('>','')
									info[i] = info[i].split(', ')

							# Handle Arrays and Vectors
							elif (info[0] != '<') and (info[0] != '('):
								info = info.replace('\n', '').split()

								# Fill in repeats (Bruker added @#*(#) to save space)
								temp_info = []
								for i in range(len(info)):
									if info[i][0] == '@':
										repeats = info[i].split('@')[1].split('*')

										for j in range(int(repeats[0])):
											temp_info.append(repeats[1][1:-1])

									elif info[i][0] != '@':
										temp_info.append(info[i])

								info = temp_info
								info = np.reshape(info, field_size)

								# Get rid of duplicate fields
								if len(info) > 10:
									is_repeated = True

									for i in range(len(info)):
										if (np.any(info[i] != info[0])):
											is_repeated = False

									if is_repeated:
										info = info[0]

							header[attribute[0]] = info
							# Rewind to beginning of next attribute and/or comment
							f.seek(prev_spot)

			line = f.readline()

	f.close()

	return header

def get_affine(visu_pars):

	if (len(visu_pars['VisuCoreOrientation'].shape) > 1):
		affine_orientation = visu_pars['VisuCoreOrientation'][0]
	else:
		affine_orientation = visu_pars['VisuCoreOrientation']

	affine_orientation = np.reshape(affine_orientation, (3,3))
	for i in range(3):
		for j in range(3):
			affine_orientation[i,j] = float(affine_orientation[i,j])

	resolution = get_resolution(visu_pars)
	offset = get_offset(visu_pars)

	affine = np.eye(4, dtype=np.float32)
	affine[0:3,0:3] = affine_orientation
	affine[0:3,3] = offset

	affine = np.linalg.inv(affine)

	affine[0:3,0:3] = affine[0:3,0:3].dot(np.diag(resolution))

	return affine

def get_offset(visu_pars):

	offset = np.zeros(3)
	offset_values = np.zeros(visu_pars['VisuCorePosition'].shape)

	# For 2D
	if visu_pars['VisuCoreDim'] == '2':
		for i in range(visu_pars['VisuCorePosition'].shape[0]):

			offset_values[i,0] = float(visu_pars['VisuCorePosition'][i][0])
			offset_values[i,1] = float(visu_pars['VisuCorePosition'][i][1])
			offset_values[i,2] = float(visu_pars['VisuCorePosition'][i][2])

		offset[0] = -np.amin(offset_values[:,0])
		offset[1] = -np.amin(offset_values[:,1]) / 2.0
		offset[2] = np.amin(offset_values[:,2])

	# For 3D
	if visu_pars['VisuCoreDim'] == '3':
		offset[0] = -float(visu_pars['VisuCorePosition'][0][0])
		offset[1] = float(visu_pars['VisuCorePosition'][0][1])
		offset[2] = float(visu_pars['VisuCorePosition'][0][2])

	return offset

def get_resolution(visu_pars):

	resolution = np.zeros(3)

	# 2D case
	if visu_pars['VisuCoreDim'] == '2':
		resolution[0] = float(visu_pars['VisuCoreExtent'][0]) / float(visu_pars['VisuCoreSize'][0])
		resolution[1] = float(visu_pars['VisuCoreExtent'][1]) / float(visu_pars['VisuCoreSize'][1])

		depth = np.zeros(3)
		for i in range(len(depth)):
			depth[i] = (float(visu_pars['VisuCorePosition'][-1][i]) - float(visu_pars['VisuCorePosition'][0][i]))

		depth = np.linalg.norm(depth)

		for i in range(len(visu_pars['VisuFGOrderDesc'])):
			if visu_pars['VisuFGOrderDesc'][i][1] == 'FG_SLICE':
				num_slices = float(visu_pars['VisuFGOrderDesc'][i][0])

		resolution[2] = np.minimum(depth / (num_slices - 1), float(visu_pars['VisuCoreFrameThickness'][0]))

	# 3D case
	if visu_pars['VisuCoreDim'] == '3':
		resolution[0] = float(visu_pars['VisuCoreExtent'][0]) / float(visu_pars['VisuCoreSize'][0])
		resolution[1] = float(visu_pars['VisuCoreExtent'][1]) / float(visu_pars['VisuCoreSize'][1])
		resolution[2] = float(visu_pars['VisuCoreExtent'][2]) / float(visu_pars['VisuCoreSize'][2])

	return resolution

def get_matrix(visu_pars):

	matrix = []

	# 2D case
	if visu_pars['VisuCoreDim'] == '2':

		matrix.append(int(visu_pars['VisuCoreSize'][0]))
		matrix.append(int(visu_pars['VisuCoreSize'][1]))

		for i in range(len(visu_pars['VisuFGOrderDesc'])):
			if visu_pars['VisuFGOrderDesc'][i][1] == 'FG_SLICE':
				matrix.append(int(visu_pars['VisuFGOrderDesc'][i][0]))

			if visu_pars['VisuFGOrderDesc'][i][1] != 'FG_SLICE':
				matrix.append(int(visu_pars['VisuFGOrderDesc'][i][0]))

	# 3D case
	if visu_pars['VisuCoreDim'] == '3':

		matrix.append(int(visu_pars['VisuCoreSize'][0]))
		matrix.append(int(visu_pars['VisuCoreSize'][1]))
		matrix.append(int(visu_pars['VisuCoreSize'][2]))

		if 'VisuFGOrderDesc' in visu_pars.keys():
			matrix.append(int(visu_pars['VisuFGOrderDesc'][0][0]))

	return tuple([int(i) for i in matrix])

def convert_to_nifti(recon_data, visu_pars, out_path):

	# Get Word Size
	word_size = visu_pars['VisuCoreWordType']

	if word_size == '_16BIT_SGN_INT':
		word_size = np.int16

	elif word_size == '_32BIT_SGN_INT':
		word_size = np.int32

	# Get Matrix Size
	matrix = get_matrix(visu_pars)

	# Get affine transformation
	affine = get_affine(visu_pars)

	img_data = np.fromfile(recon_data, word_size)
	img_data = np.reshape(img_data, matrix, order='F')

	# Convert to Nifti and Save
	img = nib.Nifti2Image(img_data, affine)
	nib.save(img, out_path)
