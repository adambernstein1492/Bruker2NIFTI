import numpy as np
import nibabel as nib
import os

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
