import os

def scan_for_files(filename, start_directory, files):
    """
    scan_for_files is a recursive function that scans a given directory
    (start_directory) for all files with the name 'filename', and returns
    the absolute path of each file it finds.

    """
    
    file_and_dir_names = os.listdir(start_directory)

    if (start_directory[-1] != '/'):
        start_directory += '/'

    for i in file_and_dir_names:
        # Look for matching files in current directory
        if ((i == filename) and os.path.isfile(start_directory + i)):
            files.append(os.path.abspath(start_directory + filename))

        # Look in sub-directories
        if os.path.isdir(start_directory + i):
            files = scan_for_files(filename, os.path.abspath(start_directory + i), files)

    return files
