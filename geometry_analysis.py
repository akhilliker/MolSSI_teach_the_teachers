"""
geometry_analysis.py
This module contains the geometry analysis project.
To run in terminal: python geometry_analysis.py {file_to_analyze.xyz}
"""

# Best practice: do all imports at the top
import os.path
import numpy
import sys

def open_xyz(filepath):
    """
    Function to open an xyz file, pull the atom names and the x, y, z coordinates
    Input:
        Filename parameter is required and contains the file path for the xyz file
        Assumes that the .xyz file has two header rows, which will be skipped
        Assumes that the data array includes columns in the following order:
            Index 0: atom name
            Index 1: x coordinate
            Index 2: y coordinate
            Index 3: z coordinate
        Coordinates are in Angstroms
    """
    xyz_file = numpy.genfromtxt(fname = filepath, skip_header = 2, dtype= 'unicode')
    symbols = xyz_file[:, 0] # all rows, column 1 to get name
    coordinates = xyz_file[:, 1:] # all rows, column 2 and on to get numbers
    coordinates = coordinates.astype(numpy.float) # change data type from string to float
    # if len(coordinates) != 3:
    #     raise ValueError('Input .xyz file requires at least three inputs for x, y, z coordinates')
    # if len(coordinates) > 3:
    #     print("Note that your xyz file has more than 3 columns. Program assumes 1st 3 columns are x, y, z coordinates.  All later columns are ignored")
    return symbols, coordinates

def calculate_distance(atom1_coord, atom2_coord):
    """
    The triple quotes allows multi line strings and can be used anywhere.
    This is a doc string. This is where you write documentation for how your function works.
    There are standards for how to do this, depending on the library. MolSSI uses numpy standards.
    By putting the doc string here before your code, it will be printed when the help(your_function_name) function is called.
    Inputs: x, y, z coordinates (in that order) for atom1 and atom2.
    Returns the bond distance between those two atoms in angstroms.
    """
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_distance = numpy.sqrt(x_distance**2 + y_distance**2 + z_distance**2)
    return bond_distance

def bond_check(bond_distance, minimum_length = 0, maximum_length = 1.5):
    """
    Input: Three numeric parameters:
        1) a distance between two atoms, in angstroms (required parameter)
        2) minimum acceptable bond length (optional parameter; default is 0)
        3) maximum acceptable bond length (optional parameter; default is 1.5)
    Evaluates whether the bond is between min and max angtroms
    Output: True if bond is between 0 and 1.5 angstroms; else, False
    """
    if bond_distance > minimum_length and bond_distance <= maximum_length:
        return True
    elif bond_distance < 0:
        raise ValueError(f'Error: a bond length calculation is negative ({bond_distance}) and must be positive')
    else:
        return False

# makes all the code below an importable function that can be called in other files
if __name__ == "__main__":
    #file_path = os.path.join('data', 'water.xyz')
    # making a flexible input for file names
    # sys argv makes anything typed after python in the command line accessible
    # putting [1] means the file name will follow the script file
    print(f"Running {sys.argv[0]} to analyze {sys.argv[1]}")

    if len(sys.argv) < 2:
        # will raise an error that stops the program and prints something
        raise NameError('Incorrect input! Please specify an xzy file to be analyzed.')

    file_path = sys.argv[1]

    # To check errors quickly:
    # bond_check(-5)

    # Call the functions
    datafile = open('H2Obondlengths4.txt','w+')

    symbols, coord = open_xyz(file_path)
    num_atoms =  len(symbols)
    for i in range(0, num_atoms):
        for j in range(0, num_atoms):
            if i < j:
                D = calculate_distance(coord[i], coord[j])
                if bond_check(D) is True:
                    print(f'{symbols[i]} to {symbols[j]} bond distance: {D: .3f} angstroms.')
                    datafile.write(f'{symbols[i]} to {symbols[j]} bond distance: {D: .3f} angstroms. \n')

    datafile.close()
