#!/usr/bin/env python3

import numpy
import sys
import os

# This is a very simple script. It reads a MNI xfm file, takes the 3x3 matrix 
# encoding the rotations, scaling etc. and calculates the determinant.
# The main objective is to find out whether there is a reflection present
# in the transformation. If the determinant is negative then there is indeed
# a reflection present in the transformation.

# we only expect a single input file:
if len(sys.argv) != 2:
    sys.exit("\nError: please specify one input xfm file.\n")

# next, verify that it is an xfm file in a fairly crude manner
filename, extension = os.path.splitext(sys.argv[1])
if extension != ".xfm":
    sys.exit("\nError: please specify an MNI transformation file (.xfm) as input\n")

xfmfile = open(sys.argv[1], "r")
xfmfilelines = xfmfile.readlines()

# one more check to verify that this is an XFM file:
if xfmfilelines[0].rstrip() != "MNI Transform File":
    sys.exit("\nError: the input file does not appear to be an xfm file. Does not start with the line \"MNI Transform File\".\n")

# get index of the line where the matrix starts:
found = False
index = 0
while not found:
    if "Linear_Transform =" in xfmfilelines[index].rstrip():
        # the next line is the one we want!
        index = index + 1
        found = True
    else:
        # go to the next line
        index = index + 1

matrix_line_1 = (xfmfilelines[index].strip()).rstrip()
matrix_line_2 = (xfmfilelines[index + 1].strip()).rstrip()
matrix_line_3 = (xfmfilelines[index + 2].strip()).rstrip()

matrix_line_1_parts = matrix_line_1.split(" ")
matrix_line_2_parts = matrix_line_2.split(" ")
matrix_line_3_parts = matrix_line_3.split(" ")


matrix = [[float(matrix_line_1_parts[0]), float(matrix_line_1_parts[1]), float(matrix_line_1_parts[2])], [float(matrix_line_2_parts[0]), float(matrix_line_2_parts[1]), float(matrix_line_2_parts[2])], [float(matrix_line_3_parts[0]), float(matrix_line_3_parts[1]), float(matrix_line_3_parts[2])]]

determinant = numpy.linalg.det(matrix)

print(determinant, " is the determinant for transform: ", sys.argv[1])

sys.exit()

