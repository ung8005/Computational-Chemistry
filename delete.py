#!/usr/bin/env python3

# This code reads a density file in text mode, and parses the water density

import numpy as np
import matplotlib.pyplot as plt 
import scipy 
import os
from scipy import optimize 

# Absolute path to folder in which the files of analysis reside 
folderPath = "/Users/ashleyung/Desktop/Testing/DensityFiles/DDC100"
filePaths = [os.path.join (folderPath, name) for name in os.listdir (folderPath)]
allFiles = []

# Iterate through all files so now allFiles[0] holds the first file loaded, allFiles[1] holds the second, etc..
for path in filePaths:
	with open (path, 'r') as f: 
		densityFile = f.readlines ()
		allFiles.append (densityFile)
# checking 
print (allFiles[1])
print ("secondfile")
print (allFiles[2])
