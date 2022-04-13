import numpy as np
import matplotlib.pyplot as plt 
import scipy 
import os
from scipy import optimize 

# Absolute path to folder in which the files reside 
folderPath= "/Users/ashleyung/Desktop/Testing/DensityFiles/"
# Change the directory 
os.chdir (folderPath)
# Read text file
def readDensityData (file_path):
	with open (file_path, 'r') as f:
		print (f.read ())
# Iterate through all files 
for file in os.listdir ():
	if file.endswith (".txt"):
		file_path = f"{path}\{file}"
		read_text_file (file_path)
