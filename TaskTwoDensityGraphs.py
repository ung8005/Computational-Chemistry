#!/usr/bin/env python3

# This code reads in multiple density file in text mode, and parses the water density

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
	with open (path, 'r+') as density_File: 
		with open ("outFile.txt", "w+") as outFile: 
			# Read first line with column headings 
			header_Line = density_File.readline ()
			allFiles.append (density_File)
			#print (allFiles[0])
			# Put coloumn headings into a list for future indexing of the data 
			data_Items=header_Line[2:].split ()
			#print (data_Items)
			# Load the density data into a 2D array 
			data_Array = np.loadtxt (density_File)
			# Extract x and y arrays to plot the water density 
			depth_plot = data_Array[:,data_Items.index ('nm')]
			water_density_plot = data_Array[:, data_Items.index ('H2O')]
			# Extract subset water density array to fit to the interface width function 
			# First define the depth in nm for start and end, then determine the array indices for start and stop of data 
			analysis_start=-1.0
			analysis_start_index=np.where(data_Array[:,data_Items.index ('nm'):] == analysis_start)[0][0]
			analysis_end = 1.0
			analysis_end_index=np.where(data_Array[:,data_Items.index ('nm'):] == analysis_end)[0][0]+1
			# Extract the fitting data into x and y arrays
			Fit_Depth_Array=data_Array[analysis_start_index:analysis_end_index,data_Items.index ('nm')]
			Fit_Water_Interface = data_Array[analysis_start_index:analysis_end_index,data_Items.index ('H2O')]
			# Fitting the water density at the interface, this will give the width of the interface 
			def hyperbolicTanFunc (m, a, b, c):
				return (a/2)*(1-np.tanh ((m-b)/(c)))
			maxHeight = 1
			inflection = 0
			width = 0.4
			popt, cov = scipy.optimize.curve_fit (hyperbolicTanFunc, Fit_Depth_Array, Fit_Water_Interface)
			perr = np.sqrt(np.diag (cov))
			print ("Error:")
			print (perr)
			maxHeight, inflection, width = popt 
			# Verifying the hyperbolic function parameters 
			print ("(Maximum height, inflection point, width)")
			print (maxHeight, inflection, width)
			yNew = hyperbolicTanFunc (Fit_Depth_Array, maxHeight, inflection, width)
			# Plotting the curve fit on the same graph as the density plot of water. 
			plt.plot (Fit_Depth_Array, yNew, color = "red", label='Fitted Water Density')
			plt.plot (depth_plot, water_density_plot, color = "green", label = 'Not-fitted Water Density')
			plt.title ('Water Density at the Interface')
			plt.xlabel ('Depth (nm)')
			plt.ylabel ('Density of Water (g/ml')
			plt.show ()
			print ('DDC File Complete')
			
			
			# This code currently works for two files (DDC100_01 and DDC100_02), however, it does not 
			# read in the third, fourth, or fifth file. Need to debug and figure out why 
