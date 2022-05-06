#!/usr/bin/env python3

# This code reads a density file in text mode, and parses the water density

import numpy as np
import matplotlib.pyplot as plt 
import scipy 
from scipy import optimize 

# Open data file
density_File= open("DDC100_01_mass_density.txt")

#Read first line with column headings.
headerLine=density_File.readline()

# Put column headings into a list for future indexing of the data
data_Items=headerLine[2:].split()

print (data_Items)

# Load the density data into a 2-D array
data_Array=np.loadtxt(density_File)

# Extract x and y arrays to plot the water density
depth_plot=data_Array[:,data_Items.index('nm')] 
water_density_plot=data_Array[:,data_Items.index('H2O')]

# Extract subset water density array to fit to the interface width function
# First defind the depth in nm for start and end, then dtermine the array indicies for start and stop of data
analysis_start=-1.0
analysis_start_index=np.where(data_Array[:,data_Items.index('nm'):] == analysis_start)[0][0]
analysis_end = 1.0
analysis_end_index=np.where(data_Array[:,data_Items.index('nm'):] == analysis_end)[0][0]+1

# Extract the fitting data into x and y arrays
Fit_Depth_Array=data_Array[analysis_start_index:analysis_end_index,data_Items.index('nm')]
Fit_Water_Interface = data_Array[analysis_start_index:analysis_end_index,data_Items.index('H2O')]

# Density plot for water 
#plt.plot (depth_plot, water_density_plot)
#plt.title ('Density Plot for Water')
#plt.xlabel ('Depth (nm)')
#plt.ylabel ('Density of Water (g/ml')

# Fitting the water density at the interface, this will give the width of the interface 
def hyperbolicTanFunc (m, a, b, c):
	return (a/2)*(1-np.tanh((m-b)/(c)))
maxHeight = 1
inflection = 0
width = 0.4
Fit_Water_Interface = hyperbolicTanFunc (Fit_Depth_Array, maxHeight, inflection, width)
popt, cov = scipy.optimize.curve_fit (hyperbolicTanFunc, Fit_Depth_Array, Fit_Water_Interface)
maxHeight#!/usr/bin/env python3

# This code reads a density file in text mode, and parses the water density

import numpy as np
import matplotlib.pyplot as plt 
import scipy 
from scipy import optimize 

# Open data file
density_File= open("DDC100_01_mass_density.txt")

#Read first line with column headings.
headerLine=density_File.readline()

# Put column headings into a list for future indexing of the data
data_Items=headerLine[2:].split()

print (data_Items)

# Load the density data into a 2-D array
data_Array=np.loadtxt(density_File)

# Extract x and y arrays to plot the water density
depth_plot=data_Array[:,data_Items.index('nm')] 
water_density_plot=data_Array[:,data_Items.index('H2O')]

# Extract subset water density array to fit to the interface width function
# First defind the depth in nm for start and end, then dtermine the array indicies for start and stop of data
analysis_start=-1.0
analysis_start_index=np.where(data_Array[:,data_Items.index('nm'):] == analysis_start)[0][0]
analysis_end = 1.0
analysis_end_index=np.where(data_Array[:,data_Items.index('nm'):] == analysis_end)[0][0]+1

# Extract the fitting data into x and y arrays
Fit_Depth_Array=data_Array[analysis_start_index:analysis_end_index,data_Items.index('nm')]
Fit_Water_Interface = data_Array[analysis_start_index:analysis_end_index,data_Items.index('H2O')]

# Density plot for water 
#plt.plot (depth_plot, water_density_plot)
#plt.title ('Density Plot for Water')
#plt.xlabel ('Depth (nm)')
#plt.ylabel ('Density of Water (g/ml')

# Fitting the water density at the interface, this will give the width of the interface 
def hyperbolicTanFunc (m, a, b, c):
	return (a/2)*(1-np.tanh((m-b)/(c)))
maxHeight = 1
inflection = 0
width = 0.4
popt, cov = scipy.optimize.curve_fit (hyperbolicTanFunc, Fit_Depth_Array, Fit_Water_Interface)
perr = np.sqrt(np.diag(cov))
print (perr)
maxHeight, inflection, width = popt 
# Verifying the hyperbolic function parameters 
print (maxHeight, inflection, width)
yNew = hyperbolicTanFunc (Fit_Depth_Array, maxHeight, inflection, width)

# Curve fit for the fitted water density 
#plt.plot (Fit_Depth_Array, yNew, color = "red")
#plt.title ('Water Density at the Interface')
#plt.xlabel ('Depth (nm)')
#plt.ylabel ('Density of Water (g/ml)')

# Plotting the curve fit on the same graph as the density plot of water. 
plt.plot (Fit_Depth_Array, yNew, color = "red", label='Fitted Water Density')
plt.plot (depth_plot, water_density_plot, color = "green", label='Not-fitted Water Density')
plt.title ('Water Density at the Interface')
plt.xlabel ('Depth (nm)')
plt.ylabel ('Density of Water (g/ml')
plt.show ()

