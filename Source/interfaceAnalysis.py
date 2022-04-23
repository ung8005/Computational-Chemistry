This Program will read in all of the atom data from the simulations using FrameDataClass.py and will create a 2-d of the uppermost water atoms (i.e. the interface) for each frame and save the z height and the atom type count. The resulting interface location data is analyzed with 2d FFT and fractal analysis.  Output three files:
		nn_##.interface_comp.txt including percenrt atomic composition of interface of O, H and holes
		nn_##.interface_FFT.txt averaged trace of 2-d fft with x axis and stdev
		nn_##.interface_FD.txt  averaged fractal area ratio, x axis and stdev
"""
import sys

import os.path
# String regular expressions
import re
# import numpy for extensive array and math functions
import numpy as np

# For multiprocessing
import multiprocessing as mp

import FrameDataClass as fd
import time


import math
from decimal import Decimal
#This is the size of your 2-d array for a surface z position and element type. Must be 2^n, and bigger is more often better
global size
size = 512

# Analysis is too slow to analyze every frame in a simulation.  This is the target number of frames to analyze.
global Analysis_Target
Analysis_Target = 100

# First check for the presence of expected data files		
		
def main():
	print("Welcome To Surface Analysis")
# Check for presence of two parameters in input line
	if len(sys.argv)!=3:
		print("Two input arguments are expected: The input pdb file for the simulations, and the integer number of the simulation.")
		sys.exit(1)
		
	basefile = sys.argv[1]
	
	simnumber= int(sys.argv[2])
		
# First check for presence of epected data files
	if not os.path.isfile(basefile+".siminfo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".siminfo"))
		sys.exit(1)
	if not os.path.isfile(basefile+".simtopo"):
		print("Simulation info file {} does not exist.  Exiting!".format(basefile+".simtopo"))
		sys.exit(1)
	
	global filename
	filename="{}_{:02d}".format(basefile,simnumber)
	
	if not os.path.isfile(filename+"_traj.framebin"):
		print("Simulation info file "+filename+"_traj.framebin does not exist.  Exiting!")
		sys.exit(1)

	
	global sim_data
#Read the data for the simulation
	sim_data=fd.FrameData(basefile,simnumber)
	
	print('Simulation dimensions {} x {} x {} {}'.format(sim_data.simdimensions.x, sim_data.simdimensions.y, sim_data.simdimensions.z, sim_data.simdimensions.units))
	
	if (sim_data.simdimensions.x != sim_data.simdimensions.y):
		print("Analysis assumes square dimensional simulation. Exiting!")
		sys.exit(1)

	global space_dim
	space_dim = sim_data.simdimensions.x


		
	print('Frame count ', sim_data.framecount)
	print('Atom count ', sim_data.atomcount)


# Interface Locating constants
		

#minimum z height for simulation
	global min_z
	min_z = -4.0 #in nanometers(-40 angstroms)

# Variable to convert real position to array index and inverse.  These constants reduce calculations in time-critical core function
	global cell_dimension
	cell_dimension = space_dim/size
	global inverse_cell_dimension
	inverse_cell_dimension = size / space_dim
	

#The number of atoms in the simulation
	global atomCount
	atomCount = sim_data.atomcount

#the number of frames to Analyze
	global Analysis_Count
	Analysis_Count = min(Analysis_Target,sim_data.framecount)
	Frame_Skip = sim_data.framecount // Analysis_Count
	
#The specific molecule we are looking for on our surface. The b signifies that this is a binary string.  Assumes water w/ van der Waals Radius as specified.
	global interface_molecule
	interface_molecule = b'HOH'
	
	inflate = 1.58
	global vdwRadius
	vdwRadius = {
	1 : 0.12*inflate,
	8 : 0.15*inflate
	}	
# Defines the range above and below the interface to test for surface atoms.  in nm
	global interfaceBound
	interfaceBound = 2.0
	
	
	
#2-d arrays for surface atom type and z height.  Defined as 1-d array for purposes of multiprocessing, which only uses 1-d indexing	
	global z_type
	mpArray_type = mp.Array('l',size*size)
		
	z_type = np.frombuffer(mpArray_type.get_obj())
		
	global z_height
	mpArray_height = mp.Array('d',size*size)
		
	z_height = np.frombuffer(mpArray_height.get_obj(),dtype='float64')

# Array to keep atom type -- from the simulation topography.  Used for interface composition	
	global atomLookup
	atomLookup = np.zeros(atomCount)

# Array for composition
	global interface_comp
	interface_comp = np.empty([Analysis_Count,3],int)

# Array and constants for 2-d fft

# data array

	global interface_fft
	interface_fft = np.empty([Analysis_Count,size//2],float)
	
	global hanning2d
	han= np.hanning(size)
	hanning2d = np.sqrt(np.outer(han,han))
	
	global hanning_correct
	hanning_correct=12.0
	
# Fractal Analysis constants

# box increment
	global box_increment
	box_increment = 2  # MUST BE EVEN
	global fd_size
	fd_size = size//box_increment
	
# fractal dimension data array 
	global interface_fd
	interface_fd = np.empty((Analysis_Count,fd_size))

# Constants used for time display			
	global intervals
	intervals = (
		('weeks', 604800),  # 60 * 60 * 24 * 7
		('days', 86400),    # 60 * 60 * 24
		('hours', 3600),    # 60 * 60
		('minutes', 60),
		('seconds', 1),
		)
	

#Fills Dictionary with atomic number from atom index
	print("Creating Atom Type Dictionary")
	global waterAtoms
	waterAtoms = []
	for atom in range(atomCount):
		#This if statement ensures that we only grab the information from the water molecule atoms
		if sim_data.atominfo[atom]['molecule'] == interface_molecule:
			waterAtoms.append(atom)
			atomLookup[atom] = sim_data.atominfo[atom]['atomicno']
		
#This portion fills the z_height and z_type arrays with the max z heights at every (x,y) and what atom type they are	
	print("Starting Interface Analysis.")

	start_time = time.time()
	Frames_Remaining = Analysis_Count
	
#frame loop	
	for frame in range(0,sim_data.framecount+1-Frame_Skip,Frame_Skip):
		print('Analysis of frame {} of {}'.format(frame + 1,sim_data.framecount))
		
		frame_start_time = time.time()

# Load Frame Data --  skip over
		for i in range(Frame_Skip):
			sim_data.next()

# Clear Result Arrays	
		z_height.fill(-1.0*interfaceBound)
		z_type.fill(0)
		

# Set up the multiprocessing tu use 1/2 of available cores.		
		print("Determining Interface composition and location")
		pool = mp.Pool(mp.cpu_count()//2)
# The multiprocessing is called here and closed when all threads are finished.
		results = pool.map(Test_Atom,[atom for atom in waterAtoms])
		pool.close()
		pool.join()

# Record surface composition				

		
		append_composition(frame//Frame_Skip)
		
		
		print("Begin FFT analysis")
		Analysis_FFT(frame//Frame_Skip)
		print("Begin Fractal Dimension Analysis")
		Analysis_FD(frame//Frame_Skip)
		
		#This portion estimates the remaining time left in the data retrieval
		Frames_Remaining -= 1
		remainingSeconds = (time.time() - frame_start_time) * (Frames_Remaining)
		
		print("Frame analysis Complete\n\tEstimated Remaining Time = {} ".format(display_time(remainingSeconds)))	
	
	sim_data.done()
	print("Writing Data Files")
	write_composition()
	write_FFT()
	write_FD()
	
	print("\nExecution time = {} ".format(display_time(time.time() - start_time)))
	
	
	print("Surface Analysis Complete for {}_{:02d}".format(basefile,simnumber))	

#######	
# This function assumes atom is index of atom in water, and tests if the sphere of the atom is above the current interface. Does not need to be called sequentially, and therefore is the core multiprocess function.
def Test_Atom(atom):
	
# First test if atom is within +/- interfacebound (8 angstrom defined above) of interface 0	
	z_center = sim_data.positions[atom,2]
	atom_type = int(atomLookup[atom])
	atom_radius = vdwRadius[atom_type]
	if z_center - atom_radius <= interfaceBound and z_center + atom_radius >= -interfaceBound:	
		
		x_pos = sim_data.positions[atom,0]
		y_pos = sim_data.positions[atom,1]
		
# Calculate the surface array index bounds for the atom radius	
		x_minIndex = math.ceil((x_pos - atom_radius)*inverse_cell_dimension)
		x_maxIndex = math.ceil((x_pos + atom_radius)*inverse_cell_dimension)
		
		y_minIndex = math.ceil((y_pos - atom_radius)*inverse_cell_dimension)
		y_maxIndex = math.ceil((y_pos + atom_radius)*inverse_cell_dimension)

	

		for x_index in range(x_minIndex,x_maxIndex):
			for y_index in range(y_minIndex,y_maxIndex):
#This calculates the intersection of the atom sphere with the line in the z direction at the surface index point. Uses always the larger intersection point
				calculation = atom_radius**2 - (x_index*cell_dimension - x_pos)**2 - (y_index*cell_dimension - y_pos)**2
				if calculation >= 0.0:
					z = z_center + np.sqrt(calculation)
					
					
# Test if the calculated height of atom at array position (x,y) is greater that current. If yes, then replace then height ant type.  Note: the arrays are 1-d for multiprocessing, therefore the index must be calculated.  Also, to accout for periodic boundaries calculate the index with a %size -- i.e. the remainder.
					matrix_index=x_index%size + ((y_index%size)*size)
					if z > z_height[matrix_index]:
						z_height[matrix_index] = z
						z_type[matrix_index] = atom_type
					
						
	return True
# End of the multiprocess core function


#######
# appends composition count to array
def append_composition(frame):
	holesCount = (z_type == 0).sum()
	hydrogenCount = (z_type == 1).sum()
	oxygenCount = (z_type == 8).sum()

	print("O atom count = {:d}\nH atom count = {:d}\nHole count = {:d}".format(oxygenCount, hydrogenCount, holesCount))

	interface_comp[frame]=(hydrogenCount,oxygenCount,holesCount)

	return True

#####	
# Completes 2-D fft, extracts diagonal and appends to array	
# calculates a power spectrum with a Hanning nwindow to reduce artifacts, then does approximate amplitude correction & normalization.
#Data saved for each frame in array.  File write  average, std calculation
def Analysis_FFT(frame):
	
	z_square = np.reshape(z_height,(size,size))
	ps = np.abs(np.fft.rfft2(z_square * hanning2d))*hanning_correct/size**2
	
	for i in range(size//2):
		interface_fft[frame, i] = ps[i+1,i+1]
	
	return True
	
#####
# Do the fractal dimension analysis using pyramid analysis	This is a revision of the Clark method:  The eight panel method of W. Sun
def Analysis_FD(frame):


# Some variables precalculated for efficiency	
	global x0, xm, x1, l2, l4, box_len
#loop linearly through multiple of starting box length

	for box_len in range(box_increment, size +1, box_increment):
		Sum_FD_Area=0.0   # Collects calculated area in a row -- from MP function 
		count=len(range (0, size, box_len))  # need a counter to know the number of cells in a column -- for normalization
		
# The base of the real pyramid is invarient for a given box length.  Center of pyramic is at x=0,y=0. Corners +-box_len/2
		corner = box_len//2*cell_dimension
		l2=corner*corner
		l4=l2*l2

# Loop through columns spaced by the length of box  0 and 1 are coreners  m is midpoint  data wrapping requires %size
		for x0 in range (0, size, box_len):
			xm=(x0+box_len//2)%size
			x1=(x0+box_len)%size
			
			for y0 in range(0,size,box_len):
				Sum_FD_Area +=Sum_Cell(y0)
			

#Save the sum and normalize by the area of the 'flat' base beneath it.
		interface_fd[frame,(box_len//box_increment)-1]= Sum_FD_Area/((box_len*cell_dimension*count)**2)
	
	return True

def Sum_Cell (ystart):
	y0=ystart*size
	ym=((ystart+box_len//2)%size)*size
	y1=((ystart+box_len)%size)*size
	
# calculate the coordinates of the 4 corners and four edge points going around the box of the current analysis box	Note:  center point is 0,0,0 and all z heights adjusted to this value.  indices precalculated for efficiency
#This method uses the algebraic calculation of the 8 prism area, simplified and made efficient.  math.sqrt is 5x faster than np.sqrt!  Sage was used to calculate the return expression		
	center= z_height[xm + ym]
	z00= z_height[x0 + y0]-center
	z01= z_height[x0 + ym]-center
	z02= z_height[x0 + y1]-center
	z12= z_height[xm + y1]-center
	z22= z_height[x1 + y1]-center
	z21= z_height[x1 + ym]-center
	z20= z_height[x1 + y0]-center
	z10= z_height[xm + y0]-center
		

	return 0.5*(math.sqrt(l4 + l2*(2*(z01*z01 - z00*z01) + z00*z00)) + math.sqrt(l4 + l2*(2*(z01*z01 - z01*z02) + z02*z02)) + math.sqrt(l4 + l2*(2*(z10*z10 - z00*z10) + z00*z00)) + math.sqrt(l4 + l2*(2*(z10*z10 - z10*z20) + z20*z20)) + math.sqrt(l4 + l2*(2*(z12*z12 - z02*z12) + z02*z02)) + math.sqrt(l4 + l2*(2*(z12*z12 - z12*z22) + z22*z22)) + math.sqrt(l4 + l2*(2*(z21*z21 - z20*z21) + z20*z20)) + math.sqrt(l4 + l2*(2*(z21*z21 - z21*z22) + z22*z22)))
	


#######
# Output functions

def write_composition():
	
	f=open(filename+".interface_comp.txt",'w')
	
	comp_ave = np.average(interface_comp, axis = 0)
	comp_std = np.std(interface_comp, axis = 0)

	f.write('Interface Composition\n')
	f.write('\tH atoms\tO atoms\tholes')
	f.write('\naverage')
	for i in range(3):
		f.write('\t{:.1f}'.format(comp_ave[i]))
	
	f.write('\nst dev')
	for i in range(3):
		f.write('\t{:.1f}'.format(comp_std[i]))

	f.write('\nSample Count = {:d}\n'.format(Analysis_Count))
	
	f.close()
	
	return True
	
	
def write_FFT():

	f=open(filename+".interface_FFT.txt",'w')

	FFT_ave = np.average(interface_fft, axis = 0)
	FFT_std = np.std(interface_fft, axis = 0)
	

	f.write('Interface FFT\n')
	f.write('x values (nm^-1)')
	for i in range(size//2-1):
		f.write('\t{:.3e}'.format(space_dim/(i+1)))
	
	f.write('\naverage')
	for i in range(size//2-1):
		f.write('\t{:.3e}'.format(FFT_ave[i]))
	
	f.write('\nst dev')
	for i in range(size//2-1):
		f.write('\t{:.3e}'.format(FFT_std[i]))
	
	f.write('\nSample Count = {:d}\n'.format(Analysis_Count))
	
	f.close()
	
	return True
	
def write_FD():

	f=open(filename+".interface_FD.txt",'w')
	
	FD_x_values = np.arange(box_increment,size+1,box_increment)* cell_dimension
	FD_ave = np.average(interface_fd, axis = 0)
	FD_std = np.std(interface_fd, axis = 0)
	
	f.write('Interface Fractal Dimension')
	f.write('\nx values (nm)')
	for i in range(size//box_increment):
		f.write('\t{:.3e}'.format(FD_x_values[i]))
		
	f.write('\naverage')
	for i in range(size//box_increment):
		f.write('\t{:.3e}'.format(FD_ave[i]))
		
	f.write('\nst dev')
	for i in range(size//box_increment):
		f.write('\t{:.3e}'.format(FD_std[i]))
		
	f.write('\nSample Count = {:d}\n'.format(Analysis_Count))
	
	f.close()
	
	return True




########
		
def display_time(seconds, granularity=2):
	result = []

	for name, count in intervals:
		value = seconds // count
		if value:
			seconds -= value * count
			if value == 1:
				name = name.rstrip('s')
			result.append("{} {}".format(value, name))
	return ', '.join(result[:granularity])
	
main()
