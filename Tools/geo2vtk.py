"""
geo2vtk.py
Rupert Deese
November 2013

Script to translate geoFEST output files into VTK files, one for each timestep.

Usage: python geo2vtk.py inputDirectory

Input:

The directory containing the output from geoFEST (relative or absolute path)

Output: 

A vtk file for each timestep of the run, of the form

dispFile_ts1.vtk
dispFile_ts2.vtk
...

.png images of graphs providing displacement and velocity information for the nodes in
nodeset 4 (right now it is assumed to be surface nodes)

TODO:

- Debug/improve robustness. Question as to what kind of bad input it makes sense to
	check for in this case.
- Make sure the best possible portfolio of graphs is being produced for this purpose.
- Add some kind of control feature to choose graphs/vtk/both

POSSIBLE IMPROVEMENTS:

- condense vtk files into one file with time information.

"""

import sys
import os
import subprocess
import math
from matplotlib import pyplot as plt

secsInYear = 31557600


##################################################################

## ACTUAL PROGRAM EXECUTION STARTS HERE

##################################################################
## File handling

# Check if directory is valid.
if not os.path.exists(sys.argv[1]):
	raise OSError('Invalid directory provided.\nUsage: python geo2vtk.py inputDirectory')

toolsDir = os.getcwd()
print 'Hi,' + toolsDir

# Go into the directory provided as an argument
if sys.argv[1] != '':
	try:
		os.chdir(sys.argv[1])
	except OSError:
		raise IOError('Invalid directory provided.')

# Get name of dispFile and number of nodepoints
basic = open('basic.dat', 'r')
dispName = basic.readline().split()[0]
for line in basic:
	data = line.split()
	if len(data) > 0 and data[0] == 'NUMNP':
		numNodes = data[1]

# Write the new coordinate file
coordIn = open('coord.dat', 'r')
coordOut = open('coord_vtk.dat', 'w')
coordOut.write(numNodes+'\n')
for line in coordIn:
	coordOut.write(line)
coordOut.close()
coordIn.close()

# Write the new element data file
elIn = open('eldata.dat', 'r')
elOut = open('eldata_vtk.dat', 'w')
numEl = elIn.readline().split()[1]
for i in range(11):
	elIn.readline()
elOut.write(numEl+'\n')
for line in elIn:
	elOut.write(line)
elOut.close()
elIn.close()

# Parse the plotNodes file to find the points of interest
pointsFile = open('plotNodes.dat', 'r')
plotPoints = []
for line in pointsFile:
	plotPoints += line.split()
pointsFile.close()


# Parse the dispFile to find all timesteps, and make a dispFile for each step.

# Stores all of the plot data. CURRENTLY UNUSED
timeSteps = []
# Keeps track of the current timestep
vtkIndex = 0

dispIn = open(dispName, 'r')
line = dispIn.readline()
while line != '':
	# Marks the beginning of a timestep.
	if line == " Global coordinates & displacements & delt displacements \n":

		# Reset plot data
		plotData = []
		for i in range(6):
			plotData.append([])

		dispIn.readline()		# skip whitespace
		timeData = dispIn.readline().split()
		timePoint = float(timeData[3])
		timeStep = float(timeData[-1])
		dispIn.readline()
		dispIn.readline()

		dispOut = open('disp'+str(vtkIndex)+'.dat', 'w')
		dispOut.write(str(timeStep)+'\n')

		for i in range(int(numNodes)):
			line = dispIn.readline()
			data = line.split()

			if data[1] in plotPoints:
				# Get position
				x = float(data[2])/1000
				y = float(data[3])/1000
				z = float(data[4])/1000
				r = math.sqrt(math.pow(x,2) + math.pow(y,2))

				# Get displacements in this timestep
				ux = float(data[5])
				uy = float(data[6])
				uz = float(data[7])
				if r != 0:
					ur = ux*x/r + uy*y/r
				else:
					ur = 0

				# Get velocities in this timestep if they are valid
				# (i.e. just not the first step in the model)
				if timeStep == 0:
					vx = vy = vz = 0
				else:
					#print "Timestep is: ", str(timeStep/secsInYear), " years."
					vx = float(data[8])*1000/(timeStep/secsInYear)
					vy = float(data[9])*1000/(timeStep/secsInYear)
					vz = float(data[10])*1000/(timeStep/secsInYear)
				if r != 0:
					vr = vx*x/r + vy*y/r
				else:
					vr = 0
	
				# Add relevant quantities to the plot data.
				plotData[0].append(timePoint)
				plotData[1].append(r)
				plotData[2].append(ur)
				plotData[3].append(uz)
				plotData[4].append(vr)
				plotData[5].append(vz)

			dispOut.write(line)

		dispOut.close()

		timeSteps.append(plotData)

		# Plot surface displacements
		plt.scatter(plotData[1], plotData[2], c='r', hold=False, label='Radial component')
		plt.scatter(plotData[1], plotData[3], c='b', hold=True, label='Vertical component')
		plt.title('t='+str(int(round(plotData[0][0]/secsInYear, 1))) +' yrs Surface Displacements')
		plt.xlabel('Radial Distance (km)')
		plt.ylabel('Surface Displacement (m)')
		plt.legend(loc='lower right', frameon=False)
		plt.axis([0, 3000, -700, 100]);
		plt.savefig('ts'+str(vtkIndex)+'_displacements.png')

		# Plot surface velocities
		plt.scatter(plotData[1], plotData[4], c='r', hold=False, label='Radial component')
		plt.scatter(plotData[1], plotData[5], c='b', hold = True, label='Vertical component')
		plt.title('t='+str(int(round(plotData[0][0]/secsInYear, 1))) +' yrs Surface Velocities')
		plt.xlabel('Radial Distance (km)')
		plt.ylabel('Surface Velocity (mm/yr)')
		plt.legend(loc='lower right', frameon=False)
		plt.axis([0, 3000, -5, 100]);
		plt.savefig('ts'+str(vtkIndex)+'_velocities.png')

		vtkIndex += 1
	line = dispIn.readline()
dispIn.close()

# Run Greg's script to produce vtk files
vtkPrefix = dispName.split('.')[0]
for i in range(vtkIndex):
	vtkName = vtkPrefix+'_ts'+str(i)+'.vtk'
	print [toolsDir+'/gft2vtk', coordOut.name, 'disp'+str(i)+'.dat', 'stresses.dat', vtkName]
	try:
		subprocess.Popen([toolsDir+'/gft2vtk', coordOut.name, elOut.name, 'disp'+str(i)+'.dat', 'stresses.dat', vtkName]).wait()
	except OSError:
		raise OSError('\'./gft2vtk\' executable not found. geo2vtk requires Greg\'s parser, gft2vtk.c')





