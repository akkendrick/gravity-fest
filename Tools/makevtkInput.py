import sys
import os
import subprocess

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

newDir = os.getcwd()
print 'Hi,' + newDir

# Get name of dispFile and number of nodepoints
basic = open('basic.dat', 'r')
dispName = basic.readline().split()[0]
for line in basic:
	data = line.split()
	if len(data) > 0 and data[0] == 'NUMNP':
		numNodes = data[1]

# Get name of dispFile and number of nodepoints
time = open('time.dat', 'r')
timeData = time.readlines()
numTimeSteps = int(timeData[2]) - 1


# Go into the directory provided as an argument
print numNodes

filename = sys.argv[2]

print filename

#vtkData = open(filename, 'r+')
vtkFile = open(filename, 'r+')
vtkData = vtkFile.readlines()

i = 0
for x in range(0,numTimeSteps):
        outname = 'simout_'+str(x)+'.txt'
	print outname
	outFile = open(outname, 'w')
	
	for line in vtkData[i+1:]:
                header = line[0:4]
                nodeTag = line[5:10]
                #print nodeTag
                if header == 'node':
       			outFile.write(line)
                        if nodeTag == numNodes:
                                i = i + 1
                                break
                i = i+1
                                
        j= 0 
#	print vtkData[i+2]
	for line in vtkData[i+1:]:
                lineList = line.split()
   		if line[0:4] == 'node':
   			break
   		if len(lineList) == 11:
   			if lineList[0] != 'X':
                                outFile.write(line)
                i = i + 1
