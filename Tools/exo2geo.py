"""
exo2geo.py
Rupert Deese
October 2013

Script to translate Exodus II files (so far only produced in CUBIT) into geoFEST format.

Usage: python exo2geo.py inputFile.exo

Input:

An Exodus II file (produced in CUBIT is all that I have tested). 

NODESETS:
1: nodes fixed in x
2: nodes fixed in y
3: nodes fixed in z
4: nodes to plot (gft2vtk.py)
5: split nodes

SIDESETS:
1: surface traction (can choose arbitrary side and set to zero for no surface forces)
2-n: boundaries for buoyancy conditions

Output: 
Incomplete versions of the main files required for a geoFEST run:

coord.dat
bcc.dat
eldata.dat
surfdata.dat
buoydata.dat

TODO:

- Improve program robustness: currently program does basic checking for malformed input.
- Test with no sidesets: should behave correctly now or with minimal debugging. 
	Ideal result would be empty buoydata.dat and surfdata.dat files.
- Finish commenting code.

"""

import sys
import os
import subprocess
from argparse import ArgumentParser 

secsInYear = 31557600.0

##################################################################
## Helper functions for read/write

def readInputToList(inputFile, outputList):
	endLine = outputList[-1] if len(outputList) != 0 else ''
	line = inputFile.readline()
	while endLine != ';':
		outputList += map(lambda X: X[:-1] if X[-1] == ',' else X, line.split())
		line = inputFile.readline()
		endLine = outputList[-1]
	del outputList[-1]
def getGFESTSide(side):
	if side == 1:
		return 3
	elif side == 2:
		return 1
	elif side == 3:
		return 2
	elif side == 4:
		return 4
	else:
		sys.stderr.write("Exodus file parsing error: side # " + side + " of element " + \
			els(i) + " out of range!")
		return None

def writeCoords(numNodes):
	print 'Writing node coordinate data...'
	crdFile = open('coord.dat', 'w')
	for i in range(numNodes):
		crdFile.write(str(i+1) + ' ' + str(1000*float(coords[i])) + ' ' + str(1000*float(coords[i+numNodes])) 
			+ ' ' + str(1000*float(coords[i+2*numNodes])) + '\n')
	crdFile.write('0 0\n')
	crdFile.close()

def writeBCC(bcc):
	print 'Writing boundary condition data...'
	bccFile = open('bcc.dat', 'w')
	for line in bcc:
		bccFile.write(' '.join(map(str, line)) + '\n')
	bccFile.write("0 0\n")
	bccFile.close()

def writeElData(elems, numat, numsuf, numbuoy):
	print 'Writing element connectivity data...'
	elemFile = open('eldata.dat', 'w')
	elemFile.write('NUMEL ' + str(len(elems)) + '\n' +
					'TYPE 4\n' +
					'NUMAT ' + str(numat) +'\n' +
					'NUMSUF ' + str(numsuf) + '\n' +
					'NUMBUOY ' + str(numbuoy) + '\n' +
					'NSPLIT ' + str(numsplit) + '\n\n0\n\n')

	for i in range(1, numat+1):
		vals = raw_input("Enter space separated float values for lambda mu eta n " + \
							" gx gy gz for the material corresponding to block #" + \
							str(i) +": ")
		try:
			v = map(float, vals.split())
			if len(v) != 7:
				raise IOError('Wrong number of values.')
		except (ValueError, IOError):
			print "Non-float or wrong number of values entered!"
			elemFile.close()
			sys.exit(1)
		elemFile.write(vals+'\n')

	elemFile.write('\n')

	for line in elems:
		elemFile.write(' '.join(map(str, line)) + '\n')
	elemFile.write("0 0\n")
	elemFile.close()

def writeSurfData(elsSides, totalSimTime):
	print 'Writing surface traction data...'
	sideList = open('surfdata.dat', 'w')
	num_epochs = raw_input("Enter the number of epochs in this run: ")
	try:
		num_epochs = int(num_epochs)
	except ValueError:
		print "Not an integer!"
		sideList.close()
		sys.exit(1)

	for i in range(1, num_epochs+1):
		forceVector = raw_input("Enter the force of surface traction on sideset 1" + \
								" in epoch number " + str(i) + \
								" as a floating point vector \"X Y Z\": ")
		try:
			v = map(float, forceVector.split())
			if len(v) != 3:
				raise IOError("Wrong number of values.")
		except (ValueError, IOError):
			print "Not a valid vector of three floats."
			sideList.close()
			sys.exit(1)

		if i > 1:
			startTime = raw_input("Enter a time in years for the start of this epoch: ")
			try:
				startTime = float(startTime)*secsInYear
				if startTime > totalSimTime or startTime < 0:
					raise IOError("Bad start time.")
			except (ValueError, IOError):
				print "Not a valid start time."
				sideList.close()
				sys.exit(1)

			sideList.write('0 1\n\n'+str(startTime)+'\n')

		for line in elsSides:
			sideList.write(' '.join(map(str, line)) + ' ' + forceVector + "\n")

	sideList.write("0 0\n")
	sideList.close()

def writeBuoyData(buoyData):
	print 'Writing bouyancy data (in order of increasing sideset #)...'
	buoyFile = open('buoydata.dat', 'w')
	ssN = 2
	for sideSet in buoyData:
		forceVector = raw_input("Enter the direction of buoyancy on sideset " + str(ssN) + \
			" as a floating point vector \"X Y Z\": ")
		if len(forceVector.split()) != 3:
			sys.stderr.write("User input error: wrong dimension of vector entered.")
		buoyancy = raw_input("Enter a buoyancy value (rho*g): ")
		ssLen = len(sideSet)
		if (args.radial):
			ssLen = -ssLen
		buoyFile.write(str(ssLen) + ' ' + forceVector + ' ' + buoyancy + '\n')
		for line in sideSet:
			buoyFile.write(' '.join(map(str, line)) + "\n")
		buoyFile.write('0 0\n')
		ssN += 1
	buoyFile.close()


# FIXME PRESUMES NODE MEMBERSHIP IN ONLY ONE FAULT
def writeSplitData(splitData):
	print 'Writing split node file...'
	splitFile = open('fltdata.dat', 'w')
	print 'The Cubit file has nodesets specifying ' + str(len(splitData)) + ' faults.'
	faultNum = 1
	for faultNodes in splitData:

		bVector = raw_input("Enter a B vector for fault " + str(faultNum) + \
			" as a floating point vector \"X Y Z\": ")
		if len(bVector.split()) != 3:
			sys.stderr.write("User input error: wrong dimension of vector entered.")

		sVector = raw_input("Enter an S vector for fault " + str(faultNum) + \
			" as a floating point vector \"X Y Z\": ")
		if len(sVector.split()) != 3:
			sys.stderr.write("User input error: wrong dimension of vector entered.")

		slip = raw_input("Enter a slip amount for all nodes on fault " + str(faultNum) + \
			" as a single floating point value:")

		for line in faultNodes:
			splitFile.write(str(line) + ' 1 ' + bVector + ' ' + sVector + ' ' + str(faultNum) + ' ' + slip + '\n')
			numsplit += 1

		faultNum += 1

	splitFile.close()

def writeUndersideData(elsSides):
	print 'Writing surface traction data...'
	sideList = open('underside.dat', 'w')
	forceVector = raw_input("Enter a floating point {X,Y,Z} as \"X Y Z\" vector for surface traction: ")
	if len(forceVector.split()) != 3:
		sys.stderr.write("User input error: wrong dimension of vector entered.")
	for line in elsSides:
		sideList.write(' '.join(map(str, line)) + ' ' + forceVector + "\n")
	sideList.write("0 0\n")
	sideList.close()

def writeBasicData(filename, numNodes):
	print "Writing basic.dat..."
	basic = open('basic.dat', 'w')
	basic.write(filename+'.out\n')
	comment1 = raw_input("Write two lines of comments to identify this run. Line 1: ")
	comment2 = raw_input("Line 2: ")
	if len(comment1) == 0:
		comment1 = "No comment."
	if len(comment2) == 0:
		comment2 = "No comment."
	basic.write(comment1+'*\n')
	basic.write(comment2+'*\n\n')
	basic.write('ELASTIC1	1\n' +
				'ELAS_OUT1	1\n' +
				'REFINE		0\n' +
				'REFINE_OUT	0\n' +
				'ELASTIC2	0\n' +
				'ELAS_OUT2	0\n' +
				'VISCO		1\n\n' +
				'0\n\n' +
				'NUMNP ' + str(numNodes) + '\n' +
				'NSD 3\n' +
				'NDOF 3\n' +
				'NRATES 0\n' +
				'SAVE_SHAPE 1\n' +
				'SOLVER 2\n' +
				'NUMGROUPS 1\n' +
				'NPRINTNODES -1\n' +
				'NPRINTELS 0\n' +
				'NTIMEGROUPS 1\n' +
				'NREFORM 1\n' +
				'NBACKUP 5000\n' +
				'NFLTGROUPS ' + str(len(fltGroups)) + '\n\n' +
				'0\n\n' +
				'NO_RESTART\n\n' +
				'NO_SAVE\n')
	basic.close()

def writeBCVData():
	bcv = open('bcv.dat', 'w')
	bcv.write('0 0\n')
	bcv.close()

def writePrintData():
	pD = open('print.dat', 'w')
	pD.write('0\n')
	pD.close()

def writeTimeData():
	timeData = open('time.dat', 'w')
	print "Writing to time.dat..."

	totalSimTime = raw_input("Enter the total length of the run in years: ")
	try:
		totalSimTime = float(totalSimTime)*secsInYear
	except ValueError:
		print "Not a number!"
		bcv.close()
		sys.exit(1)

	timeStep = raw_input("Enter the length of one time step in years: ")
	try:
		timeStep = float(timeStep)*secsInYear
	except ValueError:
		print "Not a number!"
		bcv.close()
		sys.exit(1)
	
	timeData.write(str(totalSimTime)+' 1.0 '+str(timeStep)+'\n\n')

	numOuts = raw_input("Number of times in the run to produce output: ")
	try:
		numOuts = int(numOuts)
	except ValueError:
		print "Not an integer!"
		bcv.close()
		sys.exit(1)
	timeData.write(str(numOuts)+'\n')

	for i in range (1, numOuts+1):
		printTime = raw_input("Enter a time in years for output #"+str(i)+": ")
		try:
			printTime = float(printTime)*secsInYear
		except ValueError:
			print "Not a number!"
			bcv.close()
			sys.exit(1)
		if printTime > totalSimTime or printTime < 0:
			print "Not a valid print time!"
			bcv.close()
			sys.exit(1)
		timeData.write(str(printTime)+'\n')

	timeData.close()
	return totalSimTime

##################################################################

## ACTUAL PROGRAM EXECUTION STARTS HERE

##################################################################
## File handling

parser = ArgumentParser()
parser.add_argument('-ut', '-undersidetraction', default=False,
	help="Use second side set for underside surface traction.")
parser.add_argument('file', type=str, help="exo file to convert to geoFEST format.")
parser.add_argument('-split', action='store_true', help="Nodesets numbered > 4 specify split nodes");
parser.add_argument('-radial', action='store_true', help="Buoyancies are radial and their directions specify an origin.");

args = parser.parse_args()

# Get path, cd into directory of exo file.
if not os.path.exists(args.file):
	raise OSError('Invalid directory and/or file provided.\nUsage: python exo2geo.py inputFile.exo')

path = args.file.split('/')
filename, ext = path[-1].split('.')

if ext != 'exo':
	raise IOError('Non-Exodus filetype provided as input (or wrong extension used).' + \
				  '\nUsage: python exo2geo.py inputFile.exo')

# If the CWD is not the same as the exo file's, change it to the folder containing the exo file.
if len(path) > 1:
	try:
		os.chdir(os.getcwd()+'/'+reduce(lambda X,Y: X+'/'+Y,path[:-1]))
	except OSError:
		raise IOError('Invalid directory provided.')

# Convert to txt file using ncdump utility (part of netCDF)
exoFile = open(filename+'.txt', 'w')

try:
	subprocess.Popen(['ncdump', path[-1]], stdout=exoFile).wait()
except OSError:
	raise OSError('\'ncdump\' command not found. exo2geo requires netCDF to run.')

exoFile.close()

exoFile = open(filename+'.txt', 'r')


##################################################################
## Parsing text dump of exo file.

###############
# Case switching:
# 0 = nothing
# 1 = coordinate input
# 2 = nodeset input (BCs)
# 4 = Element connectivity
# 5 & 6 = Side set 1 (surface traction)
# 7 & 8 = Side sets 2 and up (buoyancy boundaries)
###############

# Starting values for key program data.
buoyData = []
elems = []
fltGroups = []
elemNum = 1
case = 0
totalSimTime = 0
numat = 0
numsuf = 0
numbuoy = 0
numsplit = 0
numNodes = 0


# Start reading the file.
line = exoFile.readline()
while len(line) != 0:
	data = line.split()

	# For every line that isn't just a newline:
	if len(data) != 0:
		# CASE 0: Parse headers to redirect to other cases.
		if case == 0:
			if data[0] == 'num_nodes':
				numNodes = int(data[2])
				bcc = map(lambda X: [X+1, 0], range(numNodes))
				totalSimTime = writeTimeData()
			elif data[0] == 'coord':
				case = 1
			elif data[0][:-1] == 'node_ns':
				case = 2
			elif data[0][:-1] == 'connect':
				numat += 1
				elemType = data[0][-1]
				#print 'elemType = ', elemType
				case = 4
			elif data[0] == 'elem_ss1':
				case = 5
			elif data[0] == 'side_ss1':
				case = 6
			elif data[0] == 'elem_ss2' and args.ut:
				case = 9
			elif data[0] == 'side_ss2'and args.ut:
				case = 10
			elif data[0][:-1] == 'elem_ss':
				numbuoy += 1
				sideSetNum = data[0][-1]
				case = 7
			elif data[0][:-1] == 'side_ss':
				case = 8

		# GET NODE COORDINATES
		if case == 1:
			print 'Getting node coordinates...'
			coords = []
			readInputToList(exoFile, coords)
			writeCoords(numNodes)
			case = 0

		# GET BOUNDARY CONDITIONS
		if case == 2:
			print 'Getting boundary conditions...'
			nodes = map(lambda X: X[:-1] if X[-1] == ',' else X, data[2:])
			readInputToList(exoFile, nodes)
			if int(data[0][-1]) > 3:
				if int(data[0][-1]) > 4 and args.split:
					print 'Getting split nodes'
					fltGroups.append(nodes)
				else:
					if int(data[0][-1]) > 4:
						print 'More than 5 nodesets are present, nothing is handled above nodeset 5'
						# raise(IOError, "too many nodesets for command line options")
					else:
						print 'Getting nodes to plot'
						plotNodes = nodes
			else:
				for i in range(numNodes):
					bcc[i].append(0) if str(i+1) in nodes else bcc[i].append(1)
			## By convention, nodes_ns1 is fixed in x, node_ns2 is fixed in y, and node_ns3 is fixed in z
			case = 0

		# GET ELEMENTS
		if case == 4:
			print 'Getting element connectivity...'
			endline = ''
			line = exoFile.readline()
			while endline != ';':
				elems.append([elemNum, 0, elemType] + map(lambda X: X[:-1] if X[-1] == ',' else X, line.split()))
				endline = elems[-1][-1]
				elemNum += 1
				line = exoFile.readline()
			del(elems[-1][-1])
			case = 0

		# GET SURFACE TRACTION DATA
		if case == 5:
			print 'Getting surface traction data...'
			els = map(lambda X: X[:-1] if X[-1] == ',' else X, data[2:])
			readInputToList(exoFile, els)
			numsuf = len(els)
			case = 0
		# GET SURFACE TRACTION DATA
		if case == 6:
			print 'Getting surface traction data...'
			sides = map(lambda X: X[:-1] if X[-1] == ',' else X, data[2:])
			readInputToList(exoFile, sides)

			elsSides = []
			for i in range(len(els)):
				elsSides.append([int(els[i])])
				elsSides[i].append(getGFESTSide(int(sides[i])))
			writeSurfData(elsSides, totalSimTime)
			case = 0

		# GET BUOYANCY DATA
		if case == 7:
			print 'Getting buoyancy data...'
			els = map(lambda X: X[:-1] if X[-1] == ',' else X, data[2:])
			readInputToList(exoFile, els)
			case = 0
		# GET BUOYANCY DATA
		if case == 8:
			print 'Getting buoyancy data...'
			sides = map(lambda X: X[:-1] if X[-1] == ',' else X, data[2:])
			readInputToList(exoFile, sides)

			elsSides = []
			for i in range(len(els)):
				elsSides.append([int(els[i])])
				elsSides[i].append(getGFESTSide(int(sides[i])))
			buoyData.append(elsSides)
			case = 0

		# GET UNDERSIDE TRACTION DATA
		if case == 9:
			print 'Getting underside traction data...'
			els = map(lambda X: X[:-1] if X[-1] == ',' else X, data[2:])
			readInputToList(exoFile, els)
			case = 0
		# GET UNDERSIDE TRACTION DATA
		if case == 10:
			print 'Getting underside traction data...'
			sides = map(lambda X: X[:-1] if X[-1] == ',' else X, data[2:])
			readInputToList(exoFile, sides)

			elsSides = []
			for i in range(len(els)):
				elsSides.append([int(els[i])])
				elsSides[i].append(getGFESTSide(int(sides[i])))
			writeUndersideData(elsSides)
			case = 0

	# Read next line
	line = exoFile.readline()

# Create dummy stress file.
out = open('stresses.dat','w')
for i in range(elemNum):
	out.write(' '.join(['0']*10) + '\n')
out.close()

# Create file giving nodes to plot
out = open('plotNodes.dat', 'w')
out.write(' '.join(plotNodes)+'\n')
out.close

# have to do this at the end since BCC and element data requires reads from multiple node lists.
writeElData(elems, numat, numsuf, numbuoy)
writeBCC(bcc)
writeBuoyData(buoyData)
writeBCVData()
writePrintData()
writeSplitData(fltGroups)
writeBasicData(filename, numNodes)
exoFile.close()




