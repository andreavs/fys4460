import numpy as np
import os
import numpy.random
import matplotlib.pyplot as plt
import glob

class GenerateMaterial:
	def __init__(self, fromexperiment, toexperiment):
		self.fromState = "experiments/" + fromexperiment + "/results/lastState.xyz"
		self.toState = "experiments/" + toexperiment + "/results/lastState.xyz"
		fromFile = open(self.fromState)
		self.totalAtoms = int(fromFile.readline())
		fromFile.readline()
		self.contents = range(self.totalAtoms)
		counter = 0
		for line in fromFile:
			self.contents[counter] = line.split(' ')
			counter += 1
		#print self.contents
		#print self.totalAtoms

	def makecylinder(self, radius, orientation='x'):
		#makes a cylinder in the middle of the guy
		if orientation == 'x':
			n = 1
		pos = np.zeros(self.totalAtoms)	
		for i in range(self.totalAtoms):
			pos[i] = float(self.contents[i][n])
		size = max(pos)
		origin = size/2.
		r2 = radius*radius
		for i in range(self.totalAtoms):
			if (float(self.contents[i][2])-origin)**2 + (float(self.contents[i][3])-origin)**2 < r2:
				self.contents[i][0] = "Ne"
			else:	
				self.contents[i][10] = "1"


	def makespheres(self, nSpheres, minradius, maxradius):
		origins = np.zeros([nSpheres,3])
		origins[:,0] = numpy.random.random(nSpheres)
		origins[:,1] = numpy.random.random(nSpheres)
		origins[:,2] = numpy.random.random(nSpheres)
		radii = minradius + (maxradius-minradius)*numpy.random.random(nSpheres)
		pos = np.zeros(self.totalAtoms)	
		for i in range(self.totalAtoms):
			pos[i] = float(self.contents[i][1])
		size = max(pos)
		origins = size*origins
		r2 = radii*radii
		counter = 0
		for i in xrange(self.totalAtoms):
			for j in xrange(nSpheres):
				if (float(self.contents[i][1])-origins[j][0])**2 + (float(self.contents[i][2])-origins[j][1])**2 + (float(self.contents[i][3])-origins[j][2])**2 < r2[j]:
					self.contents[i][0] = "Ne"
			if self.contents[i][0]  != "Ne":	
				counter+=1.

			else:
				self.contents[i][10] = "1"
				self.contents[i][4] = "0.0"
				self.contents[i][5] = "0.0"
				self.contents[i][6] = "0.0"
				#counter += 1.
		return counter/self.totalAtoms
		


	def saveAlteredState(self):
		os.system("touch " + self.toState)
		toFile = open(self.toState,'w')
		toFile.write(str(self.totalAtoms) + "\n")
		toFile.write("This line has not unintentionally been left unblank \n")
		for i in xrange(self.totalAtoms):

			toFile.write(' '.join(self.contents[i]))



	def removeSomeUnFrozenAtoms(self, eps):
		newcontent = []
		for i in range(self.totalAtoms):
			if self.contents[i][10] == "0":
				test = np.random.random()
				if test < eps: 
					newcontent.append(self.contents[i])
			else: 
				newcontent.append(self.contents[i])
		self.contents = newcontent
		self.totalAtoms = len(self.contents)


if __name__ == '__main__':
	

	gm = GenerateMaterial("thermalize", "spheres")
	# p = gm.makespheres(20, 5.87, 8.81)
	
	# gm = GenerateMaterial("thermalize", "cylinder")
	# gm.makecylinder(5.78)
	# gm.removeSomeUnFrozenAtoms(0.5)
	# gm.saveAlteredState()
	

	porosities = []
	fluxi = []
	nExperiments = 100
	for i in range(nExperiments):
		#gm.makecylinder(5.87)
		gm = GenerateMaterial("thermalize", "spheres")
		p = gm.makespheres(20, 5.87, 8.81)
		print i
		print p
		porosities.append(p)
		gm.removeSomeUnFrozenAtoms(0.5)
		gm.saveAlteredState()
		os.system("./MDnano nanopourous")
		os.system("./MDnano gravitycyl")
		folder = 'experiments/gravitycyl/results/'
		N = 2
		systemFilename = folder + 'systemresults.txt'
		myresults = open(systemFilename, 'r')
		#myresults.readline()
		contents = []
		myresults.readline()
		i = 0
		for line in myresults:
			contents.append(line.split(' '))
			print contents[i]
			i = i + 1

		print contents[N-1][7]
		fluxi.append(float(contents[-1][-2]))

	print fluxi
	print porosities
	fluxi = np.array(fluxi)
	fluxi = 0.01*fluxi
	plt.plot(porosities, fluxi, 'o')
	plt.title("Porosity vs permeability")
	plt.xlabel("porosity")
	plt.ylabel("permeability")
	plt.show()