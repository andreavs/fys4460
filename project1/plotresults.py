import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import os, sys
import glob
import subprocess
from numpy import linalg as LA

	

def mergeresults():
	mywritefile = open('results/results.xyz', 'w')
	folder = 'betterArgon-build-desktop-Qt_4_8_1_in_PATH__System__Release/'
	N = len(glob.glob(folder + 'results*.xyz'))
	for i in range(N):
		filenumber = i;
		filenumberstring = str(filenumber).zfill(4)
		filename = folder + 'results'+ (filenumberstring) + '.xyz'
		myreadfile = open(filename, 'r')
		for line in myreadfile.readlines():
			mywritefile.write(line)
		myreadfile.close()
		os.remove(filename)
	mywritefile.close()


def plotsystemproperties():
	folder = 'betterArgon-build-desktop-Qt_4_8_1_in_PATH__System__Release/'
	N = len(glob.glob(folder + 'results*.xyz'))
	systemFilename = folder + 'systemresults.txt'
	myresults = open(systemFilename, 'r')
	energy = np.zeros(N)
	kinetic = np.zeros(N)
	potential = np.zeros(N)
	time = np.zeros(N)
	pressure = np.zeros(N)
	temperature = np.zeros(N)
	displacement = np.zeros(N)

	myresults.readline()
	counter = 0
	for line in myresults.readlines():
		content = line.split(' ')
		time[counter] = float(content[0])
		energy[counter] = float(content[1])
		kinetic[counter] = float(content[2])
		potential[counter] = float(content[3])
		pressure[counter] = float(content[4])
		temperature[counter] = float(content[5])
		displacement[counter] = float(content[6])
		counter += 1

	energyfluct = sum(energy[50:]**2)/len(energy[50:]) - \
	 (sum(energy[50:]/len(energy[50:])))**2
	print "The variance of the energies is: " 
	print energyfluct
	myplot = plt.figure()
	myplot.text(.5, .95, 'Energies in liquid phase', horizontalalignment='center') 
	plt.subplot2grid((2,2),(0,0))
	plt.plot(time,kinetic)
	plt.legend(['kinetic energy'], loc='upper right')
	plt.subplot2grid((2,2),(0,1))
	plt.plot(time, potential)
	plt.legend(['potential energy'], loc='lower right')
	plt.subplot2grid((2,2),(1,0),colspan=2)
	plt.plot(time, energy)
	plt.legend(['Total energy'], loc='lower right')
	myplot2 = plt.figure()
	myplot2.text(.5, .95, 'Physical properties of the system', horizontalalignment='center') 
	plt.subplot2grid((2,1),(0,0))
	plt.plot(time, 119.74*temperature)
	plt.title('Temperature in kelvin')
	plt.subplot2grid((2,1),(1,0))
	plt.plot(time,pressure)
	plt.title('Pressure')

	plt.figure()
	plt.plot(time, displacement)
	plt.title("Average square displacement as a function of time")
	plt.ylabel("Average square displacement [MD units]")
	plt.xlabel("Time [MD units]")

	plt.figure()
	plt.plot(temperature,pressure, 'p')
	plt.xlabel('Temperature [MD units]')
	plt.ylabel('Pressure [MD units]')
	plt.title('Pressure as a function of Temperature')
	plt.show()

def plotRadialDist():
	folder = 'betterArgon-build-desktop-Qt_4_8_1_in_PATH__System__Release/'
	N = len(glob.glob(folder + 'results*.xyz'))
	positions = np.zeros((3,2048))
	Lhalf = 4*1.545
	L = 2*Lhalf
	#otherpos = np.zeros((3,2048))
	diff = np.zeros(3)
	numBins = 100
	distancebins = np.zeros(numBins)
	for i in xrange(50,70):
		print i
		filenumber = i;
		filenumberstring = str(filenumber).zfill(4)
		filename = folder + 'results'+ (filenumberstring) + '.xyz'
		myreadfile = open(filename, 'r')
		counter = 0
		myreadfile.readline()
		myreadfile.readline()
		for line in myreadfile.readlines():
			contents = line.split()
			positions[0,counter] = float(contents[1])
			positions[1,counter] = float(contents[2])
			positions[2,counter] = float(contents[3])
			counter += 1

		for k in xrange(2048):
			for l in xrange(k+1,2048):
				diff = np.copy(positions[:,k] - positions[:,l])
				for h in range(3):
					diff[h] = (abs(diff[h]) < diff[h] + L)*(abs(diff[h])<abs(diff[h] - L))*diff[h] + \
						(diff[h] - L)*(abs(diff[h])>abs(diff[h] - L))*(abs(diff[h]) + L>abs(diff[h] - L)) + \
						(diff[h] + L)*(diff[h]+L<abs(diff[h] - L))*(abs(diff[h]) > diff[h]+L)
				distance = LA.norm(positions[:,k] - positions[:,l],2)
				index = int(distance*numBins/Lhalf)
				if(index < numBins):
					distancebins[index] += 1
		myreadfile.close()
	dr = Lhalf/numBins
	for i in range(numBins):
		r = i*Lhalf/numBins
		volume = 4*np.pi*((r+dr)**3 - r**3)/3.
		distancebins[i] /= volume
	plt.bar(range(numBins),distancebins, width=1.0, bottom=0)
	plt.title("relative Radial distribution function from 0 to L/2")
	plt.xlabel("Distance")
	plt.ylabel("Relative distribution")

	plt.show()

def plotdistributions():
	folder = 'betterArgon-build-desktop-Qt_4_8_1_in_PATH__System__Release/'
	N = len(glob.glob(folder + 'results*.xyz'))
	N2 = 8*8*8*4
	vx = np.zeros((N,N2))
	v = np.zeros((N,N2))
	displacement = np.zeros((N,N2))
	for i in xrange(N):
		filenumber = i;
		filenumberstring = str(filenumber).zfill(4)
		filename = folder + 'results'+ (filenumberstring) + '.xyz'
		myreadfile = open(filename, 'r')
		myreadfile.readline()
		myreadfile.readline()
		j = 0
		for line in myreadfile:
			content = line.split()
			vx[i,j] = float(content[4])
			v[i,j] = float(content[9])
			displacement[i,j] = float(content[8])
			j = j+1
		myreadfile.close()
	vaxis = [0,5,0,400]
	animatedistribution(v,vaxis,'velocity\ distribution')
	vxaxis = [-5,5,0,400]
	animatedistribution(vx,vxaxis,'x-velocity\ component')
	displacementaxis = [0,15,0,2000]
	animatedistribution(displacement, displacementaxis, 'square\ displacement\ distribution')

def animatedistribution(data, axis, myfilename):
	N = len(data[:,1])

	for n in xrange(N):
		t = 0.005*n
		plt.hist(data[n,:],20,hold=True)
		plt.xlabel('x-position')
		plt.ylabel('y-position')
		#
		# Notice the use of LaTeX-like markup.
		#
		plt.title(myfilename + r' $t=%6f$' % t, fontsize=20)
		plt.axis(axis)
		filename = str('%06d' % n) + '.png'
		plt.savefig(filename, dpi=100)
		
		#
		# Let the user know what's happening.
		#
		print 'Wrote file', filename
		
		#
		# Clear the figure to make way for the next image.
		#
		plt.clf()


	command = ('mencoder'+ ' ' +
	'mf://*.png'+ ' ' +
	'-mf'+ ' ' +
	'type=png:w=800:h=600:fps=25'+ ' ' +
	'-ovc'+ ' ' +
	'lavc'+ ' ' +
	'-lavcopts'+ ' ' +
	'vcodec=mpeg4'+ ' ' +
	'-oac'+ ' ' +
	'copy'+ ' ' +
	'-o'+ ' ' +
	myfilename + '.avi')
	
	print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)

	os.system(command)

	print "\n\n The movie was written to '" + myfilename + "'"

	os.system('rm *.png')


if __name__ == '__main__':
	#plotRadialDist()
	#plotdistributions()
	plotsystemproperties()
	#mergeresults()