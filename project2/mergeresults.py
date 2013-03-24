import glob

def mergeresults():
	mywritefile = open('results/results.xyz', 'w')
	folder = 'MDnano-build-desktop-Qt_4_8_1_in_PATH__System__Release/experiments/thermalize/results/'
	N = len(glob.glob(folder + 'results*.xyz'))
	for i in range(N):
		filenumber = i;
		filenumberstring = str(filenumber).zfill(4)
		filename = folder + 'results'+ (filenumberstring) + '.xyz'
		myreadfile = open(filename, 'r')
		for line in myreadfile.readlines():
			mywritefile.write(line)
		myreadfile.close()
	mywritefile.close()

if __name__ == '__main__':
	mergeresults()