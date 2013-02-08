mywritefile = open('results.xyz', 'w')
folder = 'betterArgon-build-desktop-Qt_4_8_1_in_PATH__System__Release/'
N = 499
for i in range(N+1):
	filenumber = i;
	filenumberstring = str(filenumber).zfill(3)
	filename = folder + 'results'+ (filenumberstring) + '.xyz'
	myreadfile = open(filename, 'r')
	for line in myreadfile.readlines():
		mywritefile.write(line)
	myreadfile.close()
mywritefile.close()

