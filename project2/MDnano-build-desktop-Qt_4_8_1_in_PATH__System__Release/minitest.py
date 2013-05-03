import generateMaterial
import os
import time
a = time.time()

#os.system("./MDnano thermalize")
os.system("./MDnano nanopourous")
#os.system("./MDnano gravitycyl")
#c = os.system("./MDnano aftergrav")
b = time.time()
print b-a
#print c
