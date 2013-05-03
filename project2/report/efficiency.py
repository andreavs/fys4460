import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

t1 = [4.00549912453, 2.20352911949,2.71974420547,2.54189491272]
c1 = [1,2,3,4]

plt.plot(c1,t1)
t1 = [t1[0]/1, t1[0]/2, t1[0]/3, t1[0]/4]
plt.hold('on')
plt.plot(c1,t1)

t2 = [42.3537809849, 24.1625180244, 19.3942358494, 17.8031280041]
plt.plot(c1,t2)
t2 = [t2[0]/1, t2[0]/2, t2[0]/3, t2[0]/4]
plt.plot(c1,t2)

rc('text', usetex=True)
#rc('font', family='serif')
plt.title('Time spent', fontsize=18)
plt.ylabel('time [s]', fontsize=16)
plt.xlabel('number of coures', fontsize=16)
plt.legend(['actual time, ' + r'N=2048', 'ideal time, ' + r'N=2048', 'actual time, ' + r'N=32000', 'ideal time, ' + r'N=32000'], loc='upper right')
plt.show()






plt.show()