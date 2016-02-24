import matplotlib.pyplot as plt
import numpy as np
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
freq_1  =np.array([.2061,.4163,.627,.8285,1.0324,1.2462,1.4482,1.6594,1.88514,2.057810,2.277810,
    2.466310, 2.67500,2.887,3.0967,3.2990,3.4960,3.704,3.902,4.116,4.3298])
freq_2= np.array([.3096,.6243,.940,1.244,1.550,1.8716,2.173,2.491,2.780,3.090,3.418,3.703,4.017,
    4.335,4.650,4.955,5.236,5.561,5.870,6.19,6.5168])
current=np.arange(0,2.1, .1)

m1, yint1 = np.polyfit(current, freq_1,1)
m2, yint2=np.polyfit(current, freq_2,1)
best1 = current*m1+yint1
best2=current*m2 +yint2
plt.plot(current, freq_1,'o',ms=10,label='freq1')
plt.plot(current, freq_2,'*', ms=10,label='freq2')
plt.plot(current, best1,'r', label='bestfit1')
plt.plot(current, best2,'g', label='bestfit2')
plt.legend(loc=2, fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=15)
plt.xlabel('Current [A]',size=40)
plt.ylabel('Resonant Frequency [MHz]', size=40)

plt.show()

plt.plot(current, freq_1-best1,'*',label='fit1 resid')
plt.plot(current, freq_2-best2,'o', label='fit2 resid')
plt.legend(loc=2, fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=15)
plt.xlabel('Current [A]',size=40)
plt.ylabel(r'$\Delta$ in Resonant Frequency [MHz]', size=40)
plt.show()
print np.std(freq_1-best1,ddof=1), np.mean(freq_1-best1)
print np.std(freq_2-best2,ddof=1), np.mean(freq_2-best2)


