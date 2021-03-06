import matplotlib.pyplot as plt
import numpy as np
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
freq_1  =np.array([-3.88, -3.67, -3.48, -3.2769, -3.07,
       -2.87, -2.66, -2.455, -2.25, -2.044,
       -1.839, -1.634, -1.43, -1.22 , -1.017,
       -0.81, -0.60651513, -0.401 , -0.196,  0.00975,.2061,.4163,.627,.8285,1.0324,1.2462,1.4482,1.6594,1.88514,2.057810,2.277810,
    2.466310, 2.67500,2.887,3.0967,3.2990,3.4960,3.704,3.902,4.116,4.3298])

freq_2= np.array([-5.862, -5.553, -5.244, -4.935, -4.62,
       -4.3175, -4.01, -3.6996, -3.391, -3.08,
       -2.773, -2.46385632, -2.15491827, -1.846, -1.54,
       -1.228, -0.9192, -0.610, -0.301,  0.0076481 ,.3096,.6243,.940,1.244,1.550,1.8716,2.173,2.491,2.780,3.090,3.418,3.703,4.017,
    4.335,4.650,4.955,5.236,5.561,5.870,6.19,6.5168])
current=np.arange(-2,2.1, .1)

m1, yint1 = np.polyfit(current, freq_1,1)
m2, yint2=np.polyfit(current, freq_2,1)
best1 = current*m1+yint1
best2=current*m2 +yint2
def data():
	plt.plot(current, freq_1,'o',ms=20,label=r'$^{85}Rb$')
	plt.plot(current, freq_2,'*', ms=20,label=r'$^{87}Rb$')
	plt.plot(current, best1,'r')#, label='bestfit1')
	plt.plot(current, best2,'g')#, label='bestfit2')
	plt.legend(loc=2, fontsize=65)
	plt.tick_params(axis='both', which='major', labelsize=45)
	plt.tick_params(axis='both', which='minor', labelsize=15)
	plt.xlabel('Current [A]',size=50)
	plt.ylabel('Resonant Frequency [MHz]', size=50)
	plt.xlim([-2.25,2.25])

	plt.show()
def resid():
	plt.plot(current, freq_1-best1,'*',label=' resid')
	plt.plot(current, freq_2-best2,'o', label='fit2 resid')
	plt.legend(loc=2, fontsize=25)
	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.tick_params(axis='both', which='minor', labelsize=15)
	plt.xlabel('Current [A]',size=40)
	plt.ylabel(r'$\Delta$ in Resonant Frequency [MHz]', size=40)
	plt.show()
print np.std(freq_1-best1,ddof=1), np.mean(freq_1-best1)
print np.std(freq_2-best2,ddof=1), np.mean(freq_2-best2)
def nucspin(m):
	i = 1./2 * (2.5191e-2*135/(m*.275) - 1)
	return i
def earthfield(b,i):
	return b*(2*i+1)/(2.799)
def earthfield2(nuplus, numinus, I):
	return  1./2 *( 2*I+1)*(nuplus-numinus)/27.99
	
print nucspin(m1)
print nucspin(m2)
earth1 =  earthfield(yint1, 5./2) 
earth2 = earthfield(yint2, 3./2)
earthf1 = earthfield2(freq_1[21:],freq_1[0:20], 5./2)
earthf2 = earthfield2(freq_2[21:],freq_2[0:20], 3./2)
print np.mean([earthf1,earthf2])
print np.std([earthf1, earthf2],ddof=2)
print np.std([earthf1, earthf2],ddof=2)/np.sqrt(40)
# forty measurements, 2 degrees of freedom
# mean of each

