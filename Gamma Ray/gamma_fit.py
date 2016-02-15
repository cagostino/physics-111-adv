import numpy as np
import matplotlib.pyplot as plt
import os
plt.rc('font',family='serif')
plt.rc('text',usetex=True)
lstfiles = os.listdir('.')
lstfiles.sort()

spec_Na =np.transpose(np.loadtxt('02_03_2016_15_19_04.dat',skiprows=2,usecols=(0,1)))
spec_Cs =np.transpose(np.loadtxt('02_03_2016_15_25_15.dat',skiprows=2,usecols=(0,1)))
spec_Co = np.transpose(np.loadtxt('02_03_2016_15_29_31.dat',skiprows=2, usecols=(0,1)))
def standdev(data):
	return np.sqrt( np.sum( (data-np.mean(data) )**2)/np.float(len(data)-1.))

def linleastsquares_general(data,order): #many orders!
	n  = len(data[0])
	x= data[0]
	y=data[1]
	nmatleft= np.empty( (order+1,order+1))
	
	nmatright = np.empty( (order+1,1))
	for i in range(order+1):
		for j in range(order+1):
			if i == 0 and j==0:
				nmatleft[0][0] = n
			else:
				nmatleft[i][j] = np.sum(x**(i+j))
		nmatright[i][0] = np.sum(x**i * y)
	inv_left = np.linalg.inv(nmatleft)
	final=np.dot(inv_left,nmatright)
	
	return final
def findpeaks(spectrum,limit,width,typ,cent): #spectrum is 1 two-dimensional array so spectrum[0] is the x-array, spectrum[1] y array
	spec_x=spectrum[0]
	spec_y = spectrum[1]
	peaks_x = []
	peaks_y = []
	for i in range(len(spec_x)-1): # -2 because itll b
		if all(spec_y[i] > spec_y[i-width:i]) and all(spec_y[i] >spec_y[i+1:i+width]) and spec_y[i] > limit:
			peaks_x.append(np.where(spec_x == spec_x[i])[0])	
			peaks_y.append(spec_y[i])
	centroids= []	
	errors=[]
	for peak in peaks_x: #finds centroids of peaks
		xran= spec_x[peak-5:peak+5]
		
		intens_ran= spec_y[peak-5:peak+5]
		numer=[]		
		centroid = np.sum(np.array(xran*intens_ran))/np.sum(intens_ran)
		
		centroids.append(centroid)
		#errors.append( standdev(xran)**2/(np.sum(intens_ran) ))
	if cent and centroids != []:
		for centroid in centroids: #plots centroids 
			plt.axvline(x=centroid,color='r',linestyle='--')
		plt.axvline(centroids[0],color='r',linestyle='--',label=typ+' Peak Centroid')
	
	
		plt.plot(spec_x,spec_y,linewidth=2,label=typ)

		plt.tick_params(axis='both', which='major', labelsize=20)
		plt.tick_params(axis='both', which='minor', labelsize=15)
		plt.show()

	return np.array(centroids)


cent_Na = findpeaks(spec_Na,150,10,'Neon', False)
cent_Cs = findpeaks(spec_Cs,500,10,'Cesium', False)
cent_Co = findpeaks(spec_Co,1000,10,'Cobalt', False)


arr_cent = np.array([3.0613289748953973, 4.0025605026929982, 3.4341010115135413,
       3.9516175791641639, 4.0530974946344882])
def three_spec():
	fig = plt.figure()
	ax1 = fig.add_subplot(131)
	plt.axvline(x=cent_Na[0], color='r',linestyle='--',label='Sodium-22 Peaks')
	plt.axvline(x=cent_Na[1], color='r',linestyle='--')
	plt.plot(spec_Na[0],spec_Na[1], linewidth=2, label='Sodium-22')
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.ylabel('Intensity[Counts]',size=40)
	plt.tick_params(axis='both', which='minor', labelsize=25)
	plt.legend(loc=2, fontsize=19)
	ax2 = fig.add_subplot(132)
	plt.axvline(x=cent_Cs[0], color='r',linestyle='--',label='Cesium-137 Peak')
	plt.plot(spec_Cs[0],spec_Cs[1], linewidth=2, label='Cesium-137')
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)
	plt.xlabel('Energy',size=40)
	plt.legend(fontsize=20)
	ax3= fig.add_subplot(133)
	plt.axvline(x=cent_Co[0], color='r',linestyle='--',label='Cobalt-60 Peaks')
	plt.axvline(x=cent_Co[1], color='r',linestyle='--')
	plt.plot(spec_Co[0],spec_Co[1], linewidth=2, label='Cobalt-60')
	plt.legend(loc=2, fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)
	

	plt.show()
#three_spec()
arr_cent.sort()
arr_energies = np.array([0.511,1.28,0.6616,1.17,1.33])
arr_energies.sort()
cfit, mfit = linleastsquares_general([arr_cent, arr_energies], 1)
quad_0, quad_1, quad_2= linleastsquares_general([arr_cent,arr_energies],2)
print cfit, mfit, 'linear'
print quad_0, quad_1, quad_2,'quadratic'
ran = np.linspace(min(arr_cent), max(arr_cent), 100)
linfit=mfit*arr_cent+cfit
quadfit = quad_0+quad_1*arr_cent+quad_2*arr_cent**2
linfitpl =  mfit*ran + cfit
quadfitpl = quad_0 + quad_1*ran + quad_2 * ran**2
def plot_fits():
	fig =plt.figure()
	ax1= fig.add_subplot(121)
	plt.plot(arr_cent, arr_energies,'go',ms=10,label='Data')
	plt.ylabel('Energy(MeV)', size=40)
	plt.xlabel('Channel', size=40)
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)
	plt.plot(ran, linfitpl,'r', linewidth=3,label='Linear Fit')
	plt.plot(ran, quadfitpl,'b', linewidth=3,label='Quad Fit')
	plt.legend(loc=2, fontsize=20)
	ax2=fig.add_subplot(122)
	plt.plot(arr_cent, linfit-arr_energies,'ro', ms=10,label='Linear Residuals')
	plt.plot(arr_cent, quadfit-arr_energies,'b*', ms=10, label='Quadratic Residuals')
	plt.ylabel(r'$\Delta$ Energy(MeV)', size=40)
	plt.xlabel('Channel', size=40)
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)

	plt.legend(fontsize=20)	
	plt.show()
#plot_fits()
def ideal_quadfn(x,quad_2,quad_1,quad_0):
	return quad_2*x**2+quad_1*x+quad_0
def ideal_lin(x, m, c):
	return x*m + c

#spec_Na[0] = ideal_quadfn(spec_Na[0],quad_2,quad_1, quad_0)

spec_Na_use = np.where(spec_Na[0] > 2.95)[0]
spec_Na_x = spec_Na[0][spec_Na_use]
spec_Na_x = ideal_quadfn(spec_Na_x, quad_2,quad_1,quad_0)
cent_Na = ideal_quadfn(cent_Na,quad_2,quad_1, quad_0)

spec_co_x = spec_Co[0][np.where(spec_Co[0] > 3 ) ]
spec_co_y = spec_Co[1][np.where(spec_Co[0] > 3 ) ]
spec_co_x = ideal_quadfn(spec_co_x, quad_2,quad_1,quad_0)
cent_co = ideal_quadfn(cent_Co,quad_2,quad_1, quad_0)

spec_cs_x = spec_Cs[0][np.where(spec_Cs[0] >3 )]
spec_cs_y = spec_Cs[1][np.where(spec_Cs[0] >3)]

spec_cs_x = ideal_quadfn(spec_cs_x,quad_2,quad_1,quad_0)
cent_Cs = ideal_quadfn(cent_Cs,quad_2,quad_1, quad_0)

#plotting
def plot_cal():
	fig = plt.figure()
	ax1= fig.add_subplot(131)
	plt.plot(spec_Na_x,spec_Na[1][spec_Na_use],label='Sodium')
	plt.axvline(x=cent_Na[0], color='r',linestyle='--',label='Sodium-22 Peaks')
	plt.axvline(x=cent_Na[1], color='r',linestyle='--')
	plt.legend(loc=2, fontsize=20)
	plt.ylabel('Intensity', size=40)
	#plt.xlabel('Energy[MeV]', size=40)
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)

	plt.xlim([0,1.5])
	ax2 = fig.add_subplot(132)
	plt.plot(spec_cs_x,spec_cs_y,label='Cesium')
	#plt.xlabel('Energy [MeV]', size=40)
	#plt.ylabel('Intensity', size=40)
	plt.axvline(x=cent_Cs[0], color='r',linestyle='--',label='Cesium-137 Peak')

	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)
	plt.xlim([.5,1])

	plt.legend(loc=1, fontsize=20)

	
	ax3 = fig.add_subplot(133)
	plt.plot(spec_co_x,spec_co_y,label='Cobalt')
	plt.xlabel('Energy [MeV]', size=40)
	#plt.ylabel('Intensity', size=40)
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)
	plt.axvline(x=cent_co[0], color='r',linestyle='--',label='Cobalt-60 Peaks')
	plt.axvline(x=cent_co[1], color='r',linestyle='--')

	plt.legend(loc=1, fontsize=20)

	plt.xlim([.5, 2])

	plt.show()
#plot_cal()


def plot_inv_sq():
	lstfil = ['02_03_2016_15_40_37.dat','02_03_2016_15_48_35.dat',
	'02_03_2016_16_04_25.dat']#,'02_03_2016_16_24_24.dat']
	specs= [np.transpose(np.loadtxt(i,skiprows=2, usecols=(0,1))) for i in lstfil ]
	cent = [findpeaks(specs[i],100, 10,'Sodium-22',True) for i in range(len(lstfil))]
	xr=[26.9, 53.8,  79.3]
	return cent
cen = plot_inv_sq()
def plot_gain():
	varied_voltages= lstfiles[31:73]
	specs = [np.transpose(np.loadtxt(varied_voltages[i]
		,skiprows=2,usecols=(0,1))) for i in range(len(varied_voltages)-1)]
	peaks = np.array([])
	for i in range(len(varied_voltages)-1):
		print i
		peak_spot = np.where(specs[i][1] == max(specs[i][1]))[0][0]
		peak_x = specs[i][0][peak_spot]
		peaks = np.append(peaks, peak_x)

	xr= np.arange(400,810,10 )
	xr =np.delete(xr, 6)
	peaks = np.delete(peaks,6)
	plt.plot(xr, peaks,'o', ms=10,label='Peak Position of Cesium-137')
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=25)
	plt.xlabel('Voltage',size=40)
	plt.ylabel('Peak Position',size=40)
	plt.legend(loc=2,fontsize=25)

	quad_0,quad_1, quad_2=linleastsquares_general([xr[15:],peaks[15:]],2)
	peak_fit = ideal_quadfn(xr, quad_2,quad_1,quad_0)
	plt.plot(xr,peak_fit,'--', label='fit')
	plt.show()
	return peaks

#peakdat = plot_gain()
