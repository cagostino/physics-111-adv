import numpy as np
import matplotlib.pyplot as plt
def ivcurv():
    press10_i =np.array([5,6,7,8])
    
    press10_v1=np.array([56,56,54,53])
    press10_v2 =np.array([10,12,11,10.5])
    
    press12_v1=np.array([60,63,64,58])
    press12_v2 =np.array([12,13,14,12])
    
    press14_v1=np.array([67,67,65,69.5])
    press14_v2 =np.array([14,14,15,18])
    
    press16_v1=np.array([71,67,71,70])
    press16_v2 =np.array([18,18,20,21])
    
    press18_v1=np.array([73,72,71,75])
    press18_v2 =np.array([24.5,24,21,22])
    plt.plot(press10_i, press10_v2,'o', label='10 mm',ms=20)
    plt.plot(press10_i, press12_v2, '*',label='12 mm',ms=20)
    
    plt.plot(press10_i, press14_v2,'>',label='14 mm',ms=20)
    
    plt.plot(press10_i, press16_v2,'<',label='16 mm',ms=20)
    
    plt.plot(press10_i, press18_v2,'s', label='18 mm',ms=20)
    for i in range(len(press14_v1)):
        plt.errorbar(press10_i[i], press10_v2[i],xerr=.1, yerr=.5)
        plt.errorbar(press10_i[i], press12_v2[i],xerr=.1, yerr=.5)
        plt.errorbar(press10_i[i], press14_v2[i],xerr=.1, yerr=.5)
        plt.errorbar(press10_i[i], press16_v2[i],xerr=.1, yerr=.5)
        plt.errorbar(press10_i[i], press18_v2[i],xerr=.1, yerr=.5)

    plt.xlabel('Current [mA]',size=40)
    plt.ylabel('Voltage [kV]',size=40)
    plt.title('I-V Curve and Gas Pressure',size=40)
    
    plt.tick_params(axis='both', which='major', labelsize=40)
    plt.tick_params(axis='both', which='minor', labelsize=30)
    plt.xlim([4.5, 8.5])
    plt.legend()
    plt.show()
def vary():
    variability= np.array([95,92,91,96,98,92,96,92,97,94,91,97,92])/100.
    time = np.arange(0, 10*len(variability),10)
    plt.tick_params(axis='both', which='major', labelsize=40)
    plt.tick_params(axis='both', which='minor', labelsize=30)
    plt.axhline(y=np.mean(variability),ls='--')
    print np.mean(variability)
    print np.std(variability)
    for i in range(len(variability)):
        plt.errorbar(time[i], variability[i], yerr=.005,xerr=.001)
    plt.xlabel('Time [s]',size=40)
    plt.ylabel('Output Power [W]',size=40)
    plt.title('Output Variability',size=40)
    plt.plot(time, variability,':.')
    plt.show()
vary()