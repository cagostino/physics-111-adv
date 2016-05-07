import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
def noise_floor():
    img_low = np.genfromtxt('AFM_scan_2016-04-20_15.12.56, 0um, 512 pixels, Left.wsf', skiprows =31,skip_footer=1)
    img_hi = np.genfromtxt('AFM_scan_2016-04-20_15.23.23, 0um, 512 pixels, Left.wsf', skiprows =31,skip_footer=1)
    fig=plt.figure()
    tim = np.linspace(0,1,512.)
    freq = 1./tim
    ax1 = fig.add_subplot(121)
    ax1.imshow(img_low,origin='lower', cmap='gray')
    plt.xlabel('X [pixels]',size=40)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=15)
    plt.ylabel('Y [pixels]', size=40)
    ax2=fig.add_subplot(122)
    ax2.imshow(img_hi,origin='lower', cmap='gray')
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=15)
    plt.xlabel('X [pixels]',size=40)
   # plt.colorbar(orientation ='vertical', fraction=0.07 )
    plt.show()
    noise_pick_low = img_low[256][:]
    noise_pick_hi = img_hi[256][:]
    fft_low =abs(np.fft.fft(noise_pick_low))
    fft_hi =abs(np.fft.fft(noise_pick_hi))
    fig=plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(tim,noise_pick_low,label='Low Noise')
    ax1.plot(tim,noise_pick_hi,label='High Noise')
    plt.xlabel('Time [s]',size=40)
    plt.ylabel('Z height', size=40)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=15)
    plt.legend(fontsize=20)
    plt.show()
    ax2 = fig.add_subplot(122)
    ax2.plot(freq,fft_low,label='FFT Low')
    ax2.plot(freq,fft_hi,label='FFT High')
    plt.xlabel('Frequency (Hz)',size=40)
    plt.ylabel('Amplitude',size=40)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=15)
    plt.xlim([0,300])
    plt.ylim([0,400])
    plt.legend(fontsize=20)
    plt.show()
def force_d():
    lstfils = os.listdir(".")
    up_curve =np.zeros(485)
    down_curve= np.zeros(485)
    for i in lstfils:
        if i.endswith(".csv"):   
            data = np.transpose(np.genfromtxt(i,skiprows =4,skip_footer=1,delimiter=','))
            if data.shape[1]==485:
                up_curve +=data[0]
                down_curve+=data[1]
    
    up_curve/= 7.
    down_curve /=7.
    dist = np.linspace(0, 1000, 485.)
    lin_up = up_curve[:65]
    lin_down = down_curve[:106]
    m1,y1 = np.polyfit(dist[:65], lin_up,1)
    m2,y2=np.polyfit(dist[:106], lin_down,1)
    plt.plot(dist,up_curve,label='Up')
    plt.plot(dist,down_curve,label='Down')
    plt.xlabel('Distance (nm)',size=40)
    print np.median(down_curve)-min(down_curve),'adhesive'
    plt.ylabel('Force (nN)',size=40)
    plt.legend(fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=15)
    plt.show()
    fig = plt.figure()
    plt.plot(dist[:65],up_curve[:65],'o',label='Up')
    plt.plot(dist[:106],down_curve[:106],'*',label='Down')
    plt.plot(dist[:65],m1*dist[:65] +y1,':', label='Best Fit Up')
    plt.plot(dist[:106],m2*dist[:106] +y2,':', label='Best Fit Down')
    print m1,y1,'m1,y1'
    print m2,y2,'my,y2'
    plt.xlabel('Distance (nm)',size=40)
    print np.median(down_curve)-min(down_curve),'adhesive'
    plt.ylabel('Force (nN)',size=40)
    plt.legend(fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=15)
    plt.show()
def thermal_noise(T):
    dat = np.loadtxt('thermal_noise.txt')
    zsquare = np.mean(dat**2)
    return 1.38e-23*T/(zsquare)
    