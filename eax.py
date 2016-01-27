import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy.stats import norm
def linleastsquares_general(data,order): #many orders!
    n = len(data[0])
    x = data[0]
    y=data[1]
    nmatleft=np.empty( ( order+1,order+1))
    nmatright = np.empty( ( order+1,1))
    for i in range(order+1):
        for j in range(order+1):
            if i ==0 and j==0:
                nmatleft[0][0]=n
            else:
                nmatleft[i][j] = np.sum(x**(i+j))
        nmatright[i][0]=np.sum(x**i*y)
    inv_left=np.linalg.inv(nmatleft)
    final=np.dot(inv_left,nmatright)
    return final

def func3b(n):
    
    arr = np.random.normal(0,1,n)
    mean = np.mean(arr)
    standdev = np.std(arr,ddof=1)
    error = standdev/np.sqrt(n)
    return [arr,mean,standdev,error]
def func3c(n):
    arrs = []
    means = []
    standdevs=[]
    errors=[]
    for i in range(1000):
        arr, mean, standdev, error = func3b(n)
        arrs.append(arr)
        means.append(mean)
        standdevs.append(standdev)
        errors.append(error)
    #plt.hist(arrs,20)
    #plt.show()
    meanofmeans = np.mean(means)
    stdmeans=np.std(means,ddof=1)
    fig = plt.figure()
    ax1 = fig.add_subplot(131)
    plt.hist(means,30,label=r'$\mu$ = '+str(meanofmeans)[0:6] +r' $\sigma$ = '+ str(stdmeans)[0:6])
    plt.xlabel('Means',size=40)
    plt.ylabel('Counts', size=40)

    plt.tick_params(axis='both', which='major', labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30) 
    #plt.legend(fontsize=15)
    #plt.show()
    print meanofmeans, stdmeans
    meanofstd= np.mean(standdevs)
    stdstd = np.std(standdevs,ddof=1)
    ax2=fig.add_subplot(132)
    plt.hist(standdevs,30,label=r'$\mu$ = ' +str(meanofstd)[0:6]+r' $\sigma$= ' +str(stdmeans)[0:6])
    plt.xlabel('Standard Deviations', size=40)
    #plt.ylabel('Coutns', size=40)
    #plt.legend(fontsize=15)
    plt.tick_params(axis='both',which='major',labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30)
    #plt.show()
    print meanofstd, stdstd
    ax3=fig.add_subplot(133) 
    stderr = np.std(errors,ddof=1)
    meanoferr=np.mean(errors)
    plt.hist(errors,30,label=r'$\mu$ = ' +str(meanoferr) +r' $\sigma$= '+str(stderr)[0:6])
    plt.xlabel('Errors',size=40)
    #plt.ylabel('Counts',size=40)
    plt.tick_params(axis='both', which='major', labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30)
    #plt.legend(fontsize=15)
    print meanoferr, stderr
    plt.show()
def func4(n):
    arr = np.random.exponential(1,n)
    mean = np.mean(arr)
    standdev = np.std(arr,ddof=1)
    error = standdev/np.sqrt(n)
    return [arr,mean,standdev,error]
def func4b(n):
    arrs = []
    means = []
    standdevs=[]
    errors=[]
    for i in range(1000):
        arr, mean, standdev, error = func4(n)
        arrs.append(arr)
        means.append(mean)
        standdevs.append(standdev)
        errors.append(error)
    #plt.hist(arrs,20)
    #plt.show()
    meanofmeans = np.mean(means)
    stdmeans=np.std(means,ddof=1)
    fig = plt.figure()
    ax1 = fig.add_subplot(131)
    plt.hist(means,30,label=r'$\mu$ = '+str(meanofmeans)[0:6] +r' $\sigma$ = '+ str(stdmeans)[0:6])
    plt.xlabel('Means',size=40)
    plt.ylabel('Counts', size=40)

    plt.tick_params(axis='both', which='major', labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30) 
    #plt.legend(fontsize=15)
    #plt.show()
    print meanofmeans, stdmeans
    meanofstd= np.mean(standdevs)
    stdstd = np.std(standdevs,ddof=1)
    ax2=fig.add_subplot(132)
    plt.hist(standdevs,30,label=r'$\mu$ = ' +str(meanofstd)[0:6]+r' $\sigma$= ' +str(stdmeans)[0:6])
    plt.xlabel('Standard Deviations', size=40)
    #plt.ylabel('Coutns', size=40)
    #plt.legend(fontsize=15)
    plt.tick_params(axis='both',which='major',labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30)
    #plt.show()
    print meanofstd, stdstd
    ax3=fig.add_subplot(133) 
    stderr = np.std(errors,ddof=1)
    meanoferr=np.mean(errors)
    plt.hist(errors,30,label=r'$\mu$ = ' +str(meanoferr) +r' $\sigma$= '+str(stderr)[0:6])
    plt.xlabel('Errors',size=40)
    #plt.ylabel('Counts',size=40)
    plt.tick_params(axis='both', which='major', labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30)
    #plt.legend(fontsize=15)
    print meanoferr, stderr
    plt.show()
data5 = np.loadtxt('peak.dat')
def func5():
    mu,sigma = norm.fit(data5)
    n,bins,patches=plt.hist(data5,60,normed=1, facecolor='green', alpha=0.75)
    y = mlab.normpdf(bins,mu,sigma)
    l= plt.plot(bins,y,'r--', linewidth=2)
    plt.tick_params(axis='both', which='major', labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30)
    plt.show()
def func6():
    current=np.array([0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2])
    freq=np.array([.14,0.60,1.21,1.94,2.47,3.07,3.83,4.16,4.68,5.60,6.31,6.78])
    plt.plot(current, freq,'o')
    b, m = linleastsquares_general([current,freq],1)
      
    linfit = m*current + b
    plt.tick_params(axis='both', which='major', labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=30)

    plt.plot(current, linfit,'b-')
    plt.xlabel('Current',size=40)
    plt.ylabel('Frequency',size=40)
    plt.show()
