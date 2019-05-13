from scipy.stats import uniform
import matplotlib.pyplot as plt
import os
import numpy as np
import peakutils
from numpy import log,exp,pi
from numpy import array
from numpy import cov
from numpy import linalg as LA
from numpy.linalg import inv
from numpy.linalg import det
from mpl_toolkits import mplot3d
import pandas as pd
import corner
import h5py
from chainconsumer import ChainConsumer






def pxd(samples,means,covs,include_prior_rej):
        sigma_i=inv(covs)
#        x=samples-means
        x2=[]
        for m in range(0,len(samples),1):
            if m in [2,5,6] and include_prior_rej == 'T':
#                print(m,samples[m][0])
                if m == 2:
                    tmp_x=min([np.absolute(samples[m][0]-means[m][0]),np.absolute(samples[m][0]-2*np.arccos(-1)-means[m][0]),np.absolute(samples[m][0]+2*np.arccos(-1)-means[m][0])])
                    if tmp_x==np.absolute(samples[m][0]-means[m][0]):
                        x2.append([samples[m][0]-means[m][0]])
                    elif tmp_x== np.absolute(samples[m][0]-2*np.arccos(-1)-means[m][0]):
                        x2.append([samples[m][0]-2*np.arccos(-1)-means[m][0]])
                    else:
                        x2.append([samples[m][0]+2*np.arccos(-1)-means[m][0]])
                if m == 5:
                    tmp_x=min([np.absolute(samples[m][0]-means[m][0]),np.absolute(samples[m][0]-np.arccos(-1)-means[m][0]),np.absolute(samples[m][0]+np.arccos(-1)-means[m][0])])
                    if tmp_x==np.absolute(samples[m][0]-means[m][0]):
                        x2.append([samples[m][0]-means[m][0]])
                    elif tmp_x== np.absolute(samples[m][0]-np.arccos(-1)-means[m][0]):
                            x2.append([samples[m][0]-np.arccos(-1)-means[m][0]])
                    else:
                        x2.append([samples[m][0]+np.arccos(-1)-means[m][0]])
                if m == 6:
                    tmp_x=min([np.absolute(samples[m][0]-means[m][0]),np.absolute(samples[m][0]-2*np.arccos(-1)-means[m][0]),np.absolute(samples[m][0]+2*np.arccos(-1)-means[m][0])])
                    if tmp_x==np.absolute(samples[m][0]-means[m][0]):
                        x2.append([samples[m][0]-means[m][0]])
                    elif tmp_x== np.absolute(samples[m][0]-2*np.arccos(-1)-means[m][0]):
                        x2.append([samples[m][0]-2*np.arccos(-1)-means[m][0]])
                    else:
                        x2.append([samples[m][0]+2*np.arccos(-1)-means[m][0]])
            else:
                x2.append([samples[m][0]-means[m][0]])
#        print(x2)
        x=np.asarray(x2)
        p=[]
        for k in range(0,len(x[0]),1):
                sample_k=[x[0][k],x[1][k],x[2][k],x[3][k],x[4][k],x[5][k],x[6][k],x[7][k]]
                v1=np.dot(sigma_i,sample_k)
                v=np.dot(sample_k,v1)
                p.append(np.exp(-0.5*v)/(np.sqrt((2*np.pi)**2*det(covs))))
#        print(v,x2[0],x2[1],x2[2],x2[3],x2[4],x2[5],x2[6],x2[7])
#        print(sigma_i[0][0],sigma_i[0][1],sigma_i[0][2],sigma_i[0][3],sigma_i[0][4],sigma_i[0][5],sigma_i[0][6],sigma_i[0][7])
        return p

#def pxd(samples,means,covs):
#        sigma_i=inv(covs)
#        x=np.asarray(samples)-np.asarray(mean1)
#        sample_k=x
#        v1=np.dot(sigma_i,sample_k)
#        v=np.dot(sample_k.T,v1)
#        p=np.exp(-0.5*v)/(np.sqrt((2*np.pi)**2*det(covs)))
#        return p


def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

f = h5py.File('/Users/klackeos/data/LDC1-3_VGB_v1_FD_noiseless.hdf5', 'r') 
amp_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['Amplitude']: 
        amp_ldc.append(i)

ecl_lat_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['EclipticLatitude']: 
        ecl_lat_ldc.append(np.cos(np.pi/2.0-i))

ecl_long_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['EclipticLongitude']: 
        ecl_long_ldc.append(i)
        for kk in range(0,len(ecl_long_ldc),1):
                if ecl_long_ldc[kk] < 0:
                        ecl_long_ldc[kk]+=2*np.pi

freq_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['Frequency']: 
        freq_ldc.append(i)

dot_freq_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['FrequencyDerivative']: 
        dot_freq_ldc.append(i)

inc_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['Inclination']: 
        inc_ldc.append(np.cos(i))
print(inc_ldc)
phi0_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['InitialPhase']: 
        phi0_ldc.append(i)

pol_ldc=[] 
for i in f['H5LISA']['GWSources']['GalBinaries']['Polarization']: 
        pol_ldc.append(i)
        for kk in range(0,len(pol_ldc),1):
                if pol_ldc[kk] > np.pi:
                        pol_ldc[kk]-=np.pi
f2=open('/Users/klackeos/data/LDC1-3_VGB_v1_3_months.dat','r')
data=np.loadtxt(f2)
freq=data[:,0]


print("Enter observation time (1.5 mo = 3932160, 3 mo=7864320, 6 mo=15728640, 1 yr=31457280, 2 yr=62914560):")
T=input()
T=int(T)
print("Sample from prior? Enter T or F:")
prior_tf=input()
include_prior_rej="T"
if prior_tf == "F":
    print("Include prior rejection? Enter T or F:")
    include_prior_rej=input()

if T==3932160:
    scale=32
elif T==7864320:
    scale=64
elif T==15728640:
    scale=128
elif T==31457280:
    scale=256
else:
    print("Enter scale factor:")
    scale=input()
    scale=int(scale)


#read in chains
#f2=open('/Users/klackeos/Desktop/dimension_chain_new1.dat.1','r')
f2=open('/Users/klackeos/ldasoft/master/bin/ldc_vgb_noise/dft/source4/18mo_test_cov/chains/dimension_chain.dat.1','r')
samples=np.loadtxt(f2)
f0=samples[:,0]*T
dfdt=samples[:,1]*T**2
amp=np.log(samples[:,2])
phi=samples[:,3]
costheta=samples[:,4]
cosi=samples[:,5]
psi=samples[:,6]
phi0=samples[:,7]

print("Enter fmin:")
fmin=input()
prior00 = np.int(float(fmin)*T)
prior01 = prior00+int(scale)

#//colatitude
prior10 = -1.0
prior11 =  1.0

#//longitude
prior20 = 0.0
prior21 = 2*np.arccos(-1)

#//log amplitude
prior30 = -55
prior31 = -45

#//cos inclination
prior40 = -1.0
prior41 =  1.0

#//polarization
prior50 = 0.0
prior51 = np.arccos(-1)

#//phase
prior60 = 0.0
prior61 = 2*np.arccos(-1)

#//fdot

fdotmin = -0.000005*(prior00/T)**(13./3.)
fdotmax = 0.0000008*(prior01/T)**(11./3.)
prior70 = (fdotmin*T**2)
prior71 = (fdotmax*T**2)

prior=np.array([[prior00, prior01], [prior10, prior11], [prior20, prior21], [prior30, prior31], [prior40, prior41], [prior50, prior51], [prior60, prior61], [prior70, prior71]])
#prior=np.array([[prior00, prior01],[prior70, prior71],[prior30, prior31],[prior20, prior21],[prior10, prior11],[prior40, prior41],[prior50, prior51],[prior60, prior61]])
phi0_hist=plt.hist(phi0[phi0<np.pi],bins=200);
plt.xlabel('phi0')
plt.show()
phi0_hist_freq=list(phi0_hist[0])
max_freq = max(phi0_hist_freq)
max_index = phi0_hist_freq.index(max_freq)
phi0_hist_phi0=list(phi0_hist[1])   
phi0_peak_loc1=phi0_hist_phi0[max_index]  
plt.hist(phi0[phi0<np.pi],bins=200);
plt.axvline(x=phi0_peak_loc1,color='k')
plt.xlabel('phi0')
plt.show()

phi0_peak_loc2=phi0_peak_loc1+np.pi

plt.hist(phi0,bins=200);
plt.axvline(x=phi0_peak_loc1,color='k')
plt.axvline(x=phi0_peak_loc2,color='k')
plt.xlabel('phi0')
plt.show()

psi_hist=plt.hist(psi[psi<np.pi/2],bins=200)
plt.xlabel('psi')
plt.show()
psi_hist_freq=list(psi_hist[0])
max_freq1 = max(psi_hist_freq)
max_index = psi_hist_freq.index(max_freq1)
psi_hist_psi=list(psi_hist[1])
psi_peak_loc1a=psi_hist_psi[max_index]
plt.hist(psi[psi<np.pi],bins=200);
plt.axvline(x=psi_peak_loc1a,color='k')
plt.xlabel('psi')
plt.show()

psi_hist=plt.hist(psi[psi>np.pi/2],bins=200)
plt.xlabel('psi')
plt.show()
psi_hist_freq=list(psi_hist[0])
max_freq2 = max(psi_hist_freq)
max_index = psi_hist_freq.index(max_freq2)
psi_hist_psi=list(psi_hist[1])
psi_peak_loc1b=psi_hist_psi[max_index]
plt.hist(psi[psi<np.pi],bins=200);
plt.axvline(x=psi_peak_loc1b,color='k')
plt.xlabel('psi')
plt.show()
if max_freq2 > max_freq:
        psi_peak_loc2=psi_peak_loc1b-np.pi/2
        psi_peak_loc1=psi_peak_loc1b
else:
        psi_peak_loc2=psi_peak_loc1a+np.pi/2
        psi_peak_loc1=psi_peak_loc1a

plt.hist(psi,bins=200);
plt.axvline(x=psi_peak_loc1,color='k')
plt.axvline(x=psi_peak_loc2,color='k')
plt.xlabel('psi')
plt.show()

if phi0_peak_loc1-np.pi/2 < 0:  
        phi0_min_angle_A=phi0_peak_loc1+2*np.pi-np.pi/2
        phi0_max_angle_A=phi0_peak_loc1+2*np.pi+np.pi/2
        phi0_peak_locA=phi0_peak_loc1+2*np.pi
else:
        phi0_min_angle_A=phi0_peak_loc1-np.pi/2
        phi0_max_angle_A=phi0_peak_loc1+np.pi/2
        phi0_peak_locA=phi0_peak_loc1


if phi0_peak_loc2-np.pi/2 < 0:  
        phi0_min_angle_B=phi0_peak_loc2+2*np.pi-np.pi/2
        phi0_max_angle_B=phi0_peak_loc2+2*np.pi+np.pi/2
        phi0_peak_locB=phi0_peak_loc2+2*np.pi
else:
        phi0_min_angle_B=phi0_peak_loc2-np.pi/2
        phi0_max_angle_B=phi0_peak_loc2+np.pi/2
        phi0_peak_locB=phi0_peak_loc2



if phi0_peak_locA > phi0_peak_locB:
        phi0_peak_loc=phi0_peak_locA
        phi0_min_angle=phi0_min_angle_A
        phi0_max_angle=phi0_max_angle_A
else:
        phi0_peak_loc=phi0_peak_locB
        phi0_min_angle=phi0_min_angle_B
        phi0_max_angle=phi0_max_angle_B




if psi_peak_loc1-np.pi/4 < 0:
        psi_min_angle_A=psi_peak_loc1+np.pi-np.pi/4
        psi_max_angle_A=psi_peak_loc1+np.pi+np.pi/4
        psi_peak_locA=psi_peak_loc1+np.pi
else:
        psi_min_angle_A=psi_peak_loc1-np.pi/4
        psi_max_angle_A=psi_peak_loc1+np.pi/4
        psi_peak_locA=psi_peak_loc1


if psi_peak_loc2-np.pi/4 < 0:
        psi_min_angle_B=psi_peak_loc2+np.pi-np.pi/4
        psi_max_angle_B=psi_peak_loc2+np.pi+np.pi/4
        psi_peak_locB=psi_peak_loc2+np.pi
else:
        psi_min_angle_B=psi_peak_loc2-np.pi/4
        psi_max_angle_B=psi_peak_loc2+np.pi/4
        psi_peak_locB=psi_peak_loc2


if psi_peak_locA > psi_peak_locB:
        psi_peak_loc=psi_peak_locA
        psi_min_angle=psi_min_angle_A
        psi_max_angle=psi_max_angle_A
else:
        psi_peak_loc=psi_peak_locB
        psi_min_angle=psi_min_angle_B
        psi_max_angle=psi_max_angle_B

angles=[0,np.pi/2,np.pi,3*np.pi/2,2*np.pi]
nearest_angle_index=(np.abs(phi0_peak_loc-angles)).argmin()
angles_nearest=angles[nearest_angle_index]


f0_a=[]
f0_b=[]
dfdt_a=[]
dfdt_b=[]
amp_a=[]
amp_b=[]
phi_a=[]
phi_b=[]
costheta_a=[]
costheta_b=[]
cosi_a=[]
cosi_b=[]
phi0_a=[]
phi0_b=[]
psi_a=[]
psi_b=[]


phi0_c=[]



if phi0_peak_loc < np.pi/6:
        phi0_peak_loc = 2*np.pi
        angles_nearest=2*np.pi
for i in range(0,len(phi0),1):
        min_angle_A=phi0_peak_loc-np.pi/2
        max_angle_A=phi0_peak_loc+np.pi/2
        if min_angle_A < 0:
                min_angle_A+=2*np.pi
                max_angle_A+=2*np.pi
                if (phi0[i] + 2*np.pi  >= min_angle_A and phi0[i] + 2*np.pi <= max_angle_A):
                        phi0_a.append(phi0[i]+2*np.pi)
                        phi0_c.append(phi0[i])
                        psi_a.append(psi[i])
                        f0_a.append(f0[i])
                        dfdt_a.append(dfdt[i])
                        amp_a.append(amp[i])
                        phi_a.append(phi[i])
                        costheta_a.append(costheta[i])
                        cosi_a.append(cosi[i])
                elif (phi0[i] >= min_angle_A and phi0[i] <= max_angle_A):
                        phi0_a.append(phi0[i])
                        phi0_c.append(phi0[i])
                        psi_a.append(psi[i])
                        f0_a.append(f0[i])
                        dfdt_a.append(dfdt[i])
                        amp_a.append(amp[i])
                        phi_a.append(phi[i])
                        costheta_a.append(costheta[i])
                        cosi_a.append(cosi[i])
                else:
                        if np.absolute(phi0[i]+2*np.pi-np.mean(phi0_a)) < np.absolute(phi0[i]+2*np.pi-np.mean(phi0_b)):
                                phi0_a.append(phi0[i]+2*np.pi)
                                phi0_c.append(phi0[i])
                                psi_a.append(psi[i])
                                f0_a.append(f0[i])
                                dfdt_a.append(dfdt[i])
                                amp_a.append(amp[i])
                                phi_a.append(phi[i])
                                costheta_a.append(costheta[i])
                                cosi_a.append(cosi[i])
                        else:
                                phi0_b.append(phi0[i])
                                psi_b.append(psi[i])
                                f0_b.append(f0[i])
                                dfdt_b.append(dfdt[i])
                                amp_b.append(amp[i])
                                phi_b.append(phi[i])
                                costheta_b.append(costheta[i])
                                cosi_b.append(cosi[i])
        elif np.round(angles_nearest,2) == 6.28 or np.round(angles_nearest,2) == 4.71:
                if (phi0[i] + 2*np.pi  >= min_angle_A and phi0[i] + 2*np.pi <= max_angle_A) :
                        phi0_a.append(phi0[i]+2*np.pi)
                        phi0_c.append(phi0[i])
                        psi_a.append(psi[i])
                        f0_a.append(f0[i])
                        dfdt_a.append(dfdt[i])
                        amp_a.append(amp[i])
                        phi_a.append(phi[i])
                        costheta_a.append(costheta[i])
                        cosi_a.append(cosi[i])
                elif (phi0[i] >= min_angle_A and phi0[i] <= max_angle_A):
                        phi0_a.append(phi0[i])
                        phi0_c.append(phi0[i])
                        psi_a.append(psi[i])
                        f0_a.append(f0[i])
                        dfdt_a.append(dfdt[i])
                        amp_a.append(amp[i])
                        phi_a.append(phi[i])
                        costheta_a.append(costheta[i])
                        cosi_a.append(cosi[i])
                else:
                        phi0_b.append(phi0[i])
                        psi_b.append(psi[i])
                        f0_b.append(f0[i])
                        dfdt_b.append(dfdt[i])
                        amp_b.append(amp[i])
                        phi_b.append(phi[i])
                        costheta_b.append(costheta[i])
                        cosi_b.append(cosi[i])
        else:
                if (phi0[i] >= min_angle_A and phi0[i] <= max_angle_A):
                        phi0_a.append(phi0[i])
                        phi0_c.append(phi0[i])
                        psi_a.append(psi[i])
                        f0_a.append(f0[i])
                        dfdt_a.append(dfdt[i])
                        amp_a.append(amp[i])
                        phi_a.append(phi[i])
                        costheta_a.append(costheta[i])
                        cosi_a.append(cosi[i])
                else:
                        phi0_b.append(phi0[i])
                        psi_b.append(psi[i])
                        f0_b.append(f0[i])
                        dfdt_b.append(dfdt[i])
                        amp_b.append(amp[i])
                        phi_b.append(phi[i])
                        costheta_b.append(costheta[i])
                        cosi_b.append(cosi[i])




if len(phi0_b)>0:
    if np.round(max(phi0_b)-min(phi0_b),2)==6.28:
            for k in range(0,len(phi0_b),1):
                    if phi0_b[k]<phi0_peak_loc-np.pi:
                            phi0_b[k]+=np.pi

if len(psi_b)>0:
    if np.round(max(psi_b)-min(psi_b),2)==3.14:
        for k in range(0,len(psi_b),1):
            if psi_b[k]<psi_peak_loc-np.pi/4:
                if psi_b[k] + np.pi < psi_peak_loc+np.pi/4:
                    psi_b[k]+=np.pi

    if np.round(max(psi_b)-min(psi_b),2) == 3.14:
        for k in range(0,len(psi_b),1):
            if psi_b[k] > np.pi/2 and psi_b[k]-np.pi>0:
                psi_b[k]-=np.pi

    if np.round(max(psi_b)-min(psi_b),2) == 3.14:
        for k in range(0,len(psi_b),1):
            if psi_b[k] < np.pi/2:
                psi_b[k]+=np.pi
    for i in range(0,len(phi_b),1):
        if np.absolute(psi_b[i]-np.mean(psi_b))> np.pi/2:
            if psi_b[i] <np.mean(psi_b):
                psi_b[i]+=np.pi

if len(phi0_a)>0:
    if len(phi0_a)>0:
            if np.round(max(phi0_a)-min(phi0_a),2)==6.28:
                    for k in range(0,len(phi0_a),1):
                            if phi0_a[k]<phi0_peak_loc-np.pi:
                                    phi0_a[k]+=np.pi





if len(psi_a)>0:
        if np.round(max(psi_a)-min(psi_a),2)==3.14:
                for k in range(0,len(psi_a),1):
                        if psi_a[k]<psi_peak_loc-np.pi/4:
                                if psi_a[k]<psi_peak_loc-np.pi/4:
                                        if psi_a[k] + np.pi < psi_peak_loc+np.pi/4:
                                                psi_a[k]+=np.pi

        if np.round(max(psi_a)-min(psi_a),2) == 3.14:
                for k in range(0,len(psi_a),1):
                        if psi_a[k] > np.pi/2 and psi_a[k]-np.pi>0:
                                psi_a[k]-=np.pi

        if np.round(max(psi_a)-min(psi_a),2) == 3.14:
                for k in range(0,len(psi_a),1):
                        if psi_a[k] < np.pi/2 :
                                psi_a[k]+=np.pi



if len(phi0_a)>0 and len(phi0_b)>0:
    f, (ax1)=plt.subplots(1,1)
    ax1.scatter(phi0,psi,s=15,c='k',marker='o')
    plt.xlabel('phi0',Fontsize=15)
    plt.ylabel('psi',Fontsize=15)
    ax1.scatter(phi0_a,psi_a,s=5,c='r',marker='o')
    plt.xlabel('phi0',Fontsize=15)
    plt.ylabel('psi',Fontsize=15)
    ax1.scatter(phi0_b,psi_b,s=5,c='g',marker='o')
    plt.xlabel('phi0',Fontsize=15)
    plt.ylabel('psi',Fontsize=15)
    plt.show()

    f, (ax1)=plt.subplots(1,1)
    ax1.scatter(phi0,f0,s=15,c='k',marker='o')
    plt.xlabel('phi0',Fontsize=15)
    plt.ylabel('f0',Fontsize=15)
    ax1.scatter(phi0_a,f0_a,s=5,c='r',marker='o')
    plt.xlabel('phi0',Fontsize=15)
    plt.ylabel('f0',Fontsize=15)
    ax1.scatter(phi0_b,f0_b,s=5,c='g',marker='o')
    plt.ylim(min(f0_b),max(f0_b))
    plt.xlabel('phi0',Fontsize=15)
    plt.ylabel('f0',Fontsize=15)
    plt.show()


print("Try again? (y/n)")
ans=input()

if ans=="y":
    phi0_hist=plt.hist(phi0[phi0<np.pi],bins=200);
    plt.xlabel('phi0')
    plt.show()
    phi0_hist_freq=list(phi0_hist[0])
    max_freq = max(phi0_hist_freq)
    max_index = phi0_hist_freq.index(max_freq)
    phi0_hist_phi0=list(phi0_hist[1])
    phi0_peak_loc1=phi0_hist_phi0[max_index]
    plt.hist(phi0[phi0<np.pi],bins=200);
    plt.axvline(x=phi0_peak_loc1,color='k')
    plt.xlabel('phi0')
    plt.show()

    phi0_peak_loc2=phi0_peak_loc1+np.pi

    plt.hist(phi0,bins=200);
    plt.axvline(x=phi0_peak_loc1,color='k')
    plt.axvline(x=phi0_peak_loc2,color='k')
    plt.xlabel('phi0')
    plt.show()

    psi_hist=plt.hist(psi[psi<np.pi/2],bins=200)
    plt.xlabel('psi')
    plt.show()
    psi_hist_freq=list(psi_hist[0])
    max_freq1 = max(psi_hist_freq)
    max_index = psi_hist_freq.index(max_freq1)
    psi_hist_psi=list(psi_hist[1])
    psi_peak_loc1a=psi_hist_psi[max_index]
    plt.hist(psi[psi<np.pi],bins=200);
    plt.axvline(x=psi_peak_loc1a,color='k')
    plt.xlabel('psi')
    plt.show()

    psi_hist=plt.hist(psi[psi>np.pi/2],bins=200)
    plt.xlabel('psi')
    plt.show()
    psi_hist_freq=list(psi_hist[0])
    max_freq2 = max(psi_hist_freq)
    max_index = psi_hist_freq.index(max_freq2)
    psi_hist_psi=list(psi_hist[1])
    psi_peak_loc1b=psi_hist_psi[max_index]
    plt.hist(psi[psi<np.pi],bins=200);
    plt.axvline(x=psi_peak_loc1b,color='k')
    plt.xlabel('psi')
    plt.show()
    if max_freq2 > max_freq:
        psi_peak_loc2=psi_peak_loc1b-np.pi/2
        psi_peak_loc1=psi_peak_loc1b
    else:
        psi_peak_loc2=psi_peak_loc1a+np.pi/2
        psi_peak_loc1=psi_peak_loc1a

    plt.hist(psi,bins=200);
    plt.axvline(x=psi_peak_loc1,color='k')
    plt.axvline(x=psi_peak_loc2,color='k')
    plt.xlabel('psi')
    plt.show()

    if phi0_peak_loc1-np.pi/2 < 0:
        phi0_min_angle_A=phi0_peak_loc1+2*np.pi-np.pi/2
        phi0_max_angle_A=phi0_peak_loc1+2*np.pi+np.pi/2
        phi0_peak_locA=phi0_peak_loc1+2*np.pi
    else:
        phi0_min_angle_A=phi0_peak_loc1-np.pi/2
        phi0_max_angle_A=phi0_peak_loc1+np.pi/2
        phi0_peak_locA=phi0_peak_loc1


    if phi0_peak_loc2-np.pi/2 < 0:
        phi0_min_angle_B=phi0_peak_loc2+2*np.pi-np.pi/2
        phi0_max_angle_B=phi0_peak_loc2+2*np.pi+np.pi/2
        phi0_peak_locB=phi0_peak_loc2+2*np.pi
    else:
        phi0_min_angle_B=phi0_peak_loc2-np.pi/2
        phi0_max_angle_B=phi0_peak_loc2+np.pi/2
        phi0_peak_locB=phi0_peak_loc2



    if phi0_peak_locA > phi0_peak_locB:
        phi0_peak_loc=phi0_peak_locA
        phi0_min_angle=phi0_min_angle_A
        phi0_max_angle=phi0_max_angle_A
    else:
        phi0_peak_loc=phi0_peak_locB
        phi0_min_angle=phi0_min_angle_B
        phi0_max_angle=phi0_max_angle_B




    if psi_peak_loc1-np.pi/4 < 0:
        psi_min_angle_A=psi_peak_loc1+np.pi-np.pi/4
        psi_max_angle_A=psi_peak_loc1+np.pi+np.pi/4
        psi_peak_locA=psi_peak_loc1+np.pi
    else:
        psi_min_angle_A=psi_peak_loc1-np.pi/4
        psi_max_angle_A=psi_peak_loc1+np.pi/4
        psi_peak_locA=psi_peak_loc1


    if psi_peak_loc2-np.pi/4 < 0:
        psi_min_angle_B=psi_peak_loc2+np.pi-np.pi/4
        psi_max_angle_B=psi_peak_loc2+np.pi+np.pi/4
        psi_peak_locB=psi_peak_loc2+np.pi
    else:
        psi_min_angle_B=psi_peak_loc2-np.pi/4
        psi_max_angle_B=psi_peak_loc2+np.pi/4
        psi_peak_locB=psi_peak_loc2


    if psi_peak_locA > psi_peak_locB:
        psi_peak_loc=psi_peak_locA
        psi_min_angle=psi_min_angle_A
        psi_max_angle=psi_max_angle_A
    else:
        psi_peak_loc=psi_peak_locB
        psi_min_angle=psi_min_angle_B
        psi_max_angle=psi_max_angle_B

    angles=[0,np.pi/2,np.pi]
    nearest_angle_index=(np.abs(psi_peak_loc-angles)).argmin()
    angles_nearest=angles[nearest_angle_index]


    f0_a=[]
    f0_b=[]
    dfdt_a=[]
    dfdt_b=[]
    amp_a=[]
    amp_b=[]
    phi_a=[]
    phi_b=[]
    costheta_a=[]
    costheta_b=[]
    cosi_a=[]
    cosi_b=[]
    phi0_a=[]
    phi0_b=[]
    psi_a=[]
    psi_b=[]


    phi0_c=[]

    if phi0_peak_loc < np.pi/6:
        phi0_peak_loc = 2*np.pi
    for i in range(0,len(phi0),1):
        min_angle_phi0A=phi0_peak_loc-np.pi/2
        max_angle_phi0A=phi0_peak_loc+np.pi/2
        if min_angle_phi0A < 0:
            min_angle_phi0A+=2*np.pi
            max_angle_phi0A+=2*np.pi


    if psi_peak_loc < np.pi/6:
        psi_peak_loc = np.pi
        angles_nearest=np.pi
    for i in range(0,len(psi),1):
        min_angle_A=psi_peak_loc-np.pi/4
        max_angle_A=psi_peak_loc+np.pi/4
        if min_angle_A < 0:
            min_angle_A+=np.pi
            max_angle_A+=np.pi
            if (psi[i] + np.pi  >= min_angle_A and psi[i] + np.pi <= max_angle_A):
                phi0_a.append(phi0[i])
                phi0_c.append(phi0[i])
                psi_a.append(psi[i]+np.pi)
                f0_a.append(f0[i])
                dfdt_a.append(dfdt[i])
                amp_a.append(amp[i])
                phi_a.append(phi[i])
                costheta_a.append(costheta[i])
                cosi_a.append(cosi[i])
            elif (psi[i] >= min_angle_A and psi[i] <= max_angle_A):
                phi0_a.append(phi0[i])
                phi0_c.append(phi0[i])
                psi_a.append(psi[i])
                f0_a.append(f0[i])
                dfdt_a.append(dfdt[i])
                amp_a.append(amp[i])
                phi_a.append(phi[i])
                costheta_a.append(costheta[i])
                cosi_a.append(cosi[i])
            else:
                if np.absolute(psi[i]+np.pi-np.mean(psi_a)) < np.absolute(psi_b[i]+np.pi-np.mean(psi_b)):
                    phi0_a.append(phi0[i])
                    phi0_c.append(phi0[i])
                    psi_a.append(psi[i]+np.pi)
                    f0_a.append(f0[i])
                    dfdt_a.append(dfdt[i])
                    amp_a.append(amp[i])
                    phi_a.append(phi[i])
                    costheta_a.append(costheta[i])
                    cosi_a.append(cosi[i])
                else:
                    if (phi0[i]+np.pi*2 >= min_angle_phi0A and phi0[i]+np.pi*2 <= max_angle_phi0A+np.pi/4):
                        phi0_b.append(phi0[i]+np.pi*2)
                    else:
                        phi0_b.append(phi0[i])
                    psi_b.append(psi[i])
                    f0_b.append(f0[i])
                    dfdt_b.append(dfdt[i])
                    amp_b.append(amp[i])
                    phi_b.append(phi[i])
                    costheta_b.append(costheta[i])
                    cosi_b.append(cosi[i])
        elif np.round(angles_nearest,2) == 3.14:
            if (psi[i] + np.pi  >= min_angle_A and psi[i] + np.pi <= max_angle_A) :
                phi0_a.append(phi0[i])
                phi0_c.append(phi0[i])
                psi_a.append(psi[i]+ np.pi)
                f0_a.append(f0[i])
                dfdt_a.append(dfdt[i])
                amp_a.append(amp[i])
                phi_a.append(phi[i])
                costheta_a.append(costheta[i])
                cosi_a.append(cosi[i])
            elif (psi[i] >= min_angle_A and psi[i] <= max_angle_A):
                phi0_a.append(phi0[i])
                phi0_c.append(phi0[i])
                psi_a.append(psi[i])
                f0_a.append(f0[i])
                dfdt_a.append(dfdt[i])
                amp_a.append(amp[i])
                phi_a.append(phi[i])
                costheta_a.append(costheta[i])
                cosi_a.append(cosi[i])
            else:
                if (phi0[i]+np.pi*2 >= min_angle_phi0A and phi0[i]+np.pi*2 <= max_angle_phi0A+np.pi/4):
                    phi0_b.append(phi0[i]+np.pi*2)
                else:
                    phi0_b.append(phi0[i])
                psi_b.append(psi[i])
                f0_b.append(f0[i])
                dfdt_b.append(dfdt[i])
                amp_b.append(amp[i])
                phi_b.append(phi[i])
                costheta_b.append(costheta[i])
                cosi_b.append(cosi[i])
        else:
            if (psi[i] >= min_angle_A and psi[i] <= max_angle_A):
                phi0_a.append(phi0[i])
                phi0_c.append(phi0[i])
                psi_a.append(psi[i])
                f0_a.append(f0[i])
                dfdt_a.append(dfdt[i])
                amp_a.append(amp[i])
                phi_a.append(phi[i])
                costheta_a.append(costheta[i])
                cosi_a.append(cosi[i])
            else:
                if (phi0[i]+np.pi*2 >= min_angle_phi0A and phi0[i]+np.pi*2 <= max_angle_phi0A+np.pi/4):
                    phi0_b.append(phi0[i]+np.pi*2)
                else:
                    phi0_b.append(phi0[i])
                psi_b.append(psi[i])
                f0_b.append(f0[i])
                dfdt_b.append(dfdt[i])
                amp_b.append(amp[i])
                phi_b.append(phi[i])
                costheta_b.append(costheta[i])
                cosi_b.append(cosi[i])




    if len(phi0_b)>0:
        if np.round(max(phi0_b)-min(phi0_b),2)==6.28:
            for k in range(0,len(phi0_b),1):
                if phi0_b[k]<phi0_peak_loc-np.pi:
                    phi0_b[k]+=np.pi

    if len(psi_b)>0:
        if np.round(max(psi_b)-min(psi_b),2)==3.14:
            for k in range(0,len(psi_b),1):
                if psi_b[k]<psi_peak_loc-np.pi/4:
                    if psi_b[k] + np.pi < psi_peak_loc+np.pi/4:
                        psi_b[k]+=np.pi

    if np.round(max(psi_b)-min(psi_b),2) == 3.14:
        for k in range(0,len(psi_b),1):
            if psi_b[k] > np.pi/2 and psi_b[k]-np.pi>0:
                psi_b[k]-=np.pi

    if np.round(max(psi_b)-min(psi_b),2) == 3.14:
        for k in range(0,len(psi_b),1):
            if psi_b[k] < np.pi/2:
                psi_b[k]+=np.pi
        for i in range(0,len(phi_b),1):
            if np.absolute(psi_b[i]-np.mean(psi_b))> np.pi/2:
                if psi_b[i] <np.mean(psi_b):
                    psi_b[i]+=np.pi

    if len(phi0_a)>0:
        if len(phi0_a)>0:
            if np.round(max(phi0_a)-min(phi0_a),2)==6.28:
                for k in range(0,len(phi0_a),1):
                    if phi0_a[k]<phi0_peak_loc-np.pi:
                        phi0_a[k]+=np.pi





    if len(psi_a)>0:
        if np.round(max(psi_a)-min(psi_a),2)==3.14:
            for k in range(0,len(psi_a),1):
                if psi_a[k]<psi_peak_loc-np.pi/4:
                    if psi_a[k]<psi_peak_loc-np.pi/4:
                        if psi_a[k] + np.pi < psi_peak_loc+np.pi/4:
                            psi_a[k]+=np.pi
        
            if np.round(max(psi_a)-min(psi_a),2) == 3.14:
                for k in range(0,len(psi_a),1):
                    if psi_a[k] > np.pi/2 and psi_a[k]-np.pi>0:
                        psi_a[k]-=np.pi
                            
            if np.round(max(psi_a)-min(psi_a),2) == 3.14:
                for k in range(0,len(psi_a),1):
                    if psi_a[k] < np.pi/2 :
                        psi_a[k]+=np.pi



    if len(phi0_a)>0 and len(phi0_b)>0:
        f, (ax1)=plt.subplots(1,1)
        ax1.scatter(phi0,psi,s=15,c='k',marker='o')
        plt.xlabel('phi0',Fontsize=15)
        plt.ylabel('psi',Fontsize=15)
        ax1.scatter(phi0_a,psi_a,s=5,c='r',marker='o')
        plt.xlabel('phi0',Fontsize=15)
        plt.ylabel('psi',Fontsize=15)
        ax1.scatter(phi0_b,psi_b,s=5,c='g',marker='o')
        plt.xlabel('phi0',Fontsize=15)
        plt.ylabel('psi',Fontsize=15)
        plt.show()
        
        f, (ax1)=plt.subplots(1,1)
        ax1.scatter(phi0,f0,s=15,c='k',marker='o')
        plt.xlabel('phi0',Fontsize=15)
        plt.ylabel('f0',Fontsize=15)
        ax1.scatter(phi0_a,f0_a,s=5,c='r',marker='o')
        plt.xlabel('phi0',Fontsize=15)
        plt.ylabel('f0',Fontsize=15)
        ax1.scatter(phi0_b,f0_b,s=5,c='g',marker='o')
        plt.ylim(min(f0_b),max(f0_b))
        plt.xlabel('phi0',Fontsize=15)
        plt.ylabel('f0',Fontsize=15)
        plt.show()


else:
    pass


#if prior_tf=="T":
#    if len(phi0_a) == len(phi0) or len(phi0_b) == len(phi0):
#        if len(phi0_b) == len(phi0):
#            phi0_a=phi0_b
#            psi_a=psi_b
#            f0_a=f0_b
#            dfdt_a=dfdt_b
#            amp_a=amp_b
#            phi_a=phi_b
#            costheta_a=costheta_b
#            cosi_a=cosi_b
#        samplesA=np.vstack([f0_a,costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,dfdt_a])
#        SigmaA=cov(samplesA)
#        for k in range(0,8,1):
#            #compute midpoint of prior
#            mid_pt=(prior[k][0]+prior[k][1])/2
#            #adjust sample means to center of prior range
#            samplesA[k]*=mid_pt/np.mean(samplesA[k])
#
#        for k in range(0,2,1):
#            for g in range(0,8,1):
#                SigmaA[k][g]*=100
#
#        phi0_a=samplesA[6]
#        psi_a=samplesA[5]
#        f0_a=samplesA[0]
#        dfdt_a=samplesA[7]
#        amp_a=samplesA[3]
#        phi_a=samplesA[2]
#        costheta_a=samplesA[1]
#        cosi_a=samplesA[4]
#
#    else:
        #dist A
#        samplesA=np.vstack([f0_a,costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,dfdt_a])
#        SigmaA=cov(samplesA)
#        for k in range(0,8,1):
#            #compute midpoint of prior
#            mid_pt=(prior[k][0]+prior[k][1])/2
#            #adjust sample means to center of prior range
#            samplesA[k]*=mid_pt/np.mean(samplesA[k])
#
#
#        #dist B
#        samplesB=np.vstack([f0_b,costheta_b,phi_b,amp_b,cosi_b,psi_b,phi0_b,dfdt_b])
#        SigmaB=cov(samplesB)
#        for k in range(0,8,1):
#            #compute midpoint of prior
#            mid_pt=(prior[k][0]+prior[k][1])/2
#            #adjust sample means to center of prior range
#            samplesB[k]*=mid_pt/np.mean(samplesB[k])
#
#        for k in range(0,8,1):
#            for g in range(0,8,1):
#                SigmaB[k][g]*=np.absolute(prior[k][0]-prior[k][1])*np.absolute(prior[g][0]-prior[g][1])/9.0/(np.absolute(SigmaB[k][g]))
#                SigmaA[k][g]*=np.absolute(prior[k][0]-prior[k][1])*np.absolute(prior[g][0]-prior[g][1])/9.0/(np.absolute(SigmaA[k][g]))
#
#        phi0_a=samplesA[6]
#        psi_a=samplesA[5]
#        f0_a=samplesA[0]
#        dfdt_a=samplesA[7]
#        amp_a=samplesA[3]
#        phi_a=samplesA[2]
#        costheta_a=samplesA[1]
#        cosi_a=samplesA[4]
#
#        phi0_b=samplesB[6]
#        psi_b=samplesB[5]
#        f0_b=samplesB[0]
#        dfdt_b=samplesB[7]
#        amp_b=samplesB[3]
#        phi_b=samplesB[2]
#        costheta_b=samplesB[1]
#        cosi_b=samplesB[4]
if prior_tf=="T":
    samplesx=[]
    samples=[]
    var=f0_a
    for k in range(0,8,1):
        mean=(prior[k][0]+prior[k][1])/2
#        mean=np.mean(var)
        sigma=np.absolute(prior[k][1]-prior[k][0])/6
#        sigma=np.absolute(max(var)-min(var))*2
        for i in range(0,len(var),1):
            s=np.random.normal(mean, sigma, 1)
            while (s<prior[k][0] or s>prior[k][1]):
                s=np.random.normal(mean, sigma, 1)
            samplesx.append(s[0])
        samples.append(samplesx)
        samplesx=[]
    phi0_a=samples[6]
    psi_a=samples[5]
    f0_a=samples[0]
    dfdt_a=samples[7]
    amp_a=samples[3]
    phi_a=samples[2]
    costheta_a=samples[1]
    cosi_a=samples[4]

    phi0=phi0_a
    psi=psi_a
    f0=f0_a
    dfdt=dfdt_a
    amp=amp_a
    phi=phi_a
    costheta=costheta_a
    cosi=cosi_a
    #f0,dfdt,amp,phi,costheta,cosi,psi,phi0
    target=open(os.path.join('/Users/klackeos/Desktop/','inject.dat'),'w')
    target.write(str(np.mean(np.asarray(f0)/T)))
    target.write(" ")
    target.write(str(np.mean(np.asarray(dfdt)/(T**2))))
    target.write(" ")
    target.write(str(np.mean(np.pi/2-np.arccos(np.asarray(costheta)))))
    target.write(" ")
    target.write(str(np.mean(phi)))
    target.write(" ")
    target.write(str(np.mean(np.exp(np.asarray(amp)))))
    target.write(" ")
    target.write(str(np.mean(np.arccos(np.asarray(cosi)))))
    target.write(" ")
    target.write(str(np.mean(psi)))
    target.write(" ")
    target.write(str(np.mean(phi0)))
    target.close()


if len(phi0_a) == len(phi0) or len(phi0_b) == len(phi0):
        print("len(phi0_a) == len(phi0) or len(phi0_b) == len(phi0)")
        alpha=1.0
        if len(phi0_b) == len(phi0):
                phi0_a=phi0_b
                psi_a=psi_b
                f0_a=f0_b
                dfdt_a=dfdt_b
                amp_a=amp_b
                phi_a=phi_b
                costheta_a=costheta_b
                cosi_a=cosi_b
        elif prior_tf=="T":
                phi0_a=phi0
                psi_a=psi
                f0_a=f0
                dfdt_a=dfdt
                amp_a=amp
                phi_a=phi
                costheta_a=costheta
                cosi_a=cosi
        mcmc_sample_means_b=np.vstack([np.mean(f0_a),np.mean(costheta_a),np.mean(phi_a),np.mean(amp_a),np.mean(cosi_a),np.mean(psi_a),np.mean(phi0_a),np.mean(dfdt_a)])
        mcmc_samples_b=np.vstack([np.array(f0_a),costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,np.array(dfdt_a)])
        Sigma_mcmc_b=cov(mcmc_samples_b)
        target=open(os.path.join('/Users/klackeos/Desktop/','covariance.txt'),'w')
        target.write(str(alpha))
        target.write("        ")
        target.write(str(det(Sigma_mcmc_b)))
        target.write(" 0 0 0 0 0 0 ")
        target.write("\n")
        for l in range (0,8,1):
                target.write(str(mcmc_sample_means_b[l][0]))
                target.write(str('   '))
        target.write("\n")
        for k in range(0,8,1):
                for g in range(0,8,1):
                        target.write(str(Sigma_mcmc_b[k][g]))
                        target.write(str('   '))
                target.write("\n")
        target.close()
else:
        print('multi-modal')
        alpha=len(phi0_a)/len(phi0)
        
        mcmc_sample_means_a=np.vstack([np.mean(f0_a),np.mean(costheta_a),np.mean(phi_a),np.mean(amp_a),np.mean(cosi_a),np.mean(psi_a),np.mean(phi0_a),np.mean(dfdt_a)])
        mcmc_samples_a=np.vstack([np.array(f0_a),costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,np.array(dfdt_a)])
        Sigma_mcmc_a=cov(mcmc_samples_a)
        aaaa=LA.cholesky(Sigma_mcmc_a)
        print('here it is .......................................................................................................................................',aaaa[0][0])
        mcmc_sample_means_b=np.vstack([np.mean(f0_b),np.mean(costheta_b),np.mean(phi_b),np.mean(amp_b),np.mean(cosi_b),np.mean(psi_b),np.mean(phi0_b),np.mean(dfdt_b)])
        mcmc_samples_b=np.vstack([np.array(f0_b),costheta_b,phi_b,amp_b,cosi_b,psi_b,phi0_b,np.array(dfdt_b)])
        Sigma_mcmc_b=cov(mcmc_samples_b)
        
        bbbb=LA.cholesky(Sigma_mcmc_b)
        print('here it is .......................................................................................................................................',bbbb[0][0])
        target=open(os.path.join('/Users/klackeos/Desktop/','covariance.txt'),'w')
        #        target.write(str('alpha = '))
        target.write(str(alpha))
        target.write(" ")
        target.write(str(det(Sigma_mcmc_a)))
        target.write(" 0 0 0 0 0 0 ")

#        target.write("\n")
        target.write("\n")
        for l in range (0,8,1):
                target.write(str(mcmc_sample_means_a[l][0]))
                target.write(str('   '))
#        target.write("\n")
        target.write("\n")
        for k in range(0,8,1):
                for g in range(0,8,1):
                        target.write(str(Sigma_mcmc_a[k][g]))
                        target.write(str('   '))
                target.write("\n")
#        target.write("\n")
#        target.write("\n")
        target.write(str(1-alpha))
        target.write(" ")
        target.write(str(det(Sigma_mcmc_b)))
        target.write(" 0 0 0 0 0 0 ")
#        target.write("\n")
        target.write("\n")
        for l in range (0,8,1):
                target.write(str(mcmc_sample_means_b[l][0]))
                target.write(str('   '))
#        target.write("\n")
        target.write("\n")
        for k in range(0,8,1):
                for g in range(0,8,1):
                        target.write(str(Sigma_mcmc_b[k][g]))
                        target.write(str('   '))
                target.write("\n")
        target.close()













###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
if len(phi0_a) == len(phi0):
        print('phi0 is unimodal or you are drawing from prior')

        alpha=1.0
        d1=np.vstack([f0_a,costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,dfdt_a])
        Sigma1=cov(d1)
        mean1=[np.mean(f0_a),np.mean(costheta_a),np.mean(phi_a),np.mean(amp_a),np.mean(cosi_a),np.mean(psi_a),np.mean(phi0_a),np.mean(dfdt_a)]
        mean1_v=np.vstack([mean1[0],mean1[1],mean1[2],mean1[3],mean1[4],mean1[5],mean1[6],mean1[7]])

        #Create one Covariance matrix ###################################################################################################
        Sigma_a = Sigma1
        mean_a=mean1



        a1=LA.cholesky(Sigma_a)


        print("Enter no. of iterations:")
        E=input()
        E=int(E)
        print("Enter no. of burn-in steps:")
        BURN_IN=input()
        BURN_IN=int(BURN_IN)

        # Initialize the chain.
        prior_rej_no=0
        # Store the samples
        acc=[]
        p=np.array([0.]*(E-BURN_IN))
        chain_f0=np.array([0.]*(E-BURN_IN))
        chain_amp=np.array([0.]*(E-BURN_IN))
        chain_phi0=np.array([0.]*(E-BURN_IN))
        chain_dfdt=np.array([0.]*(E-BURN_IN))
        chain_cosi=np.array([0.]*(E-BURN_IN))
        chain_costheta=np.array([0.]*(E-BURN_IN))
        chain_phi=np.array([0.]*(E-BURN_IN))
        chain_psi=np.array([0.]*(E-BURN_IN))
        accepted_number=0.
        if prior_tf == "True":
            d0=np.vstack([uniform.rvs(prior[0][0],prior[0][1]),uniform.rvs(prior[1][0],prior[1][1]),uniform.rvs(prior[2][0],prior[2][1]),-uniform.rvs(45,55),uniform.rvs(prior[4][0],prior[4][1]),uniform.rvs(prior[5][0],prior[5][1]),uniform.rvs(prior[6][0],prior[6][1]),uniform.rvs(prior[7][0],prior[7][1])])
        else:
            dd=np.dot(a1,np.random.normal(0, 1, 8))+mean1
            d0=np.vstack([dd[0],dd[1],dd[2],dd[3],dd[4],dd[5],dd[6],dd[7]])
        pr0=1.0
        mean1_v=np.vstack([mean1[0],mean1[1],mean1[2],mean1[3],mean1[4],mean1[5],mean1[6],mean1[7]])
        q0=np.asarray(pxd(d0, mean1_v, Sigma1,include_prior_rej))
        L0=np.asarray(pxd(d0, mean1_v, Sigma1,include_prior_rej))
        for e in range(0,E,1):
                if e % 200000 == 0:
                        print("******************************************************************************************************************************************")
                        print("At iteration "+str(e))
                        print("******************************************************************************************************************************************")
                # (1) Generate a proposal (or a candidate) sample (x_cand,y_cand) from the proposal distribution q=alpha*q1+(alpha-1)*q2
                #70% of time candidate is drawn from q1 and 30% of time from q2, alpha = 0.7
                
                a1=LA.cholesky(Sigma1)
                N1=np.random.normal(0, 1, 8)
                cand_f0,cand_costheta,cand_phi,cand_amp,cand_cosi,cand_psi,cand_phi0,cand_dfdt=np.dot(a1,N1)+mean1

                
                if prior_tf=="T" or include_prior_rej=="T":
                    cand_phi_tmp=cand_phi
                    cand_phi0_tmp=cand_phi0
                    cand_psi_tmp=cand_psi

                    if cand_phi<0.0:
                            cand_phi_tmp+=2*np.arccos(-1.0)
                    elif cand_phi>2*np.arccos(-1.0):
                            cand_phi_tmp-=2*np.arccos(-1.0)
                    if cand_phi0<0.0:
                            cand_phi0_tmp+=2*np.arccos(-1.0)
                    elif cand_phi0>2*np.arccos(-1.0):
                            cand_phi0_tmp-=2*np.arccos(-1.0)
                    if cand_psi<0.0:
                            cand_psi_tmp+=np.arccos(-1.0)
                    elif cand_psi>np.arccos(-1.0):
                            cand_psi_tmp-=np.arccos(-1.0)
                    d=np.vstack([cand_f0,cand_costheta,cand_phi_tmp,cand_amp,cand_cosi,cand_psi_tmp,cand_phi0_tmp,cand_dfdt])
                else:
                    d=np.vstack([cand_f0,cand_costheta,cand_phi,cand_amp,cand_cosi,cand_psi,cand_phi0,cand_dfdt])

                q=np.asarray(pxd(d, mean1_v, Sigma1,include_prior_rej))
                L=np.asarray(pxd(d, mean1_v, Sigma1,include_prior_rej))
                print (q-q0,q,q0)
                if prior_tf =="T" or include_prior_rej=="T":
                        for g in range(0,8,1):
        #                        print(d[g],prior[g][0],prior[g][1])
                                if (d[g]<prior[g][0] or d[g]>prior[g][1]):
                                    pr=0.0
                                    prior_rej_no+=1
                                    print("REJECTED",g,d[g],prior[g][0],prior[g][1])
                                    if include_prior_rej=="T":
                                        q=0.0
                                    break
                                else:
                                    pr=1.0
                candidate=d
                # Compute the acceptance probability, Equation 8 and Equation 6.
                # We will do both equations in log domain here to avoid underflow.
                if prior_tf=="T":
                        accept=log(pr)+log(q0)
                        accept=accept-log(pr0)-log(q)
                else:
                        accept=np.log(q0*L)
                        accept=accept-np.log(q*L0)
                        print(accept)
#                        print(accept,q,L0,q0,L,d)
                #         accept=q/q0
                if prior_tf == "F":
                    if q != 0.0:
                        accept=min([1,accept])
                    else:
                        accept=-np.inf
    #                    print(accept)
                elif prior_tf == "T":
                    if pr == 0.0:
                        accept=-np.inf
#                    elif pr == 0.0 and pr0 == 0.0:
#                        accept=log(q0)-log(q)
#                    elif pr0 == 0:
#                        accept=0.0

                accept=exp(accept)
                acc.append(accept)
                # Accept candidate with probability accept.
                u=uniform.rvs(0,1)
                if u<accept:
                        d0=candidate
                        accepted_number=accepted_number+1
                        q0=q
                        L0=L
#                        pr0=0.0
                else:
#                        pr0=1.0
                        pass

                if (e>=BURN_IN):
                    chain_f0[e-BURN_IN]=d0[0][0]
                    chain_costheta[e-BURN_IN]=d0[1][0]
                    chain_phi[e-BURN_IN]=d0[2][0]
                    chain_amp[e-BURN_IN]=d0[3][0]
                    chain_cosi[e-BURN_IN]=d0[4][0]
                    chain_psi[e-BURN_IN]=d0[5][0]
                    chain_phi0[e-BURN_IN]=d0[6][0]
                    chain_dfdt[e-BURN_IN]=d0[7][0]
                    p[e-BURN_IN]=q
        
        print("...Summary...")
        print("Acceptance ratio is "+str(accepted_number/(E)))

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        xdata=chain_phi0
        ydata=chain_psi
        zdata=p
        ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Blues');
        plt.title('mcmc chain data')
        plt.show()
        chainxx=chain_phi0
        
        
        chainxx=[]
        phi0_aa=[]
        for k in range(0,len(chain_phi0),1):
                if chain_phi0[k]>max(phi0):
                        chainxx.append(chain_phi0[k]-2*np.pi)
                else:
                        chainxx.append(chain_phi0[k])



        for k in range(0,len(phi0_a),1):
                if phi0_a[k]>max(phi0):
                        phi0_aa.append(phi0_a[k]-2*np.pi)
                else:
                        phi0_aa.append(phi0_a[k])

#         chainyy=[]

#         for k in range(0,len(chain_psi),1):
#                 if chain_psi[k]>max(psi):
#                         chainyy.append(chain_psi[k]-np.pi)
#                 else:
#                         chainyy.append(chain_psi[k])
#
        psi_aa=[]
        for k in range(0,len(psi_a),1):
                if psi_a[k]>max(psi):
                        psi_aa.append(psi_a[k]-np.pi)
                else:
                        psi_aa.append(psi_a[k])


#        f, (ax1)=plt.subplots(1,1)
#        ax1.scatter(chainxx,chain_psi,s=10,c='k',marker='o')
#        ax1.scatter(phi0_aa,psi_aa,s=12,c='g',marker='o')
##         ax1.scatter(phi0,psi,s=10,c='b',marker='+')
#        # ax1.scatter(phi0_b2,psi_b2,s=10,c='m',marker='o')
#        plt.xlabel('phi0',Fontsize=15)
#        plt.ylabel('psi',Fontsize=15)
#        plt.show()


#        plt.hist(phi0,bins=200);
#        plt.hist(chainxx,bins=200);
#        plt.axvline(x=mean1[7],color='k')
#        plt.show()

#        plt.plot(chainxx,'r*')
#        plt.plot(mean1[7]*np.ones(len(chain_phi0)),'g-')
#        plt.show()
#
#        plt.hist(chain_psi,bins=200);
#        plt.hist(psi,bins=200);
#        plt.axvline(x=mean1[6],color='k')
#        plt.show()
#
#        plt.plot(chain_psi,'r*')
#        plt.plot(mean1[6]*np.ones(len(chain_phi0)),'g-')
#        plt.show()
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
else:
    
        alpha=len(phi0_a)/len(phi0)
        d1=[f0_a,costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,dfdt_a]
        Sigma1=cov(d1)
        mean1=[np.mean(f0_a),np.mean(costheta_a),np.mean(phi_a),np.mean(amp_a),np.mean(cosi_a),np.mean(psi_a),np.mean(phi0_a),np.mean(dfdt_a)]


        d2=[f0_b,costheta_b,phi_b,amp_b,cosi_b,psi_b,phi0_b,dfdt_b]
        Sigma2=cov(d2)
        mean2=[np.mean(f0_b),np.mean(costheta_b),np.mean(phi_b),np.mean(amp_b),np.mean(cosi_b),np.mean(psi_b),np.mean(phi0_b),np.mean(dfdt_b)]



        #Create two Covariance matrices ###################################################################################################
        Sigma_a = Sigma1
        mean_a=mean1

        Sigma_b = Sigma2
        mean_b=mean2


        a1=LA.cholesky(Sigma_a)

        a2=LA.cholesky(Sigma_b)
        mean1_v=np.vstack([mean1[0],mean1[1],mean1[2],mean1[3],mean1[4],mean1[5],mean1[6],mean1[7]])
        mean2_v=np.vstack([mean2[0],mean2[1],mean2[2],mean2[3],mean2[4],mean2[5],mean2[6],mean2[7]])
        #Draw samples from Covariance matrices  ###################################################################################################


        # MH sampler
        print("Enter no. of iterations:")
        E=input()
        E=int(E)
        print("Enter no. of burn-in steps:")
        BURN_IN=input()
        BURN_IN=int(BURN_IN)

        # Initialize the chain.

        # Store the samples
        prior_rej_no=0
        acc=[]
        p=np.array([0.]*(E-BURN_IN))
        chain_f0=np.array([0.]*(E-BURN_IN))
        chain_amp=np.array([0.]*(E-BURN_IN))
        chain_phi0=np.array([0.]*(E-BURN_IN))
        chain_dfdt=np.array([0.]*(E-BURN_IN))
        chain_cosi=np.array([0.]*(E-BURN_IN))
        chain_costheta=np.array([0.]*(E-BURN_IN))
        chain_phi=np.array([0.]*(E-BURN_IN))
        chain_psi=np.array([0.]*(E-BURN_IN))
        accepted_number=0.
        
        if prior_tf == "True":
            d0=np.vstack([uniform.rvs(prior[0][0],prior[0][1]),uniform.rvs(prior[1][0],prior[1][1]),uniform.rvs(prior[2][0],prior[2][1]),-uniform.rvs(45,55),uniform.rvs(prior[4][0],prior[4][1]),uniform.rvs(prior[5][0],prior[5][1]),uniform.rvs(prior[6][0],prior[6][1]),uniform.rvs(prior[7][0],prior[7][1])])
        else:
            N=np.random.uniform(0, 1, 1)
            if (N >= (1-alpha)):
                dd=np.dot(a1,np.random.normal(0, 1, 8))+mean1
                d0=np.vstack([dd[0],dd[1],dd[2],dd[3],dd[4],dd[5],dd[6],dd[7]])
            else:
                dd=np.dot(a1,np.random.normal(0, 1, 8))+mean2
                d0=np.vstack([dd[0],dd[1],dd[2],dd[3],dd[4],dd[5],dd[6],dd[7]])
        q0=np.asarray(pxd(d0, mean1_v, Sigma1,include_prior_rej))+np.asarray(pxd(d0, mean2_v, Sigma2,include_prior_rej))
        L0=np.asarray(pxd(d0, mean1_v, Sigma1,include_prior_rej))+np.asarray(pxd(d0, mean2_v, Sigma2,include_prior_rej))
        pr0=1.0
        for e in range(0,E,1):
            N=np.random.uniform(0, 1, 1)
            if e % 1000 == 0:
                    print("******************************************************************************************************************************************")
                    print("At iteration "+str(e))
                    print("******************************************************************************************************************************************")
            if N >= (1-alpha):
                    #q1
                    a1=LA.cholesky(Sigma1)
    
                    N1=np.random.normal(0, 1, 8)
                    cand_f0,cand_costheta,cand_phi,cand_amp,cand_cosi,cand_psi,cand_phi0,cand_dfdt=np.dot(a1,N1)+mean1
                    d=np.vstack([cand_f0,cand_costheta,cand_phi,cand_amp,cand_cosi,cand_psi,cand_phi0,cand_dfdt])
                    sigma=Sigma1
            else:
                    #q2
                    a2=LA.cholesky(Sigma2)

                    N2=np.random.normal(0, 1, 8)
                    cand_f0,cand_costheta,cand_phi,cand_amp,cand_cosi,cand_psi,cand_phi0,cand_dfdt=np.dot(a2,N2)+mean2
                    d=np.vstack([cand_f0,cand_costheta,cand_phi,cand_amp,cand_cosi,cand_psi,cand_phi0,cand_dfdt])
                    sigma=Sigma2
#            print(cand_costheta,cand_cosi)
            if prior_tf=="T" or include_prior_rej=="T":
                cand_phi_tmp=cand_phi
                cand_phi0_tmp=cand_phi0
                cand_psi_tmp=cand_psi
                
                if cand_phi<0.0:
                    cand_phi_tmp+=2*np.arccos(-1.0)
                elif cand_phi>2*np.arccos(-1.0):
                    cand_phi_tmp-=2*np.arccos(-1.0)
                if cand_phi0<0.0:
                    cand_phi0_tmp+=2*np.arccos(-1.0)
                elif cand_phi0>2*np.arccos(-1.0):
                    cand_phi0_tmp-=2*np.arccos(-1.0)
                if cand_psi<0.0:
                    cand_psi_tmp+=np.arccos(-1.0)
                elif cand_psi>np.arccos(-1.0):
                    cand_psi_tmp-=np.arccos(-1.0)
                d=np.vstack([cand_f0,cand_costheta,cand_phi_tmp,cand_amp,cand_cosi,cand_psi_tmp,cand_phi0_tmp,cand_dfdt])

            q=np.asarray(pxd(d, mean1_v, Sigma1, include_prior_rej))+np.asarray(pxd(d, mean2_v, Sigma2, include_prior_rej))
            L=np.asarray(pxd(d, mean1_v, Sigma1, include_prior_rej))+np.asarray(pxd(d, mean2_v, Sigma2, include_prior_rej))
            if prior_tf=="T" or include_prior_rej=="T":
                for g in range(0,8,1):
                    if (d[g]<prior[g][0] or d[g]>prior[g][1]):
                        pr=0.0
                        prior_rej_no+=1
                        print("REJECTED",g,d[g],prior[g][0],prior[g][1])
                        if include_prior_rej=="T":
                            q=0.0
                            break
                        else:
                            pr=1.0

            candidate=d

            # Compute the acceptance probability, Equation 8 and Equation 6.
            # We will do both equations in log domain here to avoid underflow.
            if prior_tf=="T":
                accept=log(pr)+log(q0)
                accept=accept-log(pr0)-log(q)
            else:
                accept=np.log(q*L0)
                accept=accept-np.log(q0*L)
#                print(accept,q,L)
            if prior_tf == "F":
                if q != 0.0:
                    accept=min([1,accept])
                else:
                    print('q=0')
                    accept=-np.inf
            elif prior_tf == "T":
                if pr == 0.0:
                     accept=-np.inf
            accept=exp(accept)
            acc.append(accept)
            u=uniform.rvs(0,1)
            if u<accept:
                d0=candidate
                accepted_number=accepted_number+1
                q0=q
                L0=L
            else:
                pass

            # store
            if e>=BURN_IN:
                    chain_f0[e-BURN_IN]=d0[0][0]
                    chain_costheta[e-BURN_IN]=d0[1][0]
                    chain_phi[e-BURN_IN]=d0[2][0]
                    chain_amp[e-BURN_IN]=d0[3][0]
                    chain_cosi[e-BURN_IN]=d0[4][0]
                    chain_psi[e-BURN_IN]=d0[5][0]
                    chain_phi0[e-BURN_IN]=d0[6][0]
                    chain_dfdt[e-BURN_IN]=d0[7][0]
                    p[e-BURN_IN]=q

#        fig = plt.figure()
#        ax = plt.axes(projection='3d')
#        xdata=chain_phi0
#        ydata=chain_psi
#        zdata=p
#        ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Blues');
#        plt.title('mcmc chain data')
#        plt.show()

        print("...Summary...")
        print("Acceptance ratio is "+str(accepted_number/(E)))



        # phi0_a2=[]
        # phi0_b2=[]
        # psi_a2=[]
        # psi_b2=[]
        # for i in range(0,len(chainx),1):
        #         if (chainx[i] >= 3*np.pi/2 and chainx[i] <= 5*np.pi/2):
        #                 phi0_a2.append(chainx[i])
        #                 psi_a2.append(chainy[i])
        #         else:
        #                 if chainx[i] >= np.pi/2 and chainx[i] <= 3*np.pi/2:
        #                         phi0_b2.append(chainx[i])
        #                         psi_b2.append(chainy[i])

        chainxx=[]
        phi0_aa=[]
        phi0_bb=[]
        for k in range(0,len(chain_phi0),1):
                if chain_phi0[k]>max(phi0):
                        chainxx.append(chain_phi0[k]-2*np.pi)
                else:
                        chainxx.append(chain_phi0[k])

#        for k in range(0,len(phi0_a),1):
#                if phi0_a[k]>max(phi0):
#                        phi0_aa.append(phi0_a[k]-2*np.pi)
#                else:
#                        phi0_aa.append(phi0_a[k])
#
#        for k in range(0,len(phi0_b),1):
#                if phi0_b[k]>max(phi0):
#                        phi0_bb.append(phi0_b[k]-2*np.pi)
#                else:
#                        phi0_bb.append(phi0_b[k])




        chainyy=[]
        psi_aa=[]
        psi_bb=[]
        for k in range(0,len(chain_psi),1):
                if chain_psi[k]>max(psi):
                        chainyy.append(chain_psi[k]-np.pi)
                else:
                        chainyy.append(chain_psi[k])

#        for k in range(0,len(psi_a),1):
#                if psi_a[k]>max(psi):
#                        psi_aa.append(psi_a[k]-np.pi)
#                else:
#                        psi_aa.append(psi_a[k])
#
#        for k in range(0,len(psi_b),1):
#                if psi_b[k]>max(psi):
#                        psi_bb.append(psi_b[k]-np.pi)
#                else:
#                        psi_bb.append(psi_b[k])





        f, (ax1)=plt.subplots(1,1)
        ax1.scatter(chainxx,chainyy,s=10,c='k',marker='o')
        ax1.scatter(phi0_aa,psi_aa,s=12,c='g',marker='o')
        ax1.scatter(phi0_bb,psi_bb,s=12,c='r',marker='o')
#         ax1.scatter(phi0,psi,s=10,c='b',marker='+')
        # ax1.scatter(phi0_b2,psi_b2,s=10,c='m',marker='o')
        plt.xlabel('phi0',Fontsize=15)
        plt.ylabel('psi',Fontsize=15)
        plt.show()


        plt.hist(phi0,bins=200);
        plt.hist(chainxx,bins=200);
        plt.axvline(x=mean2[7],color='k')
        plt.axvline(x=mean1[7],color='k')
        plt.show()

        plt.plot(chainxx,'r*')
        plt.plot(mean2[7]*np.ones(len(chain_phi0)),'g-')
        plt.plot(mean1[7]*np.ones(len(chain_phi0)),'g-')
        plt.show()

        plt.hist(chainyy,bins=200);
        plt.hist(psi,bins=200);
        plt.axvline(x=mean2[6],color='k')
        plt.axvline(x=mean1[6],color='k')
        plt.show()

        plt.plot(chain_psi,'r*')
        plt.plot(mean2[6]*np.ones(len(chain_phi0)),'g-')
        plt.plot(mean1[6]*np.ones(len(chain_phi0)),'g-')
        plt.show()
#
#
#
if len(phi0_a) == len(phi0) or len(phi0_b) == len(phi0) or prior_tf=="T":
        data=np.vstack([chain_f0/T,chain_dfdt/(T**2),exp(chain_amp),chain_phi,chain_costheta,chain_cosi,chain_psi,chainxx])
else:
        data=np.vstack([chain_f0/T,chain_dfdt/(T**2),exp(chain_amp),chain_phi,chain_costheta,chain_cosi,chainyy,chainxx])
        chain_psi=chainyy
target=open(os.path.join('/Users/klackeos/Desktop/','dimension_chain_new1.dat.1'),'w')
for g in range(0,len(data[0]),1):
        for k in range(0,8,1):
                target.write(str(data[k][g]))
                target.write(str('   '))
        target.write("\n")

target.close()

data0=np.vstack([f0/T,dfdt/T/T,exp(amp),phi,costheta,cosi,psi,phi0])
datanew0=np.transpose(data0)


#    chain_f0[e-BURN_IN]=d0[0][0]
#    chain_costheta[e-BURN_IN]=d0[1][0]
#    chain_phi[e-BURN_IN]=d0[2][0]
#    chain_amp[e-BURN_IN]=d0[3][0]
#    chain_cosi[e-BURN_IN]=d0[4][0]
#    chain_psi[e-BURN_IN]=d0[5][0]
#    chain_phi0[e-BURN_IN]=d0[6][0]
#    chain_dfdt[e-BURN_IN]=d0[7][0]


rng=[(min(chain_f0)/T,max(chain_f0)/T), (min(chain_dfdt)/T/T,max(chain_dfdt)/T/T), (min(exp(chain_amp)),max(exp(chain_amp))), (min(chain_phi),max(chain_phi)), (min(chain_costheta),max(chain_costheta)), (min(chain_cosi),max(chain_cosi)), (min(chain_psi),max(chain_psi)), (min(chainxx),max(chainxx))]
datanew=np.transpose(data)
figure = corner.corner(datanew,plot_datapoints=False,color='r',plot_contours=True,fill_contours=True,plot_density=True,quantiles=(0.16, 0.84), range=rng, labels=[r"$f0 (Hz)$", "dfdt", r"$amp$", r"$\phi$",r"$cos\theta$",r"$cos\iota$", r"$\psi$", r"$\phi_0$"],levels=(0.39,0.86,0.99),show_titles=True,title_fmt='01.02e')
corner.corner(datanew0, fig=figure,plot_datapoints=False,plot_contours=True,fill_contours=True,plot_density=True,range=rng,quantiles=(0.16, 0.84),levels=(0.39,0.86,0.99), weights=np.ones(len(datanew0))*len(datanew)/len(datanew0))
ndim=8
values=[]
index=int(find_nearest(freq_ldc,np.mean(f0)))
values.append(freq_ldc[index])
values.append(dot_freq_ldc[index])
values.append(amp_ldc[index])
values.append(ecl_long_ldc[index])
values.append(ecl_lat_ldc[index])
values.append(inc_ldc[index])
values.append(pol_ldc[index])
values.append(phi0_ldc[index])
print(values)

axes = np.array(figure.axes).reshape((ndim, ndim))
for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(values[i], color="g")
print(np.absolute(values[0]-np.mean(f0)))

for yi in range(ndim):
        for xi in range(yi):
                ax = axes[yi, xi]
                ax.axvline(values[xi], color="g")
                ax.axhline(values[yi], color="g")
# figure.savefig('/Users/klackeos/ldasoft/master/bin/ldc_vgb_noise/source'+no+'/'+str(obs_tm)+'/plots/pdfs/pol_costheta_phi_'+no+'.png')
plt.show()
















rng=[(min(chain_f0)/T,max(chain_f0)/T), (min(chain_dfdt)/T/T,max(chain_dfdt)/T/T), (min(exp(chain_amp)),max(exp(chain_amp))), (min(chain_phi),max(chain_phi)), (min(chain_costheta),max(chain_costheta)), (min(chain_cosi),max(chain_cosi)), (min(chain_psi),max(chain_psi)), (min(chainxx),max(chainxx))]


c = ChainConsumer()
c.add_chain(datanew, parameters=[r"$f0 (Hz)$", "df/dt", r"$amp$", r"$\phi$",r"$cos\theta$",r"$cos\iota$", r"$\psi$", r"$\phi_0$"],name="covariance matrix sampling").configure(contour_labels="confidence")
c.add_chain(datanew0, parameters=[r"$f0 (Hz)$", "df/dt", r"$amp$", r"$\phi$",r"$cos\theta$",r"$cos\iota$", r"$\psi$", r"$\phi_0$"],name="gbmcmc chains")
c.configure(shade_gradient=[0.8, 1.5],bins=25,colors=["c","k"], shade_alpha=0.8, linewidths=[1.0, 2.0], tick_font_size=8, max_ticks=4,shade=[True, False],summary=[True,True],kde=[False, False])
#c.configure(shade_gradient=[0.8],bins=25,colors=["c"], shade_alpha=0.8, linewidths=[1.0], tick_font_size=8, max_ticks=4,shade=[True],summary=[True],kde=[False])

#c.configure(bins=25, rainbow=True)
fig = c.plotter.plot(extents=rng)
plt.show()



