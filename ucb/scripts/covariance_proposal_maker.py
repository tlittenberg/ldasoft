import os
import numpy as np
import matplotlib.pyplot as plt
from numpy import log,exp,pi
from numpy import array
from numpy import cov
from numpy import linalg as LA
from numpy.linalg import inv
from numpy.linalg import det



f=open('entries.dat','r')
f1=f.readlines()

files=[]
for x in f1:
    files.append(x.split()[0]+'_chain.dat')
r=0


print("Enter observation time for covariance matrix (1.5 mo = 3932160, 3 mo=7864320, 6 mo=15728640, 1 yr=31457280, 2 yr=62914560):")
T=input()
T=int(T)

print("Enter DIRECTORY where you want to save COVARIANCE MATRIX output (warning: code will OVERWRITE old files)")
cov_dir=input()

print("Enter FILENAME for COVARIANCE MATRIX output (it is an ASCII file)")
cov_fl_nm=input()

target=open(os.path.join(str(cov_dir),str(cov_fl_nm)+'_'+str(r)),'a+')
no=0
for x in files:
        print(x)
        f2=open(x,'r')
        samples=np.loadtxt(f2)
        f0=samples[:,0]*T
        dfdt=samples[:,1]*T**2
        amp=np.log(samples[:,2])
        phi=samples[:,3]
        costheta=samples[:,4]
        cosi=samples[:,5]
        psi=samples[:,6]
        phi0=samples[:,7]


        phi0_hist=plt.hist(phi0[phi0<np.pi],bins=200);
#            plt.xlabel('phi0',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of initial GW phase, phi0 $<$ pi',Fontsize=15)
#            #plt.show()
        phi0_hist_freq=list(phi0_hist[0])
        max_freq = max(phi0_hist_freq)
        max_index = phi0_hist_freq.index(max_freq)
        phi0_hist_phi0=list(phi0_hist[1])
        phi0_peak_loc1=phi0_hist_phi0[max_index]
        plt.hist(phi0[phi0<np.pi],bins=200);
#            plt.axvline(x=phi0_peak_loc1,color='k')
#            plt.xlabel('phi0',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of initial GW phase, find peak 1',Fontsize=15)
#            #plt.show()

        phi0_peak_loc2=phi0_peak_loc1+np.pi

        plt.hist(phi0,bins=200);
#            plt.axvline(x=phi0_peak_loc1,color='k')
#            plt.axvline(x=phi0_peak_loc2,color='k')
#            plt.xlabel('phi0',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of initial GW phase, both peaks',Fontsize=15)
#            #plt.show()

        psi_hist=plt.hist(psi[psi<np.pi/2],bins=200)
#            plt.xlabel('psi',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of GW polarization, psi $<$ pi/2',Fontsize=15)
#            #plt.show()
        psi_hist_freq=list(psi_hist[0])
        max_freq1 = max(psi_hist_freq)
        max_index = psi_hist_freq.index(max_freq1)
        psi_hist_psi=list(psi_hist[1])
        psi_peak_loc1a=psi_hist_psi[max_index]
#            plt.hist(psi[psi<np.pi],bins=200);
#            plt.axvline(x=psi_peak_loc1a,color='k')
#            plt.xlabel('psi',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of GW polarization, first peak',Fontsize=15)
#            #plt.show()

        psi_hist=plt.hist(psi[psi>np.pi/2],bins=200)
#            plt.xlabel('psi',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of GW polarization, psi $>$ pi/2',Fontsize=15)
#            #plt.show()
        psi_hist_freq=list(psi_hist[0])
        max_freq2 = max(psi_hist_freq)
        max_index = psi_hist_freq.index(max_freq2)
        psi_hist_psi=list(psi_hist[1])
        psi_peak_loc1b=psi_hist_psi[max_index]
#            plt.hist(psi[psi<np.pi],bins=200);
#            plt.axvline(x=psi_peak_loc1b,color='k')
#            plt.xlabel('psi',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of GW polarization, second peak',Fontsize=15)
#            #plt.show()
        if max_freq2 > max_freq:
                psi_peak_loc2=psi_peak_loc1b-np.pi/2
                psi_peak_loc1=psi_peak_loc1b
        else:
                psi_peak_loc2=psi_peak_loc1a+np.pi/2
                psi_peak_loc1=psi_peak_loc1a

#            plt.hist(psi,bins=200);
#            plt.axvline(x=psi_peak_loc1,color='k')
#            plt.axvline(x=psi_peak_loc2,color='k')
#            plt.xlabel('psi',Fontsize=15)
#            plt.ylabel('Frequency',Fontsize=15)
#            plt.title('Histogram of GW polarization, both peaks',Fontsize=15)

#            #plt.show()

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
#                    print(i/len(phi_b))
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

        if len(phi0_a) == len(phi0) or len(phi0_b) == len(phi0):
            print("len(phi0_a) == len(phi0) or len(phi0_b) == len(phi0)")
            alpha=1.0
            mcmc_sample_means_b=np.vstack([np.mean(f0_a),np.mean(costheta_a),np.mean(phi_a),np.mean(amp_a),np.mean(cosi_a),np.mean(psi_a),np.mean(phi0_a),np.mean(dfdt_a)])
            mcmc_samples_b=np.vstack([np.array(f0_a),costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,np.array(dfdt_a)])
            Sigma_mcmc_b=cov(mcmc_samples_b)
            target=open(os.path.join(str(cov_dir),str(cov_fl_nm)+'_'+str(r)),'w')
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
        #    target.close()
#            r+=1
        else:
            print('multi-modal')
            
            
            alpha=len(phi0_a)/len(phi0)
#            print(alpha)
            mcmc_sample_means_a=np.vstack([np.mean(f0_a),np.mean(costheta_a),np.mean(phi_a),np.mean(amp_a),np.mean(cosi_a),np.mean(psi_a),np.mean(phi0_a),np.mean(dfdt_a)])
            mcmc_samples_a=np.vstack([np.array(f0_a),costheta_a,phi_a,amp_a,cosi_a,psi_a,phi0_a,np.array(dfdt_a)])
            Sigma_mcmc_a=cov(mcmc_samples_a)
            cholesky_a=LA.cholesky(Sigma_mcmc_a)
            mcmc_sample_means_b=np.vstack([np.mean(f0_b),np.mean(costheta_b),np.mean(phi_b),np.mean(amp_b),np.mean(cosi_b),np.mean(psi_b),np.mean(phi0_b),np.mean(dfdt_b)])
            mcmc_samples_b=np.vstack([np.array(f0_b),costheta_b,phi_b,amp_b,cosi_b,psi_b,phi0_b,np.array(dfdt_b)])
            Sigma_mcmc_b=cov(mcmc_samples_b)
            
            cholesky_b=LA.cholesky(Sigma_mcmc_b)
            
            if no==0:
#                print(no)
                target.write(str(len(files)))
                target.write(" 0 0 0 0 0 0 0 ")
                target.write("\n")
            target.write(str(alpha))
            target.write(" ")
            target.write(str(det(Sigma_mcmc_a)))
            target.write(" 0 0 0 0 0 0 ")

            target.write("\n")
            for l in range (0,8,1):
                    target.write(str(mcmc_sample_means_a[l][0]))
                    target.write(str('   '))
            target.write("\n")
            for k in range(0,8,1):
                    for g in range(0,8,1):
                            target.write(str(Sigma_mcmc_a[k][g]))
                            target.write(str('   '))
                    target.write("\n")
            for k in range(0,8,1):
                for g in range(0,8,1):
                    target.write(str(cholesky_a[k][g]))
                    target.write(str('   '))
                target.write("\n")
            target.write(str(1-alpha))
            target.write(" ")
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
            for k in range(0,8,1):
                for g in range(0,8,1):
                    target.write(str(cholesky_b[k][g]))
                    target.write(str('   '))
                target.write("\n")
        no+=1
target.close()
        #    r+=1

#        if len(phi0_a)>0 and len(phi0_b)>0:
#                f, (ax1)=plt.subplots(1,1)
#                ax1.scatter(phi0,psi,s=15,c='k',marker='o')
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('psi',Fontsize=15)
#                ax1.scatter(phi0_a,psi_a,s=5,c='r',marker='o',alpha=0.2)
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('psi',Fontsize=15)
#                ax1.scatter(phi0_b,psi_b,s=5,c='g',marker='o',alpha=0.2)
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('psi',Fontsize=15)
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('psi',Fontsize=15)
#                plt.title('binned psi and phi0',Fontsize=15)
#                #plt.show()

#                f, (ax1)=plt.subplots(1,1)
#                ax1.scatter(phi0,f0,s=15,c='k',marker='o')
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('f0',Fontsize=15)
#                ax1.scatter(phi0_a,f0_a,s=5,c='r',marker='o',alpha=0.2)
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('f0',Fontsize=15)
#                ax1.scatter(phi0_b,f0_b,s=5,c='g',marker='o',alpha=0.2)
#                plt.ylim(min(f0_b),max(f0_b))
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('f0',Fontsize=15)
#                plt.xlabel('phi0',Fontsize=15)
#                plt.ylabel('f0',Fontsize=15)
#                plt.title('binned phi0 and f0',Fontsize=15)
#                #plt.show()


#        print("Try again? (y/n)")
#        ans=input()
#
#        if ans=="y":
#            phi0_hist=plt.hist(phi0[phi0<np.pi],bins=200);
##                plt.xlabel('phi0',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of initial GW phase, phi0 $<$ pi',Fontsize=15)
##                #plt.show()
#            phi0_hist_freq=list(phi0_hist[0])
#            max_freq = max(phi0_hist_freq)
#            max_index = phi0_hist_freq.index(max_freq)
#            phi0_hist_phi0=list(phi0_hist[1])
#            phi0_peak_loc1=phi0_hist_phi0[max_index]
##                plt.hist(phi0[phi0<np.pi],bins=200);
##                plt.axvline(x=phi0_peak_loc1,color='k')
##                plt.xlabel('phi0',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of initial GW phase, phi0 $<$ pi',Fontsize=15)
##                #plt.show()
#
#            phi0_peak_loc2=phi0_peak_loc1+np.pi
#
##                plt.hist(phi0,bins=200);
##                plt.axvline(x=phi0_peak_loc1,color='k')
##                plt.axvline(x=phi0_peak_loc2,color='k')
##                plt.xlabel('phi0')
##                plt.xlabel('phi0',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of initial GW phase, both peaks',Fontsize=15)
##                #plt.show()
#
#            psi_hist=plt.hist(psi[psi<np.pi/2],bins=200)
##                plt.xlabel('psi',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of GW polarization, psi $<$ pi/2',Fontsize=15)
##                #plt.show()
#            psi_hist_freq=list(psi_hist[0])
#            max_freq1 = max(psi_hist_freq)
#            max_index = psi_hist_freq.index(max_freq1)
#            psi_hist_psi=list(psi_hist[1])
#            psi_peak_loc1a=psi_hist_psi[max_index]
##                plt.hist(psi[psi<np.pi],bins=200);
##                plt.axvline(x=psi_peak_loc1a,color='k')
##                plt.xlabel('psi',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of GW polarization, first peak',Fontsize=15)
##                #plt.show()
#
#            psi_hist=plt.hist(psi[psi>np.pi/2],bins=200)
##                plt.xlabel('psi',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of GW polarization, psi $>$ pi/2',Fontsize=15)
##                #plt.show()
#            psi_hist_freq=list(psi_hist[0])
#            max_freq2 = max(psi_hist_freq)
#            max_index = psi_hist_freq.index(max_freq2)
#            psi_hist_psi=list(psi_hist[1])
#            psi_peak_loc1b=psi_hist_psi[max_index]
##                plt.hist(psi[psi<np.pi],bins=200);
##                plt.axvline(x=psi_peak_loc1b,color='k')
##                plt.xlabel('psi',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of GW polarization, psi $>$ pi/2, first peak',Fontsize=15)
##                #plt.show()
#            if max_freq2 > max_freq:
#                psi_peak_loc2=psi_peak_loc1b-np.pi/2
#                psi_peak_loc1=psi_peak_loc1b
#            else:
#                psi_peak_loc2=psi_peak_loc1a+np.pi/2
#                psi_peak_loc1=psi_peak_loc1a
#
##                plt.hist(psi,bins=200);
##                plt.axvline(x=psi_peak_loc1,color='k')
##                plt.axvline(x=psi_peak_loc2,color='k')
##                plt.xlabel('psi',Fontsize=15)
##                plt.ylabel('Frequency',Fontsize=15)
##                plt.title('Histogram of GW polarization, both peaks',Fontsize=15)
##                #plt.show()
#
#            if phi0_peak_loc1-np.pi/2 < 0:
#                phi0_min_angle_A=phi0_peak_loc1+2*np.pi-np.pi/2
#                phi0_max_angle_A=phi0_peak_loc1+2*np.pi+np.pi/2
#                phi0_peak_locA=phi0_peak_loc1+2*np.pi
#            else:
#                phi0_min_angle_A=phi0_peak_loc1-np.pi/2
#                phi0_max_angle_A=phi0_peak_loc1+np.pi/2
#                phi0_peak_locA=phi0_peak_loc1
#
#
#            if phi0_peak_loc2-np.pi/2 < 0:
#                phi0_min_angle_B=phi0_peak_loc2+2*np.pi-np.pi/2
#                phi0_max_angle_B=phi0_peak_loc2+2*np.pi+np.pi/2
#                phi0_peak_locB=phi0_peak_loc2+2*np.pi
#            else:
#                phi0_min_angle_B=phi0_peak_loc2-np.pi/2
#                phi0_max_angle_B=phi0_peak_loc2+np.pi/2
#                phi0_peak_locB=phi0_peak_loc2
#
#
#
#            if phi0_peak_locA > phi0_peak_locB:
#                phi0_peak_loc=phi0_peak_locA
#                phi0_min_angle=phi0_min_angle_A
#                phi0_max_angle=phi0_max_angle_A
#            else:
#                phi0_peak_loc=phi0_peak_locB
#                phi0_min_angle=phi0_min_angle_B
#                phi0_max_angle=phi0_max_angle_B
#
#
#
#
#            if psi_peak_loc1-np.pi/4 < 0:
#                psi_min_angle_A=psi_peak_loc1+np.pi-np.pi/4
#                psi_max_angle_A=psi_peak_loc1+np.pi+np.pi/4
#                psi_peak_locA=psi_peak_loc1+np.pi
#            else:
#                psi_min_angle_A=psi_peak_loc1-np.pi/4
#                psi_max_angle_A=psi_peak_loc1+np.pi/4
#                psi_peak_locA=psi_peak_loc1
#
#
#            if psi_peak_loc2-np.pi/4 < 0:
#                psi_min_angle_B=psi_peak_loc2+np.pi-np.pi/4
#                psi_max_angle_B=psi_peak_loc2+np.pi+np.pi/4
#                psi_peak_locB=psi_peak_loc2+np.pi
#            else:
#                psi_min_angle_B=psi_peak_loc2-np.pi/4
#                psi_max_angle_B=psi_peak_loc2+np.pi/4
#                psi_peak_locB=psi_peak_loc2
#
#
#            if psi_peak_locA > psi_peak_locB:
#                psi_peak_loc=psi_peak_locA
#                psi_min_angle=psi_min_angle_A
#                psi_max_angle=psi_max_angle_A
#            else:
#                psi_peak_loc=psi_peak_locB
#                psi_min_angle=psi_min_angle_B
#                psi_max_angle=psi_max_angle_B
#
#            angles=[0,np.pi/2,np.pi]
#            nearest_angle_index=(np.abs(psi_peak_loc-angles)).argmin()
#            angles_nearest=angles[nearest_angle_index]
#
#
#            f0_a=[]
#            f0_b=[]
#            dfdt_a=[]
#            dfdt_b=[]
#            amp_a=[]
#            amp_b=[]
#            phi_a=[]
#            phi_b=[]
#            costheta_a=[]
#            costheta_b=[]
#            cosi_a=[]
#            cosi_b=[]
#            phi0_a=[]
#            phi0_b=[]
#            psi_a=[]
#            psi_b=[]
#
#
#            phi0_c=[]
#
#            if phi0_peak_loc < np.pi/6:
#                phi0_peak_loc = 2*np.pi
#            for i in range(0,len(phi0),1):
#                min_angle_phi0A=phi0_peak_loc-np.pi/2
#                max_angle_phi0A=phi0_peak_loc+np.pi/2
#                if min_angle_phi0A < 0:
#                    min_angle_phi0A+=2*np.pi
#                    max_angle_phi0A+=2*np.pi
#
#
#            if psi_peak_loc < np.pi/6:
#                psi_peak_loc = np.pi
#                angles_nearest=np.pi
#            for i in range(0,len(psi),1):
#                min_angle_A=psi_peak_loc-np.pi/4
#                max_angle_A=psi_peak_loc+np.pi/4
#                if min_angle_A < 0:
#                    min_angle_A+=np.pi
#                    max_angle_A+=np.pi
#                    if (psi[i] + np.pi  >= min_angle_A and psi[i] + np.pi <= max_angle_A):
#                        phi0_a.append(phi0[i])
#                        phi0_c.append(phi0[i])
#                        psi_a.append(psi[i]+np.pi)
#                        f0_a.append(f0[i])
#                        dfdt_a.append(dfdt[i])
#                        amp_a.append(amp[i])
#                        phi_a.append(phi[i])
#                        costheta_a.append(costheta[i])
#                        cosi_a.append(cosi[i])
#                    elif (psi[i] >= min_angle_A and psi[i] <= max_angle_A):
#                        phi0_a.append(phi0[i])
#                        phi0_c.append(phi0[i])
#                        psi_a.append(psi[i])
#                        f0_a.append(f0[i])
#                        dfdt_a.append(dfdt[i])
#                        amp_a.append(amp[i])
#                        phi_a.append(phi[i])
#                        costheta_a.append(costheta[i])
#                        cosi_a.append(cosi[i])
#                    else:
#                        if np.absolute(psi[i]+np.pi-np.mean(psi_a)) < np.absolute(psi_b[i]+np.pi-np.mean(psi_b)):
#                            phi0_a.append(phi0[i])
#                            phi0_c.append(phi0[i])
#                            psi_a.append(psi[i]+np.pi)
#                            f0_a.append(f0[i])
#                            dfdt_a.append(dfdt[i])
#                            amp_a.append(amp[i])
#                            phi_a.append(phi[i])
#                            costheta_a.append(costheta[i])
#                            cosi_a.append(cosi[i])
#                        else:
#                            if (phi0[i]+np.pi*2 >= min_angle_phi0A and phi0[i]+np.pi*2 <= max_angle_phi0A+np.pi/4):
#                                phi0_b.append(phi0[i]+np.pi*2)
#                            else:
#                                phi0_b.append(phi0[i])
#                            psi_b.append(psi[i])
#                            f0_b.append(f0[i])
#                            dfdt_b.append(dfdt[i])
#                            amp_b.append(amp[i])
#                            phi_b.append(phi[i])
#                            costheta_b.append(costheta[i])
#                            cosi_b.append(cosi[i])
#                elif np.round(angles_nearest,2) == 3.14:
#                    if (psi[i] + np.pi  >= min_angle_A and psi[i] + np.pi <= max_angle_A) :
#                        phi0_a.append(phi0[i])
#                        phi0_c.append(phi0[i])
#                        psi_a.append(psi[i]+ np.pi)
#                        f0_a.append(f0[i])
#                        dfdt_a.append(dfdt[i])
#                        amp_a.append(amp[i])
#                        phi_a.append(phi[i])
#                        costheta_a.append(costheta[i])
#                        cosi_a.append(cosi[i])
#                    elif (psi[i] >= min_angle_A and psi[i] <= max_angle_A):
#                        phi0_a.append(phi0[i])
#                        phi0_c.append(phi0[i])
#                        psi_a.append(psi[i])
#                        f0_a.append(f0[i])
#                        dfdt_a.append(dfdt[i])
#                        amp_a.append(amp[i])
#                        phi_a.append(phi[i])
#                        costheta_a.append(costheta[i])
#                        cosi_a.append(cosi[i])
#                    else:
#                        if (phi0[i]+np.pi*2 >= min_angle_phi0A and phi0[i]+np.pi*2 <= max_angle_phi0A+np.pi/4):
#                            phi0_b.append(phi0[i]+np.pi*2)
#                        else:
#                            phi0_b.append(phi0[i])
#                        psi_b.append(psi[i])
#                        f0_b.append(f0[i])
#                        dfdt_b.append(dfdt[i])
#                        amp_b.append(amp[i])
#                        phi_b.append(phi[i])
#                        costheta_b.append(costheta[i])
#                        cosi_b.append(cosi[i])
#                else:
#                    if (psi[i] >= min_angle_A and psi[i] <= max_angle_A):
#                        phi0_a.append(phi0[i])
#                        phi0_c.append(phi0[i])
#                        psi_a.append(psi[i])
#                        f0_a.append(f0[i])
#                        dfdt_a.append(dfdt[i])
#                        amp_a.append(amp[i])
#                        phi_a.append(phi[i])
#                        costheta_a.append(costheta[i])
#                        cosi_a.append(cosi[i])
#                    else:
#                        if (phi0[i]+np.pi*2 >= min_angle_phi0A and phi0[i]+np.pi*2 <= max_angle_phi0A+np.pi/4):
#                            phi0_b.append(phi0[i]+np.pi*2)
#                        else:
#                            phi0_b.append(phi0[i])
#                        psi_b.append(psi[i])
#                        f0_b.append(f0[i])
#                        dfdt_b.append(dfdt[i])
#                        amp_b.append(amp[i])
#                        phi_b.append(phi[i])
#                        costheta_b.append(costheta[i])
#                        cosi_b.append(cosi[i])
#
#
#
#
#            if len(phi0_b)>0:
#                if np.round(max(phi0_b)-min(phi0_b),2)==6.28:
#                    for k in range(0,len(phi0_b),1):
#                        if phi0_b[k]<phi0_peak_loc-np.pi:
#                            phi0_b[k]+=np.pi
#
#            if len(psi_b)>0:
#                if np.round(max(psi_b)-min(psi_b),2)==3.14:
#                    for k in range(0,len(psi_b),1):
#                        if psi_b[k]<psi_peak_loc-np.pi/4:
#                            if psi_b[k] + np.pi < psi_peak_loc+np.pi/4:
#                                psi_b[k]+=np.pi
#
#            if np.round(max(psi_b)-min(psi_b),2) == 3.14:
#                for k in range(0,len(psi_b),1):
#                    if psi_b[k] > np.pi/2 and psi_b[k]-np.pi>0:
#                        psi_b[k]-=np.pi
#
#            if np.round(max(psi_b)-min(psi_b),2) == 3.14:
#                for k in range(0,len(psi_b),1):
#                    if psi_b[k] < np.pi/2:
#                        psi_b[k]+=np.pi
#                for i in range(0,len(phi_b),1):
#                    if np.absolute(psi_b[i]-np.mean(psi_b))> np.pi/2:
#                        if psi_b[i] <np.mean(psi_b):
#                            psi_b[i]+=np.pi
#
#            if len(phi0_a)>0:
#                if len(phi0_a)>0:
#                    if np.round(max(phi0_a)-min(phi0_a),2)==6.28:
#                        for k in range(0,len(phi0_a),1):
#                            if phi0_a[k]<phi0_peak_loc-np.pi:
#                                phi0_a[k]+=np.pi
#
#
#
#
#
#            if len(psi_a)>0:
#                if np.round(max(psi_a)-min(psi_a),2)==3.14:
#                    for k in range(0,len(psi_a),1):
#                        if psi_a[k]<psi_peak_loc-np.pi/4:
#                            if psi_a[k]<psi_peak_loc-np.pi/4:
#                                if psi_a[k] + np.pi < psi_peak_loc+np.pi/4:
#                                    psi_a[k]+=np.pi
#
#                    if np.round(max(psi_a)-min(psi_a),2) == 3.14:
#                        for k in range(0,len(psi_a),1):
#                            if psi_a[k] > np.pi/2 and psi_a[k]-np.pi>0:
#                                psi_a[k]-=np.pi
#
#                    if np.round(max(psi_a)-min(psi_a),2) == 3.14:
#                        for k in range(0,len(psi_a),1):
#                            if psi_a[k] < np.pi/2 :
#                                psi_a[k]+=np.pi
#
#
#        else:
#            pass


