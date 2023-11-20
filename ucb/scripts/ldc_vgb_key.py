import h5py
import numpy as np
import sys, os, math
import sys, getopt
import cmath

print("Enter VGB input file:")
inputfile=input()
    
f = h5py.File(inputfile, 'r')
freq_ldc=f['H5LISA']['GWSources']['GalBinaries']['Frequency']
while len(freq_ldc)>10:
    print(len(freq_ldc))
    print("This is not the VGB LDC file, enter VGB input file:")
    inputfile=input()
    f = h5py.File(inputfile, 'r')
    freq_ldc=f['H5LISA']['GWSources']['GalBinaries']['Frequency']

outputfile='./LDC_VGB_key.dat'
print('')
print('			Input file is "', inputfile)
print('')
print('			Output file is "', outputfile)
print('')
 
print('')
print('			Loading LDC Fourier domain data...')
print('')
f = h5py.File(inputfile, 'r')

freq_ldc=f['H5LISA']['GWSources']['GalBinaries']['Frequency']
amp_ldc=f['H5LISA']['GWSources']['GalBinaries']['Amplitude']
ecl_lat_ldc=f['H5LISA']['GWSources']['GalBinaries']['EclipticLatitude']
ecl_long_ldc=f['H5LISA']['GWSources']['GalBinaries']['EclipticLongitude']
dot_freq_ldc=f['H5LISA']['GWSources']['GalBinaries']['FrequencyDerivative']
inc_ldc=f['H5LISA']['GWSources']['GalBinaries']['Inclination']
phi0_ldc=f['H5LISA']['GWSources']['GalBinaries']['InitialPhase']
pol_ldc=f['H5LISA']['GWSources']['GalBinaries']['Polarization']

#Sort LDC sources in order of increasing GW frequency

freq_ldc=freq_ldc[:]
indx_f = freq_ldc.argsort()
freq_ldc=freq_ldc[indx_f[::1]]
amp_ldc=amp_ldc[:]
amp_ldc=amp_ldc[indx_f[::1]]
ecl_lat_key=ecl_lat_ldc[:]
ecl_lat_key=ecl_lat_key[indx_f[::1]]
ecl_long_key=ecl_long_ldc[:]
ecl_long_key=ecl_long_key[indx_f[::1]]
dot_freq_ldc=dot_freq_ldc[:]
dot_freq_ldc=dot_freq_ldc[indx_f[::1]]
iota_key=inc_ldc[:]
iota_key=iota_key[indx_f[::1]]
phi0_ldc=phi0_ldc[:]
phi0_ldc=phi0_ldc[indx_f[::1]]
psi_key=pol_ldc[:]
psi_key=psi_key[indx_f[::1]]




ecl_lat_ldc=[]
for i in ecl_lat_key:
    ecl_lat_ldc.append(np.cos(np.pi/2.0-i))
ecl_long_ldc=[]
for i in ecl_long_key:
    ecl_long_ldc.append(i)
for i in range(0,len(ecl_long_ldc),1):
    if ecl_long_ldc[i] < 0:
        ecl_long_ldc[i]+=2*np.pi
iota_ldc=[]
for i in iota_key:
    iota_ldc.append(np.cos(i))
pol_ldc=[]
for i in psi_key:
    pol_ldc.append(i)
for i in range(0,len(pol_ldc),1):
    if pol_ldc[i] > np.pi:
        pol_ldc[i]-=np.pi
        
print("            Output is being saved to ./LDC_VGB_key.dat")
        
target=open(outputfile,'w')
for i in range(0,len(pol_ldc),1):
    target.write(str(freq_ldc[i]))
    target.write("        ")
    target.write(str(amp_ldc[i]))
    target.write("        ")
    target.write(str(dot_freq_ldc[i]))
    target.write("        ")
    target.write(str(iota_ldc[i]))
    target.write("        ")
    target.write(str(phi0_ldc[i]))
    target.write("        ")
    target.write(str(ecl_lat_ldc[i]))
    target.write("        ")
    target.write(str(ecl_long_ldc[i]))
    target.write("        ")
    target.write(str(pol_ldc[i]))
    target.write("        ")
    
    target.write("\n")
target.close()
