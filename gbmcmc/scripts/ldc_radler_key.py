import h5py
import numpy as np
import sys, os, math
import sys, getopt
import cmath

print("Enter LDC Radler galaxy input file:")
inputfile=input()
    
f = h5py.File(inputfile, 'r')
freq_ldc=f['H5LISA']['GWSources']['GalBinaries']['Frequency']
while len(freq_ldc)<10:
    print(len(freq_ldc))
    print("This is not the Radler Galaxy LDC file, enter Radler input file:")
    inputfile=input()
    f = h5py.File(inputfile, 'r')
    freq_ldc=f['H5LISA']['GWSources']['GalBinaries']['Frequency']

outputfile='./LDC_Radler_galaxy_key.dat'
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

freq_ldc=freq_ldc[:]
indx_f = freq_ldc.argsort()
freq_ldc=freq_ldc[indx_f[::1]]
amp_ldc=amp_ldc[:]
amp_ldc=amp_ldc[indx_f[::1]]
ecl_lat_ldc=ecl_lat_ldc[:]
ecl_lat_ldc=ecl_lat_ldc[indx_f[::1]]
ecl_long_ldc=ecl_long_ldc[:]
ecl_long_ldc=ecl_long_ldc[indx_f[::1]]
dot_freq_ldc=dot_freq_ldc[:]
dot_freq_ldc=dot_freq_ldc[indx_f[::1]]
inc_ldc=inc_ldc[:]
inc_ldc=inc_ldc[indx_f[::1]]
phi0_ldc=phi0_ldc[:]
phi0_ldc=phi0_ldc[indx_f[::1]]
pol_ldc=pol_ldc[:]
pol_ldc=pol_ldc[indx_f[::1]]

#Specify the minimum and maximum frequency, fmin and fmax. (NB: If you have already run gb_mcmc (discussed below) on a particular frequency range and would like to know the corresponding source parameters that exist in that range, then use the fmin and fmax defined for the data analyzed.)

print("Enter minimum frequency:")
fmin=input()

print("Enter maximum frequency:")
fmax=input()

fmin=float(fmin)
fmax=float(fmax)

while fmin>fmax:
    print('Error: fmin > fmax')
    print("Enter minimum frequency:")
    fmin=input()

    print("Enter maximum frequency:")
    fmax=input()
    
    fmin=float(fmin)
    fmax=float(fmax)

print(fmin,fmax)
srcs_f0=freq_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
no_srcs_f0=len(srcs_f0)
f_ldc=srcs_f0
amp_ldc=amp_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
ecl_lat_key=ecl_lat_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
ecl_lat_ldc=[]
for i in ecl_lat_key:
    ecl_lat_ldc.append(np.cos(np.pi/2.0-i))
ecl_long_key=ecl_long_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
ecl_long_ldc=[]
for i in ecl_long_key:
    ecl_long_ldc.append(i)
for i in range(0,len(ecl_long_ldc),1):
    if ecl_long_ldc[i] < 0:
        ecl_long_ldc[i]+=2*np.pi
dfdt_ldc=dot_freq_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
iota_key=inc_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
iota_ldc=[]
for i in iota_key:
    iota_ldc.append(np.cos(i))
phi0_ldc=phi0_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
psi_key=pol_ldc[len(freq_ldc[freq_ldc<fmin]):len(freq_ldc[freq_ldc<fmax])]
pol_ldc=[]
for i in psi_key:
    pol_ldc.append(i)
for i in range(0,len(pol_ldc),1):
    if pol_ldc[i] > np.pi:
        pol_ldc[i]-=np.pi

print("            Output is being saved to ./LDC_Radler_galaxy_key.dat")
        
target=open(outputfile,'w')
for i in range(0,len(pol_ldc),1):
    target.write(str(srcs_f0[i]))
    target.write("        ")
    target.write(str(amp_ldc[i]))
    target.write("        ")
    target.write(str(dfdt_ldc[i]))
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
