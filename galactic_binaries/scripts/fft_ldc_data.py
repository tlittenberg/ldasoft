import numpy as np
import h5py
import sys, os, math
import sys, getopt
import cmath

def compute_dft_complex(input):
    n = len(input)
    output = []
    for k in range(n):  # For each output element
        print(k/n)
        s = complex(0)
        for t in range(n):  # For each input element
            angle = 2j * cmath.pi * t * k / n
            s += input[t] * cmath.exp(-angle)
        output.append(s)
    return output


def main(argv):
    inputfile = ''
    outputfile = ''
    months = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:",["ifile=","ofile=","duration="])
    except getopt.GetoptError as err:
        print(err)
        print('ldc_gbmcmc.py -i <inputfile> -o <outputfile> -m <months>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('ldc_gbmcmc.py -i <inputfile> -o <outputfile> -m <months>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
           inputfile = arg
        elif opt in ("-o", "--ofile"):
           outputfile = arg
        elif opt in ("-m", "--duration"):
           months = arg
        if months == '':
           months = str(24)
    if inputfile == '':
        print("Enter input file:")
        inputfile=input()
    if outputfile == '':
        if np.float(months) > 24:
            outputfile=inputfile[:-5]+'_24_months.dat'
        else:
            outputfile=inputfile[:-5]+'_'+str(months)+'_months.dat'
    print('')
    print('			Input file is "', inputfile)
    print('')
    print('			Output file is "', outputfile)
    print('')
    if np.float(months) > 24:
        print('			Duration is 24 months')
    else:
        print('			Duration is', months, 'months')
    print('')
    print('')

 
    print('')
    print('			Loading LDC t,X,Y,Z data and transforming to t,A,E,T data...')
    print('')
    f = h5py.File(inputfile, 'r')
 
    time=[]
    x=[]
    y=[]
    z=[]
    A=[]
    E=[]
    T=[]
    for i in f['H5LISA']['PreProcess']['TDIdata'][:]:
       time.append(i[0])
       x.append(i[1])
       y.append(i[2])
       z.append(i[3])
       A.append((2*i[1]-i[2]-i[3])/3)
       E.append((i[3]-i[2])/np.sqrt(3))
       T.append((i[1]+i[2]+i[3])/3)

    print('')
    print('')
    print('            Now doing a Fourier transform of A(t),E(t) data with fft.c to produce f,A(f),E(f) for your selected observation time...')
    print('')
    print('')

    try:
      div = 24.0/np.float(months)
      if div < 1.0:
         div = 1.0
         print('')
         print('Months cannot be greater than 24, reverting to 24 month observation')
         print('')
    except(ValueError):
      div = 1.0
      pass



    no=int(len(time)/div)
    a=A[0:no]
    e=E[0:no]
    t=T[0:no]
    t_indx=time[0:no]

    a=np.asarray(a)
    e=np.asarray(e)
    t=np.asarray(t)
    print('')
    print('')
    print('			Observation duration is ', np.round(no*15/(2617529),1), 'months, or', np.round(no*15), 'seconds' )
    print('')
    print('')


    #python FFT

    fSa_py=np.fft.fft(a,norm=None)
    fSe_py=np.fft.fft(e,norm=None)
    f_band=np.fft.fftfreq(len(fSa_py),d=15.0)
    freq_py=f_band[0:int(len(fSa_py)/2)+1]
    norm=np.sqrt(16/len(fSa_py))
    fSa_py_positive=fSa_py[0:int(len(fSa_py)/2)]
    fSe_py_positive=fSe_py[0:int(len(fSa_py)/2)]


    print('')
    print('              Saving data as ', outputfile)
    print('')

    with open(outputfile, 'w') as fileout:
       for k in range(1,len(freq_py),1):
            try:
                fileout.write(str(freq_py[k]))
                fileout.write('    ')
                fileout.write(str(norm*fSa_py_positive.real[k-1]))
                fileout.write('    ')
                fileout.write(str(norm*fSa_py_positive.imag[k-1]))
                fileout.write('    ')
                fileout.write(str(norm*fSe_py_positive.real[k-1]))
                fileout.write('    ')
                fileout.write(str(norm*fSe_py_positive.imag[k-1]))
                fileout.write('\n')
            except(IndexError,ValueError):
                fileout.close()
                pass





if __name__ == "__main__":
   main(sys.argv[1:])







#some code for `computationally expensive` (slow) dft
#    a=a.astype(complex)
#    e=e.astype(complex)
#
#
#    fSa=compute_dft_complex(a)
#    fSe=compute_dft_complex(e)
#    f_band=np.fft.fftfreq(len(fSa),d=15.0)
#    freq=f_band[0:int(len(fSa)/2)]
#    fSa_positive=fSa[0:int(len(fSa)/2)]
#    fSe_positive=fSe[0:int(len(fSe)/2)]
#
#   print('')
#   print('            Saving data as ', outputfile)
#   print('')
#
#   outputfile1='/Users/klackeos/data/a.dat'
#   with open(outputfile1, 'w') as fileout:
#     for k in range(0,len(fSa),1):
#        fileout.write(str(fSa[k].real))
#        fileout.write('   ')
#        fileout.write(str(fSa[k].imag))
#        fileout.write('\n')
#   fileout.close()
#
#   outputfile2='/Users/klackeos/data/a.dat'
#   with open(outputfile2, 'w') as fileout:
#      for k in range(0,len(fSe),1):
#        fileout.write(str(fSe[k].real))
#        fileout.write('   ')
#        fileout.write(str(fSe[k].imag))
#        fileout.write('\n')
#    fileout.close()
#
#
#    f2=open(outputfile1,'r')
#    data=np.loadtxt(f2)
#    fa_r=data[:,0]
#    fa_i=data[:,1]
#    f2=open(outputfile2,'r')
#    data=np.loadtxt(f2)
#    fe_r=data[:,0]
#    fe_i=data[:,1]
#
#    norm=np.sqrt(16/len(fe_r))
#    #norm=1
#    with open(outputfile, 'w') as fileout:
#        for k in range(0,len(fe_r[0:int(len(fe_r)/2)]),1):
#            fileout.write(str(freq[k]))
#            fileout.write('   ')
#            fileout.write(str(fa_r[k]*norm))
#            fileout.write('   ')
#            fileout.write(str(fa_i[k]*norm))
#            fileout.write('   ')
#            fileout.write(str(fe_r[k]*norm))
#            fileout.write('   ')
#            fileout.write(str(fe_i[k]*norm))
#            fileout.write('\n')
#        fileout.close()
