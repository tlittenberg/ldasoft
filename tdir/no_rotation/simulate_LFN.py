import numpy as np
import math
import sys
import time




def filter(M,f_max,fs):

	f_c = f_max/fs
	#k = 1/(2*np.pi*f_c)


	h = np.empty(M+1)
	#symmetric point set to this to avoid divide by zero
	h[int(M/2)] = 2.0*math.pi*f_c
	#then either side is equation as normal as given in dspguide.com ch 16
	h[:int(M/2)] = (0.42-0.5*np.cos(2.0*np.pi*np.arange(int(M/2))/M)+0.08*np.cos(4.0*np.pi*np.arange(int(M/2))/M))*np.sin(2.0*np.pi*f_c*(np.arange(int(M/2))-M/2))/(np.arange(int(M/2))-M/2)
	h[int(M/2+1):] = (0.42-0.5*np.cos(2.0*np.pi*np.arange(int(M/2+1),M+1)/M)+0.08*np.cos(4.0*np.pi*np.arange(int(M/2+1),M+1)/M))*np.sin(2.0*np.pi*f_c*(np.arange(int(M/2+1),M+1)-M/2))/(np.arange(int(M/2+1),M+1)-M/2)

	sum = np.sum(h)
	h = h/sum

	return h

#1 day
dur = 1*24*3600
# 1 MHz
f_samp = 1e6
#Down-sampled sampling rate
f_s = 1
#Boundaries of the band-pass
f_min = 1e-4
f_max = 1e-1


#plus 5 for these specific parameters. ((N/seg)+m-1 = total segment, which must be 2^N)
#N = int(dur*f_samp)+1
#N = int(dur*f_samp)+2
N = int(dur*f_samp)+5


#number of segments
#seg = 103
#seg=17
#seg = 2594
seg = 319

print('seg')
print(seg)
#number in each segment
n = int(N/seg)
print('n')
print(n)
#number in filter
#m = 33801232
m = 266024518
#m = 3306256
#m = 208742
#m = 800
print('m')
print(m)

total_segment = n + m -1
print('total segment')
print(total_segment)

h_low = filter(m,0.1,f_samp)
print('low-pass portion created')


h_high  = filter(m,1.0e-4,f_samp)
print('high-pass portion created')

#spectral inversion
h_high = -1*h_high
h_high[int(m/2)]+=1

#create band-reject by adding the two
h_ = h_low + h_high

#spectral inversion of band-reject to create band-pass
h_inv = -1*h_
h_inv[int(m/2)]+=1

#pad the filter with extra zeros to make total_segment
extra = total_segment-len(h_inv)
h_pad = np.pad(h_inv,(0,extra),'constant')


h = np.fft.rfft(h_pad)
print('band-pass filter created')

#.....................................	
#Shift and add zeroes to file storing band-passed data
#.....................................

#!!!!!!!!!!!!!!!!!!!!!MAKE SURE THAT SAMPLE DELAYS ARE LESS THAN M !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#simulated delay times in seconds
L_1 = 8.339098
L_2 = 8.339104


#number of samples to create the delay
sample_delay_1 = int(2.0*L_1*f_samp)
sample_delay_2 = int(2.0*L_2*f_samp)

#number between each sample for downsampling
down_sample_factor = int(f_samp/f_s)



#.....................................	
#Overlap-add processing
#.....................................


start_time = time.time()

ds_1 = np.array([])
ds_2 = np.array([])


extra = m-1
count = 0

print('begin band-passing of segments')
for i in np.arange(seg):
	print('segment number {0} out of {1} segments'.format(i,seg))

	input = np.random.normal(0,1e-13,n)
	input = np.pad(input,(0,extra),'constant')

	f_domain = np.fft.rfft(input,norm='ortho')	
	convolved = f_domain*h

	time_domain = np.fft.irfft(convolved,norm='ortho')

	if i ==0:

		for q, p in enumerate(time_domain):

			if q<n:

				if count%down_sample_factor == 0:
					
					if count-sample_delay_1 < 0:
						ds_1 = np.append(ds_1,p-0.0)
					else:
						difference = q - sample_delay_1
						if difference < 0:

							ds_1 = np.append(ds_1, p-last_saved[difference])
						else:
							ds_1 = np.append(ds_1,p-time_domain[q-sample_delay_1])
					if count-sample_delay_2 < 0:
						ds_2 = np.append(ds_2,p-0.0)
					else:
						difference = q - sample_delay_2
						if difference < 0:
							ds_2 = np.append(ds_2, p-last_saved[difference])
						else:
							ds_2 = np.append(ds_2,p-time_domain[q-sample_delay_2])
				count+=1

	elif i>0 and i<seg-1:

		for index in np.arange(extra):
			time_domain[index]+=overlap_array[index]
		for q, p in enumerate(time_domain):

			if q<n:
				if count%down_sample_factor == 0:
					if count-sample_delay_1 < 0:
						ds_1 = np.append(ds_1,p-0.0)
					else:
						difference = q - sample_delay_1
						if difference < 0:

							ds_1 = np.append(ds_1, p-last_saved[difference])
						else:
							ds_1 = np.append(ds_1,p-time_domain[q-sample_delay_1])
					if count-sample_delay_2 < 0:
						ds_2 = np.append(ds_2,p-0.0)
					else:
						difference = q - sample_delay_2
						if difference < 0:
							ds_2 = np.append(ds_2, p-last_saved[difference])
						else:
							ds_2 = np.append(ds_2,p-time_domain[q-sample_delay_2])
				count+=1

	elif i==seg-1:
		for index in np.arange(extra):
			time_domain[index]+=overlap_array[index]
		for q, p in enumerate(time_domain):

			if count%down_sample_factor == 0:
				if count-sample_delay_1 < 0:
					ds_1 = np.append(ds_1,p-0.0)
				else:
					difference = q - sample_delay_1
					if difference < 0:

						ds_1 = np.append(ds_1, p-last_saved[difference])
					else:
						ds_1 = np.append(ds_1,p-time_domain[q-sample_delay_1])
				if count-sample_delay_2 < 0:
					ds_2 = np.append(ds_2,p-0.0)
				else:
					difference = q - sample_delay_2
					if difference < 0:
						ds_2 = np.append(ds_2, p-last_saved[difference])
					else:
						ds_2 = np.append(ds_2,p-time_domain[q-sample_delay_2])
			count+=1
		
	overlap_array = time_domain[n::]
	last_saved = time_domain[:n]
print('done band-passing')
print("Took %s seconds to band-pass all data" % (time.time() - start_time))

#Only down-sampled data is saved to file so saving entire array at the end is fine.
np.savetxt('d_1.dat',ds_1)
np.savetxt('d_2.dat',ds_2)

