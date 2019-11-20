# Fourier domain analysis of the LISA Data Challenge (LDC) data (Radler Galaxy and 10 Galactic Binary (GB) sources) with GalacticBinaryMCMC.c (gb_mcmc).


# Download LDC GB data
Username, password are required to [access challenge data](https://lisa-ldc.lal.in2p3.fr/) 
    
    Contact for LISA Data Challenges (LDC):  lisa-ldc-helpdesk-l@in2p3.fr

Filename: LDC1-3_VGB_v1.hdf5
    
    Description: two year `verification' GB timeseries with noise (10 GB's)
    GB parameters included in hdf5 file (LDC VGB key)

Filename: LDC1-4_GB_v1.hdf5
    
    Description: two year Radler Galaxy timeseries, ~30 million GB's with noise (file size is 2GB)
    GB parameters included in hdf5 file (LDC Galaxy key)

# FFT LDC data with fft_ldc_data.py
The output of fft_ldc_data.py is LDC data in the `frequency domain', which is used as input to gb_mcmc
Enter fft_ldc_data.py -h, for usage

    For example, to FFT 1.5 months of LDC1-3_VGB_v1.hdf5 data and save the output as output.dat, run from terminal the following. 
    python3 fft_ldc_data.py -i LDC1-3_VGB_v1.hdf5 -o output.dat -m 1.5

Specifying an output file (with the -o flag) is optional. If the output file is not specified, the code's default (output file) nameing system is the following: The input file name is affixed with the number of months that were Fourier transformed, e.g. LDC1-3_VGB_v1_1.5mo.dat.
(If -m flag is left unspecified, the code's default is to FT 24 months of data.)

One can also FFT the LDC Galaxy data timeseries: python3 fft_ldc_data.py -i LDC1-4_GB_v1.hdf5 -o output.dat -m 1.5
In that case, the output filename is LDC1-4_GB_v1_1.5mo.dat

# Access LDC VGB key parameters with python
The user inputs ~/LDC1-3_VGB_v1.hdf5 to the code, and the output is LDC key information, in order of increasing source frequency. Output is saved to the user's working directory. 
ldc_vgb_key.py

# Access LDC Galaxy key parameters, with python
The user inputs ~/LDC1-4_GB_v1.hdf5 to the code, and the code outputs key parameters corresponding to a specific frequency range. For example, one could input the frequency range analyzed by gb_mcmc, described next. Output is saved to the user's working directory. 
ldc_radler_key.py


# Run gb_mcmc on Fourier transformed LDC data.
The following assumes the user is analyzing a frequency segment of LDC1-3_VGB_v1_{n}mo.dat, or LDC1-4_GB_v1_{n}mo.dat, 
where n = {1.5, 3, 6, 12, 24}. The input data are obtained from fft_ldc_data.py, described above.

Select a specific frequency range for analysis by setting flags --fmin (Hz) and --samples (integer).

    One set of choices for --samples are 2^5, 2^6, 2^7, 2^8, 2^9 for analyzing 1.5, 3, 6, 12, 24 months of LDC data, respectively.
    One may prefer the range a freq. factor of two larger. In that case, --samples becomes 2^6, 2^7, 2^8, 2^9, 2^10, for the time periods listed in the line above.
    
 Specify flag --duration for an integer number of seconds 

    (e.g. use 3932160, 7864320, 15728640, 31457280, and 62914560 for 1.5, 3, 6, 12, 24 months, respectively)   

For example, to analyze 1.5 months of data for the lowest frequency VGB source, the user could run the following.

    ~/gb_mcmc --data ~/LDC1-3_VGB_v1_1.5_months.dat --frac-freq --fmin 0.001249 --samples 32 --duration 3932160 --chainseed 1234 --verbose

# Use GalacticBinaryCatalog.c (gb_catalog) to produce catalog data with the gb_mcmc output 







Below we show how to use the catalog output to build covariance matrix proposals for gb_mcmc.

# Build covariance matrix proposal from catalog output




