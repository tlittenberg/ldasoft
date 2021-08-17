/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 @file GalacticBinaryIO.h
 \brief Functions handling I/O for `gbmcmc`.
 */


#ifndef GalacticBinaryIO_h
#define GalacticBinaryIO_h

/**
 \brief Show progress bar
 */
void printProgress (double percentage);

/**
 \brief Print git hash of running version
 */
void print_version(FILE *fptr);

/**
 \brief Print functional example `gb_catalog` bash script based on input args
 */
void print_gb_catalog_script(struct Flags *flags, struct Data *data, struct Orbit *orbit);

/**
 \brief Print list of command line arguments and defaults
 */
void print_usage();

/**
 \brief Parse command line
 */
void parse(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, int Nmax, int numProc, int procID);

/**
 \brief select only part of command line defining verification binary lists
*/
void parse_vb_list(int argc, char **argv, struct Flags *flags);

/**
 \brief copy `argv` string because `parse()` does not preserve order
 */
void copy_argv(int argc, char **argv, char **new_argv);

/**
 \brief Check if file exists and handle missing files cleanly
 */
int checkfile(char filename[]);

/** @name Full Sampler State Files
 Print/read to file the full chain state for checkpointing
 */
///@{
void save_chain_state(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step);
void restore_chain_state(struct Orbit *orbit, struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int *step);
///@}

/** @name Chain State File
 Print/read current state of sampler, e.g. to Chain::chainFile
 */
///@{
void print_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int step);
void scan_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int *step);
///@}

/** @name Noise Chain File
 Print/read current state of noise model, e.g. to Chain::noiseFile
 */
///@{
void print_noise_state(struct Data *data, struct Model *model, FILE *fptr, int step);
void scan_noise_state(struct Data *data, struct Model *model, FILE *fptr, int *step);
///@}

/** @name Calibration Chain File
 **WORK IN PROGRESS** Print/read current state of calibration model, e.g. to Chain::calibrationFile
 */
///@{
void print_calibration_state(struct Data *data, struct Model *model, FILE *fptr, int step);
void scan_calibration_state(struct Data *data, struct Model *model, FILE *fptr, int *step);
///@}

/** @name Galactic Binary Chain File
 Print/read current state of source model, e.g. to Chain::parameterFile
 */
///@{
void print_source_params(struct Data *data, struct Source *source, FILE *fptr);
void scan_source_params(struct Data *data, struct Source *source, FILE *fptr);
///@}

/**
 \brief Wrapper function that calls all of the chain print functions
 */
void print_chain_files(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step);

/** @name Waveform Files
 Save various representations of waveform reconstructions
 */
///@{

/// Copy current state of reconstructed waveform (strain and power) and noise model to arrays to be used to compute reconstruction quantiles
void save_waveforms(struct Data *data, struct Model *model, int mcmc);

/// Print waveform strain
void print_waveform(struct Data *data, struct Model *model, FILE *fptr);

/// Print waveform 5, 25, 50, 75, and 95 quantiles of waveform power spectrum posteriors
void print_waveforms_reconstruction(struct Data *data, struct Flags *flags);

/// Print current state of waveform and residuals during run for diagnostics. Disabled when Flags::quiet=`TRUE`.
void print_waveform_draw(struct Data *data, struct Model *model, struct Flags *flags);
///@}

/**
 \brief Print PSD 5, 25, 50, 75, and 95 quantiles of noise model posteriors
 */
void print_noise_reconstruction(struct Data *data, struct Flags *flags);


/**
 \brief Wrapper function that calls data print functions
 */
void print_data(struct Data *data, struct TDI *tdi, struct Flags *flags, int t_index);

/**
 \brief Print evidence for each dimension model
 */
void print_evidence(struct Chain *chain,struct Flags *flags);

#endif /* GalacticBinaryIO_h */
