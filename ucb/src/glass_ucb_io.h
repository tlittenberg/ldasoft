/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 @file glass_ucb_io.h
 \brief Functions handling I/O for `ucb_mcmc`.
 */


#ifndef ucb_io_h
#define ucb_io_h

/**
 \brief Print command line options for UCB module
 */
void print_ucb_usage();

/**
 \brief Parse part of command line setting UCB module args
*/
void parse_ucb_args(int argc, char **argv, struct Flags *flags);

/**
 \brief Parse part of command line defining verification binary lists
*/
void parse_vgb_args(int argc, char **argv, struct Flags *flags);

/**
 \brief Prints constants and flags for the MCMC run
 */
void print_run_settings(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, FILE *fptr);

/**
 \brief Print functional example `ucb_catalog` bash script based on input args
 */
void print_ucb_catalog_script(struct Flags *flags, struct Data *data, struct Orbit *orbit);

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
 Print/read current state of psd model, e.g. to Chain::noiseFile
 */
///@{
void print_psd_state(struct Data *data, struct Model *model, FILE *fptr, int step);
void scan_psd_state(struct Data *data, struct Model *model, FILE *fptr, int *step);
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
void print_waveform_strain(struct Data *data, struct Model *model, FILE *fptr);

/// Print waveform power
void print_waveform(struct Data *data, struct Model *model, FILE *fptr);

/// Print waveform 5, 25, 50, 75, and 95 quantiles of waveform power spectrum posteriors
void print_waveforms_reconstruction(struct Data *data, struct Flags *flags);

/// Print current state of waveform and residuals during run for diagnostics. Disabled when Flags::quiet=`TRUE`.
void print_waveform_draw(struct Data *data, struct Model *model, struct Flags *flags);
void print_psd_draw(struct Data *data, struct Model *model, struct Flags *flags);
///@}

/**
 \brief Print evidence for each dimension model
 */
void print_evidence(struct Chain *chain,struct Flags *flags);

#endif /* ucb_io_h */
