/*
 * Copyright 2021 Tyson B. Littenberg
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
 @file glass_noise_io.h
 \brief Functions handling input/output for noise sampler.

 */
#ifndef noise_io_h
#define noise_io_h

/**
 \brief Print current state of spline model to ASCII
 */
void print_spline_state(struct SplineModel *model, FILE *fptr, int step);

/**
 \brief Print current state of instrument model to ASCII
 */
void print_instrument_state(struct InstrumentModel *model, FILE *fptr);

/**
 \brief Print current state of instrument model to ASCII
 */
void print_foreground_state(struct ForegroundModel *model, FILE *fptr);

/**
 \brief Print full PSD model to file named `filename`
 */
void print_noise_model(struct Noise *noise, char filename[]);

/**
 \brief Print data whitened by modeled variance to file named `filename`
 */
void print_whitened_data(struct Data *data, struct Noise *noise, char filename[]);

/**
 \brief Print PSD 5, 25, 50, 75, and 95 quantiles of noise model posteriors
 */
void print_noise_reconstruction(struct Data *data, struct Flags *flags);

#endif /* noise_io_h */
