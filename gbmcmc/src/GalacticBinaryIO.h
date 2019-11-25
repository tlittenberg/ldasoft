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



#ifndef GalacticBinaryIO_h
#define GalacticBinaryIO_h

void printProgress (double percentage);

void print_version(FILE *fptr);

void print_usage();
void parse(int argc, char **argv, struct Data **data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, int Nmax);
int checkfile(char filename[]);

void save_chain_state(struct Data **data, struct Model ***model, struct Chain *chain, struct Flags *flags, int step);
void print_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int step);
void scan_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int *step);\
void print_noise_state(struct Data *data, struct Model *model, FILE *fptr, int step);
void scan_noise_state(struct Data *data, struct Model *model, FILE *fptr, int *step);
void print_calibration_state(struct Data *data, struct Model *model, FILE *fptr, int step);
void scan_calibration_state(struct Data *data, struct Model *model, FILE *fptr, int *step);
void print_source_params(struct Data *data, struct Source *source, FILE *fptr);
void scan_source_params(struct Data *data, struct Source *source, FILE *fptr);
void print_chain_files(struct Data *data, struct Model ***model, struct Chain *chain, struct Flags *flags, int step);

void save_waveforms(struct Data *data, struct Model *model, int mcmc);
void print_waveform(struct Data *data, struct Model *model, FILE *fptr);
void print_waveforms_reconstruction(struct Data *data, int seg);
void print_waveform_draw(struct Data **data, struct Model **model, struct Flags *flags);


#endif /* GalacticBinaryIO_h */
