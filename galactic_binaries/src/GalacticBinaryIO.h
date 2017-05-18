//
//  GalacticBinaryModel.h
//
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//

#ifndef GalacticBinaryIO_h
#define GalacticBinaryIO_h

void print_usage();
void parse(int argc, char **argv, struct Data ***data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, int Nmax);

void print_chain_files(struct Data *data, struct Model ****model, struct Chain *chain, struct Flags *flags, int step);
void print_chain_state(struct Data *data, struct Chain *chain, struct Model **model, struct Flags *flags, FILE *fptr, int step);
void print_noise_state(struct Data *data, struct Model *model, FILE *fptr, int step);
void print_source_params(struct Data *data, struct Source *source, FILE *fptr);
void scan_source_params(struct Data *data, struct Source *source, FILE *fptr);

void save_waveforms(struct Data *data, struct Model *model, int mcmc);
void print_waveform(struct Data *data, struct Model *model, FILE *fptr);
void print_waveforms_reconstruction(struct Data *data, int seg);
void print_waveform_draw(struct Data ***data, struct Model ***model, struct Flags *flags);

#endif /* GalacticBinaryIO_h */
