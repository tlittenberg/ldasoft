This directory contains data files packaged with the installation of ldasoft

global_fit needs these in order to run

Currently they are:

`ucb_frequency_spacing.dat` -> Needs to live in the working directory of the global_fit process. This specifies the frequency ranges used by the individual galactic binary samplers. If it is not provided in the current working directory, global_fit will generate ranges.
`ldc_sangria_vgb_list.dat` -> Needs to be provided on the command line as `--known-sources ldc_sangria_vgb_list.dat`
`search_sources.dat` -> The directory of this file ($INSTALL_DIR/data) needs to be provided on the command line as `--mbh-search-path $INSTALL_DIR/data`


