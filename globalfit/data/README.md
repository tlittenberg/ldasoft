This directory contains data files packaged with the installation of ldasoft
which global_fit needs in order to run

Currently they are:

`ucb_frequency_spacing.dat` -> Needs to live in the pwd of the global_fit process. (?xcxc? If omitted I think something is constructed by the program itself??)
`ldc_sangria_vgb_list.dat` -> Needs to be provided on the command line as `--known-sources ldc_sangria_vgb_list.dat`
`search_sources.dat` -> The directory of this file ($INSTALL_DIR/sampledata) needs to be provided on the command line as `--mbh-search-path $INSTALL_DIR/sampledata`


