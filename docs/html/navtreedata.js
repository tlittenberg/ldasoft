/*
@licstart  The following is the entire license notice for the
JavaScript code in this file.

Copyright (C) 1997-2019 by Dimitri van Heesch

This program is free software; you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as published by
the Free Software Foundation

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

@licend  The above is the entire license notice
for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "LDASoft", "index.html", [
    [ "Documentation", "index.html#autotoc_md0", null ],
    [ "Acknowledgment", "index.html#autotoc_md1", null ],
    [ "Installation", "index.html#autotoc_md2", null ],
    [ "Issue tracker", "index.html#autotoc_md3", null ],
    [ "Other resources", "index.html#autotoc_md4", null ],
    [ "Installing GBMCMC", "md_gbmcmc_README.html", [
      [ "Dependencies", "md_gbmcmc_README.html#autotoc_md6", null ],
      [ "Good things to have for post-production", "md_gbmcmc_README.html#autotoc_md7", null ],
      [ "Installation", "md_gbmcmc_README.html#autotoc_md8", null ],
      [ "Example use cases for GBMCMC", "md_gbmcmc_README.html#autotoc_md9", [
        [ "Analyze on-the-fly data & source simulation", "md_gbmcmc_README.html#autotoc_md10", null ],
        [ "Analyze LDC data", "md_gbmcmc_README.html#autotoc_md11", null ],
        [ "Analyze verification binary with EM priors", "md_gbmcmc_README.html#autotoc_md12", [
          [ "Examples for AM CVn", "md_gbmcmc_README.html#autotoc_md13", null ]
        ] ],
        [ "Other use cases", "md_gbmcmc_README.html#autotoc_md14", [
          [ "Use UCBs as calibration sources", "md_gbmcmc_README.html#autotoc_md15", null ]
        ] ]
      ] ],
      [ "GBMCMC output format", "md_gbmcmc_README.html#autotoc_md16", [
        [ "<tt>chains</tt> directory", "md_gbmcmc_README.html#autotoc_md17", [
          [ "<tt>dimension_chain.dat.N</tt>", "md_gbmcmc_README.html#autotoc_md18", null ],
          [ "<tt>model_chain.dat.0</tt>", "md_gbmcmc_README.html#autotoc_md19", null ],
          [ "<tt>parameter_chain.dat.0</tt>", "md_gbmcmc_README.html#autotoc_md20", null ],
          [ "<tt>log_likelihood_chain.dat</tt>", "md_gbmcmc_README.html#autotoc_md21", null ],
          [ "<tt>noise_chain.dat.0</tt>", "md_gbmcmc_README.html#autotoc_md22", null ],
          [ "<tt>temperature_chain.dat</tt>", "md_gbmcmc_README.html#autotoc_md23", null ]
        ] ],
        [ "<tt>data</tt> directory", "md_gbmcmc_README.html#autotoc_md24", [
          [ "<tt>data_0_0.dat</tt>", "md_gbmcmc_README.html#autotoc_md25", null ],
          [ "<tt>frequency_proposal.dat</tt>", "md_gbmcmc_README.html#autotoc_md26", null ],
          [ "<tt>power_data_0_0.dat</tt>", "md_gbmcmc_README.html#autotoc_md27", null ],
          [ "<tt>power_injection_0_0.dat</tt>", "md_gbmcmc_README.html#autotoc_md28", null ],
          [ "<tt>power_noise_t0_f0.dat</tt>", "md_gbmcmc_README.html#autotoc_md29", null ],
          [ "<tt>power_reconstruction_t0_f0.dat</tt>", "md_gbmcmc_README.html#autotoc_md30", null ],
          [ "<tt>power_residual_t0_f0.dat</tt>", "md_gbmcmc_README.html#autotoc_md31", null ],
          [ "<tt>variance_residual_t0_f0.dat</tt>", "md_gbmcmc_README.html#autotoc_md32", null ],
          [ "<tt>waveform_draw_0.dat</tt>", "md_gbmcmc_README.html#autotoc_md33", null ],
          [ "<tt>waveform_injection_0_0.dat</tt>", "md_gbmcmc_README.html#autotoc_md34", null ]
        ] ],
        [ "main run directory", "md_gbmcmc_README.html#autotoc_md35", [
          [ "<tt>evidence.dat</tt>", "md_gbmcmc_README.html#autotoc_md36", null ],
          [ "<tt>gb_mcmc.log</tt>", "md_gbmcmc_README.html#autotoc_md37", null ],
          [ "<tt>run.sh</tt>", "md_gbmcmc_README.html#autotoc_md38", null ],
          [ "<tt>avg_log_likelihood.dat</tt>", "md_gbmcmc_README.html#autotoc_md39", null ]
        ] ]
      ] ],
      [ "Post processing GBMCMC", "md_gbmcmc_README.html#autotoc_md40", null ]
    ] ],
    [ "Workflow for processing LDC Data", "md_gbmcmc_README_LDC.html", [
      [ "Download LDC GB data", "md_gbmcmc_README_LDC.html#autotoc_md42", null ],
      [ "FFT LDC data with fft_ldc_data.py", "md_gbmcmc_README_LDC.html#autotoc_md43", null ],
      [ "Access LDC VGB key parameters with python", "md_gbmcmc_README_LDC.html#autotoc_md44", null ],
      [ "Access LDC Galaxy key parameters, with python", "md_gbmcmc_README_LDC.html#autotoc_md45", null ],
      [ "Run gb_mcmc on Fourier transformed LDC data.", "md_gbmcmc_README_LDC.html#autotoc_md46", null ],
      [ "Use GalacticBinaryCatalog.c (gb_catalog) to produce catalog data with the gb_mcmc output", "md_gbmcmc_README_LDC.html#autotoc_md47", null ],
      [ "Build covariance matrix proposal from catalog output", "md_gbmcmc_README_LDC.html#autotoc_md48", null ],
      [ "Run gb_mcmc with draws from a proposal distribution, for example, from cummulative distribution function (CDF proposal updates) and/or from covariance matrices", "md_gbmcmc_README_LDC.html#autotoc_md49", null ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", null ],
        [ "Variables", "functions_vars.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"BayesLine_8h_source.html",
"structFilter.html#a74848a8afe59171ffd813d41afaac404"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';