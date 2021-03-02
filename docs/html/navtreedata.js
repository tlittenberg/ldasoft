/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "LDASoft", "index.html", [
    [ "Doxygen code documentation", "index.html#autotoc_md0", null ],
    [ "Acknowledgment", "index.html#autotoc_md1", null ],
    [ "Installation", "index.html#autotoc_md2", null ],
    [ "Issue tracker", "index.html#autotoc_md3", null ],
    [ "Other resources", "index.html#autotoc_md4", null ],
    [ "GBMCMC Manual", "md_gbmcmc_README.html", [
      [ "Table of contents", "md_gbmcmc_README.html#autotoc_md6", null ],
      [ "Introduction", "md_gbmcmc_README.html#autotoc_md7", [
        [ "C Dependencies", "md_gbmcmc_README.html#autotoc_md8", null ],
        [ "Python packages for post-production", "md_gbmcmc_README.html#autotoc_md9", null ],
        [ "Installation", "md_gbmcmc_README.html#autotoc_md10", null ]
      ] ],
      [ "Example use cases for GBMCMC", "md_gbmcmc_README.html#autotoc_md11", [
        [ "Analyze single high SNR injection", "md_gbmcmc_README.html#autotoc_md12", null ],
        [ "Analyze LDC data", "md_gbmcmc_README.html#autotoc_md13", null ],
        [ "Analyze verification binary with EM priors", "md_gbmcmc_README.html#autotoc_md14", [
          [ "Examples for AM CVn", "md_gbmcmc_README.html#autotoc_md15", null ]
        ] ],
        [ "Other use cases", "md_gbmcmc_README.html#autotoc_md16", [
          [ "Use UCBs as calibration sources", "md_gbmcmc_README.html#autotoc_md17", null ]
        ] ]
      ] ],
      [ "GBMCMC output format", "md_gbmcmc_README.html#autotoc_md18", [
        [ "Parameterization", "md_gbmcmc_README.html#autotoc_md19", null ],
        [ "chains directory", "md_gbmcmc_README.html#autotoc_md20", null ],
        [ "data directory", "md_gbmcmc_README.html#autotoc_md21", null ],
        [ "main run directory", "md_gbmcmc_README.html#autotoc_md22", null ]
      ] ],
      [ "Post processing GBMCMC", "md_gbmcmc_README.html#autotoc_md23", [
        [ "GB Catalog output format", "md_gbmcmc_README.html#autotoc_md24", null ]
      ] ]
    ] ],
    [ "Workflow for processing LDC Radler Data", "md_gbmcmc_README_LDC.html", [
      [ "Download LDC GB data", "md_gbmcmc_README_LDC.html#autotoc_md26", null ],
      [ "FFT LDC data with fft_ldc_data.py", "md_gbmcmc_README_LDC.html#autotoc_md27", null ],
      [ "Access LDC VGB key parameters with python", "md_gbmcmc_README_LDC.html#autotoc_md28", null ],
      [ "Access LDC Galaxy key parameters, with python", "md_gbmcmc_README_LDC.html#autotoc_md29", null ],
      [ "Run gb_mcmc on Fourier transformed LDC data.", "md_gbmcmc_README_LDC.html#autotoc_md30", null ],
      [ "Use GalacticBinaryCatalog.c (gb_catalog) to produce catalog data with the gb_mcmc output", "md_gbmcmc_README_LDC.html#autotoc_md31", null ],
      [ "Build covariance matrix proposal from catalog output", "md_gbmcmc_README_LDC.html#autotoc_md32", null ],
      [ "Run gb_mcmc with draws from a proposal distribution, for example, from cummulative distribution function (CDF proposal updates) and/or from covariance matrices", "md_gbmcmc_README_LDC.html#autotoc_md33", null ]
    ] ],
    [ "FisherGalaxy Manual", "md_gbfisher_README.html", [
      [ "Input file format", "md_gbfisher_README.html#autotoc_md35", [
        [ "Note about sky location conventions", "md_gbfisher_README.html#autotoc_md36", null ]
      ] ],
      [ "Create orbit file with spacecraft ephimeredes", "md_gbfisher_README.html#autotoc_md37", null ],
      [ "Simulate data with detector response to full galaxy and Gaussian instrument noise", "md_gbfisher_README.html#autotoc_md38", null ],
      [ "Take data file and produce fit to confusion noise", "md_gbfisher_README.html#autotoc_md39", [
        [ "Intermediate data products for checking performance", "md_gbfisher_README.html#autotoc_md40", null ]
      ] ],
      [ "Compute residual data and detectable catalog", "md_gbfisher_README.html#autotoc_md41", [
        [ "Intermediate data products for checking performance (overwriting output from <tt>ConfusionFit</tt>)", "md_gbfisher_README.html#autotoc_md42", null ]
      ] ],
      [ "Estimate errors", "md_gbfisher_README.html#autotoc_md43", null ],
      [ "Example for input file <tt>FisherGalaxy_LDC_Radler_galaxy_key.dat</tt>:", "md_gbfisher_README.html#autotoc_md44", null ]
    ] ],
    [ "Todo List", "todo.html", null ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Variables", "functions_vars.html", "functions_vars" ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"Constants_8h.html",
"LISA_8h_source.html",
"structFilter.html#ae15e122e2a439ec6d6abdb28622dece4"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';