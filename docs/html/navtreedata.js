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
      [ "Destination directory for install is supplied via command line to install script.  Example:", "md_gbmcmc_README.html#autotoc_md11", null ],
      [ "Example use cases for GBMCMC", "md_gbmcmc_README.html#autotoc_md12", [
        [ "Analyze single high SNR injection", "md_gbmcmc_README.html#autotoc_md13", null ],
        [ "Analyze LDC data", "md_gbmcmc_README.html#autotoc_md14", null ],
        [ "Analyze verification binary with EM priors", "md_gbmcmc_README.html#autotoc_md15", [
          [ "Examples for AM CVn", "md_gbmcmc_README.html#autotoc_md16", null ]
        ] ],
        [ "Other use cases", "md_gbmcmc_README.html#autotoc_md17", [
          [ "Use UCBs as calibration sources", "md_gbmcmc_README.html#autotoc_md18", null ]
        ] ]
      ] ],
      [ "GBMCMC output format", "md_gbmcmc_README.html#autotoc_md19", [
        [ "Parameterization", "md_gbmcmc_README.html#autotoc_md20", null ],
        [ "chains directory", "md_gbmcmc_README.html#autotoc_md21", null ]
      ] ],
      [ "<tt>model_chain.dat.0</tt>: Model state information for the chain sampling the target distribution (<em>temperature</em> T=0). Each row of the file is a chain sample. The columns are", "md_gbmcmc_README.html#autotoc_md22", [
        [ "<tt>chains</tt> directory", "md_gbmcmc_README.html#autotoc_md23", null ]
      ] ],
      [ "<tt>parameter_chain.dat.0</tt>: Full set of chain samples including burn-in.", "md_gbmcmc_README.html#autotoc_md24", null ],
      [ "<tt>dimension_chain.dat.N</tt>:  The post-burnin posterior samples of", "md_gbmcmc_README.html#autotoc_md25", null ],
      [ "<tt>noise_chain.dat.0</tt>: Full set of chain samples for noise model parameters. The columns are", "md_gbmcmc_README.html#autotoc_md26", null ],
      [ "<tt>log_likelihood_chain.dat</tt>: The full log likelihood (<tt>logL + logLnorm</tt>) for each of the <tt>NC</tt> parallel tempered chains, at each chain step. Columns are", "md_gbmcmc_README.html#autotoc_md27", null ],
      [ "<tt>temperature_chain.dat</tt>: Inverse temperature for each of the <tt>NC</tt> parallel tempered chains, at each chain step. Columns are", "md_gbmcmc_README.html#autotoc_md28", [
        [ "data directory", "md_gbmcmc_README.html#autotoc_md29", null ]
      ] ],
      [ "<tt>data_0_0.dat</tt>: Fourier series of input data. For 6-link data columns are", "md_gbmcmc_README.html#autotoc_md30", [
        [ "<tt>data</tt> directory", "md_gbmcmc_README.html#autotoc_md31", null ]
      ] ],
      [ "<tt>waveform_injection_0_0.dat</tt>: Fourier series of the injected signals. For 6-link data columns are", "md_gbmcmc_README.html#autotoc_md32", null ],
      [ "<tt>power_data_0_0.dat</tt>: Power spectrum of input data. For 6-link data columns are", "md_gbmcmc_README.html#autotoc_md33", null ],
      [ "<tt>power_injection_0_0.dat</tt>: Power spectrum of injected signals. For 6-link data columns are", "md_gbmcmc_README.html#autotoc_md34", null ],
      [ "<tt>frequency_proposal.dat</tt>: Smoothed and normalized power spectrum of data used for frequency proposal. Columns are", "md_gbmcmc_README.html#autotoc_md35", null ],
      [ "<tt>power_noise_t0_f0.dat</tt>: Quantiles of the posterior distribution for the reconstructed noise model. For 6-link data the columns are", "md_gbmcmc_README.html#autotoc_md36", null ],
      [ "<tt>power_reconstruction_t0_f0.dat</tt>: Quantiles of the posterior distribution for the reconstructed signal model. For 6-link data the columns are", "md_gbmcmc_README.html#autotoc_md37", null ],
      [ "<tt>power_residual_t0_f0.dat</tt>: Quantiles of the posterior distribution for the data residual. For 6-link data the columns are", "md_gbmcmc_README.html#autotoc_md38", null ],
      [ "<tt>variance_residual_t0_f0.dat</tt>: Variance of the signal model amplitude (and therefore residual) computed by summing the variance of the real and imaginary amplitudes of the joint signal model, i.e. Var(h) = Var(Re h) + Var(Im h). For 6-link data the columns are", "md_gbmcmc_README.html#autotoc_md39", null ],
      [ "<tt>waveform_draw_0.dat</tt>: Fair draw of the data, signal model, and residual printed periodically throughout MCMC analysis to check on health of the fit. For 6-link data the columns are", "md_gbmcmc_README.html#autotoc_md40", [
        [ "main run directory", "md_gbmcmc_README.html#autotoc_md41", null ]
      ] ],
      [ "<tt>evidence.dat</tt>: Posterior for number of templates used to fit the data. Columns are", "md_gbmcmc_README.html#autotoc_md42", null ],
      [ "<tt>avg_log_likelihood.dat</tt>:  Average log likelihoood for each parallel tempered chain. These data would serve as the integrand for a thermodynamic integration step to compute the overall evidence for the model (marginalized over the number of galactic binary signals in the data). Columns are", "md_gbmcmc_README.html#autotoc_md43", null ],
      [ "<tt>example_gb_catalog.sh</tt>: Example <tt>bash</tt> script for post-processing with <tt>gb_catalog</tt>. Takes as argument the size of the model you want to post-process (e.g., maximum from <tt>evidence.dat</tt>).", "md_gbmcmc_README.html#autotoc_md44", null ],
      [ "Post processing GBMCMC", "md_gbmcmc_README.html#autotoc_md45", [
        [ "GB Catalog output format", "md_gbmcmc_README.html#autotoc_md46", null ]
      ] ],
      [ "<tt>catalog_N/LDCFFFFFFFFFF_chain.dat</tt>: The Markov chain samples for parameters consistent with <tt>LDCFFFFFFFFFF</tt>. From these samples the marginalized posteriors for the source are computed. If correlations between sources are of interest the user has to go back to the full chains. Columns are the same as the raw chain files from <tt>gb_mcmc</tt> with <strong>three</strong> additional columns appended", "md_gbmcmc_README.html#autotoc_md47", null ],
      [ "<tt>catalog_N/LDCFFFFFFFFFF_power_reconstruction.dat</tt>: Reconstructed waveform represented by quantiles of the posterior distribution of the signal's power spectrum. For 6-link data the columns are", "md_gbmcmc_README.html#autotoc_md48", null ],
      [ "<tt>catalog_N/LDCFFFFFFFFFF_waveform.dat</tt>: Point estimate of the reconstructed waveform corresponding  to the chain sample containing the <strong>median</strong> of the 1D marginalized distriubtion for the GW frequency. Columns for 6-link data are", "md_gbmcmc_README.html#autotoc_md49", null ],
      [ "<tt>catalog_N/entries.dat</tt>: Summary of sources found in catalog.  Columns are", "md_gbmcmc_README.html#autotoc_md50", null ],
      [ "<tt>catalog_N/history.dat</tt>: File associating current candidate sources to previous catalog input using <tt>--catalog</tt> and <tt>--Tcatalog</tt> arguments for <tt>gb_catalog</tt>.  The columns are", "md_gbmcmc_README.html#autotoc_md51", null ]
    ] ],
    [ "Workflow for processing LDC Radler Data", "md_gbmcmc_README_LDC.html", [
      [ "Download LDC GB data", "md_gbmcmc_README_LDC.html#autotoc_md53", null ],
      [ "FFT LDC data with fft_ldc_data.py", "md_gbmcmc_README_LDC.html#autotoc_md54", null ],
      [ "Access LDC VGB key parameters with python", "md_gbmcmc_README_LDC.html#autotoc_md55", null ],
      [ "Access LDC Galaxy key parameters, with python", "md_gbmcmc_README_LDC.html#autotoc_md56", null ],
      [ "Run gb_mcmc on Fourier transformed LDC data.", "md_gbmcmc_README_LDC.html#autotoc_md57", null ],
      [ "Use GalacticBinaryCatalog.c (gb_catalog) to produce catalog data with the gb_mcmc output", "md_gbmcmc_README_LDC.html#autotoc_md58", null ],
      [ "Build covariance matrix proposal from catalog output", "md_gbmcmc_README_LDC.html#autotoc_md59", null ],
      [ "Run gb_mcmc with draws from a proposal distribution, for example, from cummulative distribution function (CDF proposal updates) and/or from covariance matrices", "md_gbmcmc_README_LDC.html#autotoc_md60", null ]
    ] ],
    [ "FisherGalaxy Manual", "md_gbfisher_README.html", [
      [ "Input file format", "md_gbfisher_README.html#autotoc_md62", [
        [ "Note about sky location conventions", "md_gbfisher_README.html#autotoc_md63", null ]
      ] ],
      [ "Create orbit file with spacecraft ephimeredes", "md_gbfisher_README.html#autotoc_md64", null ],
      [ "Simulate data with detector response to full galaxy and Gaussian instrument noise", "md_gbfisher_README.html#autotoc_md65", null ],
      [ "Take data file and produce fit to confusion noise", "md_gbfisher_README.html#autotoc_md66", [
        [ "Intermediate data products for checking performance", "md_gbfisher_README.html#autotoc_md67", null ]
      ] ],
      [ "Compute residual data and detectable catalog", "md_gbfisher_README.html#autotoc_md68", [
        [ "Intermediate data products for checking performance (overwriting output from <tt>ConfusionFit</tt>)", "md_gbfisher_README.html#autotoc_md69", null ]
      ] ],
      [ "Estimate errors", "md_gbfisher_README.html#autotoc_md70", null ],
      [ "Example for input file <tt>FisherGalaxy_LDC_Radler_galaxy_key.dat</tt>:", "md_gbfisher_README.html#autotoc_md71", null ]
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
"structFilter.html#a3f9d003b4aaece1c94ada07b2544f493",
"structlisa__orbit.html#acc2f400f2aca07364cd34f61b7e441b1"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';