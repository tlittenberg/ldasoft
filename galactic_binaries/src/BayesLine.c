/*********************************************************************************/
/*                                                                               */
/*     BayesLine fits the LIGO/Virgo power spectra using a model made up         */
/*     of N Lorentzian lines (described by central frequency f, quality          */
/*     factor Q and amplitude A) and cubic spline with M control points.         */
/*     The number of terms in each model, N, M, are free to vary via RJMCMC      */
/*     updates. The code initializes the models in a non-Markovian fashion,      */
/*     then refines the models with a full Markovian RJMCMC subroutine. This     */
/*     subroutine (LorentzMCMC) can be called by other codes to update the       */
/*     spectral fit to the residuals (data - signal model). Doing this is        */
/*     highly recommended, as it ensures that the spectral model is not          */
/*     eating any of the gravitational wave signal. Since the code is            */
/*     transdimensional (and very free in its movement between dimensions)       */
/*     it will not "over-fit" the spectral model.                                */
/*                                                                               */
/*********************************************************************************/


#include "BayesLine.h"

/*
static void system_pause()
{
  printf("Press Any Key to Continue\n");
  getchar();
}
*/

double sample(double *fprop, double pmax, dataParams *data, gsl_rng *r)
{
  int i;
  double f, alpha;

  int ssize    = data->ncut;
  double flow  = data->flow;
  double fhigh = data->fhigh;
  double Tobs  = data->Tobs;

  int counter = 0;
  int lap=0;
  do
  {
    counter++;
    f = flow+gsl_rng_uniform(r)*(fhigh-flow);
    i = (int)(floor((f-flow)*Tobs));
    if(i < 0) i=0;
    if(i > ssize-1) i=ssize-1;
    alpha = pmax*gsl_rng_uniform(r);
    if(counter>100000)
    {
      printf("BayesLine.c:51\n");
      printf("can not propose valid spline frequency\n");
      printf("sample() counter=%i, lap=%i\n",counter,lap);
      printf("   f/flow=%g/%g\n",f,flow);
      printf("   alpha=%g/%g\n",alpha,pmax);
      printf("   fprop[%i]=%g\n",i,fprop[i]);
      abort();
      counter=0;
      lap++;
    }
  } while(alpha > fprop[i]);

  return(f);
}


double lprop(double f, double *fprop, dataParams *data)
{
  int i;

  int ssize    = data->ncut;
  double flow  = data->flow;
  double Tobs  = data->Tobs;

  i = (int)((f-flow)*Tobs);
  if(i < 0) i=0;
  else if(i > ssize-1) i=ssize-1;

  return(log(fprop[i]));
}

double loglike_fit_spline(double *respow, double *Snf, int ncut)
{
  double lgl, x;
  int i;

  lgl = 0.0;
  for(i=0; i< ncut; i++)
  {
    x = (respow[i]-Snf[i])*(respow[i]-Snf[i])/0.1;
    lgl -= (x);
  }

  return(lgl);
}

//static double logprior(double *invsigma, double *mean, double *Snf, int ilow, int ihigh)
//static double logprior(double *sigma, double *mean, double *Snf, int ilow, int ihigh)
static double logprior(double *lower, double *upper, double *Snf, int ilow, int ihigh)
{
  //double x;
  double lgp;
  int i;

  //leaving off normalizations since they cancel in Hastings ratio
  lgp = 0;
  for(i=ilow; i<ihigh; i++)
  {
    /*
    x = (mean[i] - Snf[i])*invsigma[i];
    lgp -= x*x;
     */

//    if( fabs(mean[i] - Snf[i]) > 3.*sigma[i]) lgp += -1e60;

    if(Snf[i]>upper[i] || Snf[i]<lower[i])
    {
      //printf("PSD out of range:  f=%g, Sn=%g, min=%g, max=%g\n",(double)i/4., Snf[i],lower[i],upper[i]);
      lgp=-1e60;
    }
  }

  //return (0.5*lgp);
  return (lgp);
}

/*
static double logprior_gaussian_model(double *mean, double *sigma, double *Snf, double *spline_f, int spline_n, double *lines_f, int lines_n, dataParams *data)
{
  double x,f;
  double lgp;
  int i,n;

  double flow  = data->flow;
  double Tobs  = data->Tobs;

  lgp = 0.0;

  //Lorentzian model
  for(n=0; n<lines_n; n++)
  {
    f = lines_f[n];
    i = (int)((f-flow)*Tobs);
    x = (mean[i] - Snf[i])/sigma[i];
    lgp -= x*x;


//    printf("    line %i:  f=%g, PSD=%g, , P=%g, sigma=%g, logP=%g\n",n,f,Snf[i],mean[i],fabs(x),lgp);
  }

  //Spline model
  for(n=0; n<spline_n; n++)
  {
    f = spline_f[n];
    i = (int)((f-flow)*Tobs);
    x = (mean[i] - Snf[i])/sigma[i];
    lgp -= x*x;

//    printf("    spline %i:  f=%g, PSD=%g, , P=%g, sigma=%g, logP=%g\n",n,f,Snf[i],mean[i],fabs(x),lgp);
  }

//  system_pause();

  return (0.5*lgp);
}
*/

static double logprior_gaussian(double *mean, double *sigma, double *Snf, int ilow, int ihigh)
{
  double x;
  double lgp;
  int i;

  //leaving off normalizations since they cancel in Hastings ratio
  lgp = 0;
  for(i=ilow; i<ihigh; i++)
  {
     x = (mean[i] - Snf[i])/sigma[i];
     lgp -= x*x;
  }

  return (0.5*lgp);
}

static double loglike(double *respow, double *Snf, int ncut)
{
  double lgl, x;
  int i;

  // leavimng out the log(2Pi) terms since they cancel in Hastings ratio
  lgl = 0.0;
  for(i=0; i< ncut; i++)
  {
    x = respow[i]/Snf[i];
    lgl -= (x+log(Snf[i]));
  }

  return(lgl);
}

void LorentzSplineFit(BayesLinePriors *priors, int zeroLogL, int steps, dataParams *data, lorentzianParams *lines_x, splineParams *spline, double *sfreq, double *spow, gsl_rng *r)
{
  int spass, passes, psmax;
  double logLx, logLy, logH;
  int i, j, k, ki=0, ii, jj, mc;
  int check=0;
  double alpha, heat;
  double SAmaxx, SAminn, lSAmax, lSAmin;
  double lQmin, lQmax, QQ;
  double Aminn, Amaxx, lAmin, lAmax;
  int ac0, ac1, cnt;
  int cc0, cc1, cc2;
  double *Sn, *Sbase;
  double e1, e2, e3, e4;
  double x2, x3, x4;
  double s1, s2, s3, s4;
  int typ;
  double xsm, pmax, y, z;
  double mdn, baseav;
  double logpx=0.0, logpy=0.0, x, beta;
  double *fprop;

  int ncut     = data->ncut;
  int tmax     = data->tmax;
  double flow  = data->flow;
  double fhigh = data->fhigh;

  lorentzianParams *lines_y = malloc(sizeof(lorentzianParams));
  create_lorentzianParams(lines_y,tmax);

  splineParams *spline_y = malloc(sizeof(splineParams));
  create_splineParams(spline_y,spline->n);

  int    nspline   = spline->n;
  double *spoints  = spline->points;
  double *sdatax   = spline->data;
  double *sdatay   = spline_y->data;
  double *spointsy = spline_y->points;

  Sn     = malloc((size_t)(sizeof(double)*(ncut)));
  Sbase  = malloc((size_t)(sizeof(double)*(ncut)));
  fprop  = malloc((size_t)(sizeof(double)*(ncut)));


  mdn    = exp(sdatax[0]);

  SAmaxx = 1.0e2*mdn;
  SAminn = 1.0e-2*mdn;
  Aminn = mdn;
  Amaxx = mdn*1.0e6;
  if(SAmaxx > priors->SAmax) SAmaxx = priors->SAmax;
  if(SAminn < priors->SAmin) SAminn = priors->SAmin;
  if(Amaxx  > priors->LAmax) Amaxx  = priors->LAmax;
  if(Aminn  < priors->LAmin) Aminn  = priors->LAmin;
  lQmin  = log(priors->LQmin);
  lQmax  = log(priors->LQmax);
  lAmin  = log(Aminn);
  lAmax  = log(Amaxx);
  lSAmin = log(SAminn);
  lSAmax = log(SAmaxx);

  // this is the fractional error estimate on the noise level
  s1 = 1.0/sqrt((double)(ncut));

  s2 = 0.01;
  s3 = 0.5;
  s4 = 0.5;

  // set up proposal for frequency jumps
  xsm =0.0;
  for(i=0; i< ncut; i++)
	{
    x = spow[i]/mdn;
    if(x <  100.0) x = 1.0;
    if(x >= 100.0) x = 100.0;
    fprop[i] = x;
    xsm += x;
  }

  pmax = -1.0;
  for(i=0; i< ncut; i++)
	{
    fprop[i] /= xsm;
    if(fprop[i] > pmax) pmax = fprop[i];
  }

  for(i=0; i< lines_x->n; i++)
  {
    lines_x->larray[i] = i;
    lines_y->larray[i] = i;
  }

  alpha = gsl_rng_uniform(r);

  for(i=0; i<lines_x->n; i++)
  {
    lines_x->Q[i] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
    lines_x->f[i] = sample(fprop, pmax, data, r);
    lines_x->A[i] = mdn;
  }

  spectrum_spline(Sn, Sbase, sfreq, data, lines_x, spline);

  baseav = 0.0;
  for(i=0; i< ncut; i++)
	{
    baseav += log(Sbase[i]);
  }
  baseav /= (double)(ncut);


  if(!zeroLogL) logLx = loglike(spow, Sn, ncut);
  else          logLx = 1.0;

  cnt = 0;
  passes = 0;
  spass = steps/4;

  psmax = 2;

  do
  {
    ac0 = 0;
    ac1 = 0;
    cc0 = 1;
    cc1 = 1;
    cc2 = 1;

    for(mc=0; mc < spass; mc++)
    {
      heat = 1.0;
      if(mc < spass/2) heat = pow(10.0,1.0*(double)(spass-2.0*mc)/(double)(spass));

      //copy over current state
      for(i=0; i< lines_x->n; i++)
      {
        lines_y->larray[i] = lines_x->larray[i];
        lines_y->Q[i] = lines_x->Q[i];
        lines_y->f[i] = lines_x->f[i];
        lines_y->A[i] = lines_x->A[i];
      }
      for(i=0; i< nspline; i++)
      {
        sdatay[i]   = sdatax[i];
        spointsy[i] = spoints[i];
      }
      alpha = gsl_rng_uniform(r);

      if(alpha > 0.5)  // try a transdimensional move
      {

        alpha = gsl_rng_uniform(r);
        if(alpha < 0.5)  // try and add a new term
        {
          lines_y->n = lines_x->n+1;
          typ = 2;
        }
        else // try and remove term
        {
          lines_y->n = lines_x->n-1;
          typ = 3;
        }

        check = 0;
        if(lines_y->n < 0 || lines_y->n > tmax) check = 1;


        if(check == 0)
        {

          lines_y->Q[i] = lines_x->Q[i];
          lines_y->f[i] = lines_x->f[i];
          lines_y->A[i] = lines_x->A[i];

          if(lines_y->n < lines_x->n)
          {
            i=(int)(gsl_rng_uniform(r)*(double)(lines_x->n)); // pick a term to try and kill
            k = 0;
            for(j=0; j< lines_x->n; j++)
            {
              if(j != i)
              {
                lines_y->larray[k] = lines_x->larray[j];
                jj = lines_x->larray[j];
                lines_y->A[jj] = lines_x->A[jj];
                lines_y->Q[jj] = lines_x->Q[jj];
                lines_y->f[jj] = lines_x->f[jj];
                k++;
              }
              if(j == i) ki = lines_x->larray[j];  // take note of who's demise is being proposed
            }

            logpx = lprop(lines_x->f[ki], fprop, data);
            logpy = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from alines_y->nwhere

            y = log(lines_x->A[ki])-baseav;
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa_BL*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpx += log((1.0-kappa_BL)/(lAmax-lAmin) + z);
            logpy += -log((lAmax-lAmin));
          }

          if(lines_y->n > lines_x->n)
          {
            // copy over current state
            for(j=0; j< lines_x->n; j++)
            {
              k = lines_x->larray[j];
              lines_y->A[k] = lines_x->A[k];
              lines_y->f[k] = lines_x->f[k];
              lines_y->Q[k] = lines_x->Q[k];
              lines_y->larray[j] = lines_x->larray[j];
            }

            // find a label that isn't in use

            i = -1;
            do
            {
              i++;
              k = 0;
              for(j=0; j< lines_x->n; j++)
              {
                if(i==lines_x->larray[j]) k = 1;
              }
            } while(k == 1);

            lines_y->larray[lines_x->n] = i;
            ii=i;

            // draw new term
            lines_y->Q[i] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
            lines_y->f[i] = sample(fprop, pmax, data, r);
            logpy = lprop(lines_y->f[i], fprop, data);
            logpx = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from alines_y->nwhere

            alpha = gsl_rng_uniform(r);
            if(alpha < kappa_BL)
            {
              y = fabs(gsl_ran_gaussian(r, lAwidth));
              lines_y->A[i] = exp(baseav+y);
            }
            else
            {
              lines_y->A[i] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
              y = log(lines_y->A[i])-baseav;
            }
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa_BL*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpy += log((1.0-kappa_BL)/(lAmax-lAmin) + z);
            logpx += -log((lAmax-lAmin));

            if(lines_y->A[i] > Amaxx) check = 1;
            if(lines_y->A[i] < Aminn) check = 1;
          }
        }
      }
      else  // regular MCMC update
      {

        lines_y->n = lines_x->n;

        e1 = s1;

        for(i=0; i< nspline; i++) sdatay[i] = sdatax[i] + gsl_ran_gaussian(r, e1);

        if(lines_y->n > 0)
        {

          //pick a term to update
          jj=(int)(gsl_rng_uniform(r)*(double)(lines_x->n));
          //find label of who is geting updated
          ii = lines_x->larray[jj];

          // copy over current state
          for(i=0; i< lines_x->n; i++)
          {
            k = lines_x->larray[i];
            lines_y->f[k] = lines_x->f[k];
            lines_y->A[k] = lines_x->A[k];
            lines_y->Q[k] = lines_x->Q[k];
            lines_y->larray[i] = lines_x->larray[i];
          }


          alpha = gsl_rng_uniform(r);
          if(alpha > 0.8)  // 0.8
          {
            lines_y->Q[ii] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
            lines_y->f[ii] = sample(fprop, pmax, data, r);
            logpy = lprop(lines_y->f[ii], fprop, data);
            logpx = lprop(lines_x->f[ii], fprop, data);
            typ = 0;


            alpha = gsl_rng_uniform(r);
            if(alpha < kappa_BL)
            {
              y = fabs(gsl_ran_gaussian(r, lAwidth));
              lines_y->A[ii] = exp(baseav+y);
            }
            else
            {
              lines_y->A[ii] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
              y = log(lines_y->A[ii])-baseav;
            }
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa_BL*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpy += log((1.0-kappa_BL)/(lAmax-lAmin) + z);


            y = log(lines_x->A[ii])-baseav;
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa_BL*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpx += log((1.0-kappa_BL)/(lAmax-lAmin) + z);
          }
          else
          {
            typ = 1;
            alpha = gsl_rng_uniform(r);

            if     (alpha > 0.9) beta = 1.0e+1;
            else if(alpha > 0.6) beta = 1.0e+0;
            else if(alpha > 0.3) beta = 1.0e-1;
            else                 beta = 1.0e-2;

            e2 = beta*s2;
            e3 = beta*s3;
            e4 = beta*s4;
            x2 = gsl_ran_gaussian(r, e2);
            x4 = gsl_ran_gaussian(r, e4);
            x3 = gsl_ran_gaussian(r, e3);

            lines_y->A[ii] = lines_x->A[ii]*exp(x3);
            lines_y->Q[ii] = lines_x->Q[ii]*exp(x4);
            lines_y->f[ii] = lines_x->f[ii]+x2;

            logpx = logpy = 0.0;
          }

          check =0;


          for(i=0; i<nspline; i++) if(sdatay[i] > lSAmax) check = 1;
          for(i=0; i<nspline; i++) if(sdatay[i] < lSAmin) check = 1;
          if(lines_y->A[ii] > Amaxx) check = 1;
          if(lines_y->A[ii] < Aminn) check = 1;
          if(lines_y->f[ii] < flow)  check = 1;
          if(lines_y->f[ii] > fhigh) check = 1;
          if(lines_y->Q[ii] < priors->LQmin)  check = 1;
          if(lines_y->Q[ii] > priors->LQmax)  check = 1;

        }
      }

      if(check == 0)
      {

        if(typ == 0) cc0++;
        if(typ == 1) cc1++;
        if(typ == 2) cc2++;

        if(!zeroLogL)
        {
          spectrum_spline(Sn, Sbase, sfreq, data, lines_y, spline_y);
          logLy = loglike(spow, Sn, ncut);
        }
        else
        {
          logLy = 1.0;
        }

        // prior on line number e(-zeta * n).  (this is a prior, not a proposal, so opposite sign)
        // effectively an SNR cut on lines
        logpy += zeta*(double)(lines_y->n);
        logpx += zeta*(double)(lines_x->n);


        logH = (logLy - logLx)/heat - logpy + logpx;
        alpha = log(gsl_rng_uniform(r));

        if(logH > alpha)
        {
          if(mc%1000 < 10)
          {
            baseav = 0.0;
            for(i=0; i< ncut; i++)
            {
              baseav += log(Sbase[i]);
            }
            baseav /= (double)(ncut);
          }

          if(typ == 0) ac0++;
          if(typ == 1) ac1++;
          logLx = logLy;
          lines_x->n = lines_y->n;
          for(i=0; i< nspline; i++) sdatax[i] = sdatay[i];
          for(i=0; i< tmax; i++)
          {
            lines_x->larray[i] = lines_y->larray[i];
            lines_x->A[i] = lines_y->A[i];
            lines_x->f[i] = lines_y->f[i];
            lines_x->Q[i] = lines_y->Q[i];
          }
        }
      }
    }

    passes++;

    spass *= 2;

    if(lines_x->n > 20) psmax = 3;

  }while((lines_x->n > 5) && (passes < psmax));


  QQ = -1.0;
  x = 0.0;
  // re-map the array to 0..mx ordering
  for(i=0; i< lines_x->n; i++)
  {
    k = lines_x->larray[i];
    lines_y->f[i] = lines_x->f[k];
    lines_y->A[i] = lines_x->A[k];
    lines_y->Q[i] = lines_x->Q[k];
    if(lines_x->Q[k] > QQ)
    {
      QQ = lines_x->Q[k];
      x = lines_x->A[k];
    }
  }

  // return the last value of the chain
  for(i=0; i< lines_x->n; i++)
  {
    lines_x->f[i] = lines_y->f[i];
    lines_x->A[i] = lines_y->A[i];
    lines_x->Q[i] = lines_y->Q[i];
  }

  free(fprop);
  free(Sn);
  free(Sbase);
  destroy_splineParams(spline_y);
  destroy_lorentzianParams(lines_y);
}

void spectrum_spline(double *Sn, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, splineParams *spline)
{
  int i, j, k, n;
  int istart, istop;

  int nspline     = spline->n;
  double *spoints = spline->points;
  double *sdata   = spline->data;

  double *x,*Stemp;

  n = data->ncut;

  x     = malloc((size_t)(sizeof(double)*(n)));
  Stemp = malloc((size_t)(sizeof(double)*(n)));

  //Interpolate {spoints,sdata} ==> {sfreq,x}
  CubicSplineGSL(nspline,spoints,sdata,n,sfreq,x);

  for(i=0; i< n; i++)
	{
    Sbase[i] = exp(x[i]);          // spline base model
    Sn[i] = 0.0;
  }

  for(k=0; k<lines->n; k++)
  {
    j = lines->larray[k];
    for(i=0; i<n; i++) Stemp[i]=Sn[i];
    full_spectrum_add_or_subtract(Sn, Stemp, Sbase, sfreq, data, lines, j, &istart, &istop, 1);
  }

  for(i=0; i< n; i++)
  {
    Sn[i] += Sbase[i];
  }

  free(x);
  free(Stemp);
}

double loglike_pm(double *respow, double *Sn, double *Snx, int ilow, int ihigh)
{
  double lgl, x;
  int i;

  // leavimng out the log(2Pi) terms since they cancel in Hastings ratio
  lgl = 0.0;
  for(i=ilow; i< ihigh; i++)
  {
    x = respow[i]/Sn[i]-respow[i]/Snx[i];
    lgl -= (x+log(Sn[i]/Snx[i]));
  }

  return(lgl);
}



double loglike_single(double *respow, double *Sn, double *Snx, int ilowx, int ihighx, int ilowy, int ihighy)
{
  double lgl, x;
  int i;
  int ilow, ihigh;
  int imid1, imid2;

  ilow = ilowx;
  if(ilowy < ilow) ilow = ilowy;

  if(ilow == ilowx)
  {
    if(ihighx <= ilowy)  // separate regions
    {
      imid1 = ihighx;
      imid2 = ilowy;
      ihigh = ihighy;
    }

    if(ihighx > ilowy) // overlapping regions
    {
      if(ihighx < ihighy)
      {
        imid1 = ihighx;
        imid2 = ihighx;
        ihigh = ihighy;
      }
      else
      {
        imid1 = ilowy;
        imid2 = ilowy;
        ihigh = ihighx;
      }
    }
  }

  if(ilow == ilowy)
  {
    if(ihighy <= ilowx)  // separate regions
    {
      imid1 = ihighy;
      imid2 = ilowx;
      ihigh = ihighx;
    }

    if(ihighy > ilowx) // overlapping regions
    {
      if(ihighy < ihighx)
      {
        imid1 = ihighy;
        imid2 = ihighy;
        ihigh = ihighx;
      }
      else
      {
        imid1 = ilowx;
        imid2 = ilowx;
        ihigh = ihighy;
      }
    }
  }

  // leavimng out the log(2Pi) terms since they cancel in Hastings ratio
  lgl = 0.0;
  for(i=ilow; i< imid1; i++)
  {
    x = respow[i]/Sn[i]-respow[i]/Snx[i];
    lgl -= (x+log(Sn[i]/Snx[i]));
  }
  for(i=imid2; i< ihigh; i++)
  {
    x = respow[i]/Sn[i]-respow[i]/Snx[i];
    lgl -= (x+log(Sn[i]/Snx[i]));
  }

  return(lgl);
}

void full_spectrum_add_or_subtract(double *Snew, double *Sold, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, int ii, int *ilow, int *ihigh, int flag)
{
  int i;
  double dS,f2,f4;
  double deltf;
  double fsq, x, z, deltafmax, spread;
  double amplitude;
  int istart, istop, imid, idelt;

  double A = lines->A[ii];
  double Q = lines->Q[ii];
  double f = lines->f[ii];

  int    ncut = data->ncut;
  double Tobs = data->Tobs;
  double flow = data->flow;

  // copy over current model
  for(i=0; i< ncut; i++) Snew[i] = Sold[i];

  // here we figure out how many frequency bins are needed for the line profile
  imid = (int)((f-flow)*Tobs);
  spread = (1.0e-2*Q);

  if(spread < 20.0) spread = 20.0;  // maximum half-width is f_resonance/20
  deltafmax = f/spread;
  deltf = 4.0*deltafmax;
  idelt = (int)(deltf*Tobs)+1;
  if(A < 10.0*Sbase[imid]) idelt = (int)(20.0*f*Tobs/Q)+1;

  istart = imid-idelt;
  istop = imid+idelt;
  if(istart < 0) istart = 0;
  if(istop > ncut) istop = ncut;

  *ilow  = istart;
  *ihigh = istop;


  // add or remove the old line
  f2=f*f;
  f4=f2*f2;
  amplitude = A*f4/(f2*Q*Q);
  for(i=istart; i<istop; i++)
  {
    fsq = sfreq[i]*sfreq[i];
    x = fabs(f-sfreq[i]);
    z = 1.0;
    if(x > deltafmax) z = exp(-(x-deltafmax)/deltafmax);
    //dS = z*A*f4/(f2*fsq+Q*Q*(fsq-f2)*(fsq-f2));
    if(i==0)
      dS = 0.0;
    else
      dS = z*amplitude/(fsq*(fsq-f2)*(fsq-f2));
    switch(flag)
    {
      case 1: //add new line
        Snew[i] += dS;
        break;
      case -1: //remove line
        Snew[i] -= dS;
        break;
      default:
        break;
    }
  }
}

void full_spectrum_single(double *Sn, double *Snx, double *Sbasex, double *sfreq, dataParams *data, lorentzianParams *line_x, lorentzianParams *line_y, int ii,
                          int *ilowx, int *ihighx, int *ilowy, int *ihighy)
{

  double *Stemp = malloc((size_t)(sizeof(double)*(data->ncut)));

  full_spectrum_add_or_subtract(Stemp, Snx, Sbasex, sfreq, data, line_x, ii, ilowx, ihighx,-1);
  full_spectrum_add_or_subtract(Sn, Stemp,  Sbasex, sfreq, data, line_y, ii, ilowy, ihighy, 1);

  free(Stemp);
}



void full_spectrum_spline(double *Sline, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines)
{
  int i, j, k;
  int istart, istop;

  double *Stemp = malloc((size_t)(sizeof(double)*(data->ncut)));

  for(i=0; i<data->ncut; i++) Sline[i] = 0.0;
  for(k=0; k<lines->n; k++)
  {
    j = lines->larray[k];
    for(i=0; i<data->ncut; i++) Stemp[i]=Sline[i];
    full_spectrum_add_or_subtract(Sline, Stemp, Sbase, sfreq, data, lines, j, &istart, &istop, 1);
  }

  free(Stemp);
}

void SpecFitSpline(BayesLinePriors *priors, int zeroLogL, int steps, double *freq, double *power, splineParams *spline, double *Snf, int ncut, gsl_rng *r)
{
  int i, j, k, ii, ki, mc;
  int nsy, nsx;
  int check;

  int nspline     = spline->n;
  double *sdata   = spline->data;
  double *spoints = spline->points;

  double *sdatax, *sdatay;
  double *spointsx, *spointsy;
  double *Snfx;

  double alpha;
  double logLx, logLy, logH;
  double lSAmax, lSAmin;
  double e1;

  nsx = nspline;

  // maxima and minima for the noise spectal amplitudes
  lSAmin = log(priors->SAmin);
  lSAmax = log(priors->SAmax);

  sdatax   = malloc((size_t)(sizeof(double)*(nspline)));
  sdatay   = malloc((size_t)(sizeof(double)*(nspline)));
  spointsx = malloc((size_t)(sizeof(double)*(nspline)));
  spointsy = malloc((size_t)(sizeof(double)*(nspline)));
  Snfx     = malloc((size_t)(sizeof(double)*(ncut+1)));

  for(i=0; i<nspline; i++)
  {
    sdatax[i]   = sdata[i];
    spointsx[i] = spoints[i];
  }

  // check that line amplitudes are within prior ranges
  for(i=0; i<nsx; i++)
  {
    if(sdatax[i] > lSAmax || sdatax[i] < lSAmin)
    {
      fprintf(stdout,"  Warning: spectrum priors don't seem wide enough\n");
      printf("     segment %i/%i: %lg < %lg < %lg ?\n",i,nsx, exp(lSAmin), exp(sdatax[i]), exp(lSAmax));

      sdatax[i] = lSAmin;
      fprintf(stdout,"     Setting sdata[%i] to %g (minimum)\n",i,exp(sdatax[i]));
    }
  }

  //Interpolate {spointsx,sdatax} ==> {freq,Snf}
  CubicSplineGSL(nspline, spointsx, sdatax, ncut, freq, Snf);

  if(!zeroLogL) logLx = loglike_fit_spline(power, Snf, ncut);
  else          logLx = 1.0;

  for(mc=0; mc < steps; mc++)
	{
    // copy over the current state
    nsy = nsx;
    for(i=0; i<nsx; i++)
    {
      sdatay[i]   = sdatax[i];
      spointsy[i] = spointsx[i];
    }

    alpha = gsl_rng_uniform(r);

    if(alpha > 0.8)  // try a transdimensional move
    {
      alpha = gsl_rng_uniform(r);

      if(alpha > 0.5) nsy = nsx+1; // try and add a new term
      else            nsy = nsx-1; // try and remove term

      if(nsy > 2 && nsy <= nspline)
      {

        if(nsy < nsx)
        {

          ki=1+(int)(gsl_rng_uniform(r)*(double)(nsx-2)); // pick a term to try and kill - cant be first or last term
          k = 0;
          for(j=0; j<nsx; j++)
          {
            if(j != ki)
            {
              sdatay[k]   = sdatax[j];
              spointsy[k] = spointsx[j];
              k++;
            }
          }

        }

        if(nsy > nsx)
        {
          // have to randomly pick a new point that isn't already in use
          do
          {
            ki=1+(int)(gsl_rng_uniform(r)*(double)(nspline-2));  // pick a point to add
            ii = 0;
            for(j=0; j<nsx; j++)
            {
              if(fabs(spointsx[j]-spoints[ki]) < 1.0e-3) ii = 1;  // this means the point is already in use
            }
          } while (ii == 1);

          ii = 0;
          for(j=0; j<nsx; j++)
          {
            if(spointsx[j] < spoints[ki])
            {
              sdatay[j]   = sdatax[j];
              spointsy[j] = spointsx[j];
            }

            if((spointsx[j] > spoints[ki]) && ii == 0)  // found where to slot the new term in
            {
              sdatay[j]   = lSAmin +(lSAmax - lSAmin)*gsl_rng_uniform(r);
              spointsy[j] = spoints[ki];
              ii = 1;
            }

            if((spointsx[j] > spoints[ki]) && ii == 1)
            {
              sdatay[j+1]   = sdatax[j];
              spointsy[j+1] = spointsx[j];
            }

          }
        }
      }
    }//end transdimensional move

    else  // regular MCMC update
    {
      nsy = nsx;

      //pick a term to update
      ii=(int)(gsl_rng_uniform(r)*(double)(nsx));

      // use a variety of jump sizes by using a sum of gaussians of different width
      e1 = 0.0005;
      alpha = gsl_rng_uniform(r);

      if(alpha > 0.8)      e1 = 0.002;
      else if(alpha > 0.7) e1 = 0.005;
      else if(alpha > 0.6) e1 = 0.05;

      // propose new value for the selected term
      sdatay[ii] = sdatax[ii]+gsl_ran_gaussian(r, e1);
    }


    //check that priors on dimension and line amplitude are satisfied
    check = 0;
    if(nsy < 3 || nsy > nspline-1) check = 1;
    for(i=0; i<nsy; i++) if(sdatay[i] > lSAmax || sdatay[i] < lSAmin) check = 1;

    //skip interpolation & Hastings ratio if parameters are out of bounds
    if(check == 0)
    {
      //interpolate {spointsy,sdatay} ==> {freq,Snf}
      CubicSplineGSL(nsy, spointsy, sdatay, ncut, freq, Snf);

      alpha = log(gsl_rng_uniform(r));

      if(!zeroLogL) logLy = loglike_fit_spline(power, Snf, ncut);
      else          logLy = 1.0;

      logH  = logLy - logLx;

      if(logH > alpha)
      {
        logLx = logLy;
        nsx = nsy;
        for(i=0; i< ncut; i++) Snfx[i] = Snf[i];
        for(i=0; i<nsx; i++)
        {
          sdatax[i]   = sdatay[i];
          spointsx[i] = spointsy[i];
        }
      }

    }//end prior check

  }

  // return the most recent accepted estimate for the spectrum
  for(i=0; i< ncut; i++) Snf[i] = Snfx[i];

  // pass back fit

  //interpolate {spointsx,sdatax} ==> {spoints,sdata}
  CubicSplineGSL(nsx, spointsx, sdatax, nspline, spoints, sdata);

  free(Snfx);
  free(sdatax);
  free(sdatay);
  free(spointsx);
  free(spointsy);

}

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint)
{
  int n;
  double tol=1.0e-6;


  /* set up GSL cubic spline */
  gsl_spline       *cspline = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_interp_accel *acc    = gsl_interp_accel_alloc();

  /* get derivatives */
  gsl_spline_init(cspline,x,y,N);

  /* interpolate */
  for(n=0; n<Nint; n++)
  {
    /*
     GSL cubic spline throws errors if
     interpolated points are at end of
     spline control points
     */
    if     (fabs(xint[n]-x[0])<tol)
      yint[n] = y[0];

    else if(fabs(xint[n]-x[N-1])<tol)
      yint[n] = y[N-1];

    else
      yint[n]=gsl_spline_eval(cspline,xint[n],acc);
  }

  gsl_spline_free (cspline);
  gsl_interp_accel_free (acc);

}

void create_dataParams(dataParams *data, double *f, int n)
{

  // length of segment in seconds, this should be read in from the frame file
  data->Tobs = rint(1.0/(f[1]-f[0]));

  // frequency resolution
  data->df = 1.0/data->Tobs;

//  // sample cadence in Hz, this should be read in from the frame file
//  data->cadence = pow(2.0,rint(log((double)(n))/log(2.0))+1.0)/data->Tobs;
//
//  // Nyquist frequency
//  data->fny = 2.0/data->cadence;

  // size of segments in Hz
  // If frequency snippets are too large need longer initial runs to get convergence
  data->fstep = f[256]-f[0];//.0;//9.0;

  // This sets the maximum number of Lorentzian lines per segment.
  // For segements shorter than 16 Hz this always seems to be enough
  data->tmax = 40;

  // approximate number of segments
  data->sgmts = (int)((f[n-1]-f[0])/data->fstep);

  // Maximum number of lines for full data set
  data->tfull = data->tmax*data->sgmts+1;

  //minimum frequency
  data->fmin = f[0];

  //maximum frequency
  data->fmax = f[n-1];

  //minimum Fourier bin
  data->nmin = (int)(f[0]*data->Tobs);

  // the stencil separation in Hz for the spline model. Going below 2 Hz is dangerous - will fight with line model
  data->fgrid = f[32]-f[0];//20.0;//4.0;
}

void create_lorentzianParams(lorentzianParams *lines, int size)
{
  lines->n    = 0;
  lines->size = size;

  lines->larray = malloc((size_t)(sizeof(int)*size));

  lines->f = malloc((size_t)(sizeof(double)*size));
  lines->Q = malloc((size_t)(sizeof(double)*size));
  lines->A = malloc((size_t)(sizeof(double)*size));
}

void copy_lorentzianParams(lorentzianParams *origin, lorentzianParams *copy)
{
  copy->n    = origin->n;
  copy->size = origin->size;

  int n;
  for(n=0; n<origin->size; n++)
  {
    copy->larray[n] = origin->larray[n];

    copy->f[n] = origin->f[n];
    copy->Q[n] = origin->Q[n];
    copy->A[n] = origin->A[n];
  }
}

void destroy_lorentzianParams(lorentzianParams *lines)
{
  free(lines->larray);
  free(lines->f);
  free(lines->Q);
  free(lines->A);
  free(lines);
}

void create_splineParams(splineParams *spline, int size)
{
  spline->n = size;

  spline->data   = malloc((size_t)(sizeof(double)*size));
  spline->points = malloc((size_t)(sizeof(double)*size));
}

void copy_splineParams(splineParams *origin, splineParams *copy)
{
  int n;
  copy->n = origin->n;

  for(n=0; n<origin->n; n++)
  {
    copy->data[n]   = origin->data[n];
    copy->points[n] = origin->points[n];
  }
}

void destroy_splineParams(splineParams *spline)
{
  free(spline->data);
  free(spline->points);
  free(spline);
}


void copy_bayesline_params(struct BayesLineParams *origin, struct BayesLineParams *copy)
{
  copy->data->tmax = origin->data->tmax;

  int i;
  for(i=0; i<origin->data->ncut; i++)
  {
    copy->spow[i]  = origin->spow[i];
    copy->sfreq[i] = origin->sfreq[i];
  }

  int imax = (int)(floor(origin->data->fmax  * origin->data->Tobs));
  int imin = (int)(floor(origin->data->fmin  * origin->data->Tobs));
  for(i=0; i<imax-imin; i++)
  {
    copy->Snf[i]    = origin->Snf[i];
    copy->spow[i]   = origin->spow[i];
    copy->sfreq[i]  = origin->sfreq[i];
    copy->Sbase[i]  = origin->Sbase[i];
    copy->Sline[i]  = origin->Sline[i];
  }

  copy_splineParams(origin->spline, copy->spline);
  copy_splineParams(origin->spline_x,copy->spline_x);
  copy_lorentzianParams(origin->lines_x, copy->lines_x);
  copy_lorentzianParams(origin->lines_full,copy->lines_full);

}

void BayesLineNonMarkovianFit(struct BayesLineParams *bayesline, int *nj)
{

  /******************************************************************************/
  /*                                                                            */
  /*  Rapid non-Markovian fit over small bandwidth blocks of data               */
  /*                                                                            */
  /******************************************************************************/
  int i,imax,imin;
  int j=0, jj=0;
  int kk=0;

  double mdn,sm;
  double frs,deltafmax,spread,fsq;
  double x,z;

  double *wndw=NULL;
  double *sfreq=NULL;
  double *Sn=NULL;
  double *spow=NULL;
  double *Sbase=NULL;
  double *y=NULL;

  /* Make local pointers to BayesLineParams structure members */
  dataParams *data             = bayesline->data;
  BayesLinePriors *priors      = bayesline->priors;
  lorentzianParams *lines_full = bayesline->lines_full;
  double *freq                 = bayesline->freq;
  double *power                = bayesline->power;
  double *Sna                  = bayesline->Sna;
  double *fa                   = bayesline->fa;
  gsl_rng *r                   = bayesline->r;



  lorentzianParams *lines_x = malloc( (size_t)sizeof(lorentzianParams) );

  splineParams *spline = malloc( (size_t)sizeof(splineParams) );

  // storage for line model for a segment. These get updated and passed back from the fitting routine
  create_lorentzianParams(lines_x,data->tmax);

  // start with a single line. This number gets updated and passed back
  lines_x->n = 1;


  //start with ~4 Hz resolution for the spline
  int nspline = (int)(data->fstep/data->fgrid);

  create_splineParams(spline,nspline);

  double *xint = malloc((size_t)(sizeof(double)*(nspline)));
  double *yint = malloc((size_t)(sizeof(double)*(nspline)));

  // loop over the frequency segments
  kk = 0;
  jj = 0;

  data->flow = data->fmin;
  do
  {
    data->fhigh = data->flow + data->fstep;
    if(data->fhigh>data->fmax) data->fhigh=data->fmax;
    imax = (int)(floor(data->fhigh * data->Tobs));
    imin = (int)(floor(data->flow  * data->Tobs));
    data->ncut = imax-imin;

    lines_x->n = 1;

    y     = malloc((size_t)(sizeof(double)*(data->ncut)));
    Sn    = malloc((size_t)(sizeof(double)*(data->ncut)));
    Sbase = malloc((size_t)(sizeof(double)*(data->ncut)));
    spow  = malloc((size_t)(sizeof(double)*(data->ncut)));
    sfreq = malloc((size_t)(sizeof(double)*(data->ncut)));
    wndw  = malloc((size_t)(sizeof(double)*(data->ncut)));

    sm = data->ncut/2;

    for(i=0; i<data->ncut; i++)
    {
      j = i+imin-data->nmin;

      spow[i]  = power[j];
      wndw[i]  = power[j];
      sfreq[i] = freq[j];
    }
    // The 1.47 factor converts from median to mean
    // We don't care about the even-odd issue since we have so many terms
    gsl_sort(wndw,1,data->ncut);
    mdn = log( 1.47*( gsl_stats_quantile_from_sorted_data(wndw,1,data->ncut,0.5) ) );

    for(j=0; j<nspline; j++)
    {
      spline->points[j] = data->flow+(data->fhigh-data->flow)*(double)(j)/(double)(nspline-1);
      spline->data[j]   = mdn;
    }

    LorentzSplineFit(priors, bayesline->constantLogLFlag, 40000, data, lines_x, spline, sfreq, spow, r);

    //Interpoloate {spoints,sdata} ==> {sfreq,y}
    CubicSplineGSL(spline->n,spline->points,spline->data,data->ncut,sfreq,y);

    //FILE *temp=fopen("temp2.dat","w");
    for(i=0; i< data->ncut; i++)
    {
      Sbase[i] = exp(y[i]);
      Sn[i]    = Sbase[i];
      fsq      = sfreq[i]*sfreq[i];

      for(j=0; j< lines_x->n; j++)
      {
        spread = (1.0e-2*lines_x->Q[j]);
        if(spread < 20.0) spread = 20.0;  // maximum half-width is f_resonance/100
        deltafmax = lines_x->f[j]/spread;
        frs = lines_x->f[j]*lines_x->f[j];
        x = fabs(lines_x->f[j]-sfreq[i]);
        z = 1.0;
        if(x > deltafmax) z = exp(-(x-deltafmax)/deltafmax);
        Sn[i] += z*lines_x->A[j]*frs*frs/(frs*fsq+lines_x->Q[j]*lines_x->Q[j]*(fsq-frs)*(fsq-frs));
      }

      //at frequencies below fmin output SAmax (killing lower frequencies in any future inner products)
      //fprintf(temp,"%lg %lg %lg %lg\n",sfreq[i],spow[i],Sbase[i],Sn[i]);

    }
    //fclose(temp);
    //system_pause();


    free(y);
    free(Sbase);
    free(Sn);
    free(spow);
    free(sfreq);
    free(wndw);

    for(j=0; j<nspline; j++) xint[j] = data->flow+(data->fhigh-data->flow)*(double)(j)/(double)(nspline-1);
    CubicSplineGSL(nspline,spline->points,spline->data,nspline,xint,yint);

    for(j=0; j<nspline; j++)
    {
      fa[jj]  = xint[j];
      Sna[jj] = yint[j];
      jj++;
    }

    for(j=0; j<lines_x->n; j++)
    {
      lines_full->f[kk] = lines_x->f[j];
      lines_full->Q[kk] = lines_x->Q[j];
      lines_full->A[kk] = lines_x->A[j];
      kk++;
    }

    data->flow += data->fstep;
  }while(data->fhigh < data->fmax);

  destroy_lorentzianParams(lines_x);
  destroy_splineParams(spline);
  free(xint);
  free(yint);

  // current number of terms in the Lorentzian model
  lines_full->n = kk-1;
  *nj = jj;

}

void BayesLineLorentzSplineMCMC(struct BayesLineParams *bayesline, double heat, int steps, int focus, int priorFlag, double *dan)
{
  int nsy, nsx;
  double logLx, logLy=0.0, logH;
  int ilowx, ihighx, ilowy, ihighy;
  int i, j, k=0, ki=0, ii=0, jj=0, mc;
  int check=0;
  double alpha;
  double lSAmax, lSAmin;
  double lQmin, lQmax;
  double lAmin, lAmax;
  int ac0, ac1, ac2;
  int cc0, cc1, cc2;
  double *Sn, *Sbase, *Sbasex, *Sline, *Snx;
  double *xint;
  double e1, e2, e3, e4;
  double x2, x3, x4;
  double s1, s2, s3, s4;
  int typ=0;
  double xsm, pmax, fcl, fch, dff;
  double baseav;
  double logpx=0.0, logpy=0.0, x, y, z, beta;
  double logPsy,logPsx;
  double Ac;
  double *fprop;
  double *sdatay;
  double *spointsy;
  int *foc;

  /* Make local pointers to BayesLineParams structure members */
  dataParams *data           = bayesline->data;
  lorentzianParams *lines_x  = bayesline->lines_full;
  splineParams *spline       = bayesline->spline;
  splineParams *spline_x     = bayesline->spline_x;
  BayesLinePriors *priors    = bayesline->priors;
  double *Snf                = bayesline->Snf;
  double *freq               = bayesline->sfreq;
  double *power              = bayesline->spow;
  gsl_rng *r                 = bayesline->r;

  int ncut = data->ncut;
  int tmax = data->tmax;

  double flow  = data->flow;
  double fhigh = data->fhigh;

  int nspline   = spline->n;
  int *nsplinex = &spline_x->n;

  double *sdatax = spline_x->data;

  double *spoints  = spline->points;
  double *spointsx = spline_x->points;

  Snx    = malloc((size_t)(sizeof(double)*(ncut)));
  Sn     = malloc((size_t)(sizeof(double)*(ncut)));
  Sbasex = malloc((size_t)(sizeof(double)*(ncut)));
  Sbase  = malloc((size_t)(sizeof(double)*(ncut)));
  Sline  = malloc((size_t)(sizeof(double)*(ncut)));

  sdatay   = malloc((size_t)(sizeof(double)*(nspline)));
  spointsy = malloc((size_t)(sizeof(double)*(nspline)));

  foc     = malloc((size_t)(sizeof(int)*(tmax)));

  // This keeps track of whos who in the Lorentzian model
  // Necessary complication when using delta likelihood
  // calculations that only update a single line
  lorentzianParams *lines_y = malloc(sizeof(lorentzianParams));
  create_lorentzianParams(lines_y,tmax);

  fprop = malloc((size_t)(sizeof(double)*(ncut)));
  xint  = malloc((size_t)(sizeof(double)*(ncut)));

  // maxima and minima for the noise spectal slopes and amplitudes
  // uniform priors in slope and log amplitude
  lQmin = log(priors->LQmin);
  lQmax = log(priors->LQmax);
  lAmin = log(priors->LAmin);
  lAmax = log(priors->LAmax);
  lSAmin = log(priors->SAmin);
  lSAmax = log(priors->SAmax);

  nsx = *nsplinex;

  dff = 0.01;  // half-width of frequency focus region (used if focus == 1)

  // this is the fractional error estimate on the noise level
  s1 = 1.0/sqrt((double)(ncut));

  s2 = 0.01;
  s3 = 0.5;
  s4 = 0.5;

  for(i=0; i<lines_x->n; i++)
  {
    lines_x->larray[i] = i;
    lines_y->larray[i] = i;
  }

  baseav = 0.0;

  //Interpolate {spointsx,sdatax} ==> {freq,xint}
  CubicSplineGSL(nsx,spointsx,sdatax,ncut,freq,xint);

  for(i=0; i<ncut; i++)
  {
    Sbase[i] = exp(xint[i]);
    Sbasex[i] = Sbase[i];
    baseav += xint[i];
  }

  baseav /= (double)(ncut);

  full_spectrum_spline(Sline, Sbase, freq, data, lines_x);
  for(i=0; i< ncut; i++) Sn[i] = Sbase[i]+Sline[i];

  for(i=0; i<ncut; i++) Snx[i] = Sn[i];

  if(!bayesline->constantLogLFlag)
    logLx = loglike(power, Sn, ncut);
  else
    logLx = 1.0;

  if(priorFlag==1)
  {
    logPsx = logprior(priors->lower, priors->upper, Snx, 0, ncut);
  }
  if(priorFlag==2)
  {
    logPsx = logprior_gaussian(priors->mean, priors->sigma, Snx, 0, ncut);
  }



  // set up proposal for frequency jumps
  xsm =0.0;
  pmax = -1.0;
  for(i=0; i< ncut; i++)
	{
    x = power[i]/Sn[i];
    if(x > pmax)
    {
      pmax = x;
      k = i;
    }
    if(x < 10.0) x = 1.0;
    if(x >= 10.0) x = 100.0;
    fprop[i] = x;
    xsm += x;
  }

  // define the focus region (only used if focus flag = 1)
  fcl = freq[k]-dff;
  fch = freq[k]+dff;
  Ac = power[k];
  if(fcl < freq[0]) fcl = freq[0];
  if(fch > freq[ncut-1]) fch = freq[ncut-1];

  ac0 = 0;
  ac1 = 0;
  ac2 = 0;
  cc0 = 1;
  cc1 = 1;
  cc2 = 1;

  for(mc=0; mc < steps; mc++)
  {
    typ=-1;

    //copy over current state
    lines_y->n = lines_x->n;
    nsy = nsx;
    for(i=0; i< tmax; i++)
    {
      lines_y->larray[i] = lines_x->larray[i];
      lines_y->Q[i] = lines_x->Q[i];
      lines_y->f[i] = lines_x->f[i];
      lines_y->A[i] = lines_x->A[i];
    }
    for(i=0; i<nspline; i++)
    {
      sdatay[i] = sdatax[i];
      spointsy[i] = spointsx[i];
    }

    beta = gsl_rng_uniform(r);

    if(beta > 0.5)  // update the smooth part of the spectrum
    {

      alpha = gsl_rng_uniform(r);

      logpx = logpy = 0.0;

      if(alpha > 0.8)  // try a transdimensional move
      {
        alpha = gsl_rng_uniform(r);
        if(alpha > 0.5)// || nsx<3)  // try and add a new term
        {
          nsy = nsx+1;
          typ = 5;
        }
        else // try and remove term
        {
          nsy = nsx-1;
          typ = 6;
        }

        if(nsy > 0 && nsy < nspline)
        {

          if(nsy < nsx)
          {

            ki=1+(int)(gsl_rng_uniform(r)*(double)(nsx-2)); // pick a term to try and kill - cant be first or last term
            k = 0;
            for(j=0; j<nsx; j++)
            {
              if(j != ki)
              {
                sdatay[k] = sdatax[j];
                spointsy[k] = spointsx[j];
                k++;
              }
            }

          }

          if(nsy > nsx)
          {

            // have to randomly pick a new point that isn't already in use
            do
            {
              ki=1+(int)(gsl_rng_uniform(r)*(double)(nspline-2));  // pick a point to add
              ii = 0;
              for(j=0; j<nsx; j++)
              {
                if(fabs(spointsx[j]-spoints[ki]) < 1.0e-3) ii = 1;  // this means the point is already in use
              }
            } while (ii == 1);
            ii = 0;
            for(j=0; j<nsx; j++)
            {

              if(spointsx[j] < spoints[ki])
              {
                sdatay[j] = sdatax[j];
                spointsy[j] = spointsx[j];
              }

              if((spointsx[j] > spoints[ki]) && ii == 0)  // found where to slot the new term in
              {
                sdatay[j] = lSAmin +(lSAmax - lSAmin)*gsl_rng_uniform(r);
                spointsy[j] = spoints[ki];
                ii = 1;
              }

              if((spointsx[j] > spoints[ki]) && ii == 1)
              {
                sdatay[j+1] = sdatax[j];
                spointsy[j+1] = spointsx[j];
              }

            }
          }
        }
      }
      else  // regular MCMC update
      {
        typ = 4;

        nsy = nsx;

        //pick a term to update
        ii=(int)(gsl_rng_uniform(r)*(double)(nsx));

        // use a variety of jump sizes by using a sum of gaussians of different width
        e1 = 0.0005;
        alpha = gsl_rng_uniform(r);
        if(alpha > 0.8)
        {
          e1 = 0.002;
        }
        else if(alpha > 0.6)
        {
          e1 = 0.005;
        }
        else if(alpha > 0.4)
        {
          e1 = 0.05;
        }

        // propose new value for the selected term
        sdatay[ii] = sdatax[ii]+gsl_ran_gaussian(r, e1);

      }

      check = 0;

      if(nsy < 5 || nsy > nspline-1) check = 1;

      for(i=0; i<nsy; i++)
      {
        if(sdatay[i] > lSAmax) check = 1;
        if(sdatay[i] < lSAmin) check = 1;
      }
    }
    else    // update the line model
    {

      alpha = gsl_rng_uniform(r);

      if(alpha > 0.8)  // try a transdimensional move
      {

        alpha = gsl_rng_uniform(r);
        if(alpha < 0.5)  // try and add a new term
        {
          lines_y->n = lines_x->n+1;
          typ = 2;
        }
        else // try and remove term
        {
          lines_y->n = lines_x->n-1;
          typ = 3;
        }

        check = 0;
        if(lines_y->n < 0 || lines_y->n > tmax) check = 1;


        if(check == 0)
        {
          if(lines_y->n < lines_x->n)
          {
            i=(int)(gsl_rng_uniform(r)*(double)(lines_x->n)); // pick a term to try and kill
            k = 0;
            for(j=0; j< lines_x->n; j++)
            {
              if(j != i)
              {
                lines_y->larray[k] = lines_x->larray[j];
                jj = lines_x->larray[j];
                lines_y->A[jj] = lines_x->A[jj];
                lines_y->Q[jj] = lines_x->Q[jj];
                lines_y->f[jj] = lines_x->f[jj];
                k++;
              }
              if(j == i) ki = lines_x->larray[j];  // take note of who's demise is being proposed
            }

            logpx = lprop(lines_x->f[ki], fprop, data);
            logpy = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from anywhere

            y = log(lines_x->A[ki])-baseav;
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa_BL*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpx += log((1.0-kappa_BL)/(lAmax-lAmin) + z);
            logpy += -log((lAmax-lAmin));
          }

          if(lines_y->n > lines_x->n)
          {
            // copy over current state
            for(j=0; j< lines_x->n; j++)
            {
              k = lines_x->larray[j];
              lines_y->A[k] = lines_x->A[k];
              lines_y->f[k] = lines_x->f[k];
              lines_y->Q[k] = lines_x->Q[k];
              lines_y->larray[j] = lines_x->larray[j];
            }

            // find a label that isn't in use

            i = -1;
            do
            {
              i++;
              k = 0;
              for(j=0; j< lines_x->n; j++)
              {
                if(i==lines_x->larray[j]) k = 1;
              }
            } while(k == 1);

            lines_y->larray[lines_x->n] = i;
            ii=i;

            // draw new term
            lines_y->A[i] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
            lines_y->Q[i] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
            lines_y->f[i] = sample(fprop, pmax, data, r);
            logpy = lprop(lines_y->f[i], fprop, data);
            logpx = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from anywhere
            alpha = gsl_rng_uniform(r);
            if(alpha < kappa_BL)
            {
              y = fabs(gsl_ran_gaussian(r, lAwidth));
              lines_y->A[i] = exp(baseav+y);
            }
            else
            {
              lines_y->A[i] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
              y = log(lines_y->A[i])-baseav;
            }
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa_BL*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpy += log((1.0-kappa_BL)/(lAmax-lAmin) + z);
            logpx += -log((lAmax-lAmin));


            if(focus==1) // using focused region (not Markovian - not paying penalty for proposing in such a small region)
            {
              lines_y->f[i] = fcl+(fch-fcl)*gsl_rng_uniform(r);
              lines_y->Q[i] = priors->LQmax/10.0;
              lines_y->A[i] = Ac;

              logpy = 0.0;
              logpx = 0.0;
            }
            if(lines_y->A[i] > priors->LAmax) check = 1;
            if(lines_y->A[i] < priors->LAmin) check = 1;
          }
        }
      }
      else  // regular MCMC update
      {
        lines_y->n = lines_x->n;

        if(lines_y->n > 0)
        {
          typ=1;

          //pick a term to update
          jj=(int)(gsl_rng_uniform(r)*(double)(lines_x->n));
          //find label of who is geting updated
          ii = lines_x->larray[jj];

          // copy over current state
          for(i=0; i< lines_x->n; i++)
          {
            k = lines_x->larray[i];
            lines_y->f[k] = lines_x->f[k];
            lines_y->A[k] = lines_x->A[k];
            lines_y->Q[k] = lines_x->Q[k];
            lines_y->larray[i] = lines_x->larray[i];
          }

          if(focus == 1)
          {
            // find if any lines are in the focus region
            j = 0;
            for(i=0; i< lines_x->n; i++)
            {
              k = lines_x->larray[i];
              if(lines_x->f[k] > fcl && lines_x->f[k] < fch)
              {
                foc[j] = k;
                j++;
              }
            }

            x = 0.0;
            if(j > 0)  // some lines are currently in the focus region
            {
              x = 0.8;
              jj=(int)(gsl_rng_uniform(r)*(double)(j));
              //find label of who is getting updated in the focus region
              ii = foc[jj];
            }

          }
          else
          {
            x = 0.8;
          }

          alpha = gsl_rng_uniform(r);
          if(alpha > x)
          {
            // here we try and move an exisiting line to a totally new location
            if(focus != 1)
            {
              lines_y->f[ii] = sample(fprop, pmax, data, r);
              logpy = lprop(lines_y->f[ii], fprop, data);
              logpx = lprop(lines_x->f[ii], fprop, data);

              lines_y->Q[ii] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
              alpha = gsl_rng_uniform(r);
              if(alpha < kappa_BL)
              {
                y = fabs(gsl_ran_gaussian(r, lAwidth));
                lines_y->A[ii] = exp(baseav+y);
              }
              else
              {
                lines_y->A[ii] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
                y = log(lines_y->A[ii])-baseav;
              }
              if(y < 0.0)
              {
                z = 0.0;
              }
              else
              {
                z = 2.0*kappa_BL*gsl_ran_gaussian_pdf(y, lAwidth);
              }
              logpy += log((1.0-kappa_BL)/(lAmax-lAmin) + z);
              logpx += -log((lAmax-lAmin));


            }
            else  // using focused region (not Markovian - not paying penalty for proposing in such a small region)
            {
              lines_y->f[ii] = fcl+(fch-fcl)*gsl_rng_uniform(r);
              lines_y->Q[ii] = priors->LQmax/10.0;
              lines_y->A[ii] = Ac;

              logpy = 0.0;
              logpx = 0.0;
            }
            typ = 0;
          }
          else
          {
            typ = 1;
            alpha = gsl_rng_uniform(r);

            if     (alpha > 0.9) beta = 1.0e+1;
            else if(alpha > 0.6) beta = 1.0e+0;
            else if(alpha > 0.3) beta = 1.0e-1;
            else                 beta = 1.0e-2;

            e2 = beta*s2;
            e3 = beta*s3;
            e4 = beta*s4;
            x2 = gsl_ran_gaussian(r, e2);
            x4 = gsl_ran_gaussian(r, e4);
            x3 = gsl_ran_gaussian(r, e3);

            lines_y->A[ii] = lines_x->A[ii]*exp(x3);
            lines_y->Q[ii] = lines_x->Q[ii]*exp(x4);
            lines_y->f[ii] = lines_x->f[ii]+x2;

            logpx = logpy = 0.0;
          }

          check =0;

          if(lines_y->A[ii] > priors->LAmax)  check = 1;
          if(lines_y->A[ii] < priors->LAmin)  check = 1;
          if(lines_y->f[ii] < flow)  check = 1;
          if(lines_y->f[ii] > fhigh) check = 1;
          if(lines_y->Q[ii] < priors->LQmin)  check = 1;
          if(lines_y->Q[ii] > priors->LQmax)  check = 1;

        }
        else check=1;
      }


    }

    //If line parameters satisfy priors, continue with MCMC update
    if(check == 0)
    {

      if(typ == 0) cc0++;
      if(typ == 1) cc1++;
      if(typ == 4) cc2++;

      if(!bayesline->constantLogLFlag)
      {
        if(typ > 3)  // need to do a full update (slower)
        {
          if(nsy>2)
          {

            //Interpolate {spointsy,sdatay} ==> {freq,xint}
            CubicSplineGSL(nsy,spointsy,sdatay,ncut,freq,xint);

            for(i=0; i<ncut; i++)  Sbase[i] = exp(xint[i]);

            full_spectrum_spline(Sline, Sbase, freq, data, lines_y);
            for(i=0; i< ncut; i++) Sn[i] = Sbase[i]+Sline[i];
            logLy = loglike(power, Sn, ncut);
          }
          else logLy = -1e60;
        }

        if(typ == 1 || typ == 0)  // fixed dimension MCMC of line ii
        {
          full_spectrum_single(Sn, Snx, Sbasex, freq, data, lines_x, lines_y, ii, &ilowx, &ihighx, &ilowy, &ihighy);
          logLy = logLx + loglike_single(power, Sn, Snx, ilowx, ihighx, ilowy, ihighy);
        }

        if(typ == 2)  // add new line with label ii
        {
          full_spectrum_add_or_subtract(Sn, Snx, Sbasex, freq, data, lines_y, ii, &ilowy, &ihighy,  1);
          logLy = logLx + loglike_pm(power, Sn, Snx, ilowy, ihighy);
        }

        if(typ == 3)  // remove line with label ki
        {
          full_spectrum_add_or_subtract(Sn, Snx, Sbasex, freq, data, lines_x, ki, &ilowx, &ihighx, -1);
          logLy = logLx + loglike_pm(power, Sn, Snx, ilowx, ihighx);
        }

      }
      else
      {
        logLy = 1.0;
      }

      // prior on line number e(-zeta * n).  (this is a prior, not a proposal, so opposite sign)
      // effectively an SNR cut on lines
      logpy += zeta*(double)(lines_y->n);
      logpx += zeta*(double)(lines_x->n);

      //logPsy = logprior(priors->invsigma, priors->mean, Sn, 0, ncut);
      //if(priorFlag)logPsy = logprior(priors->sigma, priors->mean, Sn, 0, ncut);
      if(priorFlag==1)
      {
        logPsy = logprior(priors->lower, priors->upper, Sn, 0, ncut);
        //logPsy = logprior_gaussian_model(priors->mean, priors->sigma, Sn, spointsy, nsy, lines_y->f, lines_y->n, data);
      }
      if(priorFlag==2)
      {
        logPsy = logprior_gaussian(priors->mean, priors->sigma, Sn, 0, ncut);
      }



      logH  = (logLy - logLx)*heat - logpy + logpx;
      if(priorFlag!=0) logH += logPsy - logPsx;

      alpha = log(gsl_rng_uniform(r));

      if(logH > alpha)
      {
        if(typ == 0) ac0++;
        if(typ == 1) ac1++;
        if(typ == 4) ac2++;
        logLx = logLy;
        //if(priorFlag!=0)
        logPsx = logPsy;
        lines_x->n = lines_y->n;
        nsx = nsy;
        for(i=0; i< ncut; i++)
        {
          Snx[i] = Sn[i];
          if(typ > 3) Sbasex[i] = Sbase[i];
        }
        for(i=0; i< tmax; i++)
        {
          lines_x->larray[i] = lines_y->larray[i];
          lines_x->A[i] = lines_y->A[i];
          lines_x->f[i] = lines_y->f[i];
          lines_x->Q[i] = lines_y->Q[i];
        }
        for(i=0; i<nspline; i++)
        {
          sdatax[i] = sdatay[i];
          spointsx[i] = spointsy[i];
        }
      }

    }//end prior check

    //Every 1000 steps update focus region
    if(mc%1000 == 0)
    {
      pmax = -1.0;
      for(i=0; i< ncut; i++)
      {
        x = power[i]/Snx[i];
        if(x > pmax)
        {
          pmax = x;
          k = i;
        }
      }

      // define the focus region (only used if focus flag = 1)
      fcl = freq[k]-dff;
      fch = freq[k]+dff;
      Ac = power[k];
      if(fcl < freq[0]) fcl = freq[0];
      if(fch > freq[ncut-1]) fch = freq[ncut-1];

      //if(focus==1)fprintf(stdout,"Focusing on [%f %f] Hz...\n", fcl, fch);


      // set up proposal for frequency jumps
      xsm =0.0;
      baseav = 0.0;
      for(i=0; i< ncut; i++)
      {
        x = power[i]/Snx[i];
        if(x < 16.0) x = 1.0;
        if(x >= 16.0) x = 100.0;
        fprop[i] = x;
        xsm += x;
        baseav += log(Sbasex[i]);
      }
      baseav /= (double)(ncut);


      pmax = -1.0;
      for(i=0; i< ncut; i++)
      {
        fprop[i] /= xsm;
        if(fprop[i] > pmax) pmax = fprop[i];
      }

    }//end focus update


  }//End MCMC loop
  
  //Interpolate {spointsx,sdatax} ==> {freq,xint}
  CubicSplineGSL(nsx,spointsx,sdatax,ncut,freq,xint);
  for(i=0; i< ncut; i++) Sbase[i] = exp(xint[i]);

  full_spectrum_spline(Sline, Sbase, freq, data, lines_x);
  for(i=0; i< ncut; i++)
  {
    Sn[i] = Sbase[i]+Sline[i];
    bayesline->Sbase[i] = Sbase[i];
    bayesline->Sline[i] = Sline[i];
  }

  // return updated spectral model
  for(i=0; i< ncut; i++) Snf[i] = Snx[i];

  // re-map the array to 0..mx ordering
  for(i=0; i< lines_x->n; i++)
  {
    k = lines_x->larray[i];
    lines_y->f[i] = lines_x->f[k];
    lines_y->A[i] = lines_x->A[k];
    lines_y->Q[i] = lines_x->Q[k];
  }
  
  // return the last value of the chain
  for(i=0; i< lines_x->n; i++)
  {
    lines_x->f[i] = lines_y->f[i];
    lines_x->A[i] = lines_y->A[i];
    lines_x->Q[i] = lines_y->Q[i];
  }
  
  // check for outliers
  pmax = -1.0;
  for(i=0; i< ncut; i++)
  {
    x = power[i]/Snx[i];
    if(x > pmax) pmax = x;
  }
  
  *dan = pmax;
  *nsplinex = nsx;
  
  free(foc);
  free(xint);
  free(fprop);
  free(Snx);
  free(Sn);
  free(Sbase);
  free(Sbasex);
  free(Sline);
  free(sdatay);
  free(spointsy);
  
  destroy_lorentzianParams(lines_y);
  
}

void BayesLineMarkovianSplineOnly(struct BayesLineParams *bayesline, int nspline, int jj)
{

  /******************************************************************************/
  /*                                                                            */
  /*  Full-data spline fit (holding lines fixed from search phase)              */
  /*                                                                            */
  /******************************************************************************/

  int i,j,k;
  double mdn;
  double x,z;

  /* Make local pointers to BayesLineParams structure members */
  dataParams *data             = bayesline->data;
  splineParams *spline         = bayesline->spline;
  BayesLinePriors *priors      = bayesline->priors;
  double *Sna                  = bayesline->Sna;
  double *fa                   = bayesline->fa;
  double *Snf                  = bayesline->Snf;
  gsl_rng *r                   = bayesline->r;


  //initialize spline grid & PSD values
  k=0;
  for(j=0; j<nspline; j++)
  {
    x = data->fmin+(data->fhigh-data->fmin)*(double)(j)/(double)(nspline-1);
    if     (x <= fa[0])    z = Sna[0];
    else if(x >= fa[jj-1]) z = Sna[jj-1];
    else
    {
      i=k-10;
      mdn = 0.0;
      do
      {
        if(i>=0) mdn = fa[i];
        i++;
      } while(x > mdn);

      k = i;
      z = Sna[i];
//      printf("z=Sna[%i]=%g ",i,Sna[i]);

    }
//    printf("\n");

    spline->points[j] = x;
    spline->data[j]   = z;
  }

  // produce an initial spline fit to the smooth part of the spectrum
  SpecFitSpline(priors, bayesline->constantLogLFlag, 100000, fa, Sna, spline, Snf, jj, r);

}


void BayesLineMarkovianFocusedAnalysis(struct BayesLineParams *bayesline)
{
  /******************************************************************************/
  /*                                                                            */
  /*  Full spline/line MCMC stage                                               */
  /*                                                                            */
  /******************************************************************************/
  
  int j;
  double dan;
  
  BayesLineLorentzSplineMCMC(bayesline, 1.0, 200000, 0, 0, &dan);

  //alternate between targeting outliers and MCMCing full solution until outliers are gone or max iterations are reached
  j = 0;
  do
  {
    BayesLineLorentzSplineMCMC(bayesline, 1.0, 50000, 1, 0, &dan);

    BayesLineLorentzSplineMCMC(bayesline, 1.0, 50000, 0, 0, &dan);

    j++;
    
  } while (dan > 25.0 && j < 8);
  
}

void BayesLineFree(struct BayesLineParams *bptr)
{
  free(bptr->data);

  // storage for line model for a segment. These get updated and passed back from the fitting routine
  destroy_lorentzianParams(bptr->lines_x);
  destroy_lorentzianParams(bptr->lines_full);

  destroy_splineParams(bptr->spline);
  destroy_splineParams(bptr->spline_x);

  free(bptr->fa);
  free(bptr->Sna);
  free(bptr->freq);
  free(bptr->power);
  free(bptr->Snf);
  free(bptr->Sbase);
  free(bptr->Sline);

  free(bptr->spow);
  free(bptr->sfreq);

  /* set up GSL random number generator */
  gsl_rng_free(bptr->r);

  
  free(bptr);

}

void BayesLineSetup(struct BayesLineParams *bptr, double *freqData, double fmin, double fmax, double deltaT, double Tobs)
{
  int i;
  int n;
  int imin, imax;
  int j;
  double reFreq,imFreq;

  bptr->data = malloc(sizeof(dataParams));

  bptr->lines_x    = malloc(sizeof(lorentzianParams));
  bptr->lines_full = malloc(sizeof(lorentzianParams));

  bptr->spline   = malloc(sizeof(splineParams));
  bptr->spline_x = malloc(sizeof(splineParams));

  bptr->TwoDeltaT = deltaT*2.0;

  /* set up GSL random number generator */
  const gsl_rng_type *T = gsl_rng_default;
  bptr->r = gsl_rng_alloc (T);
  gsl_rng_env_setup();

  imin = (int)(fmin*Tobs);
  imax = (int)(fmax*Tobs)+1;
  n    = imax-imin;

  bptr->freq  = malloc((size_t)(sizeof(double)*(n)));
  bptr->power = malloc((size_t)(sizeof(double)*(n)));
  bptr->Snf   = malloc((size_t)(sizeof(double)*(n)));
  bptr->Sbase = malloc((size_t)(sizeof(double)*(n)));
  bptr->Sline = malloc((size_t)(sizeof(double)*(n)));

  for(i=0; i<n; i++)
  {
    j = i+imin;
    bptr->freq[i] = (double)(j)/Tobs;

    reFreq = freqData[2*j];
    imFreq = freqData[2*j+1];

    bptr->power[i] = (reFreq*reFreq+imFreq*imFreq);

  }

  // storage for data meta parameters
  create_dataParams(bptr->data,bptr->freq,n);

  // storage for line model for a segment. These get updated and passed back from the fitting routine
  create_lorentzianParams(bptr->lines_x,bptr->data->tmax);

  // storage for master line model
  create_lorentzianParams(bptr->lines_full,bptr->data->tfull);

  // start with a single line. This number gets updated and passed back
  bptr->lines_x->n = 1;

  //start with ~4 Hz resolution for the spline
  int nspline = (int)(bptr->data->fstep/bptr->data->fgrid)+2;
  int smodn   = nspline*bptr->data->sgmts;

  bptr->fa  = malloc((size_t)(sizeof(double)*(smodn+1)));
  bptr->Sna = malloc((size_t)(sizeof(double)*(smodn+1)));


  //start with ~4 Hz resolution for the spline
  nspline = (int)((bptr->data->fmax-bptr->data->fmin)/bptr->data->fgrid)+2;

  create_splineParams(bptr->spline,nspline);

  create_splineParams(bptr->spline_x,nspline);


  // Re-set dataParams structure to use full dataset
  bptr->data->flow  = bptr->data->fmin;
  bptr->data->fhigh = bptr->data->fmax;
  imax  = (int)(floor(bptr->data->fhigh*bptr->data->Tobs));
  imin  = (int)(floor(bptr->data->flow*bptr->data->Tobs));
  bptr->data->ncut  = imax-imin;

  /*
   The spow and sfreq arrays hold the frequency and power of the full section of data to be whitened
   The model parameters are the number of Lorentzians, nx, and their frequency ff, amplitude AA and
   quality factor QQ, as well as the number of terms in the power law fit to the smooth part of the
   spectrum, segs, and their amplitudes at 100 Hz, SA, and their spectral slopes SP. The maximum number
   of terms in the Lorentzian model is tmax (around 200-300 should cover most spectra), and the
   maximum number of terms in the power law fit is segmax (aound 12-20 should be enough). When called
   from another MCMC code, the LorentzMCMC code needs to be passed the most recent values for the
   noise model. Arrays to hold them have to be declared somewhere in advance. The updated values are
   passed back. The first entry in the call to LoretnzMCMC are the number of iterations. It probably
   makes sense to do 100 or so each time it is called.
   */

  bptr->spow  = malloc((size_t)(sizeof(double)*(bptr->data->ncut)));
  bptr->sfreq = malloc((size_t)(sizeof(double)*(bptr->data->ncut)));

}

void BayesLineSearch(struct BayesLineParams *bptr, double *freqData, double fmin, double fmax, double deltaT, double Tobs)
{

  int i;
  int imin,imax;
  int j;


  int jj = 0;

  /******************************************************************************/
  /*                                                                            */
  /*  Rapid non-Markovian fit over small bandwidth blocks of data               */
  /*                                                                            */
  /******************************************************************************/

  BayesLineNonMarkovianFit(bptr, &jj);

  // maximum number of terms in the Lorentzian model
  bptr->data->tmax = 4*bptr->lines_full->n;

  if(bptr->data->tmax    < 20) bptr->data->tmax    = 20;
  if(bptr->lines_full->n < 1 ) bptr->lines_full->n = 1;

    // Re-set dataParams structure to use full dataset
  bptr->data->flow  = bptr->data->fmin;
  imax  = (int)(floor(bptr->data->fhigh*bptr->data->Tobs));
  imin  = (int)(floor(bptr->data->flow*bptr->data->Tobs));
  bptr->data->ncut  = imax-imin;

  /******************************************************************************/
  /*                                                                            */
  /*  Full-data spline fit (holding lines fixed from search phase)              */
  /*                                                                            */
  /******************************************************************************/

  BayesLineMarkovianSplineOnly(bptr, bptr->spline->n, jj);

  for(j=0; j<bptr->spline->n; j++)
  {
    bptr->spline_x->points[j] = bptr->spline->points[j];
    bptr->spline_x->data[j]   = bptr->spline->data[j];
  }

  bptr->data->flow  = bptr->data->fmin;

  imin  = (int)(floor(bptr->data->flow*bptr->data->Tobs));

  for(i=0; i< bptr->data->ncut; i++)
  {
    j = i + imin - bptr->data->nmin;
    bptr->spow[i] = bptr->power[j];
    bptr->sfreq[i] = bptr->freq[j];
  }

  /******************************************************************************/
  /*                                                                            */
  /*  Full spline/line MCMC stage                                               */
  /*                                                                            */
  /******************************************************************************/

  BayesLineMarkovianFocusedAnalysis(bptr);


  double dan;

  BayesLineLorentzSplineMCMC(bptr, 1.0, 200000, 0, 0, &dan);

}

void BayesLineRJMCMC(struct BayesLineParams *bayesline, double *freqData, double *psd, double *invpsd, double *splinePSD, int N, int cycle, double beta, int priorFlag)
{
  int i,j;
  int imin,imax;
  double reFreq,imFreq;
  double dan;

  //at frequencies below fmin output SAmax (killing lower frequencies in any future inner products)
  imax = (int)(floor(bayesline->data->fmax * bayesline->data->Tobs));
  imin = (int)(floor(bayesline->data->fmin * bayesline->data->Tobs));

  for(i=0; i< bayesline->data->ncut; i++)
  {
    j = i + imin;

    reFreq = freqData[2*j];
    imFreq = freqData[2*j+1];

    bayesline->spow[i]  = (reFreq*reFreq+imFreq*imFreq);
    bayesline->Sbase[i] = splinePSD[i];
  }

  BayesLineLorentzSplineMCMC(bayesline, beta, cycle, 0, priorFlag, &dan);

  /******************************************************************************/
  /*                                                                            */
  /*  Output PSD in format for BayesWave                                        */
  /*                                                                            */
  /******************************************************************************/

  for(i=0; i<N/2; i++)
  {
    if(i>=imin && i<imax)
    {
      psd[i] = bayesline->Snf[i-imin];
      splinePSD[i] = bayesline->Sbase[i-imin];
    }
    else
    {
      psd[i] = bayesline->priors->SAmax;
      splinePSD[i] = bayesline->priors->SAmax;
    }
    invpsd[i] = 1./psd[i];
  }
}


void print_line_model(FILE *fptr, struct BayesLineParams *bayesline)
{
  int j;
  
  fprintf(fptr,"%i ", bayesline->lines_full->n);
  for(j=0; j< bayesline->lines_full->n; j++) fprintf(fptr,"%e %e %e ", bayesline->lines_full->f[j], bayesline->lines_full->A[j], bayesline->lines_full->Q[j]);
  fprintf(fptr,"\n");

}
void print_spline_model(FILE *fptr, struct BayesLineParams *bayesline)
{
  int j;

  fprintf(fptr,"%i ", bayesline->spline_x->n);
  for(j=0; j<bayesline->spline_x->n; j++) fprintf(fptr,"%lg %lg ",bayesline->spline_x->points[j],bayesline->spline_x->data[j]);
  fprintf(fptr,"\n");
}

void parse_line_model(FILE *fptr, struct BayesLineParams *bayesline)
{
  int j;

  fscanf(fptr,"%i",&bayesline->lines_x->n);
  for(j=0; j< bayesline->lines_full->n; j++) fscanf(fptr,"%lg %lg %lg",&bayesline->lines_full->f[j],&bayesline->lines_full->A[j],&bayesline->lines_full->Q[j]);
}
void parse_spline_model(FILE *fptr, struct BayesLineParams *bayesline)
{
  int j;

  fscanf(fptr,"%i",&bayesline->spline_x->n);
  for(j=0; j<bayesline->spline_x->n; j++)fscanf(fptr,"%lg %lg",&bayesline->spline_x->points[j],&bayesline->spline_x->data[j]);
}


