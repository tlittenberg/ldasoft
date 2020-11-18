# lisacattools.py
# Module containing functions for working with LISA (mock) catalog data

# common modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

# constants
ep_const = 23.4 # obliquity of the ecliptic (degrees)
a_G_const = 192.85948 # r.a. of galacitc North Pole (degrees)
d_G_const = 27.12825 # declination of galactic North pole (degrees)
l_NCP_const = 122.93192 # galactic longitude of North Celestial Pole (degrees)

## Sky Coordinate Funcitons

# deg2rad
def deg2rad(degrees):
    # Convert degrees to radians
    radians = degrees*np.pi/180.0
    return radians

# rad2deg
def rad2deg(radians):
    # Convert radians to degrees
    degrees = radians * 180 / np.pi
    return degrees

# eclp2eqi
def eclp2eqi(lamb,beta,ep=ep_const):  
    # Convert ecliptic to equitorial coordinates
    # the arguments of the function are ecliptic longitude lambda, ecliptic latitude beta, and the obliquity of the ecliptic (all in degrees)
    ep = deg2rad(ep)
    lamb = deg2rad(lamb)
    beta = deg2rad(beta)
    alpha = np.arctan2(np.cos(beta) * np.sin(lamb) * np.cos(ep) - np.sin(beta) * np.sin(ep), np.cos(beta) * np.cos(lamb))
    delta = np.arcsin(np.sin(beta) * np.cos(ep) + np.cos(beta) * np.sin(ep) * np.sin(lamb))

    alpha = rad2deg(alpha)
    delta = rad2deg(delta)
    
    return alpha,delta

# eqi2gal
def eqi2gal(alpha, delta, a_G = a_G_const, d_G = d_G_const, l_NCP = l_NCP_const):
    # Convert equitorial to galactic coordinates 
    # arguments are right ascension, declination, r.a. of glactic NP, dec of galactic NP, and galactic longitude of north celestial pole (all in degrees)
    a_G = deg2rad(a_G)
    d_G = deg2rad(d_G)
    l_NCP = deg2rad(l_NCP)
    alpha = deg2rad(alpha)
    delta = deg2rad(delta)
    
    b = np.arcsin(np.sin(delta) * np.sin(d_G) + np.cos(delta) * np.cos(d_G) * np.cos(alpha - a_G))
    l = l_NCP - np.arctan2(np.cos(delta) * np.sin(alpha - a_G), np.sin(delta) * np.cos(d_G) - np.cos(delta) * np.sin(d_G) * np.cos(alpha - a_G))
    
    b = rad2deg(b)
    l  = rad2deg(l)
    return l,b

# eclp2gal
def eclp2gal(lamb, beta, ep = ep_const, a_G = a_G_const, d_G = d_G_const, l_NCP = l_NCP_const):
    # Convert equitorial to galactic coordinates 
    # arguments are right ascension, declination, r.a. of glactic NP, dec of galactic NP, and galactic longitude of north celestial pole (all in degrees)
    a,d = eclp2eqi(lamb,beta)
    l,b = eqi2gal(a,d)
    l = np.where(l>180,l-360,l)
    return l,b

# getGalcoord
def getGalcoord(df, ep = ep_const, a_G = a_G_const, d_G = d_G_const, l_NCP = l_NCP_const):
    # convert Pandas dataframe (either of sources or a chain) from Ecliptic to Galactic Coords
    lamb = df['Ecliptic Longitude']
    beta = np.pi/2-np.arccos(np.array(df['coslat']))
    lamb = rad2deg(lamb)
    beta = rad2deg(beta)
    l,b = eclp2gal(lamb,beta)
    df.insert(len(df.columns),'Galactic Longitude',l,True)
    df.insert(len(df.columns),'Galactic Latitude',b,True)
    return 

# ellipse_area
def ellipse_area(df):
    # Compute 90% ellipse area from chain samples
    cov = np.array(df.cov())
    ev = np.linalg.eigvals(cov)
    ax = 2.0*np.sqrt(4.605*ev)  # the 4.605 corresponds to 90%, for 95% we'd use 5.991
    return np.pi*ax[0]*ax[1]

# confidence_ellipse
def confidence_ellipse(df, ax, n_std=1.0, facecolor='none', **kwargs):   
    # Plot an error ellipse on some axes. It takes a 2D data frame of parameters as its input
    mean = np.array(df.mean())

    cov = np.array(df.cov())
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean[0], mean[1])

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

## Luminosity Distance Functions

# get_DL
def get_DL(df):
    # Estimate luminosity distance (in kpc) from GW amplitude, frequency, and frequency derivative
    c = 2.99e8  # speed of light in m/s
    kpc2m = 3.086e19 # 1 kiloparsec in meters
    r = (5/(96*(np.pi**2)))*(c/df['Amplitude'])*(df['Frequency Derivative'])/np.power(df['Frequency'],3)*(1/kpc2m)
    df.insert(len(df.columns),'Luminosity Distance',r,True)
    return
    



## Miscellaneous Functions


# getChain
def getChain(cat,idx,path = ''):
    # load the MCMC chain corresponding to a particular source in the catalog
    chain = pd.read_hdf(os.path.join(path,cat.loc[idx]['chain file']),key = idx + '_chain')
    return chain

# getSciRD
def getSciRD(f,Tobs):
    # Compute the LISA Sensitivity curve given a set of Fourier frequencies and an observation time
    # Curve is in units of strain amplitude and is defined by SciRD available at: https://www.cosmos.esa.int/documents/678316/1700384/SciRD.pdf

    f1 = np.array(0.4*1e-3)
    f2 = np.array(25*1e-3)
    R = 1+np.power((f/f2),2)
    S2 = 2.6*1e-41
    S1 = 5.76*1e-48*(1+np.power((f1/f),2))
    Sh = (1/2)*(20/3)*(S1/(np.power((2*np.pi*f),4))+S2)*R
    Sha = np.sqrt(Sh/Tobs)
    return Sha




