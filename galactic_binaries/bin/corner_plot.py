#!/usr/bin/env python

import numpy as np
from matplotlib import rcParams
import corner

data = np.loadtxt("PE.dat", usecols=(0,1,5))

rcParams.update({'font.size': 8})

figure = corner.corner(data, quantiles=[0.1, 0.5, 0.9], show_titles=True, title_kwargs={"fontsize": 8},labels=[r"$f$ [mHz]", r"$\log_{10}\frac{df}{dt}$ [s$^{-2}$]",r"$\log_{10}\frac{d^2f}{dt^2}$ [s$^{-3}$]"])

figure.savefig("PE.png", dpi=300)
