import numpy as np
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.ticker as mticker
import matplotlib.colors
from chainconsumer import ChainConsumer
import os

if not os.path.exists('plots'):
    os.makedirs('plots')

aqua='#008E97'
orange='#FC4C02'
blue='#005778'


figsize = (20,5)

# Load files
data = 'data/power_data.dat'
res  = 'data/power_residual.dat'
wave = 'data/power_reconstruction.dat'
rec_files = glob.glob('catalog*/*power_reconstruction.dat')
chain_files = glob.glob('catalog*/*chain.dat')
entry_files = glob.glob('catalog*/entries.dat')
entries=np.genfromtxt(entry_files[0],usecols=(0),dtype='str')

# Data
fig,ax=plt.subplots(figsize=figsize)
f = np.loadtxt(data,usecols=(0))
A = np.loadtxt(data,usecols=(1))
E = np.loadtxt(data,usecols=(2))
ax.semilogy(f*1000,A,color=orange, label='data')

ax.set_xlabel(r'$f\ ({\rm mHz})$')
ax.set_ylabel(r'$|A|$')
ax.legend()

plt.savefig('plots/tdi_data.png', dpi=300)

# Residual
fig,ax=plt.subplots(figsize=figsize)

f = np.loadtxt(data,usecols=(0))
A = np.loadtxt(data,usecols=(1))
E = np.loadtxt(data,usecols=(2))
ax.semilogy(f*1000,A,color=orange,label='data')

f = np.loadtxt(res,usecols=(0))
A = np.loadtxt(res,usecols=(1))
E = np.loadtxt(res,usecols=(2))
ax.semilogy(f*1000,A,color=aqua,label='residual')

ax.set_xlabel(r'$f\ ({\rm mHz})$')
ax.set_ylabel(r'$|A|$')
ax.legend()

plt.savefig('plots/tdi_residual.png', dpi=300)

# Joint Reconstructions
fig,ax=plt.subplots(figsize=figsize)

f = np.loadtxt(data,usecols=(0))
A = np.loadtxt(data,usecols=(1))
E = np.loadtxt(data,usecols=(2))
ax.semilogy(f*1000,A,color=orange,label='data')

f = np.loadtxt(wave,usecols=(0))
A50 = np.loadtxt(wave,usecols=(1))
A25 = np.loadtxt(wave,usecols=(2))
A75 = np.loadtxt(wave,usecols=(3))
A05 = np.loadtxt(wave,usecols=(4))
A95 = np.loadtxt(wave,usecols=(5))
ax.fill_between(f*1000,A05,A95,alpha=0.25,color=aqua)
ax.fill_between(f*1000,A25,A75,alpha=0.25,color=aqua)
ax.semilogy(f*1000,A50,color=aqua,label='ucb model')

ax.set_xlabel(r'$f\ ({\rm mHz})$')
ax.set_ylabel(r'$|A|$')
ax.legend()

plt.savefig('plots/tdi_reconstructions_joint.png', dpi=300)

# Net Reconstructions
fig,ax=plt.subplots(figsize=figsize)

f = np.loadtxt(data,usecols=(0))
A = np.loadtxt(data,usecols=(1))
E = np.loadtxt(data,usecols=(2))
ax.semilogy(f*1000,A,color=orange,label='data')

for index,file in enumerate(rec_files):
    f = np.loadtxt(file,usecols=(0))
    A50 = np.loadtxt(file,usecols=(1))
    A25 = np.loadtxt(file,usecols=(2))
    A75 = np.loadtxt(file,usecols=(3))
    A05 = np.loadtxt(file,usecols=(4))
    A95 = np.loadtxt(file,usecols=(5))
    ax.fill_between(f*1000,A05,A95,alpha=0.25,color=aqua)
    ax.fill_between(f*1000,A25,A75,alpha=0.25,color=aqua)
    if index==0:
        ax.semilogy(f*1000,A50,color=aqua,label='waveform reconstructions')
    else:
        ax.semilogy(f*1000,A50,color=aqua)
    ax.legend()
    plt.savefig('plots/tdi_reconstructions.png', dpi=300)
    
# Individual Reconstructions
def last_nonzero(arr, axis, invalid_val=-1):
    mask = arr!=0
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)

def first_nonzero(arr, axis, invalid_val=-1):
    mask = arr!=0
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)

for file in rec_files:
    name = file.replace('/','_').split('_')[2]

    fig,[ax1,ax2]=plt.subplots(1,2,figsize=figsize,sharex=True,sharey=True)
    
    f = np.loadtxt(data,usecols=(0))
    A = np.loadtxt(data,usecols=(1))
    E = np.loadtxt(data,usecols=(2))
    ax1.semilogy(f*1000,A,color=orange,label=r'data')
    ax2.semilogy(f*1000,E,color=orange)
        
    f = np.loadtxt(file,usecols=(0))
    A50 = np.loadtxt(file,usecols=(1))
    A25 = np.loadtxt(file,usecols=(2))
    A75 = np.loadtxt(file,usecols=(3))
    A05 = np.loadtxt(file,usecols=(4))
    A95 = np.loadtxt(file,usecols=(5))
    E50 = np.loadtxt(file,usecols=(6))
    E25 = np.loadtxt(file,usecols=(7))
    E75 = np.loadtxt(file,usecols=(8))
    E05 = np.loadtxt(file,usecols=(9))
    E95 = np.loadtxt(file,usecols=(10))
    
    ax1.fill_between(f*1000,A05,A95,alpha=0.25,color=aqua)
    ax1.fill_between(f*1000,A25,A75,alpha=0.25,color=aqua)
    ax1.semilogy(f*1000,A50,color=aqua, label=name)

    ax2.fill_between(f*1000,E05,E95,alpha=0.25,color=aqua)
    ax2.fill_between(f*1000,E25,E75,alpha=0.25,color=aqua)
    ax2.semilogy(f*1000,E50,color=aqua)

    imin=first_nonzero(A50,axis=0)
    imax=last_nonzero(A50,axis=0)
    ax1.set_xlim(f[imin]*1000,f[imax]*1000)
    ax1.legend()

    ax1.set_xlabel(r'$f\ ({\rm mHz})$')
    ax2.set_xlabel(r'$f\ ({\rm mHz})$')
    ax1.set_ylabel(r'$|A|$')
    ax2.set_ylabel(r'$|E|$')
    ax2.yaxis.set_tick_params(labelbottom=True)
    plt.savefig('plots/'+name+'_reconstructions.png', dpi=300)
    
# Individual Corner Plots
parameters=[r'$f_0$',r'$\dot f$',r'$\log_{10} A$',r'$\phi$',r'$\cos\theta$',r'$\cos\iota$',r'$\psi$',r'$\varphi_0$']
for file in chain_files:
    c=ChainConsumer()
    chain=np.loadtxt(file,usecols=(0,1,2,3,4,5,6,7))
    chain[:,2] = np.log10(chain[:,2])
    name = file.replace('/','_').split('_')[2]
    c.add_chain(chain,parameters=parameters,color=aqua,name=name)
    c.configure(sigmas=[0,1,2,3])
    #c.plotter.plot(filename='plots/'+name+'_corner.png')
    
# Joint Sky Map
parameters=[r'$\phi$',r'$\cos\theta$']

c=ChainConsumer()
arrays = [np.loadtxt(f, usecols=(3,4)) for f in chain_files]
chain=np.concatenate(arrays)
c.add_chain(chain,parameters=parameters,name=' ',plot_contour=True, color=aqua)
c.configure(plot_hists=False, sigmas=np.linspace(0,3,20))
plot=c.plotter.plot(figsize=(10,5),extents=[(0,2*np.pi),(-1,1)],filename='plots/pe_joint_sky.png')

# Combined Sky Map
parameters=[r'$\phi$',r'$\cos\theta$']

c=ChainConsumer()
for f in chain_files:
    chain=np.loadtxt(f, usecols=(3,4))
    name = f.replace('/','_').split('_')[2]
    c.add_chain(chain,parameters=parameters,name=name,plot_contour=True)

c.configure(plot_hists=False,sigmas=[0,1,2])
plot=c.plotter.plot(figsize=(10,5),extents=[(0,2*np.pi),(-1,1)],filename='plots/pe_net_sky.png')

# Joint f-A Plane
parameters=[r'$f_0$',r'$\log_{10} A$']

c=ChainConsumer()
arrays = [np.loadtxt(f, usecols=(0,2)) for f in chain_files]
chain=np.concatenate(arrays)
chain[:,1]=np.log10(chain[:,1])
c.add_chain(chain,parameters=parameters,name=' ',plot_contour=True, color=aqua)
c.configure(plot_hists=False, sigmas=np.linspace(0,3,20),bins=10.0)
plot=c.plotter.plot(figsize=(10,5),filename='plots/pe_joint_fA.png')

# Combined f-A Plane
parameters=[r'$f_0$',r'$A$']

c=ChainConsumer()
for f in chain_files:
    chain=np.loadtxt(f, usecols=(0,2))
    name = f.replace('/','_').split('_')[2]
    c.add_chain(chain,parameters=parameters,name=name,plot_contour=True)

c.configure(plot_hists=False,sigmas=[0,1,2])
plot=c.plotter.plot(figsize=(10,5),filename='plots/pe_net_fA.png')

# Joint f-fdot Plane
parameters=[r'$f_0$',r'$\dot f$']

c=ChainConsumer()
arrays = [np.loadtxt(f, usecols=(0,1)) for f in chain_files]
chain=np.concatenate(arrays)
c.add_chain(chain,parameters=parameters,name=' ',plot_contour=True, color=aqua)
c.configure(plot_hists=False, sigmas=np.linspace(0,3,20))
plot=c.plotter.plot(figsize=(10,5),filename='plots/pe_joint_ffdot.png')

# Combined f-fdot Plane
parameters=[r'$f_0$',r'$\dot f$']

c=ChainConsumer()
for f in chain_files:
    chain=np.loadtxt(f, usecols=(0,1))
    name = f.replace('/','_').split('_')[2]
    c.add_chain(chain,parameters=parameters,name=name,plot_contour=True)

c.configure(plot_hists=False,sigmas=[0,1,2])
plot=c.plotter.plot(figsize=(10,5),filename='plots/pe_net_ffdot.png')

# Correlation map
corr_files = glob.glob('catalog*/correlation_matrix.dat')
norm = matplotlib.colors.Normalize(-1,1)
colors = [[norm(-1.0), orange],
          [norm(0), "white"],
          [norm( 1.0), aqua]]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

file = np.loadtxt(corr_files[0])
fig,ax=plt.subplots(figsize=(10,10))
ax = sns.heatmap(file,cmap=cmap, vmin=-1, vmax=1,xticklabels = 8, yticklabels = 8, linewidths=.1, cbar=False)

plt.savefig('plots/correlation.png', dpi=300)


        


