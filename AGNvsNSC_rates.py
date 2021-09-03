#! /usr/bin/env python
#AGNvsNSC_rates.py

#Makes plot of relative rate of BBH mergers in AGN vs gas-free, SMBH containing NSC
#Uses 'Drake equation' style reasoning

from pylab import *

#import sys, os, time, string, math, commands, subprocess
import sys, os, time, string, math, subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import AxesGrid

#SI units
Msun=1.99e30 #kg per solar mass
Rsun=6.95e8 #meters per solar radius
G=6.67e-11
c=3e8
sigma_SB=5.7e-8 #stefan-boltzmann const
yr=3.15e7 #seconds per year
pc=3.086e16 #meters per parsec
AU=1.496e11 #meters per AU
h=6.626e-34 #planck const
kB=1.38e-23 #boltzmann const
m_p=1.67e-27 #mass of proton
sigma_T=6.65e-29 #Thomson xsec
PI=3.1415926

def RelativeRate(time,frac,tQbin,fbin_AoverQ):
    #rate of mergers in AGN/quiescent nuclei=fAGN*(tQbin/time)*fbin(AGN/QGN)
    #where fAGN is fraction of galaxies that are active
    #---technically relation only works for fAGN smallish
    #tQbin is the average binary lifetime in a quiescent nucleus
    #---get from Antonini&Perets's rates
    #fbin(AGN/QGN) is the fraction of BH that are in binaries in an AGN/same frac in QGN
    #---fbin should always be >=1
    #time=relevant timescale:
    #---if tAbin is the average binary lifetime in an active nucleus:
    #-----time is AGN lifetime (tauAGN) for tauAGN>tAbin
    #-----time is tAbin for tauAGN=tAbin
    #-----time is tAbin/tauAGN for tauAGN<tAbin
    relrate=frac*(tQbin/time)*fbin_AoverQ

    return relrate

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    #Saavik found this on stack overflow thanks to Tom Callister
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

if __name__ == "__main__":
    #to make plot
    #AGN fraction vs AGN lifetime, color bar for Active/Quiescent rate
    #BEGIN PHYSICAL INPUTS:
    fAGN_min=0.001 #frac of galaxies that are AGN--ie that mess with the dynamics (minimum)
    fAGN_max=0.2 #frac of galaxies that are AGN--ie that mess with the dynamics (maximum)
    tauAGN_min=0.1 #AGN liftime, minimum, units of Myr
    tauAGN_max=100.0 #AGN liftime, maximum,units of Myr
    tQbin=40000.0 #average lifetime of binary in gas-free NSC, units of Myr
    #derive above from Antonini & Perets
    fbin_AoverQ=1.0 #binary fraction in AGN/binary fraction in NSC, minimum of 1
    tbin_AGN=1.0 #average lifetime of binary in AGN, units of Myr--sets turnover
    
    #BEGIN DISPLAY INPUTS:
    #pick color map and limits
    cm=plt.cm.get_cmap('RdBu')
    #format=left, bottom, width, height
    rect1=0.12,0.12,0.65,0.85
    rect2=(rect1[0]+rect1[2]+0.03),rect1[1],0.03,rect1[3]
        
    #make figure
    fig1=plt.figure(1)
    #add axes
    ax1=fig1.add_axes(rect1)
    cax=fig1.add_axes(rect2)
    #label them
    ax1.yaxis.set_label_coords(-0.1, 0.5)
    ax1.set_ylabel(r"$f_{AGN}$", fontsize=18)
    ax1.xaxis.set_label_coords(0.5, -0.06)
    ax1.set_xlabel(r"$\tau_{AGN} (Myr)$", fontsize=18)
    

    #set up range for x-axis
    #AGN lifetime
    log_tauAGN=np.arange(log10(tauAGN_min), log10(tauAGN_max), 0.01)
    tauAGN=pow(10,log_tauAGN)
    #range for y-axis
    plt.ylim(-3.0,0.0)
    #END DISPLAY INPUTS

    #set up y-axis variables
    #AGN fraction
    log_fAGN=np.arange(log10(fAGN_min),log10(fAGN_max),0.01)
    fAGN=pow(10,log_fAGN)

    #relative rates
    #governing timescale is tauAGN where avg bin lifetime is shorter,
    #fraction of avg bin lifetime where tauAGN is shorter;
    tA=np.where(tauAGN>=tbin_AGN, tauAGN, tbin_AGN/tauAGN)
    log_tA=log10(tA)
    #set up display grid
    X, Y = np.meshgrid(tA, fAGN)
    #compute relative rates and log it
    R_AoverQ=RelativeRate(X, Y, tQbin, fbin_AoverQ)
    log_R_AoverQ=log10(R_AoverQ)

    #set up colormap
    colormax=np.max(log_R_AoverQ)
    colormin=np.min(log_R_AoverQ)
    cmid = 1 - colormax / (colormax + abs(colormin))
    shifted_cmap = shiftedColorMap(cm, midpoint=cmid, name='shifted')

    #plot it
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.pcolormesh(tauAGN, fAGN, log_R_AoverQ, cmap=shifted_cmap)

    #add actual colorbar with same normalization as plotting used
    norm=mpl.colors.Normalize(vmin=colormin, vmax=colormax)
    cb=mpl.colorbar.ColorbarBase(cax,cmap=shifted_cmap,norm=norm, orientation='vertical')
    cb.set_label(r"$log \ ({\cal{R}}_{A/G})$", fontsize=18)

    #dump plots to figure files/show on screen:
    savefig('AGNvsNSC_rates.eps')
    savefig('AGNvsNSC_rates.png')
    show()


