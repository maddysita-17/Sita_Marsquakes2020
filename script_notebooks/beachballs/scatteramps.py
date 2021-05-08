import pandas as pd
import numpy as np
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from obspy.imaging.beachball import beachball
from obspy.imaging.beachball import beach
from obspy.imaging.source import plot_radiation_pattern
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/beachballs/csvs/'

def abserr(vals, errx, erry, x, y):
    if x == 1 or y == 1:
        slope = x/y
        sig1 = ((x/y)**2)*((errx/x)**2)
        sig2 = ((erry/y)**2)*((x/y)**2)
        sig = sig1 + sig2
        print('Error: ', sig)

        xarr = np.linspace(vals[0], vals[1], 11)

        e1 = slope+sig ; e2 = slope-sig
        eline1 = e1*xarr ; eline2=e2*xarr

    else:
        slope = y/x
        print('Slope: ', slope)
        sig2 = (errx/x)**2 + (erry/y)**2
        sig = np.sqrt(sig2) * slope
        print('Error: ', sig)

        xarr = np.linspace(vals[0], vals[1], 11)

        e1 = slope+sig ; e2 = slope-sig
        eline1 = e1*xarr ; eline2=e2*xarr

    return xarr,eline1, eline2

def scatplt(event, df, Pamp, SHamp, SVamp, Perr, SHerr, SVerr, axvals):
    PSV = Pamp/SVamp
    PSH = Pamp/SHamp
    SHSV = SHamp/SVamp
    print(PSH)

    fig, ax = plt.subplots(2,3,sharex=True,sharey=True)
    fig.suptitle(event)

    try:
        #-----SH/SV Angle1-------
        sns.scatterplot(ax=ax[0,0], data=df, x='SV', y='SH', hue='Angle1')
        sns.scatterplot(ax=ax[1,0], data=df, x='SV', y='SH', hue='Sum')
        # ax[0,0].plot([axvals[0],axvals[1]], [axvals[0],SHSV*axvals[1]], 'k--', alpha=0.5)
        # ax[1,0].plot([axvals[0],axvals[1]], [axvals[0],SHSV*axvals[1]], 'k--', alpha=0.5)

        x, e1SVSH, e2SVSH = abserr(axvals, SVerr, SHerr, SVamp, SHamp)
        ax[0,0].fill_between(x, e1SVSH, e2SVSH, color='b', alpha=0.05)
        ax[1,0].fill_between(x, e1SVSH, e2SVSH, color='b', alpha=0.05)

        #------P/SV Angle2-------
        sns.scatterplot(ax=ax[0,1], data=df, x='SV', y='P', hue='Angle2')
        sns.scatterplot(ax=ax[1,1], data=df, x='SV', y='P', hue='Sum')
        # ax[0,1].plot([axvals[0],axvals[1]], [axvals[0],PSV*axvals[1]], 'k--', alpha=0.5)
        # ax[1,1].plot([axvals[0],axvals[1]], [axvals[0],PSV*axvals[1]], 'k--', alpha=0.5)

        x, e1PSV, e2PSV = abserr(axvals, SVerr, Perr, SVamp, Pamp)
        ax[0,1].fill_between(x, e1PSV, e2PSV, color='b', alpha=0.05)
        ax[1,1].fill_between(x, e1PSV, e2PSV, color='b', alpha=0.05)

        #------P/SH Angle3--------
        sns.scatterplot(ax=ax[0,2], data=df, x='SH', y='P', hue='Angle3')
        sns.scatterplot(ax=ax[1,2], data=df, x='SH', y='P', hue='Sum')
        # ax[0,2].plot([axvals[0],axvals[1]], [axvals[0],PSH*axvals[1]], 'k--', alpha=0.5)
        # ax[1,2].plot([axvals[0],axvals[1]], [axvals[0],PSH*axvals[1]], 'k--', alpha=0.5)

        x, e1PSH, e2PSH = abserr(axvals, SHerr, Perr, SHamp, Pamp)
        ax[0,2].fill_between(x, e1PSH, e2PSH, color='b', alpha=0.05)
        ax[1,2].fill_between(x, e1PSH, e2PSH, color='b', alpha=0.05)

        plt.show()

    except:
        print('Angle values could not be binned')

    return


s0173a = pd.read_csv(path + 'S0173a.csv')
scatplt('S0173a', s0173a, 195, 144, -176, 18, 47, 47, [-1.5,1.5])

s0173ab = pd.read_csv(path + 'S0173ab.csv')
scatplt('S0173ab', s0173ab, 335, 558, 412, 18, 47, 47, [1,2])

s0235b = pd.read_csv(path + 'S0235b.csv')
scatplt('S0235b', s0235b, -80, -318, -234, 10, 26, 26, [-2.5,0])

s0325a = pd.read_csv(path + 'S0325a.csv')
scatplt('S0325a', s0325a, 132, 1, 257, 26, 93, 93, [0,5])

s0325ab = pd.read_csv(path + 'S0325ab.csv')
scatplt('S0325ab', s0325ab, 228, -365, 362, 26, 93, 93, [-3,4])
