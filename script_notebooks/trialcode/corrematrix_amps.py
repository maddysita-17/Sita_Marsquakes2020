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

def abserr(errx, erry, x, y, neg=False):
    slope = y/x
    print('Slope: ', slope)
    sig2 = (errx/x)**2 + (erry/y)**2
    sig = np.sqrt(sig2) * slope
    print('Error:', sig)

    if neg==False:
        xarr = np.linspace(0, 5, 6)

        e1 = slope+sig ; e2 = slope-sig
        eline1 = e1*xarr ; eline2=e2*xarr

    elif neg==True:
        xarr = np.linspace(-2,2,6)

        e1 = slope+sig ; e2 = slope-sig
        print(e1, e2)
        eline1 = e1*xarr ; eline2=e2*xarr

    return xarr,eline1, eline2

def scatplt(event, df, Pamp, SHamp, SVamp, Perr, SHerr, SVerr):
    PSV = Pamp/SVamp
    PSH = Pamp/SHamp
    SHSV = SHamp/SVamp
    print(PSH)

    fig, ax = plt.subplots(2,3)
    fig.suptitle(event)

    #-----SH/SV Ratio1-------
    sns.scatterplot(ax=ax[0,0], data=df, x='SV', y='SH', hue='Ratio1')
    sns.scatterplot(ax=ax[1,0], data=df, x='SV', y='SH', hue='Sum')
    ax[0,0].plot([0,5], [0,SHSV*5], 'k--', alpha = 0.5)
    ax[1,0].plot([0,5], [0,SHSV*5], 'k--', alpha = 0.5)

    x, e1SVSH, e2SVSH = abserr(SVerr, SHerr, SVamp, SHamp)
    ax[0,0].fill_between(x, e1SVSH, e2SVSH, color='b', alpha=0.05)
    ax[1,0].fill_between(x, e1SVSH, e2SVSH, color='b', alpha=0.05)

    #------P/SV Ratio2-------
    sns.scatterplot(ax=ax[0,1], data=df, x='SV', y='P', hue='Ratio2')
    sns.scatterplot(ax=ax[1,1], data=df, x='SV', y='P', hue='Sum')
    ax[0,1].plot([0,5], [0,PSV*5], 'k--', alpha=0.5)
    ax[1,1].plot([0,5], [0,PSV*5], 'k--', alpha=0.5)

    x, e1PSV, e2PSV = abserr(SVerr, Perr, SVamp, Pamp)
    ax[0,1].fill_between(x, e1PSV, e2PSV, color='b', alpha=0.05)
    ax[1,1].fill_between(x, e1PSV, e2PSV, color='b', alpha=0.05)

    #------P/SH Ratio3--------
    sns.scatterplot(ax=ax[0,2], data=df, x='SH', y='P', hue='Ratio3')
    sns.scatterplot(ax=ax[1,2], data=df, x='SH', y='P', hue='Sum')
    ax[0,2].plot([0,5], [0,PSH*5], 'k--', alpha=0.5)
    ax[1,2].plot([0,5], [0,PSH*5], 'k--', alpha=0.5)

    x, e1PSH, e2PSH = abserr(SHerr, Perr, SHamp, Pamp)
    ax[0,2].fill_between(x, e1PSH, e2PSH, color='b', alpha=0.05)
    ax[1,2].fill_between(x, e1PSH, e2PSH, color='b', alpha=0.05)

    plt.show()
    return


s0173a = pd.read_csv('S0173a.csv')
scatplt('S0173a', s0173a, 195, 144, -176, 18, 47, 47)

s0235b = pd.read_csv('S0235b.csv')

# sns.set_style("darkgrid", {"axes.facecolor": ".9"})
# grid = sns.pairplot(faults, kind='reg')
# plt.show()

# sns.scatterplot(data=s0173a, x="SV", y="SH", hue="Ratio1")
# plt.plot([0,-176],[0,144], '--')
# plt.title('S0173a SH/SV')
# plt.show()


s0325a = pd.read_csv('S0325a.csv')
scatplt('S0325a', s0325a, 132, 1, 257, 26, 93, 93)

plt.subplot(231)
plt.axis('equal')

sns.scatterplot(data=s0325a, x="SV", y="SH", hue='Ratio1')
plt.plot([0,257], [0,0], '--')
plt.title('S0325a SH/SV')

plt.subplot(232)
plt.axis('equal')

sns.scatterplot(data=s0325a, x="SV", y="P", hue='Ratio2')
plt.plot([0,257], [0,132], '--')
plt.title('S0325a P/SV')

plt.subplot(233)
plt.axis('equal')

sns.scatterplot(data=s0325a, x="SH", y="P", hue='Ratio3')
plt.plot([0,0], [0,132], '--')
plt.title('S0325a P/SH')

plt.subplot(234)
plt.axis('equal')

sns.scatterplot(data=s0325a, x="SV", y="SH", hue='Sum')
plt.plot([0,257], [0,0], '--')
plt.title('S0325a SH/SV')

plt.subplot(235)
plt.axis('equal')

sns.scatterplot(data=s0325a, x="SV", y="P", hue='Sum')
plt.plot([0,257], [0,132], '--')
plt.title('S0325a P/SV')

plt.subplot(236)

sns.scatterplot(data=s0325a, x="SH", y="P", hue='Sum')
plt.plot([0,0], [0,132], '--')
plt.title('S0325a P/SH')

plt.axis('equal')
plt.show()
