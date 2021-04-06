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


s0173a = pd.read_csv('S0173a.csv')
faults = s0173a[['P', 'SV', 'SH']]
faults = faults.copy()

# sns.set_style("darkgrid", {"axes.facecolor": ".9"})
# grid = sns.pairplot(faults, kind='reg')
# plt.show()


def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

def getb(m,x):
    y = m*x
    b = y - m*x
    return b


sns.regplot(x="SV", y="SH", data=faults)
b = getb(-0.818181818, -1)
abline(-0.818181818, b)
plt.title('S0173a SH/SV')
plt.show()

# sns.regplot(x="SV", y="P", data=faults)
# b = getb(-1.110367893, 1)
# abline(-1.110367893, b)
# plt.title('S0173a P/SV')
# plt.show()
#
# sns.regplot(x="SH", y="P", data=faults)
# b = getb(1.376641327, -1)
# abline(1.376641327, b)
# plt.title('S0173aP/SH')
# plt.show()

s0325a = pd.read_csv('S0325a.csv')
faults3 = s0325a[['P', 'SV', 'SH']]
faults3 = faults3.copy()

sns.regplot(x="SV", y="SH", data=faults3)
b = getb(-0.399310873, 2)
abline(-0.399310873, b)
plt.title('S0325a SH/SV')
plt.show()
