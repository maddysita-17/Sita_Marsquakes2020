import pandas as pd
import numpy as np
import math
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from obspy.imaging.beachball import beachball
from obspy.imaging.beachball import beach
from obspy.imaging.source import plot_radiation_pattern
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

def bbb(event, data):
    # faults = data.drop_duplicates(subset = ["Strike","Dip","Rake"])
    # print(len(faults))

    y_val = int(len(data)/5)

    print('Event:', event)

    fig, ax = plt.subplots(figsize=(7,7))
    ax.set(xlim=(-1, y_val+2), ylim=(0, y_val+1))
    ax.set_title(event)

    n=0
    for index, rows in data.iterrows():
        if n <= y_val:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(0,n+0.3), width=0.6)
            ax.add_collection(bball)
            fa = [rows.Strike, rows.Dip, rows.Rake, truncate(rows['Misfit'], 4)]
            plt.annotate(fa, (-0.5, n+0.6), size=4)
            n += 1
        elif n <= y_val*2:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(1,n-y_val), width=0.6)
            ax.add_collection(bball)
            fa = [rows.Strike, rows.Dip, rows.Rake, truncate(rows['Misfit'], 4)]
            plt.annotate(fa, (0.5, (n-y_val)+0.5), size=4)
            n += 1
        elif n <= y_val*3:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(2,n-(y_val*2)), width=0.6)
            ax.add_collection(bball)
            fa = [rows.Strike, rows.Dip, rows.Rake, truncate(rows['Misfit'], 4)]
            plt.annotate(fa, (1.5, n-(y_val*2)+0.5), size=4)
            n += 1
        elif n <= y_val*4:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(3,n-(y_val*3)), width=0.6)
            ax.add_collection(bball)
            fa = [rows.Strike, rows.Dip, rows.Rake, truncate(rows['Misfit'], 4)]
            plt.annotate(fa, (2.5, n-(y_val*3)+0.5), size=4)
            n += 1
        elif n <= y_val*5:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(4,n-(y_val*4)), width=0.6)
            ax.add_collection(bball)
            fa = [rows.Strike, rows.Dip, rows.Rake, truncate(rows['Misfit'], 4)]
            plt.annotate(fa, (3.5, n-(y_val*4)+0.5), size=4)
            n += 1
        elif n <= y_val*6:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(5,n-(y_val*5)), width=0.6)
            ax.add_collection(bball)
            fa = [rows.Strike, rows.Dip, rows.Rake, truncate(rows['Misfit'], 4)]
            plt.annotate(fa, (4.5, n-(y_val*5)+0.5), size=4)
            n += 1
    plt.show()

uga = pd.read_csv('sort_uganda.csv')
bbb('Uganda', uga)
#
# nm = pd.read_csv('newmexico.csv')
# bbb('NewMexico', nm)

# id = pd.read_csv('idaho_g.csv')
# bbb('Idaho - Greenland', id)
# # id = pd.read_csv('idaho_g2.csv')
# # bbb('Idaho - Greenland2', id)
#
# id = pd.read_csv('idaho.csv')
# bbb('Idaho - CostaRica', id)

# bk = pd.read_csv('backwards.csv')
# bbb('Backwards', bk)

# mf_3d = pd.read_csv('misfit3D.csv')
# bbb('Misfit 3D', mf_3d)
#
# mf_s = pd.read_csv('misfitsum.csv')
# bbb('Misfit Sum', mf_s)
#
# mf_d = pd.read_csv('misfitd.csv')
# bbb('Misfit Dist', mf_d)
