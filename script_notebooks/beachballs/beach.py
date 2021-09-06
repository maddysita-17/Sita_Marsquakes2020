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

path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/beachballs/csvs/'
oldp = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/beachballs/csvs/oldtrials/'

def bbb(event, data):
    # faults = data.drop_duplicates(subset = ["Strike","Dip","Rake"])
    # print(len(faults))

    y_val = int(len(data)/5)

    print('Event:', event)

    fig, ax = plt.subplots(figsize=(7,7))
    ax.set(xlim=(-1, 6), ylim=(0, y_val+1))
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

# S0173a = pd.read_csv(path + 'resp_S0173a.csv')
# top173a = S0173a[:20]
# bbb('S0173a', top173a)
#
# S0173ab = pd.read_csv(path + 'resp_S0173ab.csv')
# top173ab = S0173ab[:20]
# bbb('S0173ab', top173ab)
#
# S0235b = pd.read_csv(path + 'resp_S0235b.csv')
# top235b = S0235b[:40]
# bbb('S0235b', top235b)
#
# # S0235bi = pd.read_csv('S0235bi.csv')
# # bbb('S0235bi', S0235bi)
#
# # S0325a = pd.read_csv(path + 'sort_325a.csv')
# # bbb('S0325a', S0325a)
#
# S0325ab = pd.read_csv(path + 'resp_S0325ab.csv')
# top325ab = S0325ab[:20]
# bbb('S0325ab', top325ab)

S0173a = pd.read_csv(path + 'check173a.csv')
S0173ab = pd.read_csv(path + 'check173ab.csv')
S0235b = pd.read_csv(path + 'check235b.csv')
S0325ab = pd.read_csv(path + 'check325ab.csv')

bbb('S0173a', S0173a)
bbb('S0173ab', S0173ab)
bbb('S0235b', S0235b)
bbb('S0325ab', S0325ab)
