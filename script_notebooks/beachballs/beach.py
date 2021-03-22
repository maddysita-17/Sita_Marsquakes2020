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

def bbb(event, data):
    print(len(data))
    faults = data.drop_duplicates(subset = ["Strike","Dip","Rake"])
    print(len(faults))

    y_val = int(len(faults)/5)

    print('Event:', event)

    fig, ax = plt.subplots(figsize=(12,10))
    ax.set(xlim=(-1, 6), ylim=(0, y_val+1))

    n=1
    for index, rows in faults.iterrows():
        if n <= y_val:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(0,n), width=0.6)
            ax.add_collection(bball)
            plt.annotate(f, (-0.5, n+0.5), size=5)
            n += 1
        elif n <= y_val*2:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(1,n-y_val), width=0.6)
            ax.add_collection(bball)
            plt.annotate(f, (0.5, (n-y_val)+0.5), size=5)
            n += 1
        elif n <= y_val*3:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(2,n-(y_val*2)), width=0.6)
            ax.add_collection(bball)
            plt.annotate(f, (1.5, n-(y_val*2)+0.5), size=5)
            n += 1
        elif n <= y_val*4:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(3,n-(y_val*3)), width=0.6)
            ax.add_collection(bball)
            plt.annotate(f, (2.5, n-(y_val*3)+0.5), size=5)
            n += 1
        elif n <= y_val*5:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(4,n-(y_val*4)), width=0.6)
            ax.add_collection(bball)
            plt.annotate(f, (3.5, n-(y_val*4)+0.5), size=5)
            n += 1
        elif n <= y_val*6:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(5,n-(y_val*5)), width=0.6)
            ax.add_collection(bball)
            plt.annotate(f, (4.5, n-(y_val*5)+0.5), size=5)
            n += 1
    plt.show()

# S0173a = pd.read_csv('S0173a.csv')
# bbb('S0173a', S0173a)

S0173ab = pd.read_csv('S0173ab.csv')
bbb('S0173ab', S0173ab)

# S0235b = pd.read_csv('S0235b.csv')
# df1 = S0235b[S0235b['Sum']<= 0.16]
# df2 = df1[df1['Dip'] ]
# df3 = df2[df1['Strike'] <= ]
# bbb('S0235b', df1)
# bbb('S0235b', df2)

# S0235bi = pd.read_csv('S0235bi.csv')
# bbb('S0235bi', S0235bi)
#
# S0325a = pd.read_csv('S0325a.csv')
# bbb('S0325a', S0325a)
#
# S0325ab = pd.read_csv('S0325ab.csv')
# bbb('S0325ab', S0325ab)
