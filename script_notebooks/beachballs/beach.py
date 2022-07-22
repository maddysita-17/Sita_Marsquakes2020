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

def bbb(event, data, savefig=False):
    # faults = data.drop_duplicates(subset = ["Strike","Dip","Rake"])
    # print(len(faults))

    y_val = int(len(data)/5)

    print('Event:', event)

    fig, ax = plt.subplots(figsize=(10,7))
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

    if savefig == True:
        plt.savefig('S0' + str(event) +'_modelbeach.pdf')
    else:
        pass

    plt.show()

def bbb2(event,data):
    print('Event:', event)

    fig = plt.figure(figsize=(10,20))
    fig.suptitle(event)

    misfit = data['Misfit']
    print(misfit)

    max_mis = misfit.max()
    min_mis = misfit.min()

    norm_ls = []
    for value in misfit:
        newval = (max_mis-value)/(max_mis-min_mis)
        norm_ls.append(newval)

    misfit_norm = data.assign(Normalized = norm_ls)

    n = 0
    for index,rows in misfit_norm.iterrows():
        f = [rows.Strike, rows.Dip, rows.Rake]

        if rows.Extra[1:3] == '25':
            c = 'r'
        elif rows.Extra[1:3] == '35':
            c = 'g'
        elif rows.Extra[1:3] == '45':
            c = 'b'
        elif rows.Extra[1:3] == '55':
            c = 'orange'

        if n <= 20:
            fault = beachball(f, xy=(0,210*n), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 40:
            fault = beachball(f, xy=(210,210*(n-21)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 60:
            fault = beachball(f, xy=(420,210*(n-41)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 80:
            fault = beachball(f, xy=(630,210*(n-61)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 100:
            fault = beachball(f, xy=(840,210*(n-81)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 120:
            fault = beachball(f, xy=(1050,210*(n-101)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 140:
            fault = beachball(f, xy=(1260,210*(n-121)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 160:
            fault = beachball(f, xy=(1470,210*(n-141)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 180:
            fault = beachball(f, xy=(1680,210*(n-161)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 200:
            fault = beachball(f, xy=(1890,210*(n-181)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 220:
            fault = beachball(f, xy=(2100,210*(n-201)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 240:
            fault = beachball(f, xy=(2310,210*(n-221)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 260:
            fault = beachball(f, xy=(2520,210*(n-241)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 280:
            fault = beachball(f, xy=(2730,210*(n-261)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 300:
            fault = beachball(f, xy=(2940,210*(n-281)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 320:
            fault = beachball(f, xy=(3150,210*(n-301)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 340:
            fault = beachball(f, xy=(3360,210*(n-321)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 360:
            fault = beachball(f, xy=(3570,210*(n-341)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 380:
            fault = beachball(f, xy=(3780,210*(n-361)), fig=fig, facecolor=c, alpha = rows.Normalized)

        elif n <= 400:
            fault = beachball(f, xy=(3990,210*(n-381)), fig=fig, facecolor=c, alpha = rows.Normalized)

        n+=1

    # if savefig == True:
    #     plt.savefig('S0' + str(event) +'_allbeach.pdf')
    # else:
    #     pass

S0173a = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/top50_all173a.csv')
# g173a = S0173a.iloc[:200, :]
# print(g173a)
# t173a = S0173a.iloc[201:400, :]
# m173a = S0173a.iloc[401:600, :]
# c173a = S0173a.iloc[601:, :]
#
# bbb2('S0173a', g173a)
# bbb2('S0173a', t173a)


# S0173ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/check173ab.csv')
# S0235b = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/check235b.csv')
# S0325ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/check325ab.csv')

#-----sorted------
# for event in ['173a','173ab','235b','325ab']:
#     event_df = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/check'+ event + '.csv')
#     top_df = event_df.sort_values(by = ['Misfit'])
#     bbb(event, top_df, savefig=True)

# #-----sorted------
# for event in ['173a','173ab','235b','325ab']:
#     event_df = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/top50_all'+ event + '.csv')
#     g = event_df.iloc[:200, :]
#     print(g)
#     t = event_df.iloc[201:400, :]
#     m = event_df.iloc[401:600, :]
#     c = event_df.iloc[601:, :]
#
#     for modelz in [g,t,m,c]:
#         bbb2(event, modelz)

# #-----by model-----
# for event in ['173a','173ab','235b','325ab']:
#     event_df = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/check'+ event + '.csv')
#     bbb(event, event_df, savefig=True)

bbb2('S0173a',S0173a)
