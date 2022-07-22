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

def bbb(event, data, chi):
    data0 = data[data['Sum'] < chi]
    X = data0['Strike']
    Y = data0['Dip']
    Z = data0['Rake']
    C = data0['Sum']

    d = {'Strike': X, 'Dip': Y, 'Rake': Y}
    df = pd.DataFrame.from_dict(d)
    faults = df.drop_duplicates(subset = ["Strike","Dip","Rake"])
    print(len(faults))

    y_val = int(len(X)/5)

    print('Event:', event)

    fig, ax = plt.subplots(figsize=(12,10))
    ax.set(xlim=(-1, 6), ylim=(0, y_val+1))

    n=1
    for index, rows in faults.iterrows():
        if n <= y_val:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(0,n), width=0.6, nofill=True)
            ax.add_collection(bball)
            plt.annotate(f, (-0.5, n+0.5), size=5)
            n += 1
        elif n <= y_val*2:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(1,n-y_val), width=0.6, nofill=True)
            ax.add_collection(bball)
            plt.annotate(f, (0.5, (n-y_val)+0.5), size=5)
            n += 1
        elif n <= y_val*3:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(2,n-(y_val*2)), width=0.6, nofill=True)
            ax.add_collection(bball)
            plt.annotate(f, (1.5, n-(y_val*2)+0.5), size=5)
            n += 1
        elif n <= y_val*4:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(3,n-(y_val*3)), width=0.6, nofill=True)
            ax.add_collection(bball)
            plt.annotate(f, (2.5, n-(y_val*3)+0.5), size=5)
            n += 1
        elif n <= y_val*5:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(4,n-(y_val*4)), width=0.6, nofill=True)
            ax.add_collection(bball)
            plt.annotate(f, (3.5, n-(y_val*4)+0.5), size=5)
            n += 1
        elif n <= y_val*6:
            f = [rows.Strike, rows.Dip, rows.Rake]
            bball = beach(f, xy=(5,n-(y_val*5)), width=0.6, nofill=True)
            ax.add_collection(bball)
            plt.annotate(f, (4.5, n-(y_val*5)+0.5), size=5)
            n += 1
    plt.show()

def plot3d(X, Y, Z, C):
    markercolor = C

    fig1 = go.Scatter3d(x=X,
                        y=Y,
                        z=Z,
                        marker=dict(color=markercolor,
                                    opacity=1,
                                    reversescale=True,
                                    colorscale='Greens',
                                    size=2, showscale=True),
                        line=dict (width=0.2),
                        mode='markers')


    mylayout = go.Layout(scene=dict(xaxis=dict( title="Strike"),
                                    yaxis=dict( title="Dip"),
                                    zaxis=dict(title="Rake")),)

    #Plot and save html
    plotly.offline.plot({"data": [fig1],
                         "layout": mylayout},
                         auto_open=True,
                         filename=("4DPlot"))
    return

def histo(pf):
    st = pf['Strike']
    dp = pf['Dip']
    rk = pf['Rake']
    chi = pf['Sum']

    fig = make_subplots(rows=4, cols=1, subplot_titles=("Strike", "Dip", "Rake"))
    strike = go.Histogram(x=st, nbinsx=18)

    dip = go.Histogram(x=dp, nbinsx=9)

    rake = go.Histogram(x=rk, nbinsx=9)

    misfit = go.Histogram(x=chi)

    fig.append_trace(strike, 1, 1)
    fig.append_trace(dip, 2, 1)
    fig.append_trace(rake, 3, 1)
    fig.append_trace(misfit, 4, 1)
    fig.show()

def beach_plot(event, data, chi, bball=True, plot=False, hist=False):
    C = data['Sum']

    data0 = data[data['Sum'] < chi]
    X = data0['Strike']
    Y = data0['Dip']
    Z = data0['Rake']
    C = data0['Sum']

    if bball == True:
        bbb(event, data, chi)
    else:
        faults = zip(X, Y, Z)

    d = {'Strike': X, 'Dip': Y, 'Rake': Z, 'Sum': C}
    df = pd.DataFrame.from_dict(d)
    df.to_csv(event +'_possible_faults.csv', index=False)
    return df

    if plot == True:
        plot3d(X, Y, Z, C)
    else:
        pass

    if hist == True:
        histo(df)
    else:
        pass


# S0173a = pd.read_csv('S0173a_NewGudkova_35.csv')
# pf = beach_plot('S0173a', S0173a, 0.00414)
# pf['depth'] = 35
#
# bS0173a = pd.read_csv('S0173a_NewGudkova_25.csv')
# bpf = beach_plot('S0173a', bS0173a, 0.00414, bball=False)
# bpf['depth'] = 25
#
# cS0173a = pd.read_csv('S0173a_NewGudkova_45.csv')
# cpf = beach_plot('S0173a', cS0173a, 0.00414, bball=False)
# cpf['depth'] = 45
#
# faults173a = pd.concat([pf, bpf, cpf])
#
# faults173a = faults173a.drop(columns=['Sum'])
# sns.pairplot(faults173a, diag_kind = 'hist', hue='depth', palette="Set2", plot_kws=dict(edgecolor="white", linewidth=0.1))
# plt.title('S0173a')
# plt.show()

# S0173ab = pd.read_csv('S0173ab_NewGudkova_35.csv')
# pf = beach_plot('S0173ab', S0173ab, 0.03)
# faults173ab = pf.drop(columns=['Sum'])
# sns.pairplot(faults173ab)
# plt.title('S0173ab')
# plt.show()
#
# S0235b = pd.read_csv('S0235b_NewGudkova_35.csv')
# pf = beach_plot('S0235b', S0235b, 0.000170)
# faults235b = pf.drop(columns=['Sum'])
# sns.pairplot(faults235b)
# plt.title('S0235b')
# plt.show()

# S0325a = pd.read_csv('S0325a_NewGudkova_35.csv')
# pf = beach_plot('S0325a', S0325a, 0.0207, bball=False)
# faults325a = pf.drop(columns=['Sum'])
# sns.set_style("darkgrid", {"axes.facecolor": ".9", "lines.color":"indigo"})
# sns.pairplot(faults325a, plot_kws=dict(marker="+", linewidth=1),
#     diag_kws=dict(fill=True))
# plt.show()

# S0325ab = pd.read_csv('S0325ab_NewGudkova_35.csv')
# pf = beach_plot('S0325ab', S0325ab, 0.002762, bball=False)
# faults325ab = pf.drop(columns=['Sum'])
# sns.pairplot(faults325ab)
# plt.title('S0325ab')
# plt.show()

fifty = []; ten = []
# S0325ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/event-by-event/S0325ab/csvs/ng325ab_sqrt3.csv')
# S0325ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/event-by-event/S0235b/csvs/ng235b.csv')
# S0325ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/event-by-event/S0173/csvs/ng173a.csv')
# S0325ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/event-by-event/S0173/csvs/ng173ab.csv')
# S0325ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/event-by-event/VANUATU/vanuatu_faults.csv')

S0325ab = pd.read_csv('/Users/maddysita/Desktop/CIERA_REU/event-by-event/S0235b/csvs/ng235b_enhanP.csv')

faults325a = S0325ab.drop(columns=['Index','Extra'])
# print(faults325a[0:5])
rev = faults325a.iloc[::-1]

# sns.set_style("darkgrid", {"axes.facecolor": "dbe9ef", "lines.color": 'dbe9ef'})

sns.set_style("darkgrid", {"axes.facecolor": "0.9", "lines.color": '1'})
g = sns.PairGrid(rev, hue=None, corner=True)

branges = iter([(0,180),(0,90),(-180,180),(0,0.025423677579)])
# branges = iter([(0,180),(0,90),(-180,180)])
def myhist(*args, **kwargs):
    plt.hist(*args, range=next(branges), **kwargs)

g.map_diag(myhist, color="0.4", bins=12)
# g.map_diag(myhist, color="#4B6D97", bins=12)

g.map_offdiag(sns.scatterplot,linewidth=0.1, color='0.4')
# g.map_offdiag(sns.scatterplot,linewidth=0.1, color='#4B6D97')


g.axes[2,0].set_xlim(-2,182); g.axes[2,1].set_xlim(-2,92); g.axes[2,2].set_xlim(-182,182)

# for y in range(0,2):
#     g.axes[3,y].set_ylim(-0.003,0.0016)

#S0235b : 78,28,86,0.000392270514954
def single_pt(st,dp,rk,mf):
    #column 0 - strike on x axis
    g.axes[1,0].scatter(st,dp, linewidth=0.1)
    g.axes[2,0].scatter(st,rk, linewidth=0.1)
    g.axes[3,0].scatter(st,mf, linewidth=0.1)
    #column 1 - dip on x axis
    g.axes[2,1].scatter(dp,rk, linewidth=0.1)
    g.axes[3,1].scatter(dp,mf, linewidth=0.1)
    #column 2 - rake on x axis
    g.axes[3,2].scatter(rk,mf, linewidth=0.1)

# def single_pt_p(st,dp,rk,mf,c):
#     #column 0 - strike on x axis
#     g.axes[1,0].scatter(st,dp, linewidth=0.1,color=c); g.axes[2,0].scatter(st,rk, linewidth=0.1,color=c); #g.axes[3,0].scatter(st,mf, linewidth=0.1,color=c)
#     #column 1 - dip on x axis
#     g.axes[2,1].scatter(dp,rk, linewidth=0.1,color=c); #g.axes[3,1].scatter(dp,mf, linewidth=0.1,color=c)
#     #column 2 - rake on x axis
#     #g.axes[3,2].scatter(rk,mf, linewidth=0.1,color=c)
#
#     # row 0 - strike on y axis
#     g.axes[0,1].scatter(dp,st, linewidth=0.1,color=c); g.axes[0,2].scatter(rk,st, linewidth=0.1,color=c); #g.axes[0,3].scatter(mf,st, linewidth=0.1,color=c)
#     # row 1 - dip on y axis
#     g.axes[1,2].scatter(rk,dp, linewidth=0.1,color=c); #g.axes[1,3].scatter(mf,dp, linewidth=0.1,color=c)
#     # row 2 - rake on y axis
#     # g.axes[2,3].scatter(mf,rk, linewidth=0.1,color=c)



single_pt(78,28,86,0.00105289202079) #s0235b
# single_pt(70,76,-116,0.000955321846323)  #s0325ab
# single_pt(46,38,148,0.00176499374116)   #s0173a
# single_pt(46,50,-92,0.000638672817604) #VANUATU

# single_pt_p(78,28,86,0.000392270514954, '#B8D96E') #s0235b
# single_pt_p(70,76,-116,0.00169152994718, '#B8D96E') #s0325ab


g.add_legend()
plt.show()
g.savefig('/Users/maddysita/Desktop/CIERA_REU/script_notebooks/paper_figures/3dplots/235b_enhanP.png')
