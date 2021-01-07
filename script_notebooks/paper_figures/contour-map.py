import matplotlib.pyplot as plt
import numpy as np
#from obspy.taup import TauPyModel
import pandas as pd

# df = pd.read_csv('173a-35.csv')
# print(df.head())
#
# SHSV = df['SH/SV']
# PSV = df['P/SV']
# PSH = df['P/SH']
#
# ratio1_ls = [(-0.626168224-i)**2 for i in SHSV]
# ratio2_ls = [(-0.831775701-i)**2 for i in PSV]
# ratio3_ls = [(1.328358209-i)**2 for i in PSH]
#
# df['Ratio1'] = ratio1_ls
# df['Ratio2'] = ratio2_ls
# df['Ratio3'] = ratio3_ls
#
# sum = df['Ratio1'] + df['Ratio2'] + df['Ratio3']
# df['Sum'] = sum
# print(df)

def ratiosum(df, obs_SHSV, obs_PSV, obs_PSH):
    SHSV = df['SH/SV']
    PSV = df['P/SV']
    PSH = df['P/SH']

    ratio1_ls = [(obs_SHSV-i)**2 for i in SHSV]
    ratio2_ls = [(obs_PSV-i)**2 for i in PSV]
    ratio3_ls = [(obs_PSH-i)**2 for i in PSH]

    df['Ratio1'] = ratio1_ls
    df['Ratio2'] = ratio2_ls
    df['Ratio3'] = ratio3_ls

    sum = df['Ratio1'] + df['Ratio2'] + df['Ratio3']
    df['Sum'] = sum

    return df


df173a = pd.read_csv('try5_173a-35.csv')
df235b = pd.read_csv('try_235b-35.csv')
df325a = pd.read_csv('try2_325a-35.csv')
df325ab = pd.read_csv('325ab-35.csv')

#gudkova model observed amplitude ratios from seismograms
df173a = ratiosum(df173a, -0.626168224, -0.831775701, 1.328358209)
#df173a.to_csv('new_173a-sum.csv', index = False)
#print(df173a)
df235b = ratiosum(df235b, 1.035519126, 0.43715847, 0.422163588)
df325a = ratiosum(df325a, -0.4, 0.584615385, -1.461538462)
df325ab = ratiosum(df325ab, -1.005586592, 0.656424581, -0.652777778)


# Z = df173a.pivot_table(index='Azimuth', columns='Plunge', values='Sum').T.values
# X_unique = np.sort(df173a.Azimuth.unique())
# Y_unique = np.sort(df173a.Plunge.unique())
# X, Y = np.meshgrid(X_unique, Y_unique)
#
# levels = np.arange(0,8,0.25)
#
# cp = plt.contourf(X, Y, Z)
# plt.colorbar(cp)
# plt.show()


def contour_plot(event, df, az, plP, plS):
    Z = df.pivot_table(index='Azimuth', columns='Plunge', values='Sum').T.values
    X_unique = np.sort(df.Azimuth.unique())
    Y_unique = np.sort(df.Plunge.unique())
    X, Y = np.meshgrid(X_unique, Y_unique)
    levels = np.arange(0,8,0.25)
    cp = plt.contourf(X, Y, Z, levels=levels, cmap='nipy_spectral')
    plt.colorbar(cp)
    plt.plot(az, plP, 'bo', label='P')
    plt.annotate('P', (az, plP+1))
    plt.plot(az, plS, 'ro', label='S')
    plt.annotate('S', (az, plS+1))
    plt.xlabel('Azimuth')
    plt.ylabel('Plunge')
    plt.title(event)
    plt.savefig('contour_' + event + '.png')
    plt.show()

contour_plot('S0173a', df173a, 258, 28, 27)
contour_plot('S0235b', df235b, 274, 27, 26)
contour_plot('S0325a', df325a, 300, 33, 31)
contour_plot('S0325ab', df325ab, 290, 30, 28)
