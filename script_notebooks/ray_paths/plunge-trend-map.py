import matplotlib.pyplot as plt
import numpy as np
#from obspy.taup import TauPyModel
import pandas as pd

path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/ray_paths/plunge-trend-maps'

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


df173a = pd.read_csv('S0173a_[160][60][80].csv')
df235b = pd.read_csv('S0235b_[130][90][-90].csv')
df325a = pd.read_csv('S0325a_[35][60][-90].csv')

df325ab = pd.read_csv('S0325ab_[35][60][90].csv')
df173ab = pd.read_csv('S0173ab_[130][90][-100].csv')

# df325ab = pd.read_csv('S0325ab_[45][90][-90].csv')
# df173ab = pd.read_csv('S0173ab_[160][60][80].csv')

df183a = pd.read_csv('S0183a_[255][80][-100].csv')

# #gudkova model observed amplitude ratios from seismograms
# df173a = ratiosum(df173a, -0.626168224, -0.831775701, 1.328358209)
# df235b = ratiosum(df235b, 1.035519126, 0.43715847, 0.422163588)
# df325a = ratiosum(df325a, -0.4, 0.584615385, -1.461538462)
# df325ab = ratiosum(df325ab, -1.005586592, 0.656424581, -0.652777778)
# df183a = ratiosum(df183a, 0.402906209, 0.911492734, 2.262295082)

#gudkova model observed amplitude ratios from seismograms ignoring the signs of the amplitude ratios
df173a = ratiosum(df173a, 0.626168224, 0.831775701, 1.328358209)
df235b = ratiosum(df235b, 1.035519126, 0.43715847, 0.422163588)
df325a = ratiosum(df325a, 0.4, 0.584615385, 1.461538462)
df325ab = ratiosum(df325ab, 1.005586592, 0.656424581, 0.652777778)
df173ab = ratiosum(df173ab, 1.696969697, 0.742424242, 0.4375)
df183a = ratiosum(df183a, 0.402906209, 0.911492734, 2.262295082)


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


def contour_sum(df, az, plP, plS):
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
    plt.title('Sum')
    plt.show()

def contour_1r(df, az, plP, plS):
    Z = df.pivot_table(index='Azimuth', columns='Plunge', values='Ratio1').T.values
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
    plt.title('Ratio1 - SH/SV')
    plt.show()

def contour_2r(df, az, plP, plS):
    Z = df.pivot_table(index='Azimuth', columns='Plunge', values='Ratio2').T.values
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
    plt.title('Ratio2 - P/SV')
    plt.show()

def contour_3r(df, az, plP, plS):
    Z = df.pivot_table(index='Azimuth', columns='Plunge', values='Ratio3').T.values
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
    plt.title('Ratio3 - P/SH')
    plt.show()

#contour_sum(df173a, 258, 28, 27)
#contour_1r(df173a, 258, 28, 27)
#contour_2r(df173a, 258, 28, 27)
#contour_3r(df173a, 258, 28, 27)

#contour_sum(df235b, 274, 27, 26)
#contour_1r(df235b, 274, 27, 26)
#contour_2r(df235b, 274, 27, 26)
#contour_3r(df235b, 274, 27, 26)

contour_sum(df325a, 300, 33, 31)
# contour_1r(df325a, 300, 33, 31)
# contour_2r(df325a, 300, 33, 31)
# contour_3r(df325a, 300, 33, 31)

#contour_sum(df325ab, 290, 30, 28)
#contour_1r(df325ab, 290, 30, 28)
#contour_2r(df325ab, 290, 30, 28)
#contour_3r(df325ab, 290, 30, 28)

# contour_sum(df173ab, 268, 27, 26)
# contour_1r(df173ab, 268, 27, 26)
# contour_2r(df173ab, 268, 27, 26)
# contour_3r(df173ab, 268, 27, 26)

# contour_sum(df183a, 247, 36, 33)
# contour_1r(df183a, 247, 36, 33)
# contour_2r(df183a, 247, 36, 33)
# contour_3r(df183a, 247, 36, 33)
