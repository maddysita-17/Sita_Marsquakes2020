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

    ratio1_ls = [(np.arctan(obs_SHSV)-np.arctan(i))**2 for i in SHSV]
    ratio2_ls = [(np.arctan(obs_PSV)-np.arctan(i))**2 for i in PSV]
    ratio3_ls = [(np.arctan(obs_PSH)-np.arctan(i))**2 for i in PSH]

    df['Ratio1'] = ratio1_ls
    df['Ratio2'] = ratio2_ls
    df['Ratio3'] = ratio3_ls

    sum = df['Ratio1'] + df['Ratio2'] + df['Ratio3']
    df['Sum'] = sum

    return df


df173a = pd.read_csv('S0173a_[155][60][80].csv')
df235b = pd.read_csv('S0235b_[130][90][-90].csv')
df325a = pd.read_csv('S0325a_[35][60][-90].csv')

df325ab = pd.read_csv('S0325ab_[35][60][90].csv')
df173ab = pd.read_csv('S0173ab_[130][90][-100].csv')

#df183a = pd.read_csv('S0183a_[255][80][-100].csv')

#gudkova model observed amplitude ratios from seismograms ignoring the signs of the amplitude ratios
df173a = ratiosum(df173a, 0.80657748, 1.110367893, 1.376641327)
df235b = ratiosum(df235b, 1.172223922, 0.249311716, 0.212682672)
df325a = ratiosum(df325a, 0.399310873, 0.518759571, 1.299137105)
df325ab = ratiosum(df325ab, 1.009103448, 0.631448276, 0.625751777)
df173ab = ratiosum(df173ab, 1.354580708, 0.813378575, 0.600465199)
#df183a = ratiosum(df183a, 0.402906209, 0.911492734, 2.262295082)


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

# sum_df = df173a['Sum']
# min = sum_df[np.where(sum_df == 0)]
# print(min)


def contour_sum(df, az, plP, plS):
    Z = df.pivot_table(index='Azimuth', columns='Plunge', values='Sum').T.values
    X_unique = np.sort(df.Azimuth.unique())
    Y_unique = np.sort(df.Plunge.unique())
    X, Y = np.meshgrid(X_unique, Y_unique)
    levels = np.arange(0,2.5,0.1)
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

contour_sum(df173a, 273, 27.4, 26)
# contour_1r(df173a, 273, 27.4, 26)
# contour_2r(df173a, 273, 27.4, 26)
# contour_3r(df173a, 273, 27.4, 26)

contour_sum(df235b, 259, 28.6, 26.9)
#contour_1r(df235b, 259, 28.6, 26.9)
#contour_2r(df235b, 259, 28.6, 26.9)
#contour_3r(df235b, 259, 28.6, 26.9)

contour_sum(df325a, 300, 32.9, 30.4)
# contour_1r(df325a, 300, 32.9, 30.4)
# contour_2r(df325a, 300, 32.9, 30.4)
# contour_3r(df325a, 300, 32.9, 30.4)

contour_sum(df325ab, 315, 32.9, 30.4)
#contour_1r(df325ab, 315, 32.9, 30.4)
#contour_2r(df325ab, 315, 32.9, 30.4)
#contour_3r(df325ab, 315, 32.9, 30.4)

contour_sum(df173ab, 269, 27.4, 26)
# contour_1r(df173ab, 269, 27.4, 26)
# contour_2r(df173ab, 269, 27.4, 26)
# contour_3r(df173ab, 269, 27.4, 26)

# contour_sum(df183a, 247, 36, 33)
# contour_1r(df183a, 247, 36, 33)
# contour_2r(df183a, 247, 36, 33)
# contour_3r(df183a, 247, 36, 33)
