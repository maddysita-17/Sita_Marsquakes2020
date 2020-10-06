import matplotlib.pyplot as plt
import numpy as np
#from obspy.taup import TauPyModel
import pandas as pd

df = pd.read_csv('173a-35.csv')
print(df.head())

#attempted a different "ratio" calcultion to avoid /0 or getting large numbers
SHSV = (df['SH']-df['SV']) / (df['SH']+df['SV'])
PSV = (df['SV']-df['P']) / (df['P']+df['SV'])
PSH = (df['P']-df['SH']) / (df['P']+df['SH'])

ratio1_ls = [(-0.012135584-i)**2 for i in SHSV]
ratio2_ls = [(0.010290875-i)**2 for i in PSV]
ratio3_ls = [(0.001844709-i)**2 for i in PSH]

df['Ratio1'] = ratio1_ls
df['Ratio2'] = ratio2_ls
df['Ratio3'] = ratio3_ls

sum = df['Ratio1'] + df['Ratio2'] + df['Ratio3']
df['Sum'] = sum
print(df)
df.to_csv('173a-sum.csv', index = False)


# def ratiosum(df, obs_SHSV, obs_PSV, obs_PSH):
#     SHSV = (df['SH']-df['SV']) / (df['SH']+df['SV'])
#     PSV = df['P/SV']
#     PSH = df['P/SH']
#
#     ratio1_ls = [(obs_SHSV-i)**2 for i in SHSV]
#     ratio2_ls = [(obs_PSV-i)**2 for i in PSV]
#     ratio3_ls = [(obs_PSH-i)**2 for i in PSH]
#
#     df['Ratio1'] = ratio1_ls
#     df['Ratio2'] = ratio2_ls
#     df['Ratio3'] = ratio3_ls
#
#     sum = df['Ratio1'] + df['Ratio2'] + df['Ratio3']
#     df['Sum'] = sum
#
#     return df
#
#
# df173a = pd.read_csv('173a-35.csv')
# df235b = pd.read_csv('235b-35.csv')
# df325a = pd.read_csv('325a-35.csv')
#
# #gudkova model observed amplitude ratios from seismograms
# df173a = ratiosum(df173a, -0.626168224, -0.831775701, 1.328358209)
# df173a.to_csv('173a-sum.csv', index = False)
# print(df173a)
# df235b = ratiosum(df235b, 1.035519126, 0.43715847, 0.422163588)
# df325a = ratiosum(df325a, -0.4, 0.584615385, -1.461538462)


Z = df.pivot_table(index='Azimuth', columns='Plunge', values='Sum').T.values
X_unique = np.sort(df.Azimuth.unique())
Y_unique = np.sort(df.Plunge.unique())
X, Y = np.meshgrid(X_unique, Y_unique)

levels = np.arange(0,8,0.25)

cp = plt.contourf(X, Y, Z)
plt.colorbar(cp)
plt.show()


# def contour_plot(df):
#     Z = df.pivot_table(index='Azimuth', columns='Plunge', values='Sum').T.values
#     X_unique = np.sort(df.Azimuth.unique())
#     Y_unique = np.sort(df.Plunge.unique())
#     X, Y = np.meshgrid(X_unique, Y_unique)
#     cp = plt.contourf(X, Y, Z)
#     plt.colorbar(cp)
#     plt.show()
#
# contour_plot(df173a)
# contour_plot(df235b)
# contour_plot(df325a)
