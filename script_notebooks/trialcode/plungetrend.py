import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel
import pandas as pd

def eventbuild(event, dist):
    path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/faultdata/' + event + '/'
    mtimes = mars.get_travel_times(source_depth_in_km = depth, distance_in_degree = dist, phase_list=["P", "S"])

    #incident angle at the station
    Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
    Pvel = radius*np.sin(np.radians(Pa))/Pp
    try:
        Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
    except:
        Sp = 0; Sa = 0
        print('Within S-wave shadow zone')
    Svel = radius*np.sin(np.radians(Sa))/Sp

    return path, Pp, Sp, Pa, Sa

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

def contour_sum(df, az, plP, plS):
    pl = df['Plunge']
    df['plungeP'] = pl - plP
    df['plungeS'] = pl - plS


    fig = plt.figure(figsize=(10,7))

    Z = df.pivot_table(index='Azimuth', columns='plungeP', values='Sum').T.values
    X_unique = np.sort(df.Azimuth.unique())
    Y_unique = np.sort(df.plungeP.unique())
    X, Y = np.meshgrid(X_unique, Y_unique)
    levels = np.arange(0,3,0.1)
    cp = plt.contourf(X, Y, Z, levels=levels, alpha = 0.25, cmap='nipy_spectral')

    Z = df.pivot_table(index='Azimuth', columns='plungeS', values='Sum').T.values
    X_unique = np.sort(df.Azimuth.unique())
    Y_unique = np.sort(df.plungeS.unique())
    X, Y = np.meshgrid(X_unique, Y_unique)
    levels = np.arange(0,3,0.1)
    cp = plt.contourf(X, Y, Z, levels=levels, alpha = 0.25, cmap='nipy_spectral')

    plt.colorbar(cp)
    plt.plot(az, plP, 'bo', label='P')
    plt.annotate('P', (az, plP+1))
    plt.plot(az, plS, 'ro', label='S')
    plt.annotate('S', (az, plS+1))
    plt.title('Sum')
    plt.show()

path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/ray_paths/plunge-trend-maps/csvs/'
df173a = pd.read_csv(path + 'S0173a_[165][75][80].csv')

df173a = ratiosum(df173a, 0.80657748, 1.110367893, 1.376641327)
contour_sum(df173a, 273, 27.4, 26)
