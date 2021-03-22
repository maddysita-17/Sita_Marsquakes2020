import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel
import pandas as pd
from pathlib import Path

def getmt(fault):
    rad = 180./np.pi; m0 = 1
    st,dp,rk = fault

    st = st/rad; dp = dp/rad; rk = rk/rad
    sd = np.sin(dp); cd = np.cos(dp)
    sd2 = np.sin(2*dp); cd2 = np.cos(2*dp)
    ss = np.sin(st); cs = np.cos(st)
    ss2 = np.sin(2*st); cs2 = np.cos(2*st)
    sr = np.sin(rk); cr = np.cos(rk)

    # formulas from Aki & Richards Box 4.4
    # mt(1-6): Mrr, Mtt, Mff, Mrt, Mrf, Mtf
    # mt(1-6): Mzz, Mxx, Myy, Mzx, -Mzy, -Mxy
    mt = [ sr*sd2, -1.*sd*cr*ss2 - sd2*sr*ss*ss, sd*cr*ss2 - sd2*sr*cs*cs, -1.*cd*cr*cs - cd2*sr*ss, cd*cr*ss - cd2*sr*cs, -1.*sd*cr*cs2 - 0.5*sd2*sr*ss2]
    return mt


def Rpattern(fault,azimuth,incidence_angles):
    """
    Calculate predicted amplitudes of P, SV, and SH waves.
    IN: fault = [strike, dip, rake]
             = faulting mechanism, described by a list of strike, dip, and rake
             (note, strike is measure clockwise from N, dip is measured positive downwards
             (between 0 and 90) w.r.t. a horizontal that is 90 degrees clockwise from strike,
             and rake is measured positive upwards (counterclockwise)
        azimuth: azimuth with which ray path leaves source (clockwise from N)
        incidence_angles = [i, j]
              i = angle between P ray path & vertical in the source model layer
              j = angle between S ray path & vertical in the source model layer
    OUT: Amplitudes for P, SV, and SH waves
    P as measured on L (~Z) component, SV measured on Q (~R) component, and SH measured on T component.
    All input is in degrees.
    """

    strike,dip,rake = fault
    a = azimuth; rela = strike - azimuth
    sinlam = np.sin(np.radians(rake))
    coslam = np.cos(np.radians(rake))
    sind = np.sin(np.radians(dip))
    cosd = np.cos(np.radians(dip))
    cos2d = np.cos(np.radians(2*dip))
    sinrela = np.sin(np.radians(rela))
    cosrela = np.cos(np.radians(rela))
    sin2rela = np.sin(np.radians(2*rela))
    cos2rela = np.cos(np.radians(2*rela))

    sR = sinlam*sind*cosd
    qR = sinlam*cos2d*sinrela + coslam*cosd*cosrela
    pR = coslam*sind*sin2rela - sinlam*sind*cosd*cos2rela
    pL = sinlam*sind*cosd*sin2rela + coslam*sind*cos2rela
    qL = -coslam*cosd*sinrela + sinlam*cos2d*cosrela

    iP = np.radians(incidence_angles[0])
    jS = np.radians(incidence_angles[1])

    AP = np.abs(sR*(3*np.cos(iP)**2 - 1) - qR*np.sin(2*iP) - pR*np.sin(iP)**2)
    ASV = np.abs(1.5*sR*np.sin(2*jS) + qR*np.cos(2*jS) + 0.5*pR*np.sin(2*jS))
    ASH = np.abs(-qL*np.cos(jS) - pL*np.sin(jS))

    return AP,ASV,ASH


model_ls = ['NewGudkova']
path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/trialcode/'


#---get fault function--
def getfault(az, st, dp, rk):
    """
    INPUT: az = azimuth in degrees from lat-long-az.py
    st, dp, rk = strike, dip, rake from 3 separte lists
    OUTPUT: df = data frame containing model outputs for P,SV,SH amplitudes at each depth input
    ip, ij = exit angles in degrees for each model computed in the Rpattern function
    """
    strike_ls = []; dip_ls = []; rake_ls = []
    P_ls = []; SH_ls =[]; SV_ls=[]
    exit_angP = []; exit_angS = []
    for strike in st:
        for dip in dp:
            for rake in rk:
                strike_ls.append(strike)
                dip_ls.append(dip)
                rake_ls.append(rake)

                fault = [strike, dip, rake]
                mt = getmt(fault)
                azimuth = az

                iP = np.degrees(np.arcsin(Pvelz*Pp/radius))
                jS = np.degrees(np.arcsin(Svelz*Sp/radius))
                exit_angP.append(iP); exit_angS.append(jS)

                P,SV,SH = Rpattern(fault,azimuth,[iP,jS])
                scalefactor = (Pvelz/Svelz)**3
                SV,SH = SV*scalefactor, SH*scalefactor
                P_ls.append(P); SH_ls.append(SH); SV_ls.append(SV)

    ratio1 = [i / j for i, j in zip(SH_ls, SV_ls)]
    ratio2 = [i / j for i, j in zip(P_ls, SV_ls)]
    ratio3 = [i / j for i, j in zip(P_ls, SH_ls)]

    data = {'Strike': strike_ls,
            'Dip': dip_ls,
            'Rake': rake_ls,
            'P': P_ls,
            'SH': SH_ls,
            'SV': SV_ls,
            'SH/SV': ratio1,
            'P/SV': ratio2,
            'P/SH': ratio3,
            'Plunge P': exit_angP,
            'Plunge S': exit_angS}

    df = pd.DataFrame(data, columns = ['Strike', 'Dip', 'Rake', 'P', 'SV', 'SH', 'SH/SV', 'P/SV', 'P/SH', 'Plunge P', 'Plunge S'])
    return df, iP, jS

def eventbuild(event, dist):
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

def autofault(df, obs_SHSV, obs_PSV, obs_PSH):
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


    faults = df

    return faults


strike = [*range(0, 182, 2)]
dip = [*range(0,92,2)]
rake = [*range(-100,105,5)]

mod ='NewGudkova'
mars = TauPyModel(model=mod)
radius = 3389.5

Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []

Gudkova_depth = [45, 35, 25]

for depth in Gudkova_depth:
    if depth <= 50 and depth > 42:
        Pvelz = 7.12500; Svelz = 4.00300    #rounded
    elif depth <= 42 and depth > 21:
        Pvelz = 7.13900; Svelz = 4.01900
    elif depth <=21 and depth > 16:
        Pvelz = 7.14300; Svelz = 4.02300
    elif depth <= 16 and depth > 10:
        Pvelz = 7.15000; Svelz = 4.03000    #rounded
    elif depth <= 10 and depth > 9:
        Pvelz = 6.48200; Svelz = 3.59900
    elif depth <= 9 and depth > 8.586:        #barrier at 8.081km
        Pvelz = 5.94800; Svelz = 3.25500
    elif depth <= 8.586 and depth > 8.081:
        Pvelz = 5.68100; Svelz = 3.08300
#missing a lot of velocities in the middle but can add if we want to run more depths
    elif depth <= 5.051 and depth > 4.545:
        Pvelz = 3.13700; Svelz = 1.60200
    else:
        print("There is no computed velcocity at this depth")

    #----S0173a----

    path, Pp, Sp, Pa, Sa= eventbuild('173a', 28.4)

    data173a, Pe, Se = getfault(-87.86, strike, dip, rake)
    data173a = autofault(data173a,0.80657748, 1.110367893, 1.376641327)
    data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

    Pia173a.append(Pa)
    Sia173a.append(Sa)
    Pe173a.append(Pe)
    Se173a.append(Se)


    # #-----S0235b-----
    # path, Pp, Sp, Pa, Sa = eventbuild('235b', 27)
    #
    # data235b, Pe, Se = getfault(-102.31, strike, dip, rake)
    # data235b = autofault(data235b,1.172223922, 0.249311716, 0.212682672)
    # data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)
    #
    # Pia235b.append(Pa)
    # Sia235b.append(Sa)
    # Pe235b.append(Pe)
    # Se235b.append(Se)
    #
    # #---S0325a---
    # path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.4)
    #
    # data325a, Pe, Se = getfault(-60.38, strike, dip, rake)
    # data325a = autofault(data325a,0.399310873, 0.518759571, 1.299137105)
    # data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)
    #
    # Pia325a.append(Pa)
    # Sia325a.append(Sa)
    # Pe325a.append(Pe)
    # Se325a.append(Se)
    #
    # #---S0325ab---
    # path, Pp, Sp, Pa, Sa = eventbuild('325ab', 38.4)
    #
    # data325ab, Pe, Se = getfault(-45.69, strike, dip, rake)
    # data325ab = autofault(data325ab,1.009103448, 0.631448276, 0.625751777)
    # data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)
    #
    # Pia325ab.append(Pa)
    # Sia325ab.append(Sa)
    # Pe325ab.append(Pe)
    # Se325ab.append(Se)
    #
    # #----S0173ab----
    # path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.4)
    #
    # data173ab, Pe, Se = getfault(-91.37, strike, dip, rake)
    # data173ab = autofault(data173ab,1.354580708, 0.81337857, 0.600465199)
    # data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)
    #
    # Pia173ab.append(Pa)
    # Sia173ab.append(Sa)
    # Pe173ab.append(Pe)
    # Se173ab.append(Se)
