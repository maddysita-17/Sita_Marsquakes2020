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


model_ls = ['DWAK', 'EH45Tcold', 'EH45TcoldCrust1b', 'Gudkova', 'LFAK', 'MAAK', 'TAYAK']


#---get fault function--
def getfault(az, st, dp, rk):
    """
    INPUT: az = azimuth in degrees from lat-long-az.py
    st, dp, rk = strike, dip, rake from 3 separte lists
    OUTPUT: df = data frame containing model outputs for P,SV,SH amplitudes at each depth input
    ip, ij = exit angles in degrees for each model computed in the Rpattern function
    """
    strike_ls = []; dip_ls = []; rake_ls = []
    P_ls = []; SH_ls =[]; SV_ls=[];
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
            'P/SH': ratio3}

    df = pd.DataFrame(data, columns = ['Strike', 'Dip', 'Rake', 'P', 'SV', 'SH', 'SH/SV', 'P/SV', 'P/SH'])
    return df, iP, jS

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

def autofault(df, ratio1, ratio2, ratio3):
    n = 0
    for ratio in [ratio1, ratio2, ratio3]:
        rmin = ratio - (0.1 * ratio)
        rmax = ratio + (0.1 * ratio)
        if n == 0:
            ratio1df = df[abs(df['SH/SV'])>rmin]
            ratio1df = ratio1df[abs(ratio1df['SH/SV'])<rmax]
        elif n == 1:
            ratio2df = df[abs(df['P/SV'])>rmin]
            ratio2df = ratio2df[abs(ratio2df['P/SV'])<rmax]
        elif n == 2:
            ratio3df = df[abs(df['P/SH'])>rmin]
            ratio3df = ratio3df[abs(ratio3df['P/SH'])<rmax]
        n += 1

    frames = [ratio1df, ratio2df, ratio3df]
    faults = pd.concat(frames)

    return faults

strike173a = [40,45,50,110,112,115,125,130,140,150,155,160]
strike235b = [90,110,112,113,115,125,130,150,175]
strike325a = [0,1,2,3,40,45,50,90,140,150,160]
strike325ab = [0,1,2,3,40,45,50,85,90,100,140,150,160]

strike173ab = [40,45,50,110,112,115,125,130,140,150,155,160]
strike183a = [10,50,100,150,200,230,240,250,255,260,270]


dip = [0, 20, 45, 60, 70, 80, 90]
rake = [-100, -90, 45, 0, 45, 90, 100]

Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []

Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []

# Mars:
radius = 3389.5


for mod in model_ls:
    mars = TauPyModel(model=mod)

#--model velocity input---

    if mod=='DWAK':
        DWAK_depth = [65, 45, 35, 25, 15, 10, 5]
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
        Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []
        for depth in DWAK_depth:
            if depth <= 66 and depth > 10:
                Pvelz = 5.90405; Svelz = 3.30798
            elif depth <= 10 and depth > 1:
                Pvelz = 5.05377; Svelz = 2.83298
            elif depth <= 1 and depth > 0:
                Pvelz =  3.80908; Svelz = 1.79894
            else:
                print("There is no computed velocity at this depth")

            #----S0173a----

            path, Pp, Sp, Pa, Sa= eventbuild('173a', 29)

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            data173a = autofault(data173a, 0.720430108, 1.064516129, 1.47761194)
            data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            path, Pp, Sp, Pa, Sa = eventbuild('235b', 27.5)

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            data235b = autofault(data235b, 1.169753086, 0.435185185, 0.372031662)
            data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.5)

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            data325a = autofault(data325a, 0.401544402, 0.517374517, 1.288461538)
            data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            path, Pp, Sp, Pa, Sa = eventbuild('325ab', 33.6)

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            data325ab = autofault(data325ab, 0.994475138, 0.621546961, 0.625)
            data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.3)

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            data173ab = autofault(data173ab, 1.365853659, 0.709756098, 0.519642857)
            data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            path, Pp, Sp, Pa, Sa = eventbuild('183a', 43.4)

            data183a, Pe, Se = getfault(-113.9, strike183a, dip, rake)
            data183a = autofault(data183a, 0.512396694, 1.190082645, 2.322580645)
            data183a.to_csv(path + 'S0183a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia183a.append(Pa)
            Sia183a.append(Sa)
            Pe183a.append(Pe)
            Se183a.append(Se)


        incid = {'Model': 'DWAK',
                'Depth': DWAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab,
                '173ab Pa': Pia173ab,
                '173ab Sa': Sia173ab,
                '183a Pa': Pia183a,
                '183a Sa': Sia183a}

        exit = {'Model': 'DWAK',
                'Depth': DWAK_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab,
                '173ab Pe': Pe173ab,
                '173ab Se': Se173ab,
                '183a Pe': Pe183a,
                '183a Se': Se183a}

        aDWAK = pd.DataFrame.from_dict(incid)
        eDWAK = pd.DataFrame.from_dict(exit)


    elif mod=='EH45Tcold':
        EH45_depth = [45, 35, 25, 15, 10, 5]
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
        Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []
        for depth in EH45_depth:
            Pvelz = 6.78574; Svelz = 3.91775

            #----S0173a----

            path, Pp, Sp, Pa, Sa= eventbuild('173a', 29)

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            data173a = autofault(data173a, 1.472527473,1.945054945,1.320895522)
            data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            path, Pp, Sp, Pa, Sa = eventbuild('235b', 27.5)

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            data235b = autofault(data235b, 1.96373057,0.445595855,0.226912929)
            data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.5)

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            data325a = autofault(data325a,0.530612245,0.433673469,0.817307692)
            data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            path, Pp, Sp, Pa, Sa = eventbuild('325ab', 33.6)

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            data325ab = autofault(data325ab, 1.263157895,0.596491228,0.472222222)
            data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.3)

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            data173ab = autofault(data173ab, 1.302325581,0.634883721,0.4875)
            data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            path, Pp, Sp, Pa, Sa = eventbuild('183a', 43.4)

            data183a, Pe, Se = getfault(-113.9, strike183a, dip, rake)
            data183a = autofault(data183a, 0.616161616,1.212121212,1.967213115)
            data183a.to_csv(path + 'S0183a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia183a.append(Pa)
            Sia183a.append(Sa)
            Pe183a.append(Pe)
            Se183a.append(Se)


        incid = {'Model': 'EH45Tcold',
                'Depth': EH45_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab,
                '173ab Pa': Pia173ab,
                '173ab Sa': Sia173ab,
                '183a Pa': Pia183a,
                '183a Sa': Sia183a}

        exit = {'Model': 'EH45Tcold',
                'Depth': EH45_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab,
                '173ab Pe': Pe173ab,
                '173ab Se': Se173ab,
                '183a Pe': Pe183a,
                '183a Se': Se183a}

        aEH45 = pd.DataFrame.from_dict(incid)
        eEH45 = pd.DataFrame.from_dict(exit)

    elif mod=='EH45TcoldCrust1b':
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
        coldCrust_depth = [85, 75, 65, 55, 45, 35, 25, 15, 10, 5]
        for depth in coldCrust_depth:
            if depth <= 85 and depth > 80:
                Pvelz = 7.00139; Svelz = 4.04225
            elif depth <= 80 and depth > 75:
                Pvelz = 6.97778; Svelz = 4.02862
            elif depth <= 75 and depth > 70:
                Pvelz = 6.95417; Svelz = 4.01499
            elif depth <= 70 and depth > 66:
                Pvelz = 6.93056; Svelz = 4.00136
            elif depth <= 66 and depth > 61:
                Pvelz = 6.90694; Svelz = 3.98773
            elif depth <= 61 and depth > 56:
                Pvelz = 6.88333; Svelz = 3.97409
            elif depth <= 56 and depth > 51:
                Pvelz = 6.85972; Svelz = 3.96046
            elif depth <= 51 and depth > 47:
                Pvelz = 6.83611; Svelz = 3.94683
            elif depth <= 47 and depth > 42:
                Pvelz = 5.69750; Svelz = 3.28945
            elif depth <= 42 and depth > 37:
                Pvelz = 5.66444; Svelz = 3.27037
            elif depth <= 37 and depth > 33:
                Pvelz = 5.63139; Svelz = 3.25128
            elif depth <= 33 and depth > 28:
                Pvelz = 5.59833; Svelz = 3.23220
            elif depth <= 28 and depth > 23:
                Pvelz = 5.56528; Svelz = 3.21311
            elif depth <= 23 and depth > 18:
                Pvelz = 5.53222; Svelz = 3.19403
            elif depth <= 18 and depth > 14:
                Pvelz = 5.49917; Svelz = 3.17495
            elif depth <= 14 and depth > 9:
                Pvelz = 5.46611; Svelz = 3.15586
            elif depth <= 9 and depth > 4:
                Pvelz = 5.43306; Svelz = 3.13678
            else:
                print("There is no computed velocity at this depth")

            #----S0173a----

            path, Pp, Sp, Pa, Sa= eventbuild('173a', 29)

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            data173a = autofault(data173a, 0.701570681,1.015706806,1.447761194)
            data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            path, Pp, Sp, Pa, Sa = eventbuild('235b', 27.5)

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            data235b = autofault(data235b, 1.134730539,0.45508982,0.401055409)
            data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.5)

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            data325a = autofault(data325a,0.395437262,0.539923954,1.365384615)
            data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            path, Pp, Sp, Pa, Sa = eventbuild('325ab', 33.6)

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            data325ab = autofault(data325ab,0.983606557,0.636612022,0.647222222)
            data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.3)

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            data173ab = autofault(data173ab,1.393034826,0.701492537,0.503571429)
            data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            path, Pp, Sp, Pa, Sa = eventbuild('183a', 43.4)

            data183a, Pe, Se = getfault(-113.9, strike183a, dip, rake)
            data183a = autofault(data183a, 0.616161616,1.212121212,1.967213115)
            data183a.to_csv(path + 'S0183a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia183a.append(Pa)
            Sia183a.append(Sa)
            Pe183a.append(Pe)
            Se183a.append(Se)


        incid = {'Model': 'EH45TcoldCrust1b',
                'Depth': coldCrust_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab}
                #'173ab Pa': Pia173ab,
                #'173ab Sa': Sia173ab,
                #'183a Pa': Pia183a,
                #'183a Sa': Sia183a}

        exit = {'Model': 'EH45TcoldCrust1b',
                'Depth': coldCrust_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab}
                # '173ab Pe': Pe173ab,
                # '173ab Se': Se173ab,
                # '183a Pe': Pe183a,
                # '183a Se': Se183a}

        acoldCrust = pd.DataFrame.from_dict(incid)
        ecoldCrust = pd.DataFrame.from_dict(exit)

    elif mod=='Gudkova':
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
        Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []
        Gudkova_depth = [45, 35, 25, 15, 10, 5]
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

            path, Pp, Sp, Pa, Sa= eventbuild('173a', 29)

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            data173a = autofault(data173a,0.626168224,0.831775701,1.328358209)
            data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            path, Pp, Sp, Pa, Sa = eventbuild('235b', 27.5)

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            data235b = autofault(data235b,1.035519126,0.43715847,0.422163588)
            data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.5)

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            data325a = autofault(data325a,1.035519126,0.43715847,0.422163588)
            data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            path, Pp, Sp, Pa, Sa = eventbuild('325ab', 33.6)

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            data325ab = autofault(data325ab,1.005586592,0.656424581,0.652777778)
            data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.3)

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            data173ab = autofault(data173ab,1.696969697,0.742424242,0.4375)
            data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            path, Pp, Sp, Pa, Sa = eventbuild('183a', 43.4)

            data183a, Pe, Se = getfault(-113.9, strike183a, dip, rake)
            data183a = autofault(data183a, 0.402906209,0.911492734,2.262295082)
            data183a.to_csv(path + 'S0183a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia183a.append(Pa)
            Sia183a.append(Sa)
            Pe183a.append(Pe)
            Se183a.append(Se)


        incid = {'Model': 'Gudkova',
                'Depth': Gudkova_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab,
                '173ab Pa': Pia173ab,
                '173ab Sa': Sia173ab,
                '183a Pa': Pia183a,
                '183a Sa': Sia183a}

        exit = {'Model': 'Gudkova',
                'Depth': Gudkova_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab,
                '173ab Pe': Pe173ab,
                '173ab Se': Se173ab,
                '183a Pe': Pe183a,
                '183a Se': Se183a}

        aGudkova = pd.DataFrame.from_dict(incid)
        eGudkova = pd.DataFrame.from_dict(exit)

    elif mod=='LFAK':
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
        Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []
        LFAK_depth = [55, 45, 35, 25, 15, 10, 5] #issue at depth = 25 for 235b
        for depth in LFAK_depth:
            if depth <= 56 and depth > 10:
                Pvelz = 6.11300; Svelz = 3.46103
            elif depth <= 10 and depth > 1:
                Pvelz = 5.32184; Svelz = 3.01450
            else:
                print("There is no computed velocity at this depth")

            #----S0173a----

            path, Pp, Sp, Pa, Sa= eventbuild('173a', 29)

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            data173a = autofault(data173a,0.72826087,1.081521739,1.485074627)
            data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            path, Pp, Sp, Pa, Sa = eventbuild('235b', 27.5)

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            data235b = autofault(data235b,1.184375,0.434375,0.366754617)
            data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.5)

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            data325a = autofault(data325a,0.403100775,0.507751938,1.259615385)
            data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            path, Pp, Sp, Pa, Sa = eventbuild('325ab', 33.6)

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            data325ab = autofault(data325ab,1.002785515,0.621169916,0.619444444)
            data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.3)

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            data173ab = autofault(data173ab,1.330166271,0.688836105,0.517857143)
            data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            path, Pp, Sp, Pa, Sa = eventbuild('183a', 43.4)

            data183a, Pe, Se = getfault(-113.9, strike183a, dip, rake)
            data183a = autofault(data183a,0.50330033,1.179867987,2.344262295)
            data183a.to_csv(path + 'S0183a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia183a.append(Pa)
            Sia183a.append(Sa)
            Pe183a.append(Pe)
            Se183a.append(Se)


        incid = {'Model': 'LFAK',
                'Depth': LFAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab,
                '173ab Pa': Pia173ab,
                '173ab Sa': Sia173ab,
                '183a Pa': Pia183a,
                '183a Sa': Sia183a}

        exit = {'Model': 'LFAK',
                'Depth': LFAK_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab,
                '173ab Pe': Pe173ab,
                '173ab Se': Se173ab,
                '183a Pe': Pe183a,
                '183a Se': Se183a}

        aLFAK = pd.DataFrame.from_dict(incid)
        eLFAK = pd.DataFrame.from_dict(exit)

    elif mod=='MAAK':
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
        Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []
        MAAK_depth = [65, 55, 45, 35, 25, 15, 10, 5]
        for depth in MAAK_depth:
            if depth <= 68 and depth > 10:
                Pvelz = 5.94027; Svelz = 3.33676
            elif depth <= 10 and depth > 1:
                Pvelz = 5.09729; Svelz = 2.86612
            else:
                print("There is no computed velocity at this depth")

            #----S0173a----

            path, Pp, Sp, Pa, Sa= eventbuild('173a', 29)

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            data173a = autofault(data173a,0.732240437,1.081967213,1.47761194)
            data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            path, Pp, Sp, Pa, Sa = eventbuild('235b', 27.5)

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            data235b = autofault(data235b,1.169753086,0.441358025,0.377308707)
            data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.5)

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            data325a = autofault(data325a,0.393939394,0.515151515,1.307692308)
            data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            path, Pp, Sp, Pa, Sa = eventbuild('325ab', 33.6)

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            data325ab = autofault(data325ab,0.991735537,0.619834711,0.625)
            data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.3)

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            data173ab = autofault(data173ab,1.355932203,0.702179177,0.517857143)
            data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            path, Pp, Sp, Pa, Sa = eventbuild('183a', 43.4)

            data183a, Pe, Se = getfault(-113.9, strike183a, dip, rake)
            data183a = autofault(data183a,0.504918033,1.180327869,2.337662338)
            data183a.to_csv(path + 'S0183a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia183a.append(Pa)
            Sia183a.append(Sa)
            Pe183a.append(Pe)
            Se183a.append(Se)


        incid = {'Model': 'MAAK',
                'Depth': MAAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab,
                '173ab Pa': Pia173ab,
                '173ab Sa': Sia173ab,
                '183a Pa': Pia183a,
                '183a Sa': Sia183a}

        exit = {'Model': 'MAAK',
                'Depth': MAAK_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab,
                '173ab Pe': Pe173ab,
                '173ab Se': Se173ab,
                '183a Pe': Pe183a,
                '183a Se': Se183a}

        aMAAK = pd.DataFrame.from_dict(incid)
        eMAAK = pd.DataFrame.from_dict(exit)

    elif mod=='TAYAK':
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Pia173ab = []; Sia173ab = []; Pe173ab = []; Se173ab = []
        Pia183a = []; Sia183a = []; Pe183a = []; Se183a = []
        TAYAK_depth = [75, 65, 55, 45, 35, 25, 15, 10, 5]
        for depth in TAYAK_depth:
            if depth <= 77 and depth > 10:
                Pvelz = 5.84666; Svelz = 3.28116
            elif depth <= 10 and depth >1:
                Pvelz = 4.95225; Svelz = 2.78097
            else:
                print("There is no computed velocity at this depth")

            #----S0173a----

            path, Pp, Sp, Pa, Sa= eventbuild('173a', 29)

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            data173a = autofault(data173a,0.71657754,1.064171123,1.485074627)
            data173a.to_csv(path + 'S0173a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            path, Pp, Sp, Pa, Sa = eventbuild('235b', 27.5)

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            data235b = autofault(data235b,1.155487805,0.43597561,0.377308707)
            data235b.to_csv(path + 'S0235b_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.5)

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            data325a = autofault(data325a,0.398467433,0.521072797,1.307692308)
            data325a.to_csv(path + 'S0325a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            path, Pp, Sp, Pa, Sa = eventbuild('325ab', 33.6)

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            data325ab = autofault(data325ab,0.989010989,0.618131868,0.625)
            data325ab.to_csv(path + 'S0325ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.3)

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            data173ab = autofault(data173ab,1.355932203,0.697336562,0.514285714)
            data173ab.to_csv(path + 'S0173ab_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            path, Pp, Sp, Pa, Sa = eventbuild('183a', 43.4)

            data183a, Pe, Se = getfault(-113.9, strike183a, dip, rake)
            data183a = autofault(data183a,0.501644737,1.184210526,2.360655738)
            data183a.to_csv(path + 'S0183a_' + str(mod) + '_' + str(depth) + '.csv', index=False)

            Pia183a.append(Pa)
            Sia183a.append(Sa)
            Pe183a.append(Pe)
            Se183a.append(Se)


        incid = {'Model': 'TAYAK',
                'Depth': TAYAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab,
                '173ab Pa': Pia173ab,
                '173ab Sa': Sia173ab,
                '183a Pa': Pia183a,
                '183a Sa': Sia183a}

        exit = {'Model': 'TAYAK',
                'Depth': TAYAK_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab,
                '173ab Pe': Pe173ab,
                '173ab Se': Se173ab,
                '183a Pe': Pe183a,
                '183a Se': Se183a}

        aTAYAK = pd.DataFrame.from_dict(incid)
        eTAYAK = pd.DataFrame.from_dict(exit)


# dfs = [aDWAK, aEH45, acoldCrust, aGudkova, aLFAK, aMAAK, aTAYAK]
# incid_angles = pd.concat(dfs, ignore_index=True)
# incid_angles.to_csv('incident_angles.csv', index=False)
#
# edfs = [eDWAK, eEH45, ecoldCrust, eGudkova, eLFAK, eMAAK, eTAYAK]
# exit_angles = pd.concat(edfs, ignore_index=True)
# exit_angles.to_csv('exit_angles.csv', index=False)

def fault_search(event):
    path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/faultdata/' + event + '/'
    source_files = sorted(Path(path).glob('*.csv'))

    dataframes = []
    for file in source_files:
        df = pd.read_csv(file)
        df['source'] = file.name
        dataframes.append(df)

    faults = pd.concat(dataframes)
    faults.to_csv(path + "faults_" + event + '.csv', index=False)

    unique_faults = faults.drop_duplicates(subset = ["Strike","Dip","Rake"])
    path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/faultdata/'
    unique_faults.to_csv(path + "uniquefaults_" + event + '.csv', index=False)

fault_search('173a')
fault_search('235b')
fault_search('325a')
fault_search('325ab')
fault_search('173ab')
fault_search('183a')
