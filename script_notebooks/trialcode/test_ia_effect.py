import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel
import pandas as pd

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

    AP = sR*(3*np.cos(iP)**2 - 1) - qR*np.sin(2*iP) - pR*np.sin(iP)**2
    ASV = 1.5*sR*np.sin(2*jS) + qR*np.cos(2*jS) + 0.5*pR*np.sin(2*jS)
    ASH = -qL*np.cos(jS) - pL*np.sin(jS)

    return AP,ASV,ASH


model_ls = ['DWAK', 'EH45Tcold', 'EH45TcoldCrust1b', 'Gudkova', 'LFAK', 'MAAK', 'TAYAK']
#incid_angles = pd.DataFrame(columns = ['Model', 'Depth', '173a Pa', '173a Sa', '235b Pa', '235b Sa', '325a Pa', '325a Sa'])

#---get fault function--
def getfault(az, st, dp, rk, ia):
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

                P,SV,SH = Rpattern(fault,azimuth,ia)
                scalefactor = (Pvelz/Svelz)**3
                SV,SH = SV*scalefactor, SH*scalefactor
                P_ls.append(P); SH_ls.append(SH); SV_ls.append(SV)

    ratio = [i / j for i, j in zip(SH_ls, SV_ls)]

    data = {'Strike': strike_ls,
            'Dip': dip_ls,
            'Rake': rake_ls,
            'P': P_ls,
            'SH': SH_ls,
            'SV': SV_ls,
            'SH/SV': ratio}

    df = pd.DataFrame(data, columns = ['Strike', 'Dip', 'Rake', 'P', 'SV', 'SH', 'SH/SV'])
    return df

strike173a = [110,113]
strike235b = [110,112,114]
strike325a = [1,2]


dip = [60, 70, 80]
rake = [-90, 0, 90]

Pia173a = []; Sia173a = []
Pia235b = []; Sia235b = []
Pia325a = []; Sia325a = []

# Mars:
radius = 3389.5


for mod in model_ls:
    mars = TauPyModel(model=mod)

#--model velocity input---

    if mod=='DWAK':
        DWAK_depth = [45, 15, 5]
        Pia173a = []; Sia173a = []
        Pia235b = []; Sia235b = []
        Pia325a = []; Sia325a = []
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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data173a = getfault(-89.94, strike173a, dip, rake, [27.9,24.0])
            print('ia: 27.9 & 24.0')
            print(data173a)
            data173a = getfault(-89.94, strike173a, dip, rake, [28.4,24.5])
            print('ia: 28.4 & 24.5')
            print(data173a)

            #data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)

            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data235b = getfault(-102.22, strike235b, dip, rake, [28.0, 24.1])
            print('ia: 28.0 & 24.1')
            print(data235b)
            data235b = getfault(-102.22, strike235b, dip, rake, [28.5, 24.6])
            print('ia: 28.5 & 24.6')
            print(data235b)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data325a = getfault(-60.46, strike325a, dip, rake, [26.4, 23.2])
            print('ia: 26.4 & 23.2')
            print(data325a)
            data325a = getfault(-60.46, strike325a, dip, rake, [26.9, 23.7])
            print('ia: 26.9 & 23.7')
            print(data325a)

            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)


        incid = {'Model': 'DWAK',
                'Depth': DWAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a}

        aDWAK = pd.DataFrame.from_dict(incid)

    elif mod=='EH45Tcold':
        EH45_depth = [45, 15, 5]
        Pia173a = []; Sia173a = []
        Pia235b = []; Sia235b = []
        Pia325a = []; Sia325a = []
        for depth in EH45_depth:
            Pvelz = 6.78574; Svelz = 3.91775

            #----S0173a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data173a = getfault(-89.94, strike173a, dip, rake, [56.4, 57.9])
            print(data173a)
            data173a = getfault(-89.94, strike173a, dip, rake, [56.9, 58.4])
            print(data173a)
            #data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)

            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data235b = getfault(-102.22, strike235b, dip, rake, [56.7, 57.9])
            print('ia: 56.7 & 57.9')
            print(data235b)
            data235b = getfault(-102.22, strike235b, dip, rake, [57.2, 58.4])
            print('ia: 57.2 & 58.4')
            print(data235b)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data325a = getfault(-60.46, strike325a, dip, rake, [54.1, 57.0])
            print('ia: 54.1 & 57.0')
            print(data325a)
            data325a = getfault(-60.46, strike325a, dip, rake, [54.6, 57.5])
            print('ia: 54.6 & 57.5')
            print(data325a)

            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)

        incid = {'Model': 'EH45Tcold',
                'Depth': EH45_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a}

        aEH45 = pd.DataFrame.from_dict(incid)

    elif mod=='EH45TcoldCrust1b':
        Pia173a = []; Sia173a = []
        Pia235b = []; Sia235b = []
        Pia325a = []; Sia325a = []
        coldCrust_depth = [85, 45, 15, 5]
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

            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data173a = getfault(-89.94, strike173a, dip, rake, [19.3, 18.9])
            print('ia: 19.3 & 18.9')
            print(data173a)
            data173a = getfault(-89.94, strike173a, dip, rake, [19.8, 18.9])
            print('ia: 19.8 & 18.9')
            print(data173a)

            #data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)

            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data235b = getfault(-102.22, strike235b, dip, rake, [19.3, 19.6])
            print('ia: 19.3 & 19.6')
            print(data235b)
            data235b = getfault(-102.22, strike235b, dip, rake, [19.8, 20.1])
            print('ia: 19.8 & 20.1')
            print(data235b)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data325a = getfault(-60.46, strike325a, dip, rake, [18.7, 18.7])
            print('ia: 18.7 & 18.7')
            print(data325a)
            data325a = getfault(-60.46, strike325a, dip, rake, [19.2, 19.2])
            print('ia: 19.2 & 19.2')
            print(data325a)

            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)


        incid = {'Model': 'EH45TcoldCrust1b',
                'Depth': coldCrust_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a}

        acoldCrust = pd.DataFrame.from_dict(incid)

    elif mod=='Gudkova':
        Pia173a = []; Sia173a = []
        Pia235b = []; Sia235b = []
        Pia325a = []; Sia325a = []
        Gudkova_depth = [45, 15, 5]
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

            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data173a = getfault(-89.94, strike173a, dip, rake, [1.8, 1.9])
            print('ia: 1.8 & 1.9')
            print(data173a)
            data173a = getfault(-89.94, strike173a, dip, rake, [2.3, 2.4])
            print('ia: 2.3, 2.4')
            print(data173a)
            #data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)

            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data235b = getfault(-89.94, strike173a, dip, rake, [1.8, 1.9])
            print('ia: 1.8 & 1.9')
            print(data235b)
            data235b = getfault(-89.94, strike173a, dip, rake, [2.3, 2.4])
            print('ia: 2.3, 2.4')
            print(data235b)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data325a = getfault(-60.46, strike325a, dip, rake, [1.7, 1.8])
            print('ia: 1.7 & 1.8')
            print(data325a)
            data325a = getfault(-60.46, strike325a, dip, rake, [2.2, 2.3])
            print('ia: 2.2 & 2.3')
            print(data325a)

            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)

        incid = {'Model': 'Gudkova',
                'Depth': Gudkova_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a}

        aGudkova = pd.DataFrame.from_dict(incid)

    elif mod=='LFAK':
        Pia173a = []; Sia173a = []
        Pia235b = []; Sia235b = []
        Pia325a = []; Sia325a = []
        LFAK_depth = [55, 25, 15, 5]
        for depth in LFAK_depth:
            if depth <= 56 and depth > 10:
                Pvelz = 6.11300; Svelz = 3.46103
            elif depth <= 10 and depth > 1:
                Pvelz = 5.32184; Svelz = 3.01450
            else:
                print("There is no computed velocity at this depth")

            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data173a = getfault(-89.94, strike173a, dip, rake, [30.3, 26.2])
            print('ia: 30.2 & 26.2')
            print(data173a)
            data173a = getfault(-89.94, strike173a, dip, rake, [30.8, 26.7])
            print('ia: 30.8 & 26.7')
            print(data173a)
            #data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)

            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data235b = getfault(-89.94, strike173a, dip, rake, [30.5, 26.3])
            print('ia: 30.5 & 26.3')
            print(data235b)
            data235b = getfault(-89.94, strike173a, dip, rake, [31.0, 26.8])
            print('ia: 31.0 & 26.8')
            print(data235b)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data325a = getfault(-60.46, strike325a, dip, rake, [29.1, 25.8])
            print('ia: 29.1 & 25.8')
            print(data325a)
            data325a = getfault(-60.46, strike325a, dip, rake, [29.6, 26.3])
            print('ia: 29.6 & 26.3')
            print(data325a)

            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)

        incid = {'Model': 'LFAK',
                'Depth': LFAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a}

        aLFAK = pd.DataFrame.from_dict(incid)

    elif mod=='MAAK':
        Pia173a = []; Sia173a = []
        Pia235b = []; Sia235b = []
        Pia325a = []; Sia325a = []
        MAAK_depth = [65, 45, 15, 5]
        for depth in MAAK_depth:
            if depth <= 68 and depth > 10:
                Pvelz = 5.94027; Svelz = 3.33676
            elif depth <= 10 and depth > 1:
                Pvelz = 5.09729; Svelz = 2.86612
            else:
                print("There is no computed velocity at this depth")

            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data173a = getfault(-89.94, strike173a, dip, rake, [27.4, 23.9])
            print('ia: 27.4 & 23.9')
            print(data173a)
            data173a = getfault(-89.94, strike173a, dip, rake, [27.9, 24.4])
            print('ia: 27.9 & 24.4')
            print(data173a)
            #data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)

            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data235b = getfault(-89.94, strike173a, dip, rake, [27.5, 24.0])
            print('ia: 27.5 & 24.0')
            print(data235b)
            data235b = getfault(-89.94, strike173a, dip, rake, [28.0, 24.5])
            print('ia: 28.0 & 24.5')
            print(data235b)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data325a = getfault(-60.46, strike325a, dip, rake, [25.9, 23.5])
            print('ia: 25.9 & 23.5')
            print(data325a)
            data325a = getfault(-60.46, strike325a, dip, rake, [26.4, 24.0])
            print('ia: 26.4 & 24.0')
            print(data325a)

            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)

        incid = {'Model': 'MAAK',
                'Depth': MAAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a}

        aMAAK = pd.DataFrame.from_dict(incid)

    elif mod=='TAYAK':
        Pia173a = []; Sia173a = []
        Pia235b = []; Sia235b = []
        Pia325a = []; Sia325a = []
        TAYAK_depth = [75, 45, 15, 5]
        for depth in TAYAK_depth:
            if depth <= 77 and depth > 10:
                Pvelz = 5.84666; Svelz = 3.28116
            elif depth <= 10 and depth >1:
                Pvelz = 4.95225; Svelz = 2.78097
            else:
                print("There is no computed velocity at this depth")

            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data173a = getfault(-89.94, strike173a, dip, rake, [26.5, 22.8])
            print('ia: 26.5 & 22.8')
            print(data173a)
            data173a = getfault(-89.94, strike173a, dip, rake, [27.0, 23.3])
            print('ia: 27.0 & 23.3')
            print(data173a)
            #data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)

            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data235b = getfault(-89.94, strike173a, dip, rake, [26.6, 22.9])
            print('ia: 26.6 & 22.9')
            print(data235b)
            data235b = getfault(-89.94, strike173a, dip, rake, [27.1, 23.4])
            print('ia: 27.1 & 23.4')
            print(data235b)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print(Pa)
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print(Sa)

            data325a = getfault(-60.46, strike325a, dip, rake, [25.2, 22.2])
            print('ia: 25.2 & 22.2')
            print(data325a)
            data325a = getfault(-60.46, strike325a, dip, rake, [25.7, 22.7])
            print('ia: 25.7 & 22.7')
            print(data325a)

            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)


        incid = {'Model': 'TAYAK',
                'Depth': TAYAK_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a}

        aTAYAK = pd.DataFrame.from_dict(incid)
