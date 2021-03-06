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

strike173a = [40,45,50,110,112,115,125,130,140,150,155,160]
strike235b = [90,110,112,113,115,125,150,175]
strike325a = [0,1,2,3,40,45,50,90,140,150,160]
strike325ab = [0,1,2,3,40,45,50,85,90,100,140,150,160]

strike173ab = [40,45,50,110,112,115,125,130,140,150,155,160]
strike183a = [230,240,250,255,260,270]


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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            #data173a.to_csv('D_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            #data235b.to_csv('D_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            #data325a.to_csv('D_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            #data325ab.to_csv('D_325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=28.3, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            #data173ab.to_csv('D_173ab-' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=43.4, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data183a, Pe, Se = getfault(-97.89, strike183a, dip, rake)
            #data183a.to_csv('D_183a-' + str(depth) + '.csv', index=False)

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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            #data173a.to_csv('EH_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            #data235b.to_csv('EH_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            #data325a.to_csv('EH_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            #data325ab.to_csv('EH_325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=28.3, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            #data173ab.to_csv('EH_173ab-' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=43.4, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data183a, Pe, Se = getfault(-97.89, strike183a, dip, rake)
            #data183a.to_csv('EH_183a-' + str(depth) + '.csv', index=False)

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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            #data173a.to_csv('cC_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            #data235b.to_csv('cC_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            #data325a.to_csv('cC_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            #data325ab.to_csv('cC_325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=28.3, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            #data173ab.to_csv('cC_173ab-' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=43.4, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data183a, Pe, Se = getfault(-97.89, strike183a, dip, rake)
            #data183a.to_csv('cC_183a-' + str(depth) + '.csv', index=False)

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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            #data173a.to_csv('G_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            #data235b.to_csv('G_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            #data325a.to_csv('G_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            #data325ab.to_csv('G_325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=28.3, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            #data173ab.to_csv('G_173ab-' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=43.4, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data183a, Pe, Se = getfault(-97.89, strike183a, dip, rake)
            #data183a.to_csv('G_183a-' + str(depth) + '.csv', index=False)

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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            #data173a.to_csv('L_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            #data235b.to_csv('L_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            #data325a.to_csv('L_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            #data325ab.to_csv('L_325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=28.3, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            #data173ab.to_csv('L_173ab-' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=43.4, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data183a, Pe, Se = getfault(-97.89, strike183a, dip, rake)
            #data183a.to_csv('L_183a-' + str(depth) + '.csv', index=False)

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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            #data173a.to_csv('M_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            #data235b.to_csv('M_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            #data325a.to_csv('M_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            #data325ab.to_csv('M_325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=28.3, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            #data173ab.to_csv('M_173ab-' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=43.4, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data183a, Pe, Se = getfault(-97.89, strike183a, dip, rake)
            #data183a.to_csv('M_183a-' + str(depth) + '.csv', index=False)

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
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=29, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173a, Pe, Se = getfault(-89.94, strike173a, dip, rake)
            #data173a.to_csv('T_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(Pe)
            Se173a.append(Se)


            #-----S0235b-----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=27.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data235b, Pe, Se = getfault(-102.22, strike235b, dip, rake)
            #data235b.to_csv('T_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(Pe)
            Se235b.append(Se)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, Pe, Se = getfault(-60.46, strike325a, dip, rake)
            #data325a.to_csv('T_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, Pe, Se = getfault(-70.73, strike325ab, dip, rake)
            #data325ab.to_csv('T_325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(Pe)
            Se325ab.append(Se)

            #----S0173ab----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=28.3, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data173ab, Pe, Se = getfault(-91.38, strike173ab, dip, rake)
            #data173ab.to_csv('T_173ab-' + str(depth) + '.csv', index=False)

            Pia173ab.append(Pa)
            Sia173ab.append(Sa)
            Pe173ab.append(Pe)
            Se173ab.append(Se)

            #----S0183a----
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=43.4, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            try:
                Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            except:
                print('Within S-wave shadow zone')
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data183a, Pe, Se = getfault(-97.89, strike183a, dip, rake)
            #data183a.to_csv('T_183a-' + str(depth) + '.csv', index=False)

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


dfs = [aDWAK, aEH45, acoldCrust, aGudkova, aLFAK, aMAAK, aTAYAK]
incid_angles = pd.concat(dfs, ignore_index=True)
incid_angles.to_csv('incident_angles.csv', index=False)

edfs = [eDWAK, eEH45, ecoldCrust, eGudkova, eLFAK, eMAAK, eTAYAK]
exit_angles = pd.concat(edfs, ignore_index=True)
exit_angles.to_csv('exit_angles.csv', index=False)
