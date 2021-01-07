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


model_ls = ['Gudkova']

#---get fault function--
def getfault(st, dp, rk):
    """
    INPUT:
    st, dp, rk = strike, dip, rake from 3 separte lists
    OUTPUT: df = data frame containing model outputs for P,SV,SH amplitudes at each depth input
    ip, ij = exit angles in degrees for each model computed in the Rpattern function
    """
    strike_ls = []; dip_ls = []; rake_ls = []
    P_ls = []; SH_ls =[]; SV_ls=[]
    VH_ls = []; PV_ls = []; PH_ls = []
    exit_ang = []; az_ls = []
    for ex in range(0, 95, 5):
        for a in range(0, 370, 10):
            for strike in st:
                for dip in dp:
                    for rake in rk:
                        strike_ls.append(strike)
                        dip_ls.append(dip)
                        rake_ls.append(rake)
                        exit_ang.append(ex); az_ls.append(a)

                        fault = [strike, dip, rake]
                        mt = getmt(fault)
                        azimuth = a

                        iP = ex
                        jS = ex

                        P,SV,SH = Rpattern(fault,azimuth,[iP,jS])
                        scalefactor = (Pvelz/Svelz)**3
                        SV,SH = SV*scalefactor, SH*scalefactor
                        P_ls.append(P); SH_ls.append(SH); SV_ls.append(SV)


    VHratio = [i / j for i, j in zip(SH_ls, SV_ls)]
    PVratio = [i / j for i, j in zip(P_ls, SV_ls)]
    PHratio = [i / j for i, j in zip(P_ls, SH_ls)]
    plunge = [90-i for i in exit_ang]

    data = {'Exit Angle': exit_ang,
            'Plunge': plunge,
            'Azimuth': az_ls,
            'Strike': strike_ls,
            'Dip': dip_ls,
            'Rake': rake_ls,
            'P': P_ls,
            'SH': SH_ls,
            'SV': SV_ls,
            'SH/SV': VHratio,
            'P/SV': PVratio,
            'P/SH': PHratio}

    df = pd.DataFrame(data, columns = ['Exit Angle', 'Plunge', 'Azimuth','Strike', 'Dip', 'Rake', 'P', 'SV', 'SH', 'SH/SV', 'P/SV', 'P/SH'])
    return df, exit_ang

#first attempt at 173a - try5_173a-35

#first attempt: strike173a = [150]
strike173a = [45]
strike235b = [270]
strike325a = [45]
strike325ab = [90]

#fisrt attempt: dip173a = [60]
dip173a = [90]
dip235b = [70]
dip325a = [90]
dip325ab = [80]

rake = [-90]
#first attempt: rake173a = [80]
rake173a = [-100]
rake325ab = [90]

Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []

# Mars:
radius = 3389.5


for mod in model_ls:
    mars = TauPyModel(model=mod)

#--model velocity input---
    if mod=='Gudkova':
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
        Pia325ab = []; Sia325ab = []; Pe325ab = []; Se325ab = []
        Gudkova_depth = [35]
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

            data173a, ex_ang = getfault(strike173a, dip173a, rake173a)
            data173a.to_csv('jan_173a-' + str(depth) + '.csv', index=False)

            Pia173a.append(Pa)
            Sia173a.append(Sa)
            Pe173a.append(ex_ang)
            Se173a.append(ex_ang)


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

            data235b, ex_ang = getfault(strike235b, dip235b, rake)
            data235b.to_csv('jan_235b-' + str(depth) + '.csv', index=False)

            Pia235b.append(Pa)
            Sia235b.append(Sa)
            Pe235b.append(ex_ang)
            Se235b.append(ex_ang)

            #---S0325a---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=38.5, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325a, ex_ang = getfault(strike325a, dip325a, rake)
            #data325a.to_csv('try2_325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(ex_ang)
            Se325a.append(ex_ang)

            #---S0325ab---
            mtimes = mars.get_travel_times(source_depth_in_km=depth, distance_in_degree=33.6, phase_list=["P", "S"])

            #incident angle at the station
            Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
            Pvel = radius*np.sin(np.radians(Pa))/Pp
            print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
            Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
            Svel = radius*np.sin(np.radians(Sa))/Sp
            print('Check: S surface velocity = ',Svel,' ?')

            data325ab, exit_ang = getfault(strike325ab, dip325ab, rake325ab)
            #data325ab.to_csv('325ab-' + str(depth) + '.csv', index=False)

            Pia325ab.append(Pa)
            Sia325ab.append(Sa)
            Pe325ab.append(exit_ang)
            Se325ab.append(exit_ang)


        incid = {'Model': 'Gudkova',
                'Depth': Gudkova_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a,
                '325ab Pa': Pia325ab,
                '325ab Sa': Sia325ab}

        exit = {'Model': 'Gudkova',
                'Depth': Gudkova_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a,
                '325ab Pe': Pe325ab,
                '325ab Se': Se325ab}

        aGudkova = pd.DataFrame.from_dict(incid)
        eGudkova = pd.DataFrame.from_dict(exit)

    elif mod=='EH45TcoldCrust1b':
        Pia173a = []; Sia173a = []; Pe173a = []; Se173a = []
        Pia235b = []; Sia235b = []; Pe235b = []; Se235b = []
        Pia325a = []; Sia325a = []; Pe325a = []; Se325a = []
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

            data173a, Pe, Se = getfault(-89.94, strike173a, dip173a, rake)
            data173a.to_csv('173a-' + str(depth) + '.csv', index=False)

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

            data235b, Pe, Se = getfault(-102.22, strike235b, dip235b, rake)
            #data235b.to_csv('235b-' + str(depth) + '.csv', index=False)

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

            data325a, Pe, Se = getfault(-60.46, strike325a, dip325a, rake)
            #data325a.to_csv('325a-' + str(depth) + '.csv', index=False)

            Pia325a.append(Pa)
            Sia325a.append(Sa)
            Pe325a.append(Pe)
            Se325a.append(Se)


        incid = {'Model': 'EH45TcoldCrust1b',
                'Depth': coldCrust_depth,
                '173a Pa': Pia173a,
                '173a Sa': Sia173a,
                '235b Pa': Pia235b,
                '235b Sa': Sia235b,
                '325a Pa': Pia325a,
                '325a Sa': Sia325a}

        exit = {'Model': 'EH45TcoldCrust1b',
                'Depth': coldCrust_depth,
                '173a Pe': Pe173a,
                '173a Se': Se173a,
                '235b Pe': Pe235b,
                '235b Se': Se235b,
                '325a Pe': Pe325a,
                '325a Se': Se325a}

        acoldCrust = pd.DataFrame.from_dict(incid)
        ecoldCrust = pd.DataFrame.from_dict(exit)

# dfs = [aGudkova]
# incid_angles = pd.concat(dfs, ignore_index=True)
# incid_angles.to_csv('incident_angles.csv', index=False)
#
# dfe = [eGudkova]
# exit_ang = pd.concat(dfe, ignore_index = True)
# exit_ang.to_csv('exit_angles.csv', index=False)
