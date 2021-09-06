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

    # AP = np.abs(sR*(3*np.cos(iP)**2 - 1) - qR*np.sin(2*iP) - pR*np.sin(iP)**2)
    # ASV = np.abs(1.5*sR*np.sin(2*jS) + qR*np.cos(2*jS) + 0.5*pR*np.sin(2*jS))
    # ASH = np.abs(-qL*np.cos(jS) - pL*np.sin(jS))

    AP = sR*(3*np.cos(iP)**2 - 1) - qR*np.sin(2*iP) - pR*np.sin(iP)**2
    ASV = 1.5*sR*np.sin(2*jS) + qR*np.cos(2*jS) + 0.5*pR*np.sin(2*jS)
    ASH = -qL*np.cos(jS) - pL*np.sin(jS)

    #emperical finding that the SH amplitude should be rotated by 180
    ASH = ASH * -1

    return AP,ASV,ASH

def oldRpattern(fault,azimuth,incidence_angles):
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

    # AP = np.abs(sR*(3*np.cos(iP)**2 - 1) - qR*np.sin(2*iP) - pR*np.sin(iP)**2)
    # ASV = np.abs(1.5*sR*np.sin(2*jS) + qR*np.cos(2*jS) + 0.5*pR*np.sin(2*jS))
    # ASH = np.abs(-qL*np.cos(jS) - pL*np.sin(jS))

    AP = np.abs(sR*(3*np.cos(iP)**2 - 1) - qR*np.sin(2*iP) - pR*np.sin(iP)**2)
    ASV = np.abs(1.5*sR*np.sin(2*jS) + qR*np.cos(2*jS) + 0.5*pR*np.sin(2*jS))
    ASH = np.abs(-qL*np.cos(jS) - pL*np.sin(jS))

    return AP,ASV,ASH


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

                iP = np.degrees(np.arcsin(Pvelz*Pp/(radius-depth)))
                jS = np.degrees(np.arcsin(Svelz*Sp/(radius-depth)))
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
            'SV': SV_ls}

    df = pd.DataFrame(data, columns = ['Strike', 'Dip', 'Rake', 'P', 'SV', 'SH'])
    return df, iP, jS

def oldgetfault(az, st, dp, rk):
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

                P,SV,SH = oldRpattern(fault,azimuth,[iP,jS])
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

def oldautofault(df, obs_SH, obs_SV, obs_P, chi):
    SHSV = df['SH/SV']
    PSV = df['P/SV']
    PSH = df['P/SH']

    obs_SHSV = obs_SH/obs_SV
    obs_PSV = obs_P/obs_SV
    obs_PSH = obs_P/obs_SH

    ratio1_ls = [(np.arctan(obs_SHSV)-np.arctan(i))**2 for i in SHSV]
    ratio2_ls = [(np.arctan(obs_PSV)-np.arctan(i))**2 for i in PSV]
    ratio3_ls = [(np.arctan(obs_PSH)-np.arctan(i))**2 for i in PSH]

    df['Ratio1'] = ratio1_ls
    df['Ratio2'] = ratio2_ls
    df['Ratio3'] = ratio3_ls

    sum = df['Ratio1'] + df['Ratio2'] + df['Ratio3']
    df['Sum'] = sum

    faults = df[df['Sum']<= chi]

    return faults

def autofault(df, obs_P, obs_SV, obs_SH, errP, errS):
    # hypothetically observed amplitudes P, SV, SH
    vobserved = np.array([obs_P,obs_SV,obs_SH])
    vobslength = np.linalg.norm(vobserved)
    # and errors (always positive):
    eobs = np.array([errP, errS, errS])
    # normalized:
    # normal = vobserved/vobslength
    n = vobserved/vobslength


    # eca = np.arctan(np.max([eobs[0]*(1-normal[0]),eobs[1]*(1-normal[1]),eobs[2]*(1-normal[2])])/vobslength)
    eca = np.arctan(np.max([eobs[0]*(1-n[0]),eobs[1]*(1-n[1]),eobs[2]*(1-n[2])])/vobslength)


    print('cutoff value: ', eca)

    #read in modeled csv file
    xd = df['P']
    yd = df['SV']
    zd = df['SH']
    ncalc = len(zd)

    len_xd = len(xd)
    len_yd = len(yd)
    if ncalc != len_xd:
        print('ERROR')
    if ncalc != len_yd:
        print('ERROR')

    vcall = np.array([xd,yd,zd])
    vca = vcall.T

    # misfit:
    mfd = np.zeros(ncalc)
    mf3D = np.zeros(ncalc)
    mf1 = np.zeros(ncalc)
    mf2 = np.zeros(ncalc)
    mf3 = np.zeros(ncalc)
    select = []

    for i in np.arange(ncalc):
        # angle in 3 dimensions: (in radians)
        # mf3D[i] = np.arccos(np.dot(normal,vca[i])/np.linalg.norm(vca[i]))   # should be valued between 0 and pi
        mf3D[i] = np.arccos(np.dot(n,vca[i])/np.linalg.norm(vca[i]))   # should be valued between 0 and pi
        if mf3D[i] < 3*eca:
        #if mf3D[i] > -1:
            select.append(i)

    st_ls = []; dp_ls =[]; rk_ls=[]; mf_ls = []; extra = []
    for i in select:
        st = df.at[i,'Strike'] ; dp = df.at[i,'Dip'] ; rk = df.at[i,'Rake']
        st_ls.append(st); dp_ls.append(dp); rk_ls.append(rk)
        mf_ls.append(mf3D[i])
        extra.append([i,vca[i],n])

    faults = {'Strike': st_ls,
                'Dip': dp_ls,
                'Rake': rk_ls,
                'Misfit': mf_ls,
                'Extra': extra}
    posfaults = pd.DataFrame.from_dict(faults)

    return posfaults


strike173a = [*range(0, 180, 2)]
strike173ab = [*range(0, 180, 2)]

strike235b = [*range(0, 180, 2)]

strike325a = [*range(0, 180, 2)]
strike325ab = [*range(0, 180, 2)]


dip = [*range(0,90,2)]
rake = [*range(-100,100,5)]


# Mars:
radius = 3389.5

bbpath = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/beachballs/csvs/'

mod = 'NewGudkova'
mars = TauPyModel(model=mod)

#velocity at 35km depth
depth = 35
Pvelz = 7.13900; Svelz = 4.01900


#----S0173a-----
path, Pp, Sp, Pa, Sa = eventbuild('173a', 28.4)

data173a, Pe, Se = getfault(-87.86, strike173a, dip, rake)
# data173a = autofault(data173a, 144, -176, 195, 0.24)
data173a_sorted = autofault(data173a, -1.25, 0.955, -0.37, 0.024, 0.186)
data173a_sorted = data173a_sorted.sort_values(by=['Misfit'])
# slice173a = data173a_sorted[:45]
data173a_sorted.to_csv(bbpath + 'resp_S0173a.csv', index=False)

# #---strike slip solution----
# rev173a = getfault(-87.86, [152], [78], [40])
# print('Strike-Slip Soln: ', rev173a)
# # rev173a = autofault(data173a, -0.97692, 0.652866, -0.298912, 0.024, 0.186)
# # print(rev173a)
# # rev173a.to_csv(bbpath + 'rev_ss_173a.csv', index=False)
#
# #---thrust solution----
# rev173a = getfault(-87.86, [178], [66], [80])
# print('Thrust Soln: ', rev173a)
# # rev173a = autofault(data173a, -0.977194, 0.750319, -0.434096, 0.024, 0.186)
# # print(rev173a)
# # rev173a.to_csv(bbpath + 'rev_t_173a.csv', index=False)
#
# #------insight mech----
# ins173a = getfault(-87.86, [60], [70], [-90])
# print('Insight Mech: ', ins173a)


#----S0173ab----
path, Pp, Sp, Pa, Sa = eventbuild('173ab', 28.4)

data173ab, Pe, Se = getfault(-91.37, strike173ab, dip, rake)
data173ab_sorted = autofault(data173ab, 1.11, -1.21, -3.26, 0.024, 0.186)
# data173ab_sorted = data173ab_sorted.sort_values(by=['Misfit'])
data173ab_sorted.to_csv(bbpath + 'resp_S0173ab.csv', index=False)
#
# #-----strike-slip solution-----
# rev173ab = getfault(-91.37, [8], [28], [-45])
# print('Stike-Slip Soln: ', rev173ab)
# # rev173ab = autofault(data173ab, 0.842799, -1.075585, -2.49031, 0.024, 0.186)
# # print(rev173ab)
# # rev173ab.to_csv(bbpath + 'rev_ss_173ab.csv', index=False)
# #
# #----normal solution------
# rev173ab = getfault(-91.37, [148], [66], [-95])
# print('Normal Soln: ', rev173ab)
# # rev173ab = autofault(data173ab, 0.757028, -0.950093, -2.212736, 0.024, 0.186)
# # print(rev173ab)
# # rev173ab.to_csv(bbpath + 'rev_n_173ab.csv', index=False)
#
#
#-----S0235b-----
path, Pp, Sp, Pa, Sa = eventbuild('235b', 27)

data235b, Pe, Se = getfault(-102.31, strike235b, dip, rake)
# data235b = autofault(data235b,-318, -234, -80, 0.11)
data235b_sorted = autofault(data235b, 0.360, -1.81, 0, 0.35, 0.105)
# data235b_sorted = data235b_sorted.sort_values(by=['Misfit'])
data235b_sorted.to_csv(bbpath+ 'resp_S0235b.csv', index=False)
#
# #-----normal solution----
# rev235b = getfault(-102.31, [0], [36], [-75])
# print('Normal Soln: ', rev235b)
# # rev235b = autofault(data235b, 0.676262, 4.064643, -0.001297, 0.35, 0.105)
# # print(rev235b)
# # rev235b.to_csv(bbpath + 'rev_235b.csv', index=False)
#
# #-----insight mech-------
# ins235b = getfault(-102.31, [0], [40], [-80])
# print('Insight Mech: ', ins235b)

# #
# # # #-----S0235b alt-----
# # # path, Pp, Sp, Pa, Sa = eventbuild('235bi', 15)
# # #
# # # data235bi, Pe, Se = getfault(-104.34, strike235b, dip, rake)
# # # data235bi = autofault(data235bi,1, 1, -80, 2)
# # # data235bi.to_csv(bbpath + 'S0235bi.csv', index=False)
# #
# # #---S0325a---
# # path, Pp, Sp, Pa, Sa = eventbuild('325a', 38.4)
# #
# # data325a, Pe, Se = getfault(-60.38, strike325a, dip, rake)
# # data325a_sorted = autofault(data325a, -0.65, -1.48, 0.73, 1.0, 1.2)
# # data325a_sorted = data325a_sorted.sort_values(by=['Misfit'])
# # slice325a = data325a_sorted[:45]
# # slice325a.to_csv(bbpath + 'resp_S0325a.csv', index=False)
# #
#---S0325ab---
path, Pp, Sp, Pa, Sa = eventbuild('325ab', 38.4)

data325ab, Pe, Se = getfault(-45.69, strike325ab, dip, rake)
data325ab_sorted = autofault(data325ab, -1.35, -4.08, 0, 1.0, 1.2)
# data325ab_sorted = data325ab_sorted.sort_values(by=['Misfit'])
data325ab_sorted.to_csv(bbpath + 'resp_S0325ab.csv', index=False)
# #
# #----vert dip slip solution----
# rev325ab = getfault(-45.69, [2], [2], [45])
# print('Vert Dip Slip Soln: ', rev325ab)
# # rev325ab = autofault(data325ab, -0.78395, -1.887205, 1.745873, 1.0, 1.2)
# # print(rev325ab)
# # rev325ab.to_csv(bbpath + 'rev_325ab.csv', index=False)
