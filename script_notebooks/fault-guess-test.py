import matplotlib.pyplot as plot
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

def getamp(azimuth, strike, dip, rake, rayp):
    """
    INPUT:  az = azimuth in degrees from the event to the lander
            st, dp, rk = strike dip rake from 3 separate lists
            Pp, Sp = ray paramters calculated from the model in obspy
            Pvelz, Svelz = velocity @ depth from model
            radius = predefined radius of planet
    OUTPUT: df = dataframe containing synthetic amplitudes
            ip, ij = exit angles (???)
    """
    # define empty lists
    strike_ls = []; dip_ls = []; rake_ls = []
    P_ls = []; SH_ls = []; SV_ls = []

    # loop over fault plane combinations
    for st in strike:
        for dp in dip:
            for rk in rake:
                strike_ls.append(st); dip_ls.append(dp); rake_ls.append(rk)

                # define single fault for calculations
                fault = [st, dp, rk]

                # calculating exit angles using the models velocity @ depth & ray parameters
                # radius should be the radius @ depth
                iP = np.degrees(np.arcsin(Pvelz*rayp[0]/(radius-depth)))
                jS = np.degrees(np.arcsin(Svelz*rayp[1]/(radius-depth)))

                # calculating amplitudes
                P,iSV,iSH = Rpattern(fault, azimuth, [iP, jS])
                scalefactor = (Pvelz/Svelz)**3
                SV,SH = iSV*scalefactor, iSH*scalefactor
                P_ls.append(P); SH_ls.append(SH); SV_ls.append(SV)

    # creating dataframe
    data = {
            'Model': mod,
            'Depth': depth,
            'Strike': strike_ls,
            'Dip': dip_ls,
            'Rake': rake_ls,
            'P': P_ls,
            'SV': SV_ls,
            'SH': SH_ls
            }

    df = pd.DataFrame(data, columns = ['Model','Depth','Strike', 'Dip', 'Rake', 'P', 'SV', 'SH'])
    return df, iP, jS

def eventbuild(event, dist):
    # determine travel times using obspy
    mtimes = mars.get_travel_times(source_depth_in_km= depth, distance_in_degree=dist, phase_list=['P','S'])

    # ray parameters & incidence angles at the station
    Pp = mtimes[0].ray_param ; Pa = mtimes[0].incident_angle

    try:
        Sp = mtimes[1].ray_param ; Sa = mtimes[1].incident_angle
    except:
        Sp = 0 ; Sa = 0
        print('Within S-wave shadow zone')

    print('ray parameters: ', Pp, Sp)
    print('incid angles: ', Pa, Sa)

    return Pp, Sp, Pa, Sa

def autofault(df, obsP, obsSV, obsSH, errP, errS):
    # hypothetically observed amplitudes P, SV, SH:
    vobserved = np.array([obsP, obsSV, obsSH])
    vobslength = np.linalg.norm(vobserved)
    # adding errors:
    eobs = np.array([errP, errS, errS])
    # normalized:
    norm = vobserved/vobslength
    print('normalized array: ',norm)

    # defining cutoff value:
    eca_ls = [eobs[0]*(1-norm[0]),eobs[1]*(1-norm[1]),eobs[2]*(1-norm[2])]/vobslength
    eca = np.arctan(np.max(eca_ls))
    print('cutoff value: ', eca)

    xd = df['P']
    yd = df['SV']
    zd = df['SH']
    ncalc = len(zd)

    len_xd = len(xd); len_yd = len(yd)
    if ncalc != len_xd:
        print('ERROR xd LENGTH')
    if ncalc != len_yd:
        print('ERROR yd LENGHT')

    vcall = np.array([xd,yd,zd])
    vca = vcall.T

    # ------misfit-------
    # empty array
    mf3d = np.zeros(ncalc)

    # list of index values for fault planes below misfit val
    select = []

    for i in np.arange(ncalc):
        # angle in 3 dimensions: (in radians)
        mf3d[i] = np.arccos(np.dot(norm, vca[i])/np.linalg.norm(vca[i]))
        if mf3d[i] < eca:
            select.append(i)

    # pulling fault plane data associated w/ index value
    st_ls = []; dp_ls =[]; rk_ls=[]; mf_ls = []; extra = []
    for i in select:
        st = df.at[i,'Strike']
        dp = df.at[i,'Dip']
        rk = df.at[i,'Rake']
        st_ls.append(st); dp_ls.append(dp); rk_ls.append(rk)

        mf_ls.append(mf3d[i])
        extra.append([depth,i,vca[i],norm])

    faults = {'Strike': st_ls,
                'Dip': dp_ls,
                'Rake': rk_ls,
                'Misfit': mf_ls,
                'Extra': extra}
    posfaults = pd.DataFrame.from_dict(faults)
    # print('before sort: ', posfaults)

    faults_sorted = posfaults.sort_values(by=['Misfit'])
    return faults_sorted

def predictamp(data,az,rayp):
    for index, rows in data.iterrows():
        ampsdf, iP, jS = getamp(az, [rows.Strike], [rows.Dip], [rows.Rake], rayp)
        ampsdf.to_csv('amps_'+ str(az)+'.csv', mode='a', header=False)


# ---- FORWARD CALC (FROM MECH) -----

mechanism_dict = {
                    # 'S0173a': [[60,70,-90],],      #insight mech (march)
                    # 'S0235b': [[0,40,-80],],

                    'S0173a':   [[162,71,56],],       #top mechs - NewGudkova @ 35km
                    'S0173ab':  [[30,48,-12],],
                    'S0235b':   [[34,69,-36],],
                    'S0325ab':  [[16,1,76],],
                    }

def forwardcalc(depth,mechdic,alt=False):
    # ----------- S0173A ---------------
    if 'S0173a' in mechdic:
        s0173a = mechdic['S0173a'][0]
        if alt == True:
            s0173a_b = mechdic['S0173a'][1]

        print('173a')
        Pp, Sp, Pa, Sa = eventbuild('173a', 28.4)
        data173a, Pe, Se = getamp(-87.86, [s0173a[0]], [s0173a[1]], [s0173a[2]], [Pp, Sp])

        if alt == True:
            data1, Pe, Se = getamp(-87.86, [s0173a_b[0]], [s0173a_b[1]], [s0173a_b[2]], [Pp, Sp])
            data173a = pd.concat([data173a,data1],ignore_index=True)

        print(data173a)
        for i in range(len(data173a)):
            data173a.at[i,'Model'] = 'ForwardCalc'
        data173a.to_csv(path+'amps_-87.86.csv', mode='a',header=False)
    else:
        pass

    # ----------- S0173AB ---------------
    if 'S0173ab' in mechdic:
        s0173ab = mechdic['S0173ab'][0]
        if alt == True:
            s0173ab_b = mechdic['S0173ab'][1]

        print('173ab')
        Pp, Sp, Pa, Sa = eventbuild('173ab', 28.4)
        data173ab, Pe, Se = getamp(-91.37, [s0173ab[0]], [s0173ab[1]], [s0173ab[2]], [Pp, Sp])

        if alt == True:
            data1, Pe, Se = getamp(-91.37, [s0173ab_b[0]], [s0173ab_b[1]], [s0173ab_b[2]], [Pp, Sp])
            data173ab = pd.concat([data173ab,data1],ignore_index=True)

        print(data173ab)
        for i in range(len(data173ab)):
            data173ab.at[i,'Model'] = 'ForwardCalc'
        data173ab.to_csv(path+'amps_-91.37.csv', mode='a',header=False)
    else:
        pass

    # ----------- S0235B ---------------
    if 'S0235b' in mechdic:
        s0235b = mechdic['S0235b'][0]
        if alt == True:
            s0235b_b = mechdic['S0235b'][1]

        print('235b')
        Pp, Sp, Pa, Sa = eventbuild('235b', 27)
        data235b, Pe, Se = getamp(-102.31, [s0235b[0]], [s0235b[1]], [s0235b[2]], [Pp, Sp])

        if alt == True:
            data1, Pe, Se = getamp(-102.31, [s0235b_b[0]], [s0235b_b[1]], [s0235b_b[2]], [Pp, Sp])
            data235b = pd.concat([data235b,data1],ignore_index=True)

        print(data235b)
        for i in range(len(data235b)):
            data235b.at[i,'Model'] = 'ForwardCalc'
        data235b.to_csv(path+'amps_-102.31.csv', mode='a',header=False)
    else:
        pass

    # ----------- S0235B ---------------
    if 'S0325ab' in mechdic:
        s0325ab = mechdic['S0325ab'][0]
        if alt == True:
            s0325ab_b = mechdic['S0325ab'][1]

        print('325ab')
        Pp, Sp, Pa, Sa = eventbuild('325ab', 38.4)
        data325ab, Pe, Se = getamp(-60.38, [s0325ab[0]], [s0325ab[1]], [s0325ab[2]], [Pp, Sp])

        if alt == True:
            data1, Pe, Se = getamp(-60.38, [s0325ab_b[0]], [s0325ab_b[1]], [s0325ab_b[2]], [Pp, Sp])
            data325ab = pd.concat([data325ab,data1],ignore_index=True)

        print(data325ab)
        for i in range(len(data325ab)):
            data325ab.at[i,'Model'] = 'ForwardCalc'
        data325ab.to_csv(path+'amps_-60.38.csv', mode='a',header=False)
    else:
        pass

    return


# ---------- MARS -----------------
radius = 3389.5

# mod = 'NewGudkova'
# mars = TauPyModel(model=mod)

path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/'

# #Gudkova Model
# depth = 35
# Pvelz = 7.13900; Svelz = 4.01900

# ------ FAULT PLANE LISTS -------
strike_rang = [*range(0, 180, 2)]
dip_rang = [*range(0,90,1)]
rake_rang = [*range(-100,100,2)]

def eventoutput(depth,rank):
    # ----------- S0173A ---------------
    print('173a')
    Pp, Sp, Pa, Sa = eventbuild('173a', 28.4)
    data173a, Pe, Se = getamp(-87.86, strike_rang, dip_rang, rake_rang, [Pp, Sp])
    # print('Exit Angles: ', Pe, Se)
    data173a = autofault(data173a, -1.25, 0.955, -0.371, 0.024, 0.186)
    data173a[:rank].to_csv('check173a.csv', mode='a', header=False)
    predictamp(data173a[:rank],-87.86,[Pp,Sp])

    #------------ S0173AB ---------------
    print('173ab')
    Pp, Sp, Pa, Sa = eventbuild('173ab', 28.4)
    data173ab, Pe, Se = getamp(-91.37, strike_rang, dip_rang, rake_rang, [Pp, Sp])
    # print('Exit Angles: ', Pe, Se)
    data173ab = autofault(data173ab, 1.09, -1.21, -3.26, 0.024, 0.186)
    data173ab[:rank].to_csv('check173ab.csv', mode='a', header=False)
    predictamp(data173ab[:rank],-91.37,[Pp,Sp])

    #----------- S0235B ------------------
    print('235b')
    Pp, Sp, Pa, Sa = eventbuild('235b', 27)
    data235b, Pe, Se = getamp(-102.31, strike_rang, dip_rang, rake_rang, [Pp, Sp])
    # print('Exit Angles: ', Pe, Se)
    data235b = autofault(data235b, 0.360, -1.81, 0, 0.35, 0.105)
    data235b[:rank].to_csv('check235b.csv', mode = 'a', header=False)
    predictamp(data235b[:rank],-102.31,[Pp,Sp])

    #----------- S0325AB -------------------
    print('325ab')
    Pp, Sp, Pa, Sa = eventbuild('325a', 38.4)
    data325ab, Pe, Se = getamp(-60.38, strike_rang, dip_rang, rake_rang, [Pp, Sp])
    # print('Exit Angles: ', Pe, Se)
    data325ab = autofault(data325ab,  -1.35, -4.08, 0, 1.0, 1.2)
    data325ab[:rank].to_csv('check325ab.csv', mode = 'a', header=False)
    predictamp(data325ab[:rank],-60.38,[Pp,Sp])
    return

#-------EVENT OUTPUT-----------
# for mod in ['NewGudkova', 'TAYAK', 'MAAK', 'Combined']:
#     mars = TauPyModel(model=mod)
#     for depth in [25,35,45,55]:
#         #Gudkova Model
#         if mod == 'NewGudkova':
#             if depth <= 50 and depth > 42:
#                 Pvelz = 7.12500; Svelz = 4.00300    #rounded
#             elif depth <= 42 and depth > 21:
#                 Pvelz = 7.13900; Svelz = 4.01900
#             elif depth <=21 and depth > 16:
#                 Pvelz = 7.14300; Svelz = 4.02300
#             elif depth <= 16 and depth > 10:
#                 Pvelz = 7.15000; Svelz = 4.03000    #rounded
#
#         elif mod == 'TAYAK':
#             if depth <= 77 and depth > 10:
#                 Pvelz = 5.84666; Svelz = 3.28116
#             elif depth <= 10 and depth >1:
#                 Pvelz = 4.95225; Svelz = 2.78097
#
#         elif mod == 'MAAK':
#             if depth <= 68 and depth > 10:
#                 Pvelz = 5.94027; Svelz = 3.33676
#             elif depth <= 10 and depth > 1:
#                 Pvelz = 5.09729; Svelz = 2.86612
#
#         elif mod == 'Combined':
#             if depth <= 203 and depth > 50:
#                 Pvelz = 7.45400; Svelz = 4.21600
#             elif depth <= 50 and depth > 22:
#                 Pvelz = 7.12700; Svelz = 4.00200
#             elif depth <= 22 and depth > 8.6:
#                 Pvelz = 5.14700; Svelz = 2.73900
#             elif depth <= 8.6 and depth > 0:
#                 Pvelz = 3.50400; Svelz = 1.77100
#
#         eventoutput(depth,3)

#-------FORWARD CALC--------
for mod in ['NewGudkova', 'Combined']:
    mars = TauPyModel(model=mod)
    for depth in [35]:
        #Gudkova Model
        if mod == 'NewGudkova':
            if depth <= 50 and depth > 42:
                Pvelz = 7.12500; Svelz = 4.00300    #rounded
            elif depth <= 42 and depth > 21:
                Pvelz = 7.13900; Svelz = 4.01900
            elif depth <=21 and depth > 16:
                Pvelz = 7.14300; Svelz = 4.02300
            elif depth <= 16 and depth > 10:
                Pvelz = 7.15000; Svelz = 4.03000    #rounded

        elif mod == 'TAYAK':
            if depth <= 77 and depth > 10:
                Pvelz = 5.84666; Svelz = 3.28116
            elif depth <= 10 and depth >1:
                Pvelz = 4.95225; Svelz = 2.78097

        elif mod == 'MAAK':
            if depth <= 68 and depth > 10:
                Pvelz = 5.94027; Svelz = 3.33676
            elif depth <= 10 and depth > 1:
                Pvelz = 5.09729; Svelz = 2.86612

        elif mod == 'Combined':
            if depth <= 203 and depth > 50:
                Pvelz = 7.45400; Svelz = 4.21600
            elif depth <= 50 and depth > 22:
                Pvelz = 7.12700; Svelz = 4.00200
            elif depth <= 22 and depth > 8.6:
                Pvelz = 5.14700; Svelz = 2.73900
            elif depth <= 8.6 and depth > 0:
                Pvelz = 3.50400; Svelz = 1.77100

        forwardcalc(depth,mechanism_dict)
