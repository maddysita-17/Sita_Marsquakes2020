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

    print('exit angle', iP,jS)

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
    mtimes = mars.get_travel_times(source_depth_in_km= depth, distance_in_degree=dist, phase_list=['p','s'])

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

def autofault(df, obsP, obsSV, obsSH, errP, errSV, errSH):
    # # hypothetically observed amplitudes P, SV, SH:
    # vobserved = np.array([obsP, obsSV, obsSH])
    # vobslength = np.linalg.norm(vobserved)
    # # adding errors:
    # eobs = np.array([errP, errS, errS])
    # # normalized:
    # norm = vobserved/vobslength
    # # print('normalized array: ',norm)
    #
    # # defining cutoff value:
    # eca_ls = [eobs[0]*(1-norm[0]),eobs[1]*(1-norm[1]),eobs[2]*(1-norm[2])]/vobslength
    # eca = np.arctan(np.max(eca_ls))
    # print('cutoff value: ', eca)

    vobserved = np.array([obsP, obsSV, obsSH])
    vobslength = np.linalg.norm(vobserved)
    norm = vobserved/vobslength

    eobs = np.array([errP,errSV,errSH])
    eca_ls = np.sqrt(eobs[0]**2*(1-norm[0]**2) + \
                     eobs[1]**2*(1-norm[1]**2) + \
                     eobs[2]**2*(1-norm[2]**2))/vobslength
    eca = np.arctan(eca_ls)
    print('tolerance (cut-off value): ',eca, 'radians')

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
        extra.append([depth,mod,i,vca[i],norm])

    faults = {'Strike': st_ls,
                'Dip': dp_ls,
                'Rake': rk_ls,
                'Misfit': mf_ls,
                'Extra': extra}
    posfaults = pd.DataFrame.from_dict(faults)
    # print('before sort: ', posfaults)

    faults_sorted = posfaults.sort_values(by=['Misfit'])
    return faults_sorted

def predictamp(data,az,rayp,print_state=False):
    for index, rows in data.iterrows():
        ampsdf, iP, jS = getamp(az, [rows.Strike], [rows.Dip], [rows.Rake], rayp)
        # ampsdf.to_csv('amps_'+ str(az)+'.csv', mode='a', header=False)

        if print_state ==True:
            print(ampsdf)
        else:
            pass


# ------ FAULT PLANE LISTS -------      #-> 405,000 mechanisms
strike_rang = [*range(0, 180, 2)]
dip_rang = [*range(0,90,2)]
rake_rang = [*range(-180,180,2)]

#------EARTHQUAKE TEST-----------
# Earth:
radius = 6378.137
mod = 'iasp91'
mars = TauPyModel('iasp91')

# Velocity Model:
depth = 535
Pvelz = 9.6960; Svelz = 5.2820

# ----------- VANUATU ---------------
# Parameters:
# dist = 9.22; azm = 306.33; bAzm = 128.91
dist = 10.11; azm = -142.94; bAzm = 40.02

# Event Output:
print('Vanuatu')
Pp, Sp, Pa, Sa = eventbuild('Vanuatu', dist)
data_van, Pe, Se = getamp(azm, strike_rang, dip_rang, rake_rang, [Pp, Sp])
print('exit angles: ', Pe, Se)
data_van = autofault(data_van, obsP = 0, obsSV = -0.000383, obsSH = -0.000448, errP = 2.79e-07, errSV = 3.93e-06, errSH = 4.12e-06)
print(len(data_van))
print(data_van)
# data_van.to_csv('vanuatu_faults_au.csv')
# data_van.to_csv('vanuatu_faults.csv')
van_amps = predictamp(data_van,azm,[Pp,Sp],print_state=True)
print(van_amps)

# # Forward Calc:
# mechdic = {
#                     'Vanuatu': [[45,90,0],],
#                     # 'Vanuatu2': [[78,48,-68],],
#                     # 'Vanuatu3': [[80,46,-66],],
#                     }
# # ----------- UGANDA ---------------
# if 'Vanuatu' in mechdic:
#     vanu = mechdic['Vanuatu'][0]
#
#     print('Vanuatu')
#     Pp, Sp, Pa, Sa = eventbuild('Vanuatu', dist)
#     van_data, Pe, Se = getamp(azm, [vanu[0]], [vanu[1]], [vanu[2]], [Pp, Sp])
#     print(van_data)
#     # van_data.to_csv('fwdcalc-van.csv', mode='a',header=False)
#
#
# else:
#     pass
#
# if 'Vanuatu2' in mechdic:
#     vanu2 = mechdic['Vanuatu2'][0]
#
#     print('Vanuatu')
#     Pp, Sp, Pa, Sa = eventbuild('Vanuatu', dist)
#     van_data, Pe, Se = getamp(azm, [vanu2[0]], [vanu2[1]], [vanu2[2]], [Pp, Sp])
#     print(van_data)
#
# else:
#     pass
#
# if 'Vanuatu' in mechdic:
#     vanu3 = mechdic['Vanuatu3'][0]
#
#     print('Vanuatu')
#     Pp, Sp, Pa, Sa = eventbuild('Vanuatu', dist)
#     van_data, Pe, Se = getamp(azm, [vanu3[0]], [vanu3[1]], [vanu3[2]], [Pp, Sp])
#     print(van_data)
#
# else:
#     pass
