import matplotlib.pyplot as plot
import numpy as np
from obspy.taup import TauPyModel
import pandas as pd
from pathlib import Path
import math

def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

eps=1.e-6
rad=180./np.pi; halfpi = 0.5*np.pi; twopi = 2.0*np.pi

def azdp(v):
       vr=v[0]; vt=v[1]; vp=v[2]
       dparg = np.sqrt(vt*vt + vp*vp)
       if dparg>1.:
          dparg = 1.
          # print('argument error for np.arccos: ',dparg,'  Set to 1.')
       if dparg<-1.:
          dparg = -1.
          # print('argument error for np.arccos: ',dparg,'  Set to -1.')
       vdp = np.arccos(dparg)
       dp = halfpi - vdp
       vaz = halfpi + np.arctan2(vt,vp)
       if vr>0.: vaz = vaz + np.pi
       st = vaz + halfpi
       if st>=twopi: st = st - twopi
       if vaz>=twopi: vaz = vaz - twopi
       vaz = vaz*rad; st = st*rad; vdp = vdp*rad; dp = dp*rad
       return vaz,vdp,st,dp

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

def getplanes(xm):
    """
    needs function azdp. - converts MT to DC
    IN: xm = list of moment tensor elements in Harvard order (GCMT)
    OUT: strike dip rake (twice). Also: P and T vectors
    """
    xmatrix = np.array([[xm[0],xm[3],xm[4]],[xm[3],xm[1],xm[5]],[xm[4],xm[5],xm[2]]])
    tr=(xm[0]+xm[1]+xm[2])/3.
    if np.abs(tr)>eps:
       xmatrix[0,0] = xmatrix[0,0] - tr
       xmatrix[1,1] = xmatrix[1,1] - tr
       xmatrix[2,2] = xmatrix[2,2] - tr
    #print('removed isotropic component from Moment Tensor:')
    d,pt = np.linalg.eigh(xmatrix)
    jt = np.argmax(d) ; dmax = d[jt]
    jp = np.argmin(d) ; dmin = d[jp]
    for j in [0,1,2]:
        if j!=jp and j!=jt: jn=j
    if (jn+jp+jt)!=3:
        print('ERROR in axis determination')
        return 0

    p = pt[:,jp]
    t = pt[:,jt]
    n = pt[:,jn]
    if p[0] < 0.: p = -1.*p
    if t[0] < 0.: t = -1.*t
    pole1 = (t+p)/np.sqrt(2.)
    pole2 = (t-p)/np.sqrt(2.)
    if p[0] > t[0]: pole1 = -1.*pole1
    # planes' poles not part of function output, but they could be in future
    azt,dpt,st,dp = azdp(t)
    azn,dpn,st,dp = azdp(n)
    azp,dpp,st,dp = azdp(p)

    az1,dp1,st1,dip1 = azdp(pole1)
    az2,dp2,st2,dip2 = azdp(pole2)

    if -1.*d[jp]>d[jt]:
       djpt = d[jp]
    else:
       djpt = d[jt]
    clvd = d[jn]/djpt

    m0 = 0.5*(np.abs(d[jp])+np.abs(d[jt]))

    x = np.array([0.,-1*np.cos(st1/rad),np.sin(st1/rad)])
    vfin = np.dot(pole2,x)
    if vfin>1.:
       vfin = 1.
    if vfin<-1.:
       vfin = -1.
    rake1 = rad*np.arccos(vfin)
    if pole2[0]<0.: rake1 = -1.*rake1


    x = np.array([0.,-1*np.cos(st2/rad),np.sin(st2/rad)])
    vfin = np.dot(pole1,x)
    if vfin>1.:
       vfin = 1.
    if vfin<-1.:
       vfin = -1.
    rake2 = rad*np.arccos(vfin)
    if pole1[0]<0.: rake2 = -1.*rake2
    return 3*tr,clvd, m0,(azt,dpt),(azn,dpn),(azp, dpp), (st1,dip1,rake1), (st2,dip2,rake2)


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

def getamp(azimuth, strike, dip, rake, rayp, exit):
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
                Pvelz = 8.0400; Svelz = 4.4700
                # iP_1 = np.degrees(np.arcsin(Pvelz*rayp[0]/(radius-depth)))
                # jS_1 = np.degrees(np.arcsin(Svelz*rayp[1]/(radius-depth)))
                iP_2 = exit[0]
                jS_2 = exit[1]

                # print('method 1:', iP_1, jS_1)
                # print('method 2:', iP_2, jS_2)

                # calculating amplitudes
                P,iSV,iSH = Rpattern(fault, azimuth, [iP_2, jS_2])
                # scalefactor = (Pvelz/Svelz)**3
                # SV,SH = iSV*scalefactor, iSH*scalefactor
                SV,SH = iSV, iSH
                P_ls.append(P); SH_ls.append(SH); SV_ls.append(SV)

    # creating dataframe
    data = {
            'Depth': depth,
            'Strike': strike_ls,
            'Dip': dip_ls,
            'Rake': rake_ls,
            'P': P_ls,
            'SV': SV_ls,
            'SH': SH_ls
            }

    df = pd.DataFrame(data, columns = ['Model','Depth','Strike', 'Dip', 'Rake', 'P', 'SV', 'SH'])
    return df

def autofault(df, obsP, obsSV, obsSH, errP, errSV, errSH):
    # vobserved = np.array([obsP, obsSV, obsSH])
    vobserved = np.array([5*obsP, obsSV, obsSH])
    vobslength = np.linalg.norm(vobserved)
    norm = vobserved/vobslength

    # eobs = np.array([errP,errSV,errSH])
    eobs = np.array([5*errP,errSV,errSH])
    eca_ls = np.sqrt((eobs[0]**2*(1-norm[0]**2) + \
                     eobs[1]**2*(1-norm[1]**2) + \
                     eobs[2]**2*(1-norm[2]**2))/3)/vobslength
    eca = np.arctan(eca_ls)
    print('tolerance (cut-off value): ',eca,' radians')

    # xd = df['P']
    xd = 5*df['P']
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
    azt_ls=[]; dpt_ls=[]; azn_ls = []; dpn_ls = []; azp_ls = []; dpp_ls = []
    for i in select:
        st = df.at[i,'Strike']
        dp = df.at[i,'Dip']
        rk = df.at[i,'Rake']
        st_ls.append(st); dp_ls.append(dp); rk_ls.append(rk)

        mt = getmt([st,dp,rk])
        dum1,dum2,dum3,(azt,dpt),(azn,dpn),(azp, dpp), (st1,dip1,rake1), (st2,dip2,rake2) = getplanes(mt)
        azt_ls.append(truncate(azt,2)); dpt_ls.append(dpt)
        azn_ls.append(azn); dpn_ls.append(dpn)
        azp_ls.append(azp); dpp_ls.append(dpp)

        mf_ls.append(mf3d[i])
        extra.append([i,vca[i],norm])

    faults = {'AZT': azt_ls,
        'DPT': dpt_ls,
        'AZN': azn_ls,
        'DPN': dpn_ls,
        'AZP': azp_ls,
        'DPP': dpp_ls,
        'Strike': st_ls,
        'Dip': dp_ls,
        'Rake': rk_ls,
        'Misfit': mf_ls,
        'Extra': extra}
    posfaults = pd.DataFrame.from_dict(faults)
    # print('before sort: ', posfaults)

    faults_sorted = posfaults.sort_values(by=['Misfit'])
    print(len(faults_sorted))
    faults_sorted = faults_sorted.drop_duplicates(subset = ['Strike', 'Dip', 'Rake'])
    print(len(faults_sorted))
    return faults_sorted

def predictamp(data,az,rayp,print_state=False):
    for index, rows in data.iterrows():
        ampsdf, iP, jS = getamp(az, [rows.Strike], [rows.Dip], [rows.Rake], rayp)
        # ampsdf.to_csv('amps_'+ str(az)+'.csv', mode='a', header=False)

        if print_state ==True:
            print(ampsdf)
        else:
            pass


# # ---- FORWARD CALC (FROM MECH) -----
#
# mechanism_dict = {
#                     'S0173a': [[2,24,86],           # NewGud Top Mech @ 35km
#                                 [55,88,-105]],      # InSight mech from Table 2
#                     }
#
# def forwardcalc(depth,mechdic,alt=False):
#     # ----------- S0173A ---------------
#     if 'S0173a' in mechdic:
#         s0173a = mechdic['S0173a'][0]
#         if alt == True:
#             s0173a_b = mechdic['S0173a'][1]
#
#         print('173a')
#         Pp, Sp, Pa, Sa = eventbuild('173a', 28.4)
#         data173a, Pe, Se = getamp(-87.86, [s0173a[0]], [s0173a[1]], [s0173a[2]], [Pp, Sp])
#
#         print(data173a)
#         for i in range(len(data173a)):
#             data173a.at[i,'Model'] = mod
#         data173a.to_csv(path+'/fwdcalc-173a.csv', mode='a',header=False)
#
#         if alt == True:
#             data1, Pe, Se = getamp(-87.86, [s0173a_b[0]], [s0173a_b[1]], [s0173a_b[2]], [Pp, Sp])
#             # data173a = pd.concat([data173a,data1],ignore_index=True)
#
#             print(data1)
#             for i in range(len(data1)):
#                 data1.at[i,'Model'] = 'InSight'
#             data1.to_csv(path+'/fwdcalc-173a.csv', mode='a',header=False)
#     else:
#         pass
#
#     return

#EARTHQUAKE PARAMETERS
#station bAz
bAz_COR = -123
bAz_ANMO = -97
# bAz_TEIG = 76.35 (INCORRECT)
bAz_stdep = bAz_COR +180

#station azm
azm_COR = 38.93
azm_ANMO = 59.22
# azm_TEIG = 76.35

#station dist
dist_COR = 37.18
dist_ANMO = 46.13
# dist_TEIG = 62.94
dist_stdep = dist_COR

#------EARTHQUAKE TEST-----------
# Earth:
radius = 6378.137
Earth = TauPyModel(model='iasp91')
etimes = Earth.get_travel_times(source_depth_in_km = 35.1, distance_in_degree = dist_stdep, phase_list =['P','S'])
print(etimes)
Pp = etimes[0].ray_param; Sp = etimes[1].ray_param
Pa = etimes[0].incident_angle; Sa = etimes[1].incident_angle
Pe = etimes[0].takeoff_angle; Se = etimes[1].takeoff_angle

depth = 35.1

path = '/Users/maddysita/Desktop/CIERA_REU/event-by-event/EARTHQUAKE/'


# ------ FAULT PLANE LISTS -------      #-> 405,000 mechanisms
strike_rang = [*range(0, 360,2)]
dip_rang = [*range(0,90,2)]
rake_rang = [*range(-180,180,2)]


# print('IU:COR')
# data_COR = getamp(azm_COR,strike_rang, dip_rang, rake_rang, [Pp,Sp], [Pe,Se])
# data_COR = autofault(data_COR, obsP = 1.54e-06, obsSV = -2.19e-06, obsSH = -2.34e-06, errP = 1.01e-07, errSV =  2.8e-07, errSH =  2.76e-07)
# print(len(data_COR))
# data_COR.to_csv(path + 'eqCOR_march.csv')
# fwd_COR = getamp(azm_COR, [294],[68],[-34], [Pp,Sp], [Pe,Se])       #top fitting mech from march
# print(fwd_COR)
# fwd_COR = getamp(azm_COR, [138],[24],[-38], [Pp,Sp], [Pe,Se])       #second best fitting mech from march
# print(fwd_COR)


print('IU:ANMO')
# data_ANMO = getamp(azm_ANMO,strike_rang, dip_rang, rake_rang, [Pp,Sp], [Pe,Se])
# data_ANMO = autofault(data_ANMO, obsP = 1.04e-06, obsSV = -1.55e-06, obsSH = -2.74e-06, errP = 7.68e-08, errSV =  1.41e-07, errSH =  7.12e-07)
# print(len(data_ANMO))
# data_ANMO.to_csv(path + 'eqANMO_march.csv')
fwd_ANMO = getamp(azm_ANMO, [154],[24],[-28], [Pp,Sp], [Pe,Se])       #top fitting mech from march
print(fwd_ANMO)


# print('IU:TEIG')
# data_TEIG = getamp(azm_TEIG,strike_rang, dip_rang, rake_rang, [Pp,Sp], [Pe,Se])
# data_TEIG = autofault(data_TEIG, obsP = 8.98e-07, obsSV = -2.31e-06, obsSH = 2.10e-06, errP = 4.93e-08, errSV =  2.0e-07, errSH =  1.69e-07)
# print(len(data_TEIG))
# data_TEIG.to_csv(path + 'eqTEIG_wMT.csv')
