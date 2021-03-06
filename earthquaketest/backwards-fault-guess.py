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


def eventbuild(dept, dis):
    etimes = model.get_travel_times(source_depth_in_km = dept, distance_in_degree = dis, phase_list=["P", "S"])
    print(etimes)

    #incident angle at the station
    Pp = etimes[0].ray_param; Pa = etimes[0].incident_angle
    Pvel = radius*np.sin(np.radians(Pa))/Pp
    print('Pvelz check: ', Pvel)

    Sp = etimes[4].ray_param; Sa = etimes[4].incident_angle
    Svel = radius*np.sin(np.radians(Sa))/Sp
    print('Svelz check: ', Svel)

    return Pp, Sp, Pa, Sa

def autofault(df, obs_P, obs_SH, obs_SV, errP, errSV, errSH):
    # hypothetically observed amplitudes P, SV, SH
    vobserved = np.array([obs_P,obs_SV,obs_SH])
    vobslength = np.linalg.norm(vobserved)
    # and errors (always positive):
    eobs = np.array([errP, errSV, errSH])
    # normalized:
    n = vobserved/vobslength

    eca = np.arctan(np.max([eobs[0]*(1-n[0]),eobs[1]*(1-n[1]),eobs[2]*(1-n[2])])/vobslength)

    print('cutoff value: ', eca)

    #read in modeled csv file
    xd = df['P']
    yd = df['SV']
    zd = df['SH']
    ncalc = len(zd)

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
        mf3D[i] = np.arccos(np.dot(n,vca[i])/np.linalg.norm(vca[i]))   # should be valued between 0 and pi
        if mf3D[i] < eca:
            select.append(i)

    st_ls = []; dp_ls =[]; rk_ls=[]; mf_ls = []
    for i in select:
        st = df.at[i,'Strike'] ; dp = df.at[i,'Dip'] ; rk = df.at[i,'Rake']
        st_ls.append(st); dp_ls.append(dp); rk_ls.append(rk)
        mf_ls.append(mf3D[i])

    faults = {'Strike': st_ls,
                'Dip': dp_ls,
                'Rake': rk_ls,
                'Misfit': mf_ls}
    posfaults = pd.DataFrame.from_dict(faults)

    return posfaults

strike = [*range(0, 180, 5)]
dip = [*range(0,90,5)]
rake = [*range(-100,100,10)]

# Earth:
radius = 6378.137

model = TauPyModel(model="iasp91")

#------velocity at 10km depth---------
depth = 12
Pvelz = 5.8000; Svelz = 3.3600

#------event build at depth=10km-------
dist = 1.17; azm = -154; bAzm = 26

Pp, Sp, Pa, Sa = eventbuild(depth, dist)
dataf, iP, jS = getfault(azm, strike, dip, rake)
print('incid P: ', iP)
print('incid S: ', jS)
#print(dataf)

#0.637633, -2.322515, -2.959668
ugandadf = autofault(dataf, 0.637633, -2.322515, -2.959668, 0.2, 0.2, 0.2)
ugandadf.to_csv('backwards.csv', index=False)
