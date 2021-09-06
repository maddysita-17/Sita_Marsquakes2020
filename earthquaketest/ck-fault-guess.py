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


#---get fault function--
def rotamp(st, dp, rk):
    """
    INPUT: az = azimuth in degrees from lat-long-az.py
    st, dp, rk = strike, dip, rake from 3 separte lists
    OUTPUT: df = data frame containing model outputs for P,SV,SH amplitudes at each depth input
    ip, ij = exit angles in degrees for each model computed in the Rpattern function
    """
    a_ls = []; P_ls = []; SH_ls =[]; SV_ls=[]

    for a in range(0, 361, 1):
        fault = [st, dp, rk]
        mt = getmt(fault)
        azimuth = a
        a_ls.append(a)

        iP = 70
        jS = 70

        P,SV,SH = Rpattern(fault,azimuth,[iP,jS])
        scalefactor = (Pvelz/Svelz)**3
        SV,SH = SV*scalefactor, SH*scalefactor
        P_ls.append(P); SH_ls.append(SH); SV_ls.append(SV)


    data = {'Az': a_ls,
            'P': P_ls,
            'SH': SH_ls,
            'SV': SV_ls,
        }

    df = pd.DataFrame(data, columns = ['Az', 'P', 'SV', 'SH'])
    return df, iP, jS

def rotext(az, st, dp, rk):
    """
    INPUT: az = azimuth in degrees from lat-long-az.py
    st, dp, rk = strike, dip, rake from 3 separte lists
    OUTPUT: df = data frame containing model outputs for P,SV,SH amplitudes at each depth input
    ip, ij = exit angles in degrees for each model computed in the Rpattern function
    """
    e_ls = []; P_ls = []; SH_ls =[]; SV_ls=[]

    for e in range(0, 91, 1):
        fault = [st, dp, rk]
        mt = getmt(fault)
        azimuth = az
        e_ls.append(e)

        iP = e
        jS = e

        P,SV,SH = Rpattern(fault,azimuth,[iP,jS])
        scalefactor = (Pvelz/Svelz)**3
        SV,SH = SV*scalefactor, SH*scalefactor
        P_ls.append(P); SH_ls.append(SH); SV_ls.append(SV)


    data = {'Exit': e_ls,
            'P': P_ls,
            'SH': SH_ls,
            'SV': SV_ls,
        }

    df = pd.DataFrame(data, columns = ['Exit', 'P', 'SV', 'SH'])
    return df


def eventbuild(dept, dis):
    etimes = model.get_travel_times(source_depth_in_km = dept, distance_in_degree = dis, phase_list=["P", "S"])
    print(etimes)

    #incident angle at the station
    Pp = etimes[0].ray_param; Pa = etimes[0].incident_angle
    Pvel = radius*np.sin(np.radians(Pa))/Pp
    print('Pvelz check: ', Pvel)

    Sp = etimes[1].ray_param; Sa = etimes[1].incident_angle
    Svel = radius*np.sin(np.radians(Sa))/Sp
    print('Svelz check: ', Svel)

    return Pp, Sp, Pa, Sa

strike = 180
dip = 90
rake = 0

# Earth:
radius = 6378.137

model = TauPyModel(model="iasp91")

#------velocity at 10km depth---------
depth = 10
Pvelz = 5.8000; Svelz = 3.3600

dist = 42

# Pp, Sp, Pa, Sa = eventbuild(depth, dist)
# dataf, iP, jS = rotamp(strike, dip, rake)
# print('exit P: ', iP)
# print('exit S: ', jS)
# print(dataf)
#
# fig = plt.figure()
# ax = fig.add_subplot()
# ax.scatter(dataf['Az'], dataf['P'], c = np.sign(dataf.P))
# ax.set_title('P amp')
# plt.show()
#
# fig = plt.figure()
# ax = fig.add_subplot()
# ax.set_title('SH amp')
# ax.scatter(dataf['Az'], dataf['SH'], c = np.sign(dataf.SH))
# plt.show()
#
# fig = plt.figure()
# ax = fig.add_subplot()
# ax.set_title('SV amp')
# ax.scatter(dataf['Az'], dataf['SV'], c = np.sign(dataf.SV))
# plt.show()

Pp, Sp, Pa, Sa = eventbuild(depth, dist)
dataf= rotext(137, strike, dip, rake)

fig = plt.figure()
ax = fig.add_subplot()
ax.scatter(dataf['Exit'], dataf['P'], c = np.sign(dataf.P))
ax.set_title('P amp')
plt.show()

fig = plt.figure()
ax = fig.add_subplot()
ax.set_title('SH amp')
ax.scatter(dataf['Exit'], dataf['SH'], c = np.sign(dataf.SH))
plt.show()

fig = plt.figure()
ax = fig.add_subplot()
ax.set_title('SV amp')
ax.scatter(dataf['Exit'], dataf['SV'], c = np.sign(dataf.SV))
plt.show()
