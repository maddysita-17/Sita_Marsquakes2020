import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel

# run this once to allow TauPyModel(model="TAYAK")
#from obspy.taup.taup_create import build_taup_model
#mfl = 'models_nd_format/TAYAK.nd'
#build_taup_model(mfl)

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


# Earth testing
#earth = TauPyModel(model="iasp91")
#etimes = earth.get_travel_times(source_depth_in_km=55, distance_in_degree=40, phase_list=["P", "S"])
#print('EARTH')
#print(etimes)
#print(etimes[0].time, etimes[0].ray_param, etimes[0].incident_angle)

#erays = earth.get_ray_paths(source_depth_in_km=55, distance_in_degree=40, phase_list=["P", "S"])
#erays.plot_rays()

# move on to Mars:
radius = 3389.5
mars = TauPyModel(model="TAYAK")
mtimes = mars.get_travel_times(source_depth_in_km=55, distance_in_degree=40, phase_list=["P", "S"])
print('MARS')
print(mtimes)

print(mtimes[0].time, mtimes[0].ray_param, mtimes[0].incident_angle)
# testing...
#incident angle at the station
Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
Pvel = radius*np.sin(np.radians(Pa))/Pp
print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
Svel = radius*np.sin(np.radians(Sa))/Sp
print('Check: S surface velocity = ',Svel,' ?')

mrays = mars.get_ray_paths(source_depth_in_km=55, distance_in_degree=40, phase_list=["P", "S"])
mrays.plot_rays()

# source (strike, dip, rake) and P, S velocities at source depth.
# thrust fault thats roughly orientated like the faults on the CF
fault = [120, 45,90]; mt = getmt(fault)
# manual input from the models based on the depth of the marsquake (unknown so can be done later)
Pvelz = 5.84666; Svelz = 3.28116

#from obspy.imaging.beachball import beachball
#beachball(mt, size=200, linewidth=2, facecolor='b')


# configuration:
azimuth = 270 #from the marsquake to the station
# exit angles based on velocity at depth
iP = np.degrees(np.arcsin(Pvelz*Pp/radius))
jS = np.degrees(np.arcsin(Svelz*Sp/radius))
print('P exit at ',iP); print('S exit at ',jS)
P,SV,SH = Rpattern(fault,azimuth,[iP,jS])
print(P,SV,SH)

#from obspy.imaging.source import plot_radiation_pattern
#plot_radiation_pattern(mt,kind=['p_sphere', 'beachball','s_sphere', 's_quiver'], coordinate_system='RTP', p_sphere_direction='inwards', fig=None, show=True)
# loops over all possible exit angles and azimuths
