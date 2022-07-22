import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
import obspy
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel


client = Client("IRIS")


def rot2enz(st):
    if len(st) != 3:
       print('Stream does not contain 3 Traces')
       return st
    for trace in st:
        head = trace.stats
        channel = head.channel
        if channel == 'BH1': N = trace.data
        elif channel == 'BH2': E = trace.data
        elif channel == 'BHZ': Z = trace.data
        else:
            print('Trace.channel is not BH1, BH2, or BHZ')
            return st

    head.channel = 'BHE'; trE = Trace(data=E, header=head)
    head.channel = 'BHN'; trN = Trace(data=N, header=head)
    head.channel = 'BHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ

## rotation script for 'IU', 'COR', '00' ##
def rot2enz_COR(st):
    if len(st) != 3:
       print('Stream does not contain 3 Traces')
       return st
    for trace in st:
        head = trace.stats
        channel = head.channel
        if channel == 'BH1': one = trace.data
        elif channel == 'BH2': two = trace.data
        elif channel == 'BHZ': Z = trace.data
        else:
            print('Trace.channel is not BH1, BH2, or BHZ')
            return st

    E,N = rotate(two, one, -130)

    head.channel = 'BHE'; trE = Trace(data=E, header=head)
    head.channel = 'BHN'; trN = Trace(data=N, header=head)
    head.channel = 'BHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ

## rotation script for 'IU', 'ANMO', '00' ##
def rot2enz_ANMO(st):
    if len(st) != 3:
       print('Stream does not contain 3 Traces')
       return st
    for trace in st:
        head = trace.stats
        channel = head.channel
        if channel == 'BH1': one = trace.data
        elif channel == 'BH2': two = trace.data
        elif channel == 'BHZ': Z = trace.data
        else:
            print('Trace.channel is not BH1, BH2, or BHZ')
            return st

    E,N = rotate(two, one, -156)

    head.channel = 'BHE'; trE = Trace(data=E, header=head)
    head.channel = 'BHN'; trN = Trace(data=N, header=head)
    head.channel = 'BHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ


def rotate(c1,c2,a):
    """
    IN: c1,c2 (arrays) and a (angle)
    c1 c2 are the X and Y axes, respectively of a Cartesian coordinate system
    a is an angle in degrees, positive angle means a clockwise rotation of the coordinate system.
    OUT: o1, o2 (arrays)
    o1 o2 are the X and Y axes, respectively of a rotated Cartesian coordinate system
    """
    o1 = np.cos(np.radians(a))*c1 - np.sin(np.radians(a))*c2
    o2 = np.sin(np.radians(a))*c1 + np.cos(np.radians(a))*c2
    return o1, o2


#EARTHQUAKE PARAMETERS
#station bAz
bAz_COR = -123
bAz_ANMO = -97
# bAz_TEIG = 76.35 (INCORRECT)
bAz_stdep = bAz_ANMO +180

#station azm
azm_COR = 38.93
azm_ANMO = 59.22
# azm_TEIG = 76.35

#station dist
dist_COR = 37.18
dist_ANMO = 46.13
# dist_TEIG = 62.94
dist_stdep = dist_ANMO

print('EARTHQUAKE')

event_origin = UTCDateTime('2021-10-10T21:48:36.560')

# #STATION COR
# begin = event_origin + 427
# end = event_origin + 773
# arrival_diff = abs(begin-end)
# def waveforms(start, end, adjtime):
#     st_raw = client.get_waveforms('IU', 'COR', '00', "B*", start-adjtime, end+adjtime, attach_response=True)
#     st = st_raw.copy()
#     st_filt = st.filter('lowpass', freq = 0.1)
#     st_disp = st_filt.remove_response(output='DISP')
#     # st_disp.plot()
#     return st_disp
# st_z12 = waveforms(begin, end, 600)
# st_enz = rot2enz_COR(st_z12)

#STATION ANMO
begin = event_origin + 500
end = event_origin + 904
arrival_diff = abs(begin-end)
def waveforms(start, end, adjtime):
    st_raw = client.get_waveforms('IU', 'ANMO', '00', 'B*', start-adjtime, end+adjtime, attach_response=True)
    st = st_raw.copy()
    st_filt = st.filter('lowpass', freq = 0.1)
    st_disp = st_filt.remove_response(output='DISP')
    st_disp.plot()
    return st_disp
st_z12 = waveforms(begin, end, 600)
st_enz = rot2enz_ANMO(st_z12)

# #STATION TEIG
# begin = event_origin + 622
# end = event_origin + 1130
# arrival_diff = abs(begin-end)
# def waveforms(start, end, adjtime):
#     st_raw = client.get_waveforms('IU', 'TEIG', '00', 'B*', start-adjtime, end+adjtime, attach_response=True)
#     st = st_raw.copy()
#     st_filt = st.filter('bandpass', freqmin = 0.014, freqmax = 0.100)
#     st_disp = st_filt.remove_response(output='DISP')
#     st_disp.plot()
#     return st_disp
# st_z12 = waveforms(begin, end, 600)
# st_enz = rot2enz(st_z12)

stf = st_enz.copy()
stf.plot(equal_scale=True)

hhe = stf[0].data
hhn = stf[1].data
hhz = stf[2].data


# COR #
stf_P = stf.slice(UTCDateTime('2021-10-10T21:55:52.7')-5,UTCDateTime('2021-10-10T21:55:52.7')+5)
stf_P.plot()

tvall = []
alpha = np.arange(0,360,1)

# calculate Energy on channels oriented in the a direction, for all a in alpha:
# angle a is relative to orientation of channels 1 and 2.
# c1 is the x-axis, c2 is the y-axis
for a in alpha:
    hhT,hhR = rotate(stf_P[0],stf_P[1],a)
    Tenergy = np.dot(hhT,hhT)
    tvall.append(Tenergy)

tval = np.array(tvall)
mina = alpha[np.argmin(tval)]

#angle at which the energy is minimized
print('optimal angle = ',mina,' or ',mina-180)

mina_guess = tval[np.where(alpha == mina)]
mina_218 = tval[np.where(alpha == 218)]
print(mina_guess)
print(mina_218)

# bAz + 180 = rotation angle
# bAz = 90 -> az = 90+180 = 270

hhT,hhR = rotate(hhe,hhn,bAz_stdep)


streamRT = stf.copy()
streamRT[0].data = hhT
streamRT[1].data = hhR
streamRT[2].data = hhz
streamRT[0].stats.component = 'T'
streamRT[1].stats.component = 'R'
streamRT.plot(equal_scale=True);


# In[6]:


stP_og = streamRT.slice(starttime=begin-60,endtime=begin+60)
stS_og = streamRT.slice(starttime=end-60, endtime=end+60)

#stP_og.plot(equal_scale=True)
#stS_og.plot(equal_scale=True)


# In[7]:


headerP = stP_og[0].stats
headerS = stS_og[0].stats


Earth = TauPyModel(model='iasp91')
etimes = Earth.get_travel_times(source_depth_in_km = 35.1, distance_in_degree = dist_stdep, phase_list =['P','S'])
print(etimes)
Pa = etimes[0].incident_angle; Sa = etimes[1].incident_angle
Pe = etimes[0].takeoff_angle; Se = etimes[1].takeoff_angle
print('P incid angle: ', Pa)
print('S incid angle: ', Sa)
print('P exit angle: ', Pe)
print('S exit angle: ', Se)


print('P:')
stP = stP_og.copy()
hhQ,hhL = rotate(stP[1].data, stP[2].data, Pa)
t1, t2, t3 = Trace(stP[0].data, header=headerP), Trace(hhQ, header=headerP), Trace(hhL, header=headerP)
stP_LQ = Stream(traces=[t1,t2,t3])
stP_LQ[0].stats.component = 'T'
stP_LQ[1].stats.component = 'Q'
stP_LQ[2].stats.component = 'L'

stP_LQ.plot(equal_scale=True);


print('S:')
stS = stS_og.copy()
hhQ,hhL = rotate(stS[1].data, stS[2].data, Sa)
t1, t2, t3 = Trace(stS[0].data, header=headerS), Trace(hhQ, header=headerS), Trace(hhL, header=headerS)
stS_LQ = Stream(traces=[t1,t2,t3])
stS_LQ[0].stats.component = 'T'
stS_LQ[1].stats.component = 'Q'
stS_LQ[2].stats.component = 'L'

stS_LQ.plot(equal_scale=True);


sliced_P = stP_LQ.slice(begin-40,begin-20)
sliced_S = stS_LQ.slice(end-40,end-20)

std_l = sliced_P[0].std()
std_t = sliced_S[1].std()
std_q = sliced_S[2].std()

print('Noise on L (P) comp: ',std_l)
print('Noise on T (SH) comp: ', std_t)
print('Noise on Q (SV) comp: ', std_q)
Lh = plt.hist(sliced_P[0])
plt.show()
Qh = plt.hist(sliced_S[1])
plt.show()
Th =plt.hist(sliced_S[2])
plt.show()


##############################
#     IU COR OO - OREGON     #
##############################
# P (BHL) [21:55:52.7] : 1.54e-06
# SV (BHQ) [22:01:45.19] : -2.19e-06
# SH (BHT) [22:01:45.32] : -2.34e-06
# Pa : 26.23 & Sa : 27.32
# Pe : 38.04 & Se : 37.87

###################################
#     IU ANMO OO - NEW MEXICO     #
###################################
# P (BHL) [21:57:06.6] : 1.04e-06
# SV (BHQ) [22:03:59.2] : -1.55e-06
# SH (BHT) [22:04:00.2] : -2.74e-06
# Pa : 24.24 & Sa : 25.70
# Pe : 34.91 & Se : 35.45

###############################
#     IU TEIG OO - MEXICO     #
###############################
# P (BHL) [21:59:06.8] : 1.05e-07 (NOT UPDATED TO CORRECT ROTATION)
# SV (BHQ) [22:07:40.0] : -4.86e-07
# SH (BHT) [22:04:00.2] : 1.466e-06
# Pa : 20.31 & Sa : 22.23
# Pe : 28.93 & Se : 30.41
