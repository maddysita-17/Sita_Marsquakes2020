import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
import obspy
from obspy.clients.fdsn import Client

client = Client("IRIS")

def waves_cr(start, end, adjtime):
    st_raw = client.get_waveforms("II", "JTS", "10", "LH*", start-adjtime, end+adjtime, attach_response=True)
    #st.detrend(type='simple')
    st_disp = st_raw.copy()
    st_disp.remove_response(output='DISP')
    st_disp.plot()
    return st_disp

def waves_g(start, end, adjtime):
    st_raw = client.get_waveforms("IU", "SFJD", "10", "BH*", start-adjtime, end+adjtime, attach_response=True)
    #st.detrend(type='simple')
    st_disp = st_raw.copy()
    st_disp.remove_response(output='DISP')
    st_disp.plot()
    return st_disp

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

def rot_cr(st):
    if len(st) != 3:
       print('Stream does not contain 3 Traces')
       return st
    for trace in st:
        head = trace.stats
        channel = head.channel
        if channel == 'LH1': one = trace.data
        elif channel == 'LH2': two = trace.data
        elif channel == 'LHZ': Z = trace.data
        else:
            print('Trace.channel is not BH1, BH2, or BHZ')
            return st

    E,N = rotate(two, one, 62)

    head.channel = 'LHE'; trE = Trace(data=E, header=head)
    head.channel = 'LHN'; trN = Trace(data=N, header=head)
    head.channel = 'LHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ

def rot_g(st):
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


# print('CostaRica')
# arrival = UTCDateTime('2020-04-01T00:00:00')
# begin = arrival             #p-wave arrival
# end = arrival + 420        #s-wave arrival
#
# st_ugan = waves_cr(begin, end, 800)
# st_ugan_enz = rot_cr(st_ugan)
#
# stf = st_ugan_enz.copy()
# #stf.filter('bandpass', freqmin = 0.01, freqmax = 1, corners=4, zerophase=True)
#
# stf.plot()
#
#
# hhe = stf[0].data
# hhn = stf[1].data
# hhz = stf[2].data
#
# # bAz = -31.77 -> az = 180+-31.77 = 148
# hhT,hhR = rotate(hhe,hhn,148)
#
# streamRT = stf.copy()
# streamRT[0].data = hhT
# streamRT[1].data = hhR
# streamRT[2].data = hhz
# streamRT[0].stats.component = 'T'
# streamRT[1].stats.component = 'R'
#
#
# # In[6]:
#
#
# stP_og = streamRT.slice(starttime=begin-20,endtime=begin+80)
# stS_og = streamRT.slice(starttime=end-50, endtime=end+80)
#
# stP_og.plot()
#
#
# headerP = stP_og[0].stats
# headerS = stS_og[0].stats
#
#
# # model_Pangles = [62.68]
# # model_Sangles = [63.14]
#
# model_Pangles = [25]
# model_Sangles = [26]
#
# n = 0
# for a in model_Pangles:
#     print('P')
#     stP = stP_og.copy()
#     hhQ,hhL = rotate(stP[1].data, stP[2].data, a)
#     t1, t2, t3 = Trace(stP[0].data, header=headerP), Trace(hhQ, header=headerP), Trace(hhL, header=headerP)
#     stP_LQ = Stream(traces=[t1,t2,t3])
#     stP_LQ[0].stats.component = 'T'
#     stP_LQ[1].stats.component = 'Q'
#     stP_LQ[2].stats.component = 'L'
#
#     stP_LQ.plot(equal_scale=True);
#     n += 1
#
# #S-wave
# n = 0
# for a in model_Sangles:
#     print('S')
#     stS = stS_og.copy()
#     hhQ,hhL = rotate(stS[1].data, stS[2].data, a)
#     t1, t2, t3 = Trace(stS[0].data, header=headerS), Trace(hhQ, header=headerS), Trace(hhL, header=headerS)
#     stS_LQ = Stream(traces=[t1,t2,t3])
#     stS_LQ[0].stats.component = 'T'
#     stS_LQ[1].stats.component = 'Q'
#     stS_LQ[2].stats.component = 'L'
#
#     stS_LQ.plot(equal_scale=True);
#     n += 1


print('Greenland')
arrival = UTCDateTime('2020-03-31T23:52:30')
begin = arrival + 420      #p-wave arrival
end = arrival + 820        #s-wave arrival

st_ugan = waves_g(begin, end, 800)
st_ugan_enz = rot_g(st_ugan)

stf = st_ugan_enz.copy()
#stf.filter('bandpass', freqmin = 0.01, freqmax = 1, corners=4, zerophase=True)

stf.plot()


hhe = stf[0].data
hhn = stf[1].data
hhz = stf[2].data

# bAz = -90.08 -> az = 180+-90.08 = 90
hhT,hhR = rotate(hhe,hhn,90)

streamRT = stf.copy()
streamRT[0].data = hhT
streamRT[1].data = hhR
streamRT[2].data = hhz
streamRT[0].stats.component = 'T'
streamRT[1].stats.component = 'R'


# In[6]:


stP_og = streamRT.slice(starttime=begin-20,endtime=begin+80)
stS_og = streamRT.slice(starttime=end-50, endtime=end+80)

stP_og.plot()


headerP = stP_og[0].stats
headerS = stS_og[0].stats


model_Pangles = [25]
model_Sangles = [26]

n = 0
for a in model_Pangles:
    print('P')
    stP = stP_og.copy()
    hhQ,hhL = rotate(stP[1].data, stP[2].data, a)
    t1, t2, t3 = Trace(stP[0].data, header=headerP), Trace(hhQ, header=headerP), Trace(hhL, header=headerP)
    stP_LQ = Stream(traces=[t1,t2,t3])
    stP_LQ[0].stats.component = 'T'
    stP_LQ[1].stats.component = 'Q'
    stP_LQ[2].stats.component = 'L'

    stP_LQ.plot(equal_scale=True);
    n += 1

#S-wave
n = 0
for a in model_Sangles:
    print('S')
    stS = stS_og.copy()
    hhQ,hhL = rotate(stS[1].data, stS[2].data, a)
    t1, t2, t3 = Trace(stS[0].data, header=headerS), Trace(hhQ, header=headerS), Trace(hhL, header=headerS)
    stS_LQ = Stream(traces=[t1,t2,t3])
    stS_LQ[0].stats.component = 'T'
    stS_LQ[1].stats.component = 'Q'
    stS_LQ[2].stats.component = 'L'

    stS_LQ.plot(equal_scale=True);
    n += 1
