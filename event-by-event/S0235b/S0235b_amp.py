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

def waveforms(start, end, adjtime):
    st_raw = client.get_waveforms("XB", "ELYSE", "02", "B*", start-adjtime, end+adjtime, attach_response=True)
    st_disp = st_raw.copy()
    st_disp.remove_response(output='DISP')
    st_disp.plot()
    return st_disp

def uvw2enz(st):
    if len(st) != 3:
       print('Stream does not contain 3 Traces')
       return st
    for trace in st:
        head = trace.stats
        channel = head.channel
        if channel == 'BHU': U = trace.data
        elif channel == 'BHV': V = trace.data
        elif channel == 'BHW': W = trace.data
        else:
            print('Trace.channel is not BHU, BHV, or BHW')
            return st

    d = np.radians(-30)
    aU = np.radians(135)
    aV = np.radians(15)
    aW = np.radians(255)

    A = np.array([[np.cos(d)*np.sin(aU), np.cos(d)*np.cos(aU),-np.sin(d)],
                  [np.cos(d)*np.sin(aV), np.cos(d)*np.cos(aV), -np.sin(d)],
                  [np.cos(d)*np.sin(aW), np.cos(d)*np.cos(aW), -np.sin(d)]])

    B = np.linalg.inv(A)
    E,N,Z = np.dot(B,(U,V,W))

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


# S0235b

print('S0235b')

begin = UTCDateTime('2019-07-26T12:19:18')  # p-wave arrival
end = UTCDateTime('2019-07-26T12:22:05')    # s-wave arrival

st_uvw = waveforms(begin, end, 600)
st_z12 = uvw2enz(st_uvw)

stf = st_z12.copy()
stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)

hhe = stf[0].data
hhn = stf[1].data
hhz = stf[2].data

# bAz = 72 -> 72+180 = 254
hhT,hhR = rotate(hhe,hhn,252)

streamRT = stf.copy()
streamRT[0].data = hhT
streamRT[1].data = hhR
streamRT[0].stats.component = 'T'
streamRT[1].stats.component = 'R'
streamRT.plot(equal_scale=True);


# In[10]:


stP_og = streamRT.slice(starttime=begin-20,endtime=begin+8)
stS_og = streamRT.slice(starttime=end-20, endtime=end+15)


# In[11]:


headerP = stP_og[0].stats
headerS = stS_og[0].stats

# model_ls = ['NewGudkova']
# model_Pangles = [25.5]
# model_Sangles = [23.2]

#NewGudkova
Mars = TauPyModel(model='NewGudkova')
mtimes = Mars.get_travel_times(source_depth_in_km = 35, distance_in_degree = 27.1, phase_list =['P','S'])
print(mtimes)
Pa = mtimes[0].incident_angle; Sa = mtimes[1].incident_angle
print('P incid angle: ', Pa)
print('S incid angle: ', Sa)

# model incidence angles
model_Pangles = [Pa]
model_Sangles = [Sa]


n = 0
for a in model_Pangles:
    print('P:')
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
    print('S:')
    stS = stS_og.copy()
    hhQ,hhL = rotate(stS[1].data, stS[2].data, a)
    t1, t2, t3 = Trace(stS[0].data, header=headerS), Trace(hhQ, header=headerS), Trace(hhL, header=headerS)
    stS_LQ = Stream(traces=[t1,t2,t3])
    stS_LQ[0].stats.component = 'T'
    stS_LQ[1].stats.component = 'Q'
    stS_LQ[2].stats.component = 'L'

    stS_LQ.plot(equal_scale=True);
    n += 1

###################
#     S0235b     #
###################
# P (BHL) [12:19:18.70] : 3.62e-10
# SV (BHQ) [12:22:05.75] : 3.47e-9
# SH (BHT) [12:22:05.21] : -1.611e-9
# P incid angle:  26.0108155787
# S incid angle:  23.4481631286
