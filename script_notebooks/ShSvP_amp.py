import matplotlib.pyplot as plt
import numpy as np
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
import obspy
from obspy.clients.fdsn import Client

client = Client("IRIS")

def waveforms(start, end, adjtime):
    st = client.get_waveforms("XB", "ELYSE", "02", "B*", start-adjtime, end+adjtime)
    st.detrend(type='simple')
    return st

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


# In[2]:


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


# In[3]:


import pandas as pd

# ia = pd.read_csv('incident_angles.csv')
# print(ia)
#
# model_name = ia['Model']
# depth = ia['Depth']


# #S0173a
#
# print('SO173a')
#
# begin = UTCDateTime('2019-05-23T02:22:59')  # p-wave arrival
# end = UTCDateTime('2019-05-23T02:25:53')    # s-wave arrival
#
# st_uvw = waveforms(begin, end, 600)
# st_z12 = uvw2enz(st_uvw)
#
# stf = st_z12.copy()
# stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)
#
#
# # In[5]:
#
#
# hhe = stf[0].data
# hhn = stf[1].data
# hhz = stf[2].data
#
# # bAz = 90 -> az = 90+180 = 270
# hhT,hhR = rotate(hhe,hhn,270)
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
# stP_og = streamRT.slice(starttime=begin-20,endtime=begin+20)
# stS_og = streamRT.slice(starttime=end-20, endtime=end+15)
#
# #stP_og.plot(equal_scale=True)
# #stS_og.plot(equal_scale=True)
#
#
# # In[7]:
#
#
# headerP = stP_og[0].stats
# headerS = stS_og[0].stats
#
#
# # In[8]:
#
#
# #P-wave
# # model_ls = ['DWAK', 'EH45Tcold', 'EH45TcoldCrust1b', 'NewGudkova', 'LFAK', 'MAAK', 'TAYAK']
# # # Gudkova -> P = 1.8, S = 2
# # model_Pangles = [27.9,56.6,19.4,25.9,30.5,27.5,26.7]
# # model_Sangles = [24,58,19.6,23.4,26.3,24,22.9]
#
# model_ls = ['NewGudkova']
# # Gudkova -> P = 1.8, S = 2
# model_Pangles = [25.9]
# model_Sangles = [23.4]

# n = 0
# for a in model_Pangles:
#     print('P:' + model_ls[n])
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
#     noise = stP_LQ.copy()
#     noise = noise.slice(starttime=UTCDateTime('2019-05-23T02:22:38'), endtime=UTCDateTime('2019-05-23T02:22:52'))
#     a = noise[2].data
#     plt.hist(a)
#     print(np.std(a))
#     plt.show()
#
# #S-wave
# n = 0
# for a in model_Sangles:
#     print('S:' + model_ls[n])
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
#
#     noise = stS_LQ.copy()
#     noise = noise.slice(starttime=UTCDateTime('2019-05-23T02:25:33'), endtime=UTCDateTime('2019-05-23T02:25:48'))
#     a = noise[0].data
#     b = noise[1].data
#     c = a + b
#     plt.hist(c)
#     print(np.std(c))
#     plt.show()
#
# #S0235b
#
# print('S0235b')
#
# begin = UTCDateTime('2019-07-26T12:19:18')  # p-wave arrival
# end = UTCDateTime('2019-07-26T12:22:05')    # s-wave arrival
#
# st_uvw = waveforms(begin, end, 600)
# st_z12 = uvw2enz(st_uvw)
#
# stf = st_z12.copy()
# stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)
#
# hhe = stf[0].data
# hhn = stf[1].data
# hhz = stf[2].data
#
# # bAz = 74 -> 74+180 = 254
# hhT,hhR = rotate(hhe,hhn,254)
#
# streamRT = stf.copy()
# streamRT[0].data = hhT
# streamRT[1].data = hhR
# streamRT[0].stats.component = 'T'
# streamRT[1].stats.component = 'R'
#
#
# # In[10]:
#
#
# stP_og = stf.slice(starttime=begin-20,endtime=begin+8)
# stS_og = stf.slice(starttime=end-20, endtime=end+15)
#
#
# # In[11]:
#
#
# headerP = stP_og[0].stats
# headerS = stS_og[0].stats
#
#
# # In[12]:
#
#
# #P-wave
# # model_ls = ['DWAK', 'EH45Tcold', 'EH45TcoldCrust1b', 'NewGudkova', 'LFAK', 'MAAK', 'TAYAK']
# # # Gudkova -> P, 1.9, S = 2
# # model_Pangles = [27.7,56.3,19.3,25.5,30.3,27.4,26.4]
# # model_Sangles = [23.9,57.8,18.9,23.2,26.2,23.9,22.8]
#
# model_ls = ['NewGudkova']
# # Gudkova -> P, 1.9, S = 2
# model_Pangles = [25.5]
# model_Sangles = [23.2]
#
# n = 0
# for a in model_Pangles:
#     print('P:' + model_ls[n])
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
#     noise = stP_LQ.copy()
#     noise = noise.slice(starttime=UTCDateTime('2019-07-26T12:18:58'), endtime=UTCDateTime('2019-07-26T12:19:15'))
#     a = noise[2].data
#     plt.hist(a)
#     print(np.std(a))
#     plt.show()
#
# #S-wave
# n = 0
# for a in model_Sangles:
#     print('S:' + model_ls[n])
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
#
#     noise = stS_LQ.copy()
#     noise = noise.slice(starttime=UTCDateTime('2019-07-26T12:21:45'), endtime=UTCDateTime('2019-07-26T12:21:59'))
#     a = noise[0].data
#     b = noise[1].data
#     c = a + b
#     plt.hist(c)
#     print(np.std(c))
#     plt.show()
#
# # S0325a
#
# print('S0325aa')
#
# begin = UTCDateTime('2019-10-26T06:58:57')  # p-wave arrival
# end = UTCDateTime('2019-10-26T07:02:46')    # s-wave arrival
#
# st_uvw = waveforms(begin, end, 600)
# st_z12 = uvw2enz(st_uvw)
#
# stf = st_z12.copy()
# stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)
# hhe = stf[0].data
# hhn = stf[1].data
# hhz = stf[2].data
#
# #bAz = 123 -> 123+180 = 303
# hhT,hhR = rotate(hhe,hhn,303)
#
# streamRT = stf.copy()
# streamRT[0].data = hhT
# streamRT[1].data = hhR
# streamRT[0].stats.component = 'T'
# streamRT[1].stats.component = 'R'
#
#
# # In[14]:
#
#
# stP_og = stf.slice(starttime=begin-20,endtime=begin+15)
# stS_og = stf.slice(starttime=end-20, endtime=end+25)
#
#
# # In[15]:
#
#
# headerP = stP_og[0].stats
# headerS = stS_og[0].stats
#
#
# # In[16]:
#
#
# #P-wave
# # model_ls = ['DWAK', 'EH45Tcold', 'EH45TcoldCrust1b', 'NewGudkova', 'LFAK', 'MAAK', 'TAYAK']
# # #Gudkova -> P =1.8, S=1.9
# # model_Pangles = [26.6,54.4,18.8,24.4,29.3,26.1,25.3]
# # model_Sangles = [23.3,57.1,18.7,22.4,25.9,23.5,22.2]
#
# model_ls = ['NewGudkova']
# #Gudkova -> P =1.8, S=1.9
# model_Pangles = [24.4]
# model_Sangles = [22.4]
#
# n = 0
# for a in model_Pangles:
#     print('P:' + model_ls[n])
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
#     noise = stP_LQ.copy()
#     noise = noise.slice(starttime=UTCDateTime('2019-10-26T06:58:37'), endtime=UTCDateTime('2019-10-26T06:58:54'))
#     a = noise[2].data
#     plt.hist(a)
#     print(np.std(a))
#     plt.show()
#
# #S-wave
# n = 0
# for a in model_Sangles:
#     print('S:' + model_ls[n])
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
#
#     noise = stS_LQ.copy()
#     noise = noise.slice(starttime=UTCDateTime('2019-10-26T07:02:26'), endtime=UTCDateTime('2019-10-26T07:02:46'))
#     a = noise[0].data
#     b = noise[1].data
#     c = a + b
#     plt.hist(c)
#     print(np.std(c))
#     plt.show()


print('S0325ab')

begin = UTCDateTime('2019-10-26T06:59:08')  # p-wave arrival
end = UTCDateTime('2019-10-26T07:03:00')    # s-wave arrival

st_uvw = waveforms(begin, end, 600)
st_z12 = uvw2enz(st_uvw)

stf = st_z12.copy()
stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)
hhe = stf[0].data
hhn = stf[1].data
hhz = stf[2].data

#bAz = 139 -> 139+180 = 319
hhT,hhR = rotate(hhe,hhn,319)

streamRT = stf.copy()
streamRT[0].data = hhT
streamRT[1].data = hhR
streamRT[0].stats.component = 'T'
streamRT[1].stats.component = 'R'


# In[14]:


stP_og = stf.slice(starttime=begin-20,endtime=begin+15)
stS_og = stf.slice(starttime=end-20, endtime=end+25)


# In[15]:


headerP = stP_og[0].stats
headerS = stS_og[0].stats


# In[16]:


#P-wave
model_ls = ['DWAK', 'EH45Tcold', 'EH45TcoldCrust1b', 'NewGudkova', 'LFAK', 'MAAK', 'TAYAK']
# Gudkova -> P=1.8, S=1.9
model_Pangles = [26.6,54.4,18.8,24.4,29.3,26.1,25.3]
model_Sangles = [23.3,57.1,18.7,22.4,25.9,23.5,22.2]

n = 0
for a in model_Pangles:
    print('P:' + model_ls[n])
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
    print('S:' + model_ls[n])
    stS = stS_og.copy()
    hhQ,hhL = rotate(stS[1].data, stS[2].data, a)
    t1, t2, t3 = Trace(stS[0].data, header=headerS), Trace(hhQ, header=headerS), Trace(hhL, header=headerS)
    stS_LQ = Stream(traces=[t1,t2,t3])
    stS_LQ[0].stats.component = 'T'
    stS_LQ[1].stats.component = 'Q'
    stS_LQ[2].stats.component = 'L'

    stS_LQ.plot(equal_scale=True);
    n += 1

# S0173ab

print('SO173ab')

begin = UTCDateTime('2019-05-23T02:23:03')  # p-wave arrival
end = UTCDateTime('2019-05-23T02:25:57')    # s-wave arrival

st_uvw = waveforms(begin, end, 600)
st_z12 = uvw2enz(st_uvw)

stf = st_z12.copy()
stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)


# In[5]:


hhe = stf[0].data
hhn = stf[1].data
hhz = stf[2].data

#bAz = 86 -> 86+180 = 266
hhT,hhR = rotate(hhe,hhn,266)

streamRT = stf.copy()
streamRT[0].data = hhT
streamRT[1].data = hhR
streamRT[2].data = hhz
streamRT[0].stats.component = 'T'
streamRT[1].stats.component = 'R'


# In[6]:


stP_og = streamRT.slice(starttime=begin-20,endtime=begin+20)
stS_og = streamRT.slice(starttime=end-20, endtime=end+15)

# stP_og.plot(equal_scale=True)
# stS_og.plot(equal_scale=True)


# In[7]:


headerP = stP_og[0].stats
headerS = stS_og[0].stats


# In[8]:


#P-wave
#no modeled incidence angles - 'EH45Tcold', 'EH45TcoldCrust1b',
model_ls = ['DWAK', 'NewGudkova', 'LFAK', 'MAAK', 'TAYAK']
#Gudkova -> P = 1.89, S = 1.92
model_Pangles = [27.9,25.9,30.5,27.5,26.7]
model_Sangles = [24.0,23.4,26.3,24.0,22.9]

n = 0
for a in model_Pangles:
    print('P:' + model_ls[n])
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
    print('S:' + model_ls[n])
    stS = stS_og.copy()
    hhQ,hhL = rotate(stS[1].data, stS[2].data, a)
    t1, t2, t3 = Trace(stS[0].data, header=headerS), Trace(hhQ, header=headerS), Trace(hhL, header=headerS)
    stS_LQ = Stream(traces=[t1,t2,t3])
    stS_LQ[0].stats.component = 'T'
    stS_LQ[1].stats.component = 'Q'
    stS_LQ[2].stats.component = 'L'

    stS_LQ.plot(equal_scale=True);
    n += 1


# # S0183a
#
# print('SO183a')
#
# begin = UTCDateTime('2019-06-03T02:27:49')  # p-wave arrival
# end = UTCDateTime('2019-06-03T02:32:15')    # s-wave arrival
#
# st_uvw = waveforms(begin, end, 600)
# st_z12 = uvw2enz(st_uvw)
#
# stf = st_z12.copy()
# stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)
#
#
# # In[5]:
#
#
# hhe = stf[0].data
# hhn = stf[1].data
# hhz = stf[2].data
#
# #
# hhT,hhR = rotate(hhe,hhn,285)
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
# stP_og = streamRT.slice(starttime=begin-8,endtime=begin+20)
# stS_og = streamRT.slice(starttime=end-9, endtime=end+15)
#
# #stP_og.plot(equal_scale=True)
# #stS_og.plot(equal_scale=True)
#
#
# # In[7]:
#
#
# headerP = stP_og[0].stats
# headerS = stS_og[0].stats
#
#
# # In[8]:
#
#
# #P-wave
# model_ls = ['DWAK', 'EH45Tcold', 'Gudkova', 'LFAK', 'MAAK', 'TAYAK']
# model_Pangles = [25.7,52.5,1.73,27.85,25.2,24.7]
# model_Sangles = [22.85,56.4,1.8,25.55,22.5,21.8]
#
# n = 0
# for a in model_Pangles:
#     print('P:' + model_ls[n])
#     stP = stP_og.copy()
#     hhQ,hhL = rotate(stP[1].data, stP[2].data, a)
#     t1, t2, t3 = Trace(stP[0].data, header=headerP), Trace(hhQ, header=headerP), Trace(hhL, header=headerP)
#     stP_LQ = Stream(traces=[t1,t2,t3])
#     stP_LQ[0].stats.component = 'T'
#     stP_LQ[1].stats.component = 'Q'
#     stP_LQ[2].stats.component = 'L'
#     stP_LQ.plot(equal_scale=True);
#
#     n += 1
#
# #S-wave
# n = 0
# for a in model_Sangles:
#     print('S:' + model_ls[n])
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
