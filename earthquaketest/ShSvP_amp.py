import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
import obspy
from obspy.clients.fdsn import Client

client = Client("IRIS")

def waves_uganda(start, end, adjtime):
    st = client.get_waveforms("II", "MBAR", "10", "B*", start-adjtime, end+adjtime)
    st.detrend(type='simple')
    return st

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

def rot2enz(st):
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

    N,E = rotate(one, two, -14)

    head.channel = 'BHE'; trE = Trace(data=E, header=head)
    head.channel = 'BHN'; trN = Trace(data=N, header=head)
    head.channel = 'BHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ



print('Nebbi, Uganda')

# begin = UTCDateTime('2021-03-20T17:08:20')  # p-wave arrival
# end = UTCDateTime('2021-03-20T17:08:50')    # s-wave arrival

arrival = UTCDateTime('2021-03-20T17:07:27')
begin = arrival + 51.7      #p-wave arrival
end = arrival + 85.1        #s-wave arrival

st_ugan = waves_uganda(begin, end, 800)
st_ugan_enz = rot2enz(st_ugan)

stf = st_ugan_enz.copy()
stf.filter('bandpass', freqmin = 0.125, freqmax = 1.0, corners=4, zerophase=True)


# In[5]:


hhe = stf[0].data
hhn = stf[1].data
hhz = stf[2].data

# bAz = 5.97 -> az = 180+5.97 = 95.97
hhT,hhR = rotate(hhe,hhn,186)

streamRT = stf.copy()
streamRT[0].data = hhT
streamRT[1].data = hhR
streamRT[2].data = hhz
streamRT[0].stats.component = 'T'
streamRT[1].stats.component = 'R'


# In[6]:


stP_og = streamRT.slice(starttime=begin-20,endtime=begin+30)
stS_og = streamRT.slice(starttime=end-20, endtime=end+30)


st_Par = stf.slice(starttime=begin-2,endtime=begin)
scale = 1/200
hhe = scale * st_Par[0].data
hhn = scale * st_Par[1].data

tvall = []
alpha = np.arange(0,360,1)

# calculate Energy on channels oriented in the a direction, for all a in alpha:
# angle a is relative to orientation of channels 1 and 2.
# c1 is the x-axis, c2 is the y-axis
for a in alpha:
    hhT,hhR = rotate(hhe,hhn,a)
    Tenergy = np.dot(hhT,hhT)
    tvall.append(Tenergy)

tval = np.array(tvall)
mina = alpha[np.argmin(tval)]

mina_guess = tval[np.where(alpha == mina)]
mina_5 = tval[np.where(alpha == 6)]
print(mina_5, mina_guess)

#angle at which the energy is minimized
print('optimal angle = ',mina,' or ',mina-180)


# In[7]:


headerP = stP_og[0].stats
headerS = stS_og[0].stats


model_Pangles = [45.75]
model_Sangles = [48.27]

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
