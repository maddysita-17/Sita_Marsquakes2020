import matplotlib.pyplot as plt
import numpy as np
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

client = Client("IRIS")

#waveform functions to make importing data easier
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

    A = np.array([[np.cos(d)*np.sin(aU),np.cos(d)*np.cos(aU),np.sin(d)],
                  [np.cos(d)*np.sin(aV), np.cos(d)*np.cos(aV), np.sin(d)],
                  [np.cos(d)*np.sin(aW), np.cos(d)*np.cos(aW), np.sin(d)]])

    B = np.linalg.inv(A)
    E,N,Z = np.dot(B,(U,V,W))

    head.channel = 'BHE'; trE = Trace(data=E, header=head)
    head.channel = 'BHN'; trN = Trace(data=N, header=head)
    head.channel = 'BHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ

def wfplot(st,axs,iax,scale=1,sh=0,Ptime=0):
    if Ptime == 0: Ptime = st[0].stats.starttime
    time_axis = st[0].times(reftime=Ptime)
    n = min([len(time_axis)]+[len(st[i].data) for i in np.arange(len(st))])

    a = 1/200; ylim = [-2,10]
    shift = sh; dsh = a*600
    for trace in st:
        axs[iax][0].plot(time_axis[:n], a*trace.data[:n] + shift, label = trace.stats.channel)
        shift += dsh
    axs[iax][0].set_ylim(ylim)
    sENZ = uvw2enz(st)
    shift = sh
    for trace in sENZ:
        axs[iax][1].plot(time_axis[:n], a*trace.data[:n] + shift, label = trace.stats.channel)
        shift += dsh
    axs[iax][1].set_ylim(ylim)
    if iax == len(axs)-1: 
       axs[len(axs)-1][0].legend(loc='lower right')  
       axs[len(axs)-1][1].legend(loc='lower right') 


#    ax[row,column].plot(t, y, label='BHU', color = "#1f77b4", alpha = 0.7)
#    ax[row,column].plot(t, y, label='BHV', color = "#ff7f0e", alpha = 0.7)
#    ax[row,column].plot(t, y, label='BHW', color = "#2ca02c", alpha = 0.7)
#    ax[row,column].plot(t, y, label='E', color = "#77e59b", alpha = 0.8)
#    ax[row,column].plot(t, y, label='N', color = "#ffabab", alpha = 0.8)
#    ax[row,column].plot(t, y, label = 'Z', color = "#b28dff", alpha = 0.8)

#---------------------------------------------------

# set filters
bpfilters = [
[2,8],
[1.0,4],
[0.125,1.0],
[0.03,0.125]
]

# start figure
nrows = len(bpfilters)+3
fig, axs = plt.subplots(nrows, 2, figsize=(14,7))
plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97, hspace=0.3)
iax = 0

# Set times and get data:
e = 'Event S0173a'
### 2019-05-23T02:19:33
P173a = UTCDateTime('2019-05-23T02:22:59')
S173a = UTCDateTime('2019-05-23T02:25:53')
start173a = P173a - 1*60 # 1 minute early
end173a = S173a + 3*(S173a-P173a)
s0173a = waveforms(start173a, end173a, 600)

sps = s0173a[0].stats.sampling_rate; dt = 1./sps
if 2*bpfilters[0][1] > sps: 
    print('WARNING: high bandpass corner too high. f_nyq = ',0.5*sps,' Hz')
shift = 0; scale = 1

# plot raw data
stmp = s0173a.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stmp,axs,iax,scale,shift,P173a)
iax += 1

# plot high-passed data
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stm,axs,iax,scale,shift,P173a)
iax += 1

# plot band-passed data
for f in bpfilters:
    stmp = s0173a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
    stm = stmp.slice(starttime=start173a,endtime=end173a)
    #time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P173a)
    iax += 1

# plot low-passed data
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stm,axs,iax,scale,shift,P173a)

iax = 0

e = 'Event S0235b'
### 2019-07-26T12:16:15
P235b = UTCDateTime('2019-07-26T12:19:19')
S235b = UTCDateTime('2019-07-26T12:22:05')
start235b = P235b - 1*60 # 1 minute early
end235b = S235b + 3*(S235b-P235b)
s0235b = waveforms(start235b, end235b, 600)

shift = 1; scale = 2

# plot raw data
stmp = s0235b.slice(starttime=start235b,endtime=end235b)
#time_axis = np.arange(start235b-P235b,end235b-P235b+dt,dt)
wfplot(stmp,axs,iax,scale,shift,P235b)
iax += 1

# plot high-passed data
stmp = s0235b.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
stm = stmp.slice(starttime=start235b,endtime=end235b)
#time_axis = np.arange(start235b-P235b,end235b-P235b+dt,dt)
wfplot(stm,axs,iax,scale,shift,P235b)
iax += 1

# plot band-passed data
for f in bpfilters:
    stmp = s0235b.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
    stm = stmp.slice(starttime=start235b,endtime=end235b)
    #time_axis = np.arange(start235b-P235b,end235b-P235b+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P235b)
    iax += 1

# plot low-passed data
stmp = s0235b.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
stm = stmp.slice(starttime=start235b,endtime=end235b)
#time_axis = np.arange(start235b-P235b,end235b-P235b+dt,dt)
wfplot(stm,axs,iax,scale,shift,P235b)

fig.savefig('TwoQuakes.png')
plt.show()


# ....plot particle motion...
e = 'Event S0173a'
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('bandpass',freqmin=0.125, freqmax=1.0 ,corners=4, zerophase=True)
#bt = - 4; et = 12
bt = - 4; et = 8
#bt = 8; et = 12

fjg, ays = plt.subplots(2,3,figsize=(10,7))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.5)

stm = stmp.slice(starttime=P173a+bt,endtime=P173a+et)
sENZ = uvw2enz(stm)
tvec = sENZ[0].times(reftime=P173a)

x = sENZ.select(component='N')[0].data
y = sENZ.select(component='Z')[0].data
ays[0][0].scatter(x,y,c=tvec)
ays[0][0].set_ylim([-500,500]); ays[0][0].set_xlim([-500,500]); ays[0][0].set_xlabel('North'); ays[0][0].set_ylabel('Up')
ays[0][0].set_aspect('equal'); ays[0][0].set_title(e+' P-wave Particle Motion')

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='Z')[0].data
ays[0][1].scatter(x,y,c=tvec)
ays[0][1].set_ylim([-500,500]); ays[0][1].set_xlim([-500,500]); ays[0][1].set_xlabel('East'); ays[0][1].set_ylabel('Up')
ays[0][1].set_aspect('equal')

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='N')[0].data
ays[0][2].scatter(x,y,c=tvec)
ays[0][2].set_ylim([-500,500]); ays[0][2].set_xlim([-500,500]); ays[0][2].set_xlabel('East'); ays[0][2].set_ylabel('North')
ays[0][2].set_aspect('equal')


#bt = 172; et = 188
bt = 172; et = 178
#bt = 178; et = 188

stm = stmp.slice(starttime = P173a+bt, endtime = P173a+et)
sENZ = uvw2enz(stm)
tvec = sENZ[0].times(reftime=P173a)

x = sENZ.select(component='N')[0].data
y = sENZ.select(component='Z')[0].data
ays[1][0].scatter(x,y,c=tvec)
ays[1][0].set_ylim([-500,500]); ays[1][0].set_xlim([-500,500]); ays[1][0].set_xlabel('North'); ays[1][0].set_ylabel('Up')
ays[1][0].set_aspect('equal'); ays[1][0].set_title(e+' S-wave Particle Motion')

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='Z')[0].data
ays[1][1].scatter(x,y,c=tvec)
ays[1][1].set_ylim([-500,500]); ays[1][1].set_xlim([-500,500]); ays[1][1].set_xlabel('East'); ays[1][1].set_ylabel('Up')
ays[1][1].set_aspect('equal') 

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='N')[0].data
ays[1][2].scatter(x,y,c=tvec)
ays[1][2].set_ylim([-500,500]); ays[1][2].set_xlim([-500,500]); ays[1][2].set_xlabel('East'); ays[1][2].set_ylabel('North')
ays[1][2].set_aspect('equal')

fjg.savefig('ppm173a.png')
plt.show()

# ---------------------------------

e = 'Event S0235b'
stmp = s0235b.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('bandpass',freqmin=0.125, freqmax=1.0 ,corners=4, zerophase=True)
#bt = - 4; et = 12
bt = - 4; et = 8
#bt = 8; et = 12

fkg, azs = plt.subplots(2,3,figsize=(10,7))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.5)

stm = stmp.slice(starttime=P235b+bt,endtime=P235b+et)
sENZ = uvw2enz(stm)
tvec = sENZ[0].times(reftime=P235b)

x = sENZ.select(component='N')[0].data
y = sENZ.select(component='Z')[0].data
azs[0][0].scatter(x,y,c=tvec)
azs[0][0].set_ylim([-500,500]); azs[0][0].set_xlim([-500,500]); azs[0][0].set_xlabel('North'); azs[0][0].set_ylabel('Up')
azs[0][0].set_aspect('equal'); azs[0][0].set_title(e+' P-wave Particle Motion')

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='Z')[0].data
azs[0][1].scatter(x,y,c=tvec)
azs[0][1].set_ylim([-500,500]); azs[0][1].set_xlim([-500,500]); azs[0][1].set_xlabel('East'); azs[0][1].set_ylabel('Up')
azs[0][1].set_aspect('equal')

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='N')[0].data
azs[0][2].scatter(x,y,c=tvec)
azs[0][2].set_ylim([-500,500]); azs[0][2].set_xlim([-500,500]); azs[0][2].set_xlabel('East'); azs[0][2].set_ylabel('North')
azs[0][2].set_aspect('equal')

bt = 164; et = 180
stm = stmp.slice(starttime=P235b+bt,endtime=P235b+et)
sENZ = uvw2enz(stm)
tvec = np.linspace(bt,et,len(stm[0].data))

x = sENZ.select(component='N')[0].data
y = sENZ.select(component='Z')[0].data
azs[1][0].scatter(x,y,c=tvec)
azs[1][0].set_ylim([-500,500]); azs[1][0].set_xlim([-500,500]); azs[1][0].set_xlabel('North'); azs[1][0].set_ylabel('Up')
azs[1][0].set_aspect('equal'); azs[1][0].set_title(e+' S-wave Particle Motion')

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='Z')[0].data
azs[1][1].scatter(x,y,c=tvec)
azs[1][1].set_ylim([-500,500]); azs[1][1].set_xlim([-500,500]); azs[1][1].set_xlabel('East'); azs[1][1].set_ylabel('Up')
azs[1][1].set_aspect('equal')

x = sENZ.select(component='E')[0].data
y = sENZ.select(component='N')[0].data
azs[1][2].scatter(x,y,c=tvec)
azs[1][2].set_ylim([-500,500]); azs[1][2].set_xlim([-500,500]); azs[1][2].set_xlabel('East'); azs[1][2].set_ylabel('North')
azs[1][2].set_aspect('equal')

fkg.savefig('ppm235b.png')
plt.show()


