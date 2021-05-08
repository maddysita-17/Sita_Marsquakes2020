import matplotlib.pyplot as plt
import numpy as np
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

client = Client("IRIS")

#------ function definitions------
def waves_uganda(start, end, adjtime):
    st = client.get_waveforms("II", "MBAR", "10", "B*", start-adjtime, end+adjtime)
    st.detrend(type='simple')
    return st

def waves_comp(start, end, adjtime):
    st = client.get_waveforms("US", "MNTX", "00", "B*", start-adjtime, end+adjtime)
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

def rot2enz(st, comp=False):
    if len(st) != 3:
       print('Stream does not contain 3 Traces')
       return st
    for trace in st:
        head = trace.stats
        channel = head.channel
        if channel == 'BH1': one = trace.data
        elif channel == 'BH2': two = trace.data
        elif channel == 'BHZ': Z = trace.data
        elif channel == 'BHN': N = trace.data
        elif channel == 'BHE': E = trace.data
        else:
            print('Trace.channel is not BH1, BH2, or BHZ')
            return st

    if comp==False:
        N,E = rotate(one, two, -14)

    elif comp==True:
        N,E = rotate(one, two, 23)

    head.channel = 'BHE'; trE = Trace(data=E, header=head)
    head.channel = 'BHN'; trN = Trace(data=N, header=head)
    head.channel = 'BHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ


def wfplot(st,axs,iax,scale=1,sh=0,Ptime=0, comp=False):
    if Ptime == 0: Ptime = st[0].stats.starttime
    time_axis = st[0].times(reftime=Ptime)
    n = min([len(time_axis)]+[len(st[i].data) for i in np.arange(len(st))])

    a = 1/200; ylim = [-1000, 1000]
    shift = sh; dsh = a*600
    for trace in st:
        axs[iax][0].plot(time_axis[:n], a*trace.data[:n] + shift, label = trace.stats.channel)
        shift += dsh
    axs[iax][0].set_ylim(ylim)

    if comp==False:
        sENZ = rot2enz(st)
    elif comp==True:
        sENZ = rot2enz(st, comp=True)

    shift = sh
    for trace in sENZ:
        axs[iax][1].plot(time_axis[:n], a*trace.data[:n] + shift, label = trace.stats.channel)
        shift += dsh
    axs[iax][1].set_ylim(ylim)
    if iax == len(axs)-1:
       axs[len(axs)-1][0].legend(loc='lower right')
       axs[len(axs)-1][1].legend(loc='lower right')

def sigplt(event,st,axs,scale=1,sh=2,Ptime=0, comp=False):
    if Ptime == 0: Ptime = st[0].stats.starttime
    time_axis = st[0].times(reftime=Ptime)
    n = min([len(time_axis)]+[len(st[i].data) for i in np.arange(len(st))])

    a = 1/1000; ylim = [-50,50]
    shift = sh; dsh = a*10000

    if comp==False:
        sENZ = rot2enz(st)
    elif comp==True:
        sENZ = rot2enz(st, comp=True)

    shift = sh
    for trace in sENZ:
        ajs.plot(time_axis[:n], a*trace.data[:n] + shift, label = event + ": " + trace.stats.channel)
        ajs.annotate(trace.stats.channel, xy=(-25,shift+0.5), size='medium', color='k')
        shift += dsh
    ajs.set_ylim(ylim)
    ajs.title.set_text(event)

params = {"ytick.color" : "k",
          "xtick.color" : "k",
          "axes.labelcolor" : "k",
          "axes.edgecolor" : "k"}
plt.rcParams.update(params)


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

#---- Nebbi, Uganda -----
e = 'Nebbi, Uganda'
### 2021-03-20 17:07:27.065
P173a = UTCDateTime('2021-03-20T17:07:27')
S173a = UTCDateTime('2021-03-20T17:09:05')
start173a = P173a - 1*60 # 1 minute early
end173a = S173a + 3*(S173a-P173a)
s0173a = waves_uganda(start173a, end173a, 600)

sps = s0173a[0].stats.sampling_rate; dt = 1./sps
if 2*bpfilters[0][1] > sps:
    print('WARNING: high bandpass corner too high. f_nyq = ',0.5*sps,' Hz')
shift = 0; scale = 1

# plot raw data
stmp = s0173a.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stmp,axs,iax,scale,shift,P173a)
for i in [0,1]:
    axs[iax][i].set_title('Raw Data')
iax += 1

# plot high-passed data
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stm,axs,iax,scale,shift,P173a)
for i in [0,1]:
    axs[iax][i].set_title('High-Passed Data: 2-8Hz')
iax += 1

# plot band-passed data
for f in bpfilters:
    stmp = s0173a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
    stm = stmp.slice(starttime=start173a,endtime=end173a)
    #time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P173a)
    for i in [0,1]:
        axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
    iax += 1

# plot low-passed data
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stm,axs,iax,scale,shift,P173a)
for i in [0,1]:
    axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

iax = 0

plt.show()

iax = 0

fjg,ajs = plt.subplots(1,1, figsize=(14,7))
plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97, hspace=0.3)

stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('bandpass',freqmin=0.125, freqmax=1.0,corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
sigplt(e,stm,axs,scale,shift,P173a)

plt.show()

#---- New Mexico -----
e = 'Idaho'
### 2021-03-20 17:07:27.065
P173a = UTCDateTime('2020-03-31T23:52:30')
S173a = UTCDateTime('2020-03-31T23:53:00')
start173a = P173a - 1*60 # 1 minute early
end173a = S173a + 3*(S173a-P173a)
s0173a = waves_uganda(start173a, end173a, 600)

sps = s0173a[0].stats.sampling_rate; dt = 1./sps
if 2*bpfilters[0][1] > sps:
    print('WARNING: high bandpass corner too high. f_nyq = ',0.5*sps,' Hz')
shift = 0; scale = 1

# plot raw data
stmp = s0173a.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stmp,axs,iax,scale,shift,P173a, comp=True)
for i in [0,1]:
    axs[iax][i].set_title('Raw Data')
iax += 1

# plot high-passed data
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stm,axs,iax,scale,shift,P173a, comp=True)
for i in [0,1]:
    axs[iax][i].set_title('High-Passed Data: 2-8Hz')
iax += 1

# plot band-passed data
for f in bpfilters:
    stmp = s0173a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
    stm = stmp.slice(starttime=start173a,endtime=end173a)
    #time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P173a, comp=True)
    for i in [0,1]:
        axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
    iax += 1

# plot low-passed data
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
#time_axis = np.arange(start173a-P173a,end173a-P173a+dt,dt)
wfplot(stm,axs,iax,scale,shift,P173a,comp=True)
for i in [0,1]:
    axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

iax = 0

plt.show()

iax = 0

fjg,ajs = plt.subplots(1,1, figsize=(14,7))
plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97, hspace=0.3)

stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('bandpass',freqmin=0.125, freqmax=1.0,corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
sigplt(e,stm,axs,scale,shift,P173a,comp=True)

plt.show()
