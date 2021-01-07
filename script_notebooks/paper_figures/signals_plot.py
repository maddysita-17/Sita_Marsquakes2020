import matplotlib.pyplot as plt
import numpy as np
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

client = Client("IRIS")

#------ function definitions------
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

    A = np.array([[np.cos(d)*np.sin(aU),np.cos(d)*np.cos(aU),-np.sin(d)],
                  [np.cos(d)*np.sin(aV), np.cos(d)*np.cos(aV), -np.sin(d)],
                  [np.cos(d)*np.sin(aW), np.cos(d)*np.cos(aW), -np.sin(d)]])

    B = np.linalg.inv(A)
    E,N,Z = np.dot(B,(U,V,W))

    head.channel = 'BHE'; trE = Trace(data=E, header=head)
    head.channel = 'BHN'; trN = Trace(data=N, header=head)
    head.channel = 'BHZ'; trZ = Trace(data=Z, header=head)
    stENZ = Stream(traces=[trE,trN,trZ])

    return stENZ

def wfplot(event,st,axs,iax,scale=1,sh=0,Ptime=0):
    if Ptime == 0: Ptime = st[0].stats.starttime
    time_axis = st[0].times(reftime=Ptime)
    n = min([len(time_axis)]+[len(st[i].data) for i in np.arange(len(st))])

    a = 1/200; ylim = [-2,12]
    shift = sh; dsh = a*800
    sENZ = uvw2enz(st)
    shift = sh
    for trace in sENZ:
        axs[iax].plot(time_axis[:n], a*trace.data[:n] + shift, label = event + ": " + trace.stats.channel)
        axs[iax].annotate(trace.stats.channel, xy=(-25,shift+0.5), size='xx-small', color='k')
        shift += dsh
    axs[iax].set_ylim(ylim)
    axs[iax].title.set_text(event)
    #axs[iax].legend(loc='lower right', fontsize='x-small')

# params = {"ytick.color" : "w",
#           "xtick.color" : "w",
#           "axes.labelcolor" : "w",
#           "axes.edgecolor" : "w"}
# plt.rcParams.update(params)

params = {"ytick.color" : "k",
          "xtick.color" : "k",
          "axes.labelcolor" : "k",
          "axes.edgecolor" : "k"}
plt.rcParams.update(params)


#---setting up subplots-----

fig, axs = plt.subplots(3, 1, figsize=(12,10))
for j in range(3):
    axs[j].tick_params(labelsize=6)
iax = 0

#---- 235b -----

e = 'S0235b'
### 2019-07-26T12:16:15
P235b = UTCDateTime('2019-07-26T12:19:19')
S235b = UTCDateTime('2019-07-26T12:22:05')
start235b = P235b - 1*60 # 1 minute early
end235b = P235b + 600
s0235b = waveforms(start235b, end235b, 600)

shift = 1; scale = 0

#---bandpassed data----
stmp = s0235b.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('bandpass',freqmin=0.125, freqmax=1.0,corners=4, zerophase=True)
stm = stmp.slice(starttime=start235b,endtime=end235b)
#axs[iax].set_prop_cycle(color=['#3BA7BF', '#6297A3', '#898686'])
axs[iax].set_prop_cycle(color=['k', 'k', 'k'])
wfplot(e,stm,axs,iax,scale,shift,P235b)
#axs[iax].axis('off')

iax += 1


#----- 173a -----
e = 'S0173a'
P173a = UTCDateTime('2019-05-23T02:22:59')
S173a = UTCDateTime('2019-05-23T02:25:53')
start173a = P173a - 1*60 # 1 minute early
end173a = P173a + 600
s0173a = waveforms(start173a, end173a, 600)


#-----bandpassed data----
stmp = s0173a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('bandpass',freqmin=0.125, freqmax=1.0,corners=4, zerophase=True)
stm = stmp.slice(starttime=start173a,endtime=end173a)
#axs[iax].set_prop_cycle(color=['#3BA7BF', '#6297A3', '#898686'])
axs[iax].set_prop_cycle(color=['k', 'k', 'k'])
wfplot(e,stm,axs,iax,scale,shift,P173a)
#axs[iax].axis('off')
iax += 1


#----- 325a -----
e = 'S0325a'
P325a = UTCDateTime('2019-10-26T06:58:58')
S325a = UTCDateTime('2019-10-26T07:02:56')
start325a = P325a - 1*60 # 1 minute early
end325a = P325a + 600
s0325a = waveforms(start325a, end325a, 600)


#-----bandpassed data----
stmp = s0325a.copy()
stmp.taper(0.01,max_length=1)
stmp.filter('bandpass',freqmin=0.125, freqmax=1.0,corners=4, zerophase=True)
stm = stmp.slice(starttime=start325a,endtime=end325a)
#axs[iax].set_prop_cycle(color=['#BF533B', '#D8802A', '#F0AC19'])
axs[iax].set_prop_cycle(color=['k', 'k', 'k'])
wfplot(e,stm,axs,iax,scale,shift,P325a)
#axs[iax].yaxis.set_visible(False)
iax += 1

for j in range(3):
    for item in [fig, axs[j]]:
        item.patch.set_visible(False)
    # for a in ('left', 'right', 'top', 'bottom'):
    #     axs[j].spines[a].set_visible(False)
fig.patch.set_visible(False)

fig.savefig('indiv_plots.png', transparent=True)
plt.show()
