import matplotlib.pyplot as plt
import numpy as np
from obspy.core.trace import Trace, Stats
from obspy.core.stream import Stream
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

client = Client("IRIS")

params = {'font.family' : 'sans-serif', 'font.sans-serif': 'Inconsolata'}
plt.rcParams.update(params)

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

#--------grabbing waveform data-------

e = 'Event S0235b'
### 2019-07-26T12:16:15
P235b = UTCDateTime('2019-07-26T12:19:19')
S235b = UTCDateTime('2019-07-26T12:22:05')
start235b = P235b - 1*60 # 1 minute early
end235b = S235b + 3*(S235b-P235b)
s0235b = waveforms(start235b, end235b, 600)


e = 'Event S0325a'
P325a = UTCDateTime('2019-10-26T06:58:57')
S325a = UTCDateTime('2019-10-26T07:02:46')
start325a = P325a - 1*60 # 1 minute early
end325a = S325a + 3*(S325a-P325a)
s0325a = waveforms(start325a, end325a, 600)

e = 'Event S0173a'
P173a = UTCDateTime('2019-05-23T02:22:59')
S173a = UTCDateTime('2019-05-23T02:25:53')
start173a = P173a - 1*60 # 1 minute early
end173a = S173a + 3*(S173a-P173a)
s0173a = waveforms(start173a, end173a, 600)


# ....plot particle motion...
def partmot(stream, event, Pwave, Swave, begin, end):
    stmp = stream.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('bandpass',freqmin=0.125, freqmax=1.0 ,corners=4, zerophase=True)
    bt = begin; et = end

    fjg, ays = plt.subplots(2,3,figsize=(10,6))
    for i in range(2):
        for j in range(3):
            ays[i][j].tick_params(labelsize=6)
    plt.subplots_adjust(wspace = 0.3, hspace=0.2)

    stm = stmp.slice(starttime=Pwave+bt,endtime=Pwave+et)
    sENZ = uvw2enz(stm)
    tvec = sENZ[0].times(reftime=Pwave)

    x = sENZ.select(component='N')[0].data
    y = sENZ.select(component='Z')[0].data
    ays[0][0].scatter(x,y,c=tvec, cmap='summer')
    ays[0][0].set_ylim([-200,200]); ays[0][0].set_xlim([-200,200]); ays[0][0].set_xlabel('North'); ays[0][0].set_ylabel('Vertical')
    ays[0][0].set_aspect('equal'); ays[0][0].set_title('S0' + event +' P-wave Particle Motion')

    x = sENZ.select(component='E')[0].data
    y = sENZ.select(component='Z')[0].data
    ays[0][1].scatter(x,y,c=tvec, cmap='summer')
    ays[0][1].set_ylim([-200,200]); ays[0][1].set_xlim([-200,200]); ays[0][1].set_xlabel('East'); ays[0][1].set_ylabel('Vertical')
    ays[0][1].set_aspect('equal')

    x = sENZ.select(component='E')[0].data
    y = sENZ.select(component='N')[0].data
    ays[0][2].scatter(x,y,c=tvec, cmap='summer')
    ays[0][2].set_ylim([-200,200]); ays[0][2].set_xlim([-200,200]); ays[0][2].set_xlabel('East'); ays[0][2].set_ylabel('North')
    ays[0][2].set_aspect('equal')

    pc = ays[0][2].scatter(x,y,c=tvec, cmap='summer')
    fjg.colorbar(pc, ax=ays[0:, :], label = 'Seconds', location='right', shrink=0.6)

    stm = stmp.slice(starttime = Swave+bt, endtime = Swave+et)
    sENZ = uvw2enz(stm)
    tvec = sENZ[0].times(reftime=Swave)

    x = sENZ.select(component='N')[0].data
    y = sENZ.select(component='Z')[0].data
    ays[1][0].scatter(x,y,c=tvec, cmap='summer')
    ays[1][0].set_ylim([-500,500]); ays[1][0].set_xlim([-500,500]); ays[1][0].set_xlabel('North'); ays[1][0].set_ylabel('Vertical')
    ays[1][0].set_aspect('equal'); ays[1][0].set_title('S0' + event + ' S-wave Particle Motion')

    x = sENZ.select(component='E')[0].data
    y = sENZ.select(component='Z')[0].data
    ays[1][1].scatter(x,y,c=tvec, cmap='summer')
    ays[1][1].set_ylim([-500,500]); ays[1][1].set_xlim([-500,500]); ays[1][1].set_xlabel('East'); ays[1][1].set_ylabel('Vertical')
    ays[1][1].set_aspect('equal')

    x = sENZ.select(component='E')[0].data
    y = sENZ.select(component='N')[0].data
    ays[1][2].scatter(x,y,c=tvec, cmap='summer')
    ays[1][2].set_ylim([-500,500]); ays[1][2].set_xlim([-500,500]); ays[1][2].set_xlabel('East'); ays[1][2].set_ylabel('North')
    ays[1][2].set_aspect('equal')


    fjg.savefig('part2_' + event + '.png', transparent=True)

partmot(s0325a,'325a',P325a, S325a, -1, 11)

plt.show()
