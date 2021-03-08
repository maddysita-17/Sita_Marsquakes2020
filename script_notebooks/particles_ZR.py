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


def enz2trz(st, baz):
    """
    IN: st, az
    st is the reorientated ENZ data stream
    baz is the bAz
    OUT: st
    st is the rotated ZTR data stream
    """
    az = baz + 180

    hhe = st[0].data
    hhn = st[1].data
    hhz = st[2].data

    hhT,hhR = rotate(hhe,hhn,az)

    streamRT = st.copy()
    streamRT[0].data = hhT
    streamRT[1].data = hhR
    streamRT[0].stats.component = 'T'
    streamRT[1].stats.component = 'R'
    return streamRT

def diag(ang):
    """
    IN: ang
    estimated angle of orientation of the pm / estimated incidence angle
    OUT: plotted line
    estimated line on the ppm
    """
    ang_rad = np.radians(ang)
    yval = np.cos(200)*ang_rad
    xval = np.sin(200)*ang_rad

    for j in [0,1]:
        ays[j].plot((-1*xval,xval), (-1*yval,yval), 'b--')
        ays[j].annotate(str(ang), (xval-5, yval-5))

    #figure out your triangle shit, maybe make yourself a user input so you can feel cool


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

e = 'Event S0183a'
P183a = UTCDateTime('2019-06-03T02:27:49')
S183a = UTCDateTime('2019-06-03T02:32:15')
start183a = P183a - 1*60
end183a = S183a + 3*(S183a-P183a)
s0183a = waveforms(start183a, end183a, 600)


# ....plot particle motion...
def partmot(stream, event, Pwave, Swave, begin, end, bAz, ia_guess):
    stmp = stream.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('bandpass',freqmin=0.125, freqmax=1.0 ,corners=4, zerophase=True)
    bt = begin; et = end

    fjg, ays = plt.subplots(1,2,figsize=(10,7))
    for j in range(2):
        ays[j].tick_params(labelsize=6)
    plt.subplots_adjust(left=0.05, right=0.80, bottom=0.05, top=0.90, wspace = 0.3, hspace=0.2)

    stm = stmp.slice(starttime=Pwave+bt,endtime=Pwave+et)
    sENZ = uvw2enz(stm)
    sRTZ = enz2trz(sENZ,bAz)
    tvec = sENZ[0].times(reftime=Pwave)

    x = sRTZ.select(component='R')[0].data
    y = sRTZ.select(component='Z')[0].data
    ays[0].scatter(x,y,c=tvec, cmap='Set2')
    ays[0].set_ylim([-200,200]); ays[0].set_xlim([-200,200]); ays[0].set_xlabel('Radial'); ays[0].set_ylabel('Vertical')
    ays[0].set_aspect('equal'); ays[0].set_title(event +' P-wave Particle Motion')


    stm = stmp.slice(starttime = Swave+bt, endtime = Swave+et)
    sENZ = uvw2enz(stm)
    sRTZ = enz2trz(sENZ,bAz)
    tvec = sENZ[0].times(reftime=Swave)

    x = sRTZ.select(component='R')[0].data
    y = sRTZ.select(component='Z')[0].data
    ays[1].scatter(x,y,c=tvec, cmap='Set2')
    ays[1].set_ylim([-200,200]); ays[1].set_xlim([-200,200]); ays[1].set_xlabel('Radial'); ays[1].set_ylabel('Vertical')
    ays[1].set_aspect('equal'); ays[1].set_title(event+' S-wave Particle Motion')

    sc = ays[1].scatter(x,y,c=tvec, cmap='Set2', alpha=0.8)
    fjg.colorbar(sc)

    for j in [0,1]:
        ays[j].plot((-200,200), (0,0), 'k-')
        ays[j].plot((0,0), (-200,200), 'k-')

    ang_rad = np.radians(ia_guess)
    yval = np.sin((np.radians(270))-ang_rad)*200
    xval = yval/np.tan((np.radians(270))-ang_rad)

    for j in [0,1]:
        ays[j].plot((-1*xval,xval), (-1*yval,yval), 'b--')
        ays[j].annotate(str(ia_guess), (xval-5, yval-5))

    path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/paper_figures/ppm-plots/RZplots/'
    fjg.savefig(path + 'RZ_ppm' + event + '.png', transparent=False)


partmot(s0235b,'235b',P235b, S235b, -4, 4, 74, 20)

partmot(s0325a, '325a', P325a, S325a, -1, 3, 123, 20)
partmot(s0325a,'325ab', P325a, S325a, 10, 14, 139, 15)

# partmot(s0173a,'ALL 173a',P173a, S173a, -1, 8, 90, 81)
partmot(s0173a, '173a', P173a, S173a, -1, 3, 90, 30)
partmot(s0173a, '173ab', P173a, S173a, 4, 8, 86, 30)

#partmot(s0183a, '183a', P183a, S183a, -5, 4, 89, 55)

plt.show()
