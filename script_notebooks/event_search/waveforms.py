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
# [0.03,0.125]
[0.05, 0.125]
]

# start figure
nrows = len(bpfilters)+3
fig, axs = plt.subplots(nrows, 2, figsize=(14,14))
# plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97, hspace=0.3)
plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97, hspace=0.3)
iax = 0

e_input = input('Event: ')

if e_input == '173a':
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

elif e_input == '235b':
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
    for i in [0,1]:
        axs[iax][i].set_title('Raw Data')
    iax += 1

    # plot high-passed data
    stmp = s0235b.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start235b,endtime=end235b)
    #time_axis = np.arange(start235b-P235b,end235b-P235b+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P235b)
    for i in [0,1]:
        axs[iax][i].set_title('High-Passed Data: 2-8Hz')
    iax += 1

    # plot band-passed data
    for f in bpfilters:
        stmp = s0235b.copy()
        stmp.taper(0.01,max_length=1)
        stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
        stm = stmp.slice(starttime=start235b,endtime=end235b)
        #time_axis = np.arange(start235b-P235b,end235b-P235b+dt,dt)
        wfplot(stm,axs,iax,scale,shift,P235b)
        for i in [0,1]:
            axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
        iax += 1

    # plot low-passed data
    stmp = s0235b.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start235b,endtime=end235b)
    #time_axis = np.arange(start235b-P235b,end235b-P235b+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P235b)
    for i in [0,1]:
        axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

elif e_input == '325a':
    e = 'Event S0325a'
    P325a = UTCDateTime('2019-10-26T06:58:57')
    S325a = UTCDateTime('2019-10-26T07:02:49')
    start325a = P325a - 1*60 # 1 minute early
    end325a = S325a + 3*(S325a-P325a)
    s0325a = waveforms(start325a, end325a, 600)

    shift = 1; scale = 2

    # plot raw data
    stmp = s0325a.slice(starttime=start325a,endtime=end325a)
    #time_axis = np.arange(start325a-P325a,end325a-P325a+dt,dt)
    wfplot(stmp,axs,iax,scale,shift,P325a)
    for i in [0,1]:
        axs[iax][i].set_title('Raw Data')
    iax += 1

    # plot high-passed data
    stmp = s0325a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start325a,endtime=end325a)
    #time_axis = np.arange(start325a-P325a,end325a-P325a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P325a)
    for i in [0,1]:
        axs[iax][i].set_title('High-Passed Data: 2-8Hz')
    iax += 1

    # plot band-passed data
    for f in bpfilters:
        stmp = s0325a.copy()
        stmp.taper(0.01,max_length=1)
        stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
        stm = stmp.slice(starttime=start325a,endtime=end325a)
        #time_axis = np.arange(start325a-P325a,end325a-P325a+dt,dt)
        wfplot(stm,axs,iax,scale,shift,P325a)
        for i in [0,1]:
            axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
        iax += 1

    # plot low-passed data
    stmp = s0325a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start325a,endtime=end325a)
    #time_axis = np.arange(start325a-P325a,end325a-P325a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P325a)
    for i in [0,1]:
        axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

elif e_input == '183a':
    e = 'Event S0183a'
    P183a = UTCDateTime('2019-06-03T02:27:49')
    S183a = UTCDateTime('2019-06-03T02:32:15')
    start183a = P183a - 1*60 # 1 minute early
    end183a = S183a + 3*(S183a-P183a)
    s0183a = waveforms(start183a, end183a, 600)

    shift = 1; scale = 2

    # plot raw data
    stmp = s0183a.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stmp,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Raw Data')
    iax += 1

    # plot high-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('High-Passed Data: 2-8Hz')
    iax += 1

    # plot band-passed data
    for f in bpfilters:
        stmp = s0183a.copy()
        stmp.taper(0.01,max_length=1)
        stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
        stm = stmp.slice(starttime=start183a,endtime=end183a)
        #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
        wfplot(stm,axs,iax,scale,shift,P183a)
        for i in [0,1]:
            axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
        iax += 1

    # plot low-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

elif e_input == '185a':
    e = 'Event S0185a'
    P183a = UTCDateTime('2019-06-05T02:13:50')
    S183a = UTCDateTime('2019-06-05T02:19:34')
    start183a = P183a - 1*60 # 1 minute early
    end183a = S183a + 3*(S183a-P183a)
    s0183a = waveforms(start183a, end183a, 600)

    shift = 1; scale = 2

    # plot raw data
    stmp = s0183a.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stmp,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Raw Data')
    iax += 1

    # plot high-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('High-Passed Data: 2-8Hz')
    iax += 1

    # plot band-passed data
    for f in bpfilters:
        stmp = s0183a.copy()
        stmp.taper(0.01,max_length=1)
        stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
        stm = stmp.slice(starttime=start183a,endtime=end183a)
        #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
        wfplot(stm,axs,iax,scale,shift,P183a)
        for i in [0,1]:
            axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
        iax += 1

    # plot low-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

elif e_input == '167b':
    e = 'Event S0167b'
    P183a = UTCDateTime('2019-05-17T19:30:04')
    S183a = UTCDateTime('2019-05-17T19:38:33')
    start183a = P183a - 1*60 # 1 minute early
    end183a = S183a + 3*(S183a-P183a)
    s0183a = waveforms(start183a, end183a, 600)

    shift = 1; scale = 2

    # plot raw data
    stmp = s0183a.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stmp,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Raw Data')
    iax += 1

    # plot high-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('High-Passed Data: 2-8Hz')
    iax += 1

    # plot band-passed data
    for f in bpfilters:
        stmp = s0183a.copy()
        stmp.taper(0.01,max_length=1)
        stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
        stm = stmp.slice(starttime=start183a,endtime=end183a)
        #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
        wfplot(stm,axs,iax,scale,shift,P183a)
        for i in [0,1]:
            axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
        iax += 1

    # plot low-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

elif e_input == '407a':
    e = 'Event S0407a'
    P183a = UTCDateTime('2020-01-19T09:57:45')
    S183a = UTCDateTime('2020-01-19T10:00:36')
    start183a = P183a - 1*60 # 1 minute early
    end183a = S183a + 3*(S183a-P183a)
    s0183a = waveforms(start183a, end183a, 600)

    shift = 1; scale = 2

    # plot raw data
    stmp = s0183a.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stmp,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Raw Data')
    iax += 1

    # plot high-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('High-Passed Data: 2-8Hz')
    iax += 1

    # plot band-passed data
    for f in bpfilters:
        stmp = s0183a.copy()
        stmp.taper(0.01,max_length=1)
        stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
        stm = stmp.slice(starttime=start183a,endtime=end183a)
        #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
        wfplot(stm,axs,iax,scale,shift,P183a)
        for i in [0,1]:
            axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
        iax += 1

    # plot low-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

elif e_input == '484b':
    e = 'Event S0484b'
    P183a = UTCDateTime('2020-04-07T08:52:24')
    S183a = UTCDateTime('2020-04-07T08:55:29')
    start183a = P183a - 1*60 # 1 minute early
    end183a = S183a + 3*(S183a-P183a)
    s0183a = waveforms(start183a, end183a, 600)

    shift = 1; scale = 2

    # plot raw data
    stmp = s0183a.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stmp,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Raw Data')
    iax += 1

    # plot high-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('highpass',freq=bpfilters[0][1], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('High-Passed Data: 2-8Hz')
    iax += 1

    # plot band-passed data
    for f in bpfilters:
        stmp = s0183a.copy()
        stmp.taper(0.01,max_length=1)
        stmp.filter('bandpass',freqmin=f[0], freqmax=f[1],corners=4, zerophase=True)
        stm = stmp.slice(starttime=start183a,endtime=end183a)
        #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
        wfplot(stm,axs,iax,scale,shift,P183a)
        for i in [0,1]:
            axs[iax][i].set_title('Band-Passed Data: ' + str(f[0]) + '-' + str(f[1]) + 'Hz')
        iax += 1

    # plot low-passed data
    stmp = s0183a.copy()
    stmp.taper(0.01,max_length=1)
    stmp.filter('lowpass',freq=bpfilters[-1][0], corners=4, zerophase=True)
    stm = stmp.slice(starttime=start183a,endtime=end183a)
    #time_axis = np.arange(start183a-P183a,end183a-P183a+dt,dt)
    wfplot(stm,axs,iax,scale,shift,P183a)
    for i in [0,1]:
        axs[iax][i].set_title('Low-Passed Data: 0.03-0.125Hz')

fig.savefig(e_input + '.pdf')
plt.show()
