import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
client = Client("IRIS")
from matplotlib.dates import DateFormatter
import matplotlib.ticker as ticker
import numpy as np
import math
import pandas as pd

#waveform functions to make importing data easier
def waveforms(energystart, energyend, adjtime):
    s_energy = UTCDateTime(energystart)
    e_energy = UTCDateTime(energyend)
    start_t = s_energy - adjtime
    end_t = e_energy + adjtime
    event_st = client.get_waveforms("XB", "ELYSE", "02", "B*", start_t, end_t)
    return event_st

def waveform_plotter(date, filtered, event, ax, ylim = False, channels = [0,1,2]):
    d = date[0:10]

    for channel in channels:
        full_code = filtered[channel].id
        code = full_code[12:]
        offset = channel * np.full(len(filtered[channel]), fill_value=200)

        tr = filtered[channel]
        t = tr.times('matplotlib')

        if code == 'BHU':
            ax.plot(t, filtered[channel].data + offset, label=code, color = "#1f77b4", alpha = 0.7)
        elif code == 'BHV':
            ax.plot(t, filtered[channel].data + offset, label=code, color = "#ff7f0e", alpha = 0.7)
        elif code == 'BHW':
            ax.plot(t, filtered[channel].data + offset, label=code, color = "#2ca02c", alpha = 0.7)
        else:
            ax.plot(t, filtered[channel].data + offset, label="Unknown " + code, color = "black", alpha = 0.7)

    ax.xaxis_date()
    x_labels = ax.get_xticklabels()
    ax.set_xticklabels(x_labels, rotation=30, ha='right', va='center_baseline', size=9)
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(8))

    ax.axvline(pd.to_datetime(date), c='r', ls='--', alpha = 0.5)

    if ylim == True:
        ax.set_ylim(-800,1000)

    ax.set_title(event + "\n" + "UTC " + d, size=10)

def waveform_filter(stream, event_type):

    stream.detrend('linear')
    stream.taper(max_percentage=0.05, type='cosine')

    if event_type == 'lf' or 'bb':
        filtered_stream1 = stream.filter('bandpass', freqmin = 0.01, freqmax = 0.125)
        #filtered_stream1 = stream.filter('bandpass', freqmin = 0.125, freqmax = 0.5)
        return filtered_stream1
    elif event_type == 'hf':
        filtered_stream2 = stream.filter('highpass', freq = 1)
        return filtered_stream2
    elif event_type == '2.4':
        filtered_stream3 = stream.filter('bandpass', freqmin = 1, freqmax = 4)
        return filtered_stream3
    elif event_type == 'shf':
        filtered_stream4 = stream.filter('bandpass', freqmin = 8, freqmax = 15)
        return filtered_stream4
    elif event_type == 'vhf':
        filtered_stream5 = stream.filter('bandpass', freqmin == 5, freqmax = 10)
        return filtered_stream5
    else:
        text = "This isn't a valid event type"
        return text

def xyz_plotter(date, filtered, event, ax, ylim = False, channels=[0,1,2]):
    day = date[0:10]

    for channel in channels:
        full_code = filtered[channel].id
        code = full_code[12:]

        tr = filtered[channel]
        t = tr.times('matplotlib')

        if code == 'BHU':
            U = filtered[channel].data
        elif code == 'BHV':
            V = filtered[channel].data
        elif code == 'BHW':
            W = filtered[channel].data


    d = np.radians(-30)
    aU = np.radians(135)
    aV = np.radians(15)
    aW = np.radians(255)

    A = np.array([[np.cos(d)*np.sin(aU),np.cos(d)*np.cos(aU),np.sin(d)],
    [np.cos(d)*np.sin(aV), np.cos(d)*np.cos(aV), np.sin(d)],
    [np.cos(d)*np.sin(aW), np.cos(d)*np.cos(aW), np.sin(d)]])

    B = np.linalg.inv(A)

    E,N,Z = np.dot(B,(U,V,W))

    offset = np.full(len(Z), fill_value=200)

    ax.plot(t, E + 2*offset, label='East', color = "#77e59b", alpha = 0.8)
    ax.plot(t, N + offset, label='North', color = "#ffabab", alpha = 0.8)
    ax.plot(t, Z, label = 'Z', color = "#b28dff", alpha = 0.8)

    ax.xaxis_date()
    x_labels = ax.get_xticklabels()
    ax.set_xticklabels(x_labels, rotation=30, ha='right', va='center_baseline', size=9)
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(8))

    if ylim == True:
        ax.set_ylim(-1800,2000)

    ax.axvline(pd.to_datetime(date), c='r', ls='--', alpha = 0.5)

    ax.set_title(event + "\n" + "UTC " + day, size=10)
    return E,N,Z

def legend(ax, row, column):
    st = waveforms("2019-12-19T11:55:00", "2019-12-19T12:08:55", 20*60)
    tr = st[0]
    t = tr.times('matplotlib')
    y = np.zeros(len(t))

    ax[row,column].plot(t, y, label='BHU', color = "#1f77b4", alpha = 0.7)
    ax[row,column].plot(t, y, label='BHV', color = "#ff7f0e", alpha = 0.7)
    ax[row,column].plot(t, y, label='BHW', color = "#2ca02c", alpha = 0.7)
    ax[row,column].plot(t, y, label='E', color = "#77e59b", alpha = 0.8)
    ax[row,column].plot(t, y, label='N', color = "#ffabab", alpha = 0.8)
    ax[row,column].plot(t, y, label = 'Z', color = "#b28dff", alpha = 0.8)

    ax[row,column].legend(fontsize = 'x-large', loc = 'center', facecolor = 'white', mode = 'expand', framealpha = 1)
    ax[row,column].get_xaxis().set_visible(False)
    ax[row,column].get_yaxis().set_visible(False)

fig,(ax1, ax2) = plt.subplots(2,1, figsize=(12,9))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.8)

e = 'S0173a'
### 2019-05-23T02:19:33
start173a = '2019-05-23T02:22:48'
end173a = '2019-05-23T03:00:23'

P173a = '2019-05-23T02:22:59'

S173a = '2019-05-23T02:25:53'

event_st = waveforms(start173a, end173a, 0)
f_event_st = waveform_filter(event_st, 'lf')

dt = 1/(event_st[0].stats.sampling_rate)
begin = UTCDateTime(start173a) - UTCDateTime(P173a)
print(begin)
end = UTCDateTime(end173a) - UTCDateTime(P173a)
print(end)
time_axis = np.arange(begin, end, dt)
print(len(f_event_st[0].data), len(time_axis))

E173a, N173a, Z173a = xyz_plotter(start173a, f_event_st, e, ax1)

fig2, ax3 = plt.subplots(1,1)
ax3.plot(time_axis, E173a, color='#77e59b')
ax3.plot(time_axis, N173a + 150, color='#ffabab')
ax3.plot(time_axis, Z173a + 300, color='#b28dff')

e = 'S0235b'
### 2019-07-26T12:16:15
start235b = '2019-07-26T12:19:16'
end235b = '2019-07-26T13:18:08'

P235b = '2019-07-26T12:19:19'

S235b = '2019-07-26T12:22:05'

short235b = waveforms(start235b, end235b, 0)
sf235b = waveform_filter(short235b, 'bb')

E235b, N235b, Z235b = xyz_plotter(start235b, sf235b, e, ax2)

fig2, ax3 = plt.subplots(1,1)
ax3.plot(time_axis, E235b - 50, color='#77e59b')
ax3.plot(time_axis, N235b + 100, color='#ffabab')
ax3.plot(time_axis, Z235b + 250, color='#b28dff')
plt.show()
