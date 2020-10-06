#package imports
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy import read_events
import matplotlib.pyplot as plt
client = Client("IRIS")
from matplotlib.dates import DateFormatter
import matplotlib.ticker as ticker
import numpy as np
import math
import pandas as pd

#function definitions
def cat_waveforms(start, adjtime):
    start_t = start - adjtime
    end_t = start + 2*adjtime
    event_st = client.get_waveforms("XB", "ELYSE", "02", "B*", start_t, end_t)
    return event_st

def waveform_filter(stream, event_type):

    stream.detrend('linear')
    stream.taper(max_percentage=0.05, type='cosine')

    if event_type == 'lf' or 'bb':
        filtered_stream1 = stream.filter('bandpass', freqmin = 0.125, freqmax = 0.5)
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

def waveform_plotter(date, filtered, ax, ylim = False, channels = [0,1,2]):

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
    ax.set_title(date, size=10)

    ax.axvline(pd.to_datetime(date), c='r', ls='--', alpha = 0.5)

    if ylim == True:
        ax.set_ylim(-800,1000)


#reading in start times from a file already created as event times as strings
file = open('mars2_4Hzevents.txt', 'r')

event_times = []

for line in file:
    time = line.strip()
    UTCtime = UTCDateTime(time)
    event_times.append(UTCtime)

part1 = event_times[0:170]
part2 = event_times[171:]

#creating subplot figure - 343 events (rows x columns)
fig,ax = plt.subplots(19,9, figsize=(90,190))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.1)

fig1,ax1 = plt.subplots(19,9, figsize=(90,190))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.1)


#loop for plotting in each individual subplots
k = 0
for num in range(0,19):
    for numb in range(0,9):
        if k < len(part1):
            t = part1[k]
            try:
                wave = cat_waveforms(t, 3*60)
                filt = waveform_filter(wave, 'hf')
                waveform_plotter(str(t), filt, ax[num,numb])
                print(num, numb, t, k)
            except:
                print('Error occured', t)
            k += 1

m = 0
for num in range(0,19):
    for numb in range(0,9):
        if m < len(part2):
            t = part2[m]
            try:
                wave = cat_waveforms(t, 3*60)
                filt = waveform_filter(wave, 'hf')
                waveform_plotter(str(t), filt, ax1[num,numb])
                print(num, numb, t, m)
            except:
                print('Error occured', t)
            m += 1

#save figure as png
fig.savefig('2_4Hz_part1.pdf')
fig1.savefig('2_4Hz_part2.pdf')
