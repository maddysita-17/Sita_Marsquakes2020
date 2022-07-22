import matplotlib.pyplot as plt
import numpy as np
import obspy
import pandas as pd

from obspy.imaging.beachball import beachball
from obspy.imaging.beachball import beach
from obspy.imaging.source import plot_radiation_pattern

fault = []

st = input('Strike:')
fault.append(float(st))
dp = input('Dip: ')
fault.append(float(dp))
rk = input('Rake: ')
fault.append(float(rk))

# #hawaii eq
# mt = [62800000000000000,-51410000000000000,-11400000000000000,435400000000000000,27840000000000000,-277700000000000000]
# scalarmoment = 520000000000000000
# fault = [m/scalarmoment for m in mt]

fig = beachball(fault, facecolor='b')
