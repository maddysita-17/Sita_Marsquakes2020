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

fig = beachball(fault)
