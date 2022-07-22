import matplotlib
import matplotlib.pyplot as plt
plt.style.use('classic')
import numpy as np

IPython_default = plt.rcParams.copy()

from matplotlib import cycler

plt.rcParams['font.sans-serif'] = 'Helvetica Neue'

colors = cycler('color',
                ['#EE6666', '#3388BB', '#9988DD',
                 '#EECC55', '#88BB44', '#FFBBBB'])
plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
       axisbelow=True, grid=True, prop_cycle=colors)
plt.rc('grid', color='w', linestyle='solid')
plt.rc('xtick', direction='out', color='gray')
plt.rc('ytick', direction='out', color='gray')
plt.rc('patch', edgecolor='#E6E6E6')
plt.rc('lines', linewidth=2)


x = np.random.randn(1000)
plt.hist(x)
plt.show()
