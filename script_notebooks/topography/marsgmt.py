#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import struct
import pygmt as gmt
import xarray as xr


# In[1]:

#
# import pygmt
# pygmt.show_versions()


# In[2]:


def fi2(fi1,teta1,azimuth,delta):
    term1 = np.cos(delta)*np.sin(fi1)*np.cos(teta1)
    factor2 = np.sin(delta)*np.cos(azimuth)*np.sin(teta1)
    term2 = factor2*np.sin(fi1)
    factor3 = np.sin(delta)*np.sin(azimuth)
    term3 = factor3*np.cos(fi1)
    teller = term1 - term2 + term3
    term1 = np.cos(delta)*np.cos(fi1)*np.cos(teta1)
    term2 = factor2*np.cos(fi1)
    term3 = factor3*np.sin(fi1)
    rnoemer = term1 - term2 - term3
    return np.arctan2(teller,rnoemer)

def teta2(teta1,azimuth,delta):
    term1 = np.cos(delta)*np.sin(teta1)
    term2 = np.sin(delta)*np.cos(azimuth)*np.cos(teta1)
    som = term1 + term2
    return np.arcsin(som)

def azdelt(teta1d,fi1d,deltad,azimuthd):
    # IO in degrees, calculations in radians
    teta1 = teta1d*np.pi/180.; fi1 = fi1d*np.pi/180.
    delta = deltad*np.pi/180.; azimuth = azimuthd*np.pi/180.
    if teta1 > 0.5*np.pi or teta1 < -0.5*np.pi:
        print('error, non-existent latitude ', teta1)
        return 0,0
    if delta < 0.:
        print('error, non-existent delta ',delta)
        return 0,0
    f = fi2(fi1,teta1,azimuth,delta)
    t = teta2(teta1,azimuth,delta)
    return t*180/np.pi, f*180/np.pi


# * IMG format is binary, 16-bit, topo unit is m
# * 16 pixels per degree
# * got median topography from gridded data sets (randomly selected from list)
#
# DATA_SET_ID                   = "MGS-M-MOLA-5-MEGDR-L3-V1.0" <br>
# PRODUCT_ID                    = "MEGT90N000EB.IMG"           <br>
# SPACECRAFT_NAME               = "MARS GLOBAL SURVEYOR"       <br>
# INSTRUMENT_ID                 = MOLA                         <br>
# INSTRUMENT_NAME               = "MARS ORBITER LASER ALTIMETER" <br>
#
# Relevant web site: https://pds-geosciences.wustl.edu/dataserv/default.htm
# Select MGS --> MOLA --> MEGDR --> submit
# Returns a long list of .IMG files. I randomly selected one (MEGT90N000EB.IMG); it appears to be the whole planet's median topography.

# In[3]:


ppd = 16        # 16 pixels per degree
inc = 1/ppd     # increment (in degrees) between topo values
hinc = 0.5*inc
nlats = 180*ppd
nlons = 360*ppd
topom = np.zeros([nlats,nlons])
# for the quadrilateral interpretation of np.pcolormesh (shading='flat') the lat and lon array need to mark
# the boundaries of the quadrilaterals and thus each have one more element than there are pixel (topo) values
# gmt calls this grid(line) registration
lats = np.arange(90,-90-inc,-inc)
lons = np.arange(0,360+inc,inc)
# for pixel registration (needed for xarray and imshow, for example):
plats = np.arange(90-hinc,-90-hinc,-inc)
plons = np.arange(hinc,360+hinc,inc)


# In[4]:


datablock = open('megt90n000eb.img','rb')

for irow in range(nlats):
    a = datablock.read(2*nlons)      # read all bytes that represent topo values for all longitudes and one latitutde
    row = struct.unpack('>5760h',a)  # read 5760 big-endian 2-byte integers from one latitude "row" of the binary file
    topom[irow] = row                # store topo values (now integers) in a row of 2D array topom

datablock.close()
topo = 0.001*topom                   # convert topo values to km


# In[5]:


# turn topography into an xarray object for plotting with gmt
xtopo = xr.DataArray(topo,coords=(plats,plons),dims=('latitude','longitude'))


# In[6]:


centerlon = 160; centerlat = 15

distance = 20  # degrees
azs = np.arange(0,360) # degrees
y, x = azdelt(centerlat,centerlon,distance,azs)
print("small circle points calculated")


# In[ ]:


# # plot specific region,
#
# latmin = -10; latmax = 30
# lonmin = 120; lonmax = 180
#
# fig3 = gmt.Figure()
#
# gmt.makecpt(cmap="haxby", series=[-8, 8])
# fig3.grdimage(xtopo,frame=True,projection='M12',region=[lonmin,lonmax,latmin,latmax])
# fig3.colorbar(frame='+l"topography (km)"',position='JMR')
# fig3.show()


# gmt.makecpt(cmap='haxby', series=[-8, 8])
# fig3.grdimage(xtopo)
# fig3.show()

# In[ ]:


# fig5 = gmt.Figure()
# gmt.makecpt(cmap="haxby", series=[-8, 8])
# fig5.grdimage(xtopo,frame=True,projection='M12',region=[lonmin,lonmax,latmin,latmax])
#
# fig5.plot(centerlon,centerlat,style='i0.2i',color='purple',pen=True)
# fig5.plot(x,y,pen='1,black')
# fig5.show()


# To plot ellipses with pygmt.velo one needs to provide data as follows. For example, for [150,6,0,0,6,4,120,"S0111a"]:
#
# 1,2: longitude, latitude of station (150, 6)
#
# 3,4: eastward, northward velocity  (0,0)
#
# 5,6: semi-major, semi-minor axes (6,4)  # The unites of this are NOT map degrees....
#
# 7: counter-clockwise angle, in degrees, from horizontal axis to major axis of ellipse. (120)
#
# Trailing text: name of station (optional) (event label, S0111a)
#
# Use the first number in this option "r0.2/0.95" to scale the ellipse.
#

# In[ ]:


# fig6 = gmt.Figure()
# gmt.makecpt(cmap="haxby", series=[-8, 8])
# fig6.grdimage(xtopo,frame=True,projection='M12',region=[lonmin,lonmax,latmin,latmax])
#
# fig6.plot(centerlon,centerlat,style='i0.2i',color='purple',pen=True)
# fig6.plot(y,x,pen='1,black')
#
# ellipse = {"lon": [150],
#            "lat": [6],
#            "vE": [0],
#            "vN": [0],
#            "sE": [2],
#            "sN": [1],
#            "angle": [130],
#            "site": ["S0111a"],
#           }
#
# fig6.velo(data=ellipse, uncertaintycolor="darkorange", pen="black", line=True,  spec="r0.2/0.95", transparency=50)
#
# fig6.show()


# In[ ]:


# have fun


# To plot ellipses with pygmt.velo one needs to provide data as follows. For example, for [150,6,0,0,6,4,120,"S0111a"]:
#
# 1,2: longitude, latitude of station (150, 6)
#
# 3,4: eastward, northward velocity  (0,0)
#
# 5,6: semi-major, semi-minor axes (6,4)  # The unites of this are NOT map degrees....
#
# 7: counter-clockwise angle, in degrees, from horizontal axis to major axis of ellipse. (120)
#
# Trailing text: name of station (optional) (event label, S0111a)
#
# Use the first number in this option "r0.2/0.95" to scale the ellipse.
#

# In[ ]:


# s0173a ellipse
s0173a = {"lon": [164.09],
           "lat": [3.95],
           "vE": [0],
           "vN": [0],
           "sE": [1.5],
           "sN": [0.5],
           "angle": [90],
           "site": ["S0173a"],
          }


# s0235b ellipse
s0235b = {"lon": [162.04],
           "lat": [11.22],
           "vE": [0],
           "vN": [0],
           "sE": [1.5],
           "sN": [0.5],
           "angle": [95],
           "site": ["S0235b"],
          }


# s0325ab ellipse
s0325ab = {"lon": [162.103],
           "lat": [-23.94],
           "vE": [0],
           "vN": [0],
           "sE": [1.5],
           "sN": [0.5],
           "angle": [220],
           "site": ["S0325ab"],
          }


# In[ ]:


# latmin = -30; latmax = 30
latmin = -30; latmax = 18
lonmin = 120; lonmax = 180

fig7 = gmt.Figure()
gmt.makecpt(cmap="haxby", series=[-8, 8])
fig7.grdimage(xtopo,frame=True,projection='M12',region=[lonmin,lonmax,latmin,latmax])
fig7.colorbar(frame='+l"topography (km)"',position='JMR')

# lander:
centerlat=4.5024; centerlon=135.6234
fig7.plot(x=centerlon,y=centerlat,style='i0.2i',color='#836CD9',pen='0.25p,black')
fig7.text(x=centerlon,y=centerlat + 2,text='InSight')

# dist circles:
azs = np.arange(0,360) # degrees

for dist in [10,20,30,40]:
    y, x = azdelt(centerlat,centerlon,dist,azs)
    fig7.plot(x=x,y=y,pen='0.5,black,dashed')

fig7.text(x=centerlon, y=-4, text='10@.')
fig7.text(x=centerlon, y=-14, text='20@.')
fig7.text(x=centerlon, y=-24, text='30@.')
fig7.text(x=154, y=-29, text='40@.')
    
s0173a_d = [[164.09, 3.95, 90, 1.5, 0.5]]; a_1 = [[164.09, 3.95, 90, 2, 0.75]]; a_2 = [[164.09, 3.95, 90, 2.25, 1]]
s0235b_d = [[162.04, 11.22, 95, 1.5, 0.5]]; b_1 = [[162.04, 11.22, 95, 2, 0.75]]; b_2 = [[162.04, 11.22, 95, 2.25, 1]]
s0325ab_d = [[162.103, -23.94, 220, 1.7, 0.57]]; ab_1 = [[162.103, -23.94, 220, 2.28, 0.855]]; ab_2 = [[162.103, -23.94, 220, 2.56, 1.14]]
s0173ab = [[164.01, 5.95, 92.5, 1.5, 0.5]]; abb_1 = [[164.01, 5.95, 92.5, 2, 0.75]]; abb_2 = [[164.01, 5.95, 92.5, 2.25, 1]]
    
    
fig7.text(x=169, y=3.95, text='S0173a')
fig7.text(x=159, y=17, text='S0235b')
fig7.text(x=169.5, y=-23.94, text='S0325ab')
fig7.text(x=169, y=6.5, text='S0173ab')

# for ellipse in [s0173a, s0235b,s0325ab]:
#     fig7.velo(data=ellipse, uncertaintycolor="orange", pen="black", line=False,  spec="r0.2/0.95")
    
# # [[lon, lat, direction, major_axis, minor_axis]]
# data = [[164.09, 3.95, 90, 2, 1]]
# data2 = [[164.09, 3.95, 90, 2, 1.25]]
# data3 = [[164.09, 3.95, 90, 2, 1.5]]
# fig7.plot(data=data, style="e", color="orange", transparency=50)
# fig7.plot(data=data2, style="e", color="orange", transparency=70)
# fig7.plot(data=data3, style="e", color="orange", transparency=80)

# def plot_ellipse(value_ls):
#     fig7.plot(data=value_ls[2], style="e", color="darkorange", transparency=30)
#     fig7.plot(data=value_ls[1], style="e", color="orange", transparency=40)
#     fig7.plot(data=value_ls[0], style="e", color="yellow", transparency=90)


#836CD9
def plot_ellipse(value_ls):
    fig7.plot(data=value_ls[2], style="e", color="#bb3b9d", transparency=70)
    fig7.plot(data=value_ls[1], style="e", color="#9e54bc", transparency=60)
    fig7.plot(data=value_ls[0], style="e", color="#836CD9", transparency=50)
    
plot_ellipse([s0173a_d,a_1,a_2])
plot_ellipse([s0235b_d,b_1,b_2])
plot_ellipse([s0325ab_d,ab_1,ab_2])
plot_ellipse([s0173ab,abb_1,abb_2])

# for ellipse in [s0173a, s0235b,s0325ab]:
#     fig7.velo(data=ellipse, uncertaintycolor="#975ac3", pen="black", line=False,  spec="r0.2/0.95")

# fig7.velo(data=s0173a, uncertaintycolor='orange', pen="black", line=False,  spec="r0.2/0.95")

fig7.show()
