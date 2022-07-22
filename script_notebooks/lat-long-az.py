import numpy as np

def azdelt(deld,az):
    """
    Ccompute destination lat,lon, given a starting lat,lon, bearing (azimuth), and a distance to travel
    Travel is along a great circle.
    IN:  lat1, lon1: lattitude na dlongitude of starting point
     deld and az: distance and azimuth (both in degrees) (azimuth is measured clockwise from N)
     OUT: lat2, lon2 : latitude and longitude of destination point (in degrees)
    """

    teta1 = np.radians(4.5024)
    fi1 = np.radians(135.6234)
    delta = np.radians(deld)
    azimuth = np.radians(az)

    if teta1 > 0.5*np.pi or teta1 < -0.5*np.pi:
       print('error, non-existent latitude')
       return 0,0

    if delta < 0.:
       print('error, non-existent delta')
       return 0,0

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
    fi2 = np.arctan2(teller,rnoemer)

    term1 = np.cos(delta)*np.sin(teta1)
    term2 = np.sin(delta)*np.cos(azimuth)*np.cos(teta1)
    som = term1 + term2
    teta2 = np.arcsin(som)

    return np.degrees(teta2), np.degrees(fi2)


def deltaz(deta2,di2):
    '''
    IN: long1, lat1, long2, lat2
    OUT: distance in degrees, azimuth to go from point 1 to point 2 and azimuth to go from point 2 to point 1
    '''
    teta1 = np.radians(4.5024); fi1 = np.radians(135.6234)
    teta2 = np.radians(deta2); fi2 = np.radians(di2)
    c1 = np.cos(teta1); c2 = np.cos(teta2)
    s1 = np.sin(teta1); s2 = np.sin(teta2)
    c21 = np.cos(fi2-fi1); s21 = np.sin(fi2-fi1)
    som  = s1*s2 + c21*c1*c2
    delta = np.arccos(som)
    teller = c2*s21
    rnoemer =  s2*c1 - c2*s1*c21
    azimuth = np.arctan2(teller,rnoemer)
    numerator = -s21*c1
    denominator =  s1*c2 - c1*s2*c21
    baz = np.arctan2(numerator,denominator)
    return np.degrees(delta), np.degrees(azimuth), np.degrees(baz)

lat, long = azdelt(38.4, 57)
dist, bAz, az = deltaz(lat, long)
print('dist in deg: ', dist)
print('az in deg: ', az)
print('baz in deg: ', bAz)

#IU TEIG
# dist in deg:  67.25018767660207
# az in deg:  178.6120724692286
# baz in deg:  0.045865040487818594


#----location of the lander: lat=4.5024, long=135.6234-------
#S0235b should be at lat=10, long=160
#S0172a should be at lat=5, long=162

# lat235, long235 = azdelt(27.1, 72)
# print(lat235, long235)
# d235, az235, baz235 = deltaz(lat235, long235)
# print(d235, az235, baz235)
#code prints lat=11.33, long=162.5

#
# lat173, long173 = azdelt(28.4, 90)
# print(lat173, long173)
# d173, az173, baz173 = deltaz(lat173, long173)
# print(d173, az173, baz173)
# #code prints lat=3.45, long=164.7
#
# lat325aa, long325aa = azdelt(38.4, 123)
# print(lat325aa, long325aa)
# d325aa, az325aa, baz325aa = deltaz(lat325aa, long325aa)
# print(d325aa, az325aa, baz325aa)
# #lat = -16.055244971011064, long = 168.53025369688243
#
# lat325ab, long325ab = azdelt(38.4, 139)
# print(lat325ab, long325ab)
# d325ab, az325ab, baz325ab = deltaz(lat325ab, long325ab)
# print(d325ab, az325ab, baz325ab)
# #lat = -7.082705874675253, long = 167.22546622753717
#
# lat173ab, long173ab = azdelt(28.4, 86)
# print(lat173ab, long173ab)
# d173ab, az173ab, baz173ab = deltaz(lat173ab, long173ab)
# print(d173ab, az173ab, baz173ab)
# #lat = 5.859355684400125, long = 164.0099159152684
#
# print('alt 235b')
# lat235, long235 = azdelt(15, 74)
# print(lat235, long235)
# d235, az235, baz235 = deltaz(lat235, long235)
# print(d235, az235, baz235)

# lat183a, long183a = azdelt(43.4, 55)
# print(lat183a, long183a)
# d183a, az183a, baz183a = deltaz(lat183a, long183a)
# print(d183a, az183a, baz183a)
#lat = 13.551500461355511, long = 178.67695267837766

# lat, long
# dist, baz (matches baz_calc.ipynb), az (input for fault-guess.py)

#235b
# 12.135196511781126 161.92853212353742
# 27.100000000000005 71.99999999999999 -104.11873538665255

#235b alt @ dist=15
# 8.449986768393497 150.19106815416754
# 15.000000000000018 73.99999999999999 -104.34713212886746

#173a
# 3.959605944799909 164.09755205854992
# 28.39999999999999 90.0 -87.85513325237781

#325a
# -15.907168040890824 168.2588561748905
# 38.19999999999999 123.00000000000001 -60.38590124628378

#325ab
# -23.942455195258926 162.1027760883247
# 38.399999999999984 139.00000000000003 -45.6932308101358

#173ab
# 5.861760657364189 164.11041225165621
# 28.39999999999999 86.0 -91.37294978832101

#183a
# 26.738408454295232 174.689637707948
# 43.4 54.99999999999998 -113.87918160927178
