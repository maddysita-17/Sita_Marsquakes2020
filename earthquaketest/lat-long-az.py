import numpy as np

def azdelt(deld,az):
    """
    Ccompute destination lat,lon, given a starting lat,lon, bearing (azimuth), and a distance to travel
    Travel is along a great circle.
    IN:  lat1, lon1: lattitude na dlongitude of starting point
     deld and az: distance and azimuth (both in degrees) (azimuth is measured clockwise from N)
     OUT: lat2, lon2 : latitude and longitude of destination point (in degrees)
    """

    teta1 = np.radians(-0.6019); fi1 = np.radians(30.7382)
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


def deltaz(deta1,di1,deta2,di2):
    '''
    IN: long1, lat1, long2, lat2
    OUT: distance in degrees, azimuth to go from point 1 to point 2 and azimuth to go from point 2 to point 1
    '''
    teta1 = np.radians(deta1); fi1 = np.radians(di1)
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



#----location of the station - lat: -0.6019, long: 30.7382 -------
#----location of earthquake (march) - lat: 2.605, long: 31.074 -----------

dist, bAz, az = deltaz(-29.04,167.93,-21.104, 174.895)
print('dist in deg: ', dist)
print('az in deg: ', az)
print('baz in deg: ', bAz)

# -- March --
# dist in deg:  3.224423257004671
# az in deg:  -174.01942403450002
# baz in deg:  5.974704048488569

# -- April --
# dist in deg:  3.1331929170081385
# az in deg:  -171.76100668517643
# baz in deg:  8.231560867241633

# -- Vanuatu Earthquake --
# dist in deg:  9.234710930645091
# az in deg:  -53.50397811409166
# baz in deg:  128.91428225834656

# AU Station --
# dist in deg:  10.132414687404298
# az in deg:  -142.94107469309904
# baz in deg:  40.020472861169274
