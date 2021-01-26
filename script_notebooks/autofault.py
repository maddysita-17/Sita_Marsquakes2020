def autofault(df, ratio1, ratio2, ratio3):
    n = 0
    for ratio in [ratio1, ratio2, ratio3]:
        rmax = ratio - (0.1 * ratio)
        rmin = ratio + (0.1 * ratio)
        if n=1:
            ratio1df = df[df['SH/SV']>-rmin and df['SH/SV']<rmax]
        elif n = 2:
            ratio2df = df[df['P/SV']>-rmin and df['P/SV']<rmax]
        elif = 3:
            ratio3df = df[df['P/SH']>-rmin and df['P/SH']<rmax]
        n += 1

    frames = [ratio1df, ratio2df, ratio3df]
    faults = pd.concat(frames)

    return faults

def eventbuild(dist):
    mtimes = mars.get_travel_times(source_depth_in_km = depth, distance_in_degree = dist, phase_list=["P", "S"])

    #incident angle at the station
    Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
    Pvel = radius*np.sin(np.radians(Pa))/Pp
    Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
    Svel = radius*np.sin(np.radians(Sa))/Sp

    return Pa, Sa, Pe, Se

def fault_search(event):
    path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/faultdata/' + event + '/'
    source_files = sorted(Path(path).glob('*.csv'))

    dataframes = []
    for file in source_files:
        df = pd.read_csv(file) # additional arguments up to your need
        df['source'] = file.name
        dataframes.append(df)

    faults = pd.concat(dataframes)
    faults.to_csv(path + "faults_" + event + '.csv')
