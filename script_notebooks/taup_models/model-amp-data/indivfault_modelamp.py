def indivfault(event, depth, dist, eventdata, az, st, dp, rk):
    mtimes = mars.get_travel_tiems(source_depth_in_km = depth, distance_in_degree = dist, phase_list=["P", "S"])

    Pp = mtimes[0].ray_param; Pa = mtimes[0].incident_angle
    Pvel = radius*np.sin(np.radians(Pa))/Pp
    print('Check: P surface velocity = ',Pvel,' ?') #check w/ the model file of what the vel is at the surface
    Sp = mtimes[1].ray_param; Sa = mtimes[1].incident_angle
    Svel = radius*np.sin(np.radians(Sa))/Sp
    print('Check: S surface velocity = ',Svel,' ?')

    eventdata, Pe, Se, = getfault(az, st, dp, rk)
    eventdata.to_csv(event + '-' + str(depth) + '.csv', index=False)

    return Pa, Sa, Pe, Se
