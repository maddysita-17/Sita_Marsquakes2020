import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path = '/Users/maddysita/Desktop/CIERA_REU/script_notebooks/beachballs/csvs/'

def abserr(errx, erry, x, y):
    if x == 1 or y == 1:
        slope = x/y
        sig1 = ((x/y)**2)*((errx/x)**2)
        sig2 = ((erry/y)**2)*((x/y)**2)
        sig = sig1 + sig2
        #print('Error: ', sig)

        xarr = np.linspace(-3, 3, 11)

        e1 = slope+sig ; e2 = slope-sig
        eline1 = e1*xarr ; eline2=e2*xarr

    else:
        slope = y/x
        #print('Slope: ', slope)
        sig2 = (errx/x)**2 + (erry/y)**2
        sig = np.sqrt(sig2) * slope
        #print('Error: ', sig)

        xarr = np.linspace(-3, 3, 11)

        e1 = slope+sig ; e2 = slope-sig
        eline1 = e1*xarr ; eline2=e2*xarr

    return eline1, eline2

def misfit3D(eventfile, obsP, obsSV, obsSH, errP, errSV, errSH):
    # hypothetically observed amplitudes P, SV, SH
    vobserved = np.array([obsP,obsSV,obsSH])
    vobslength = np.linalg.norm(vobserved)
    # and errors (always positive):
    eobs = np.array([errP, errSV, errSH])
    # normalized:
    n = vobserved/vobslength

    # eca is the estimated error angle that forms a cone around the observed vector of amplitudes. (in radians)
    # it is calculated by projecting the error in each dimension onto the plane that is perpendicular to v_observed,
    # followed by choosing the maximum of the three projected errors.
    eca = np.arctan(np.max([eobs[0]*(1-n[0]),eobs[1]*(1-n[1]),eobs[2]*(1-n[2])])/vobslength)
    # This is the cut-off value for misfit value mf3D
    # mf3D is the best choice for misfit becasue it works in 3D and allows the error cone estimation.
    # that said, mfd should not be much different from mf3D, except it is independent of distance to origin (unlike mf3D).
    print('cutoff value: ', eca)

    #read in modeled csv file
    event = pd.read_csv(path + str(eventfile) + '.csv')
    xd = event['P']
    yd = event['SV']
    zd = event['SH']
    print(len(zd))
    ncalc = len(zd)

    vcall = np.array([xd,yd,zd])
    vca = vcall.T

    # misfit:
    mfd = np.zeros(ncalc)
    mf3D = np.zeros(ncalc)
    mf1 = np.zeros(ncalc)
    mf2 = np.zeros(ncalc)
    mf3 = np.zeros(ncalc)
    select = []

    for i in np.arange(ncalc):
        # original arctan based misfits (3 different ones):
        mf1[i] = np.sin(0.5*(np.arctan2(vca[i,0],vca[i,1]) - np.arctan2(n[0],n[1])))**2   # P/SV
        mf2[i] = np.sin(0.5*(np.arctan2(vca[i,0],vca[i,2]) - np.arctan2(n[0],n[2])))**2   # P/SH
        mf3[i] = np.sin(0.5*(np.arctan2(vca[i,1],vca[i,2]) - np.arctan2(n[1],n[2])))**2   # SV/SH
        # angle in 3 dimensions: (in radians)
        mf3D[i] = np.arccos(np.dot(n,vca[i])/np.linalg.norm(vca[i]))   # should be valued between 0 and pi
        if mf3D[i] < eca:
            select.append(i)
        # shortest distance between calculated point (P,SV,SH) to observed ratio line:
        mfd[i] = np.linalg.norm(vca[i] - np.dot(vca[i],n)*n)

    #print('selected: ',select)

    # error calculations:
    e1PSV, e2PSV = abserr(errSV, errP, obsSV, obsP)
    e1PSH, e2PSH = abserr(errSH, errP, obsSH, obsP)
    e1SVSH, e2SVSH = abserr(errSH, errSV, obsSH, obsSV)


    # plotting:
    hor = np.linspace(-3,3,11)
    selsiz = 144
    selcol = 'red'

    plt.subplot(331)
    plt.ylabel('P'); plt.xlabel('SV')
    slope = n[0]/n[1]   # P/SV
    ver = hor*slope
    plt.plot(hor,ver)
    plt.scatter(yd,xd,c=mf1)
    plt.fill_between(hor, e1PSV, e2PSV, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])

    plt.subplot(332)
    plt.ylabel('P'); plt.xlabel('SH')
    slope = n[0]/n[2]   # P/SH
    ver = hor*slope
    plt.plot(hor,ver)
    plt.scatter(zd,xd,c=mf2)
    plt.fill_between(hor, e1PSH, e2PSH, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])
    plt.title('2D angle')

    plt.subplot(333)
    plt.ylabel('SV'); plt.xlabel('SH')
    slope = n[1]/n[2]   # SV/SH
    ver = hor*slope
    plt.plot(hor,ver)
    plt.scatter(zd,yd,c=mf3)
    plt.fill_between(hor, e1SVSH, e2SVSH, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])

    plt.subplot(334)
    plt.ylabel('P'); plt.xlabel('SV')
    slope = n[0]/n[1]   # P/SV
    ver = hor*slope
    plt.plot(hor,ver)
    for i in select:
        plt.scatter(yd[i],xd[i],s=selsiz,c=selcol,marker='s')
        st = event['Strike'].values[i] ; dp = event['Dip'].values[i] ; rk = event['Rake'].values[i]
        mech = (st,dp,rk)
        plt.annotate(f'{mech}', (yd[i], xd[i]), (yd[i]+0.001, xd[i]+0.001), size=4)
    plt.scatter(yd,xd,c=mf3D)
    plt.fill_between(hor, e1PSV, e2PSV, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])

    plt.subplot(335)
    plt.ylabel('P'); plt.xlabel('SH')
    slope = n[0]/n[2]   # P/SH
    ver = hor*slope
    plt.plot(hor,ver)
    for i in select:
        plt.scatter(zd[i],xd[i],s=selsiz,c=selcol,marker='s')
        st = event['Strike'].values[i] ; dp = event['Dip'].values[i] ; rk = event['Rake'].values[i]
        mech = (st,dp,rk)
        plt.annotate(f'{mech}', (zd[i], xd[i]), (zd[i]+0.001, xd[i]+0.001), size=4)
    plt.scatter(zd,xd,c=mf3D)
    plt.fill_between(hor, e1PSH, e2PSH, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])
    plt.title('3D angle (best misfit definition)')

    plt.subplot(336)
    plt.ylabel('SV'); plt.xlabel('SH')
    slope = n[1]/n[2]   # SV/SH
    ver = hor*slope
    plt.plot(hor,ver)
    for i in select:
        plt.scatter(zd[i],yd[i],s=selsiz,c=selcol,marker='s')
        st = event['Strike'].values[i] ; dp = event['Dip'].values[i] ; rk = event['Rake'].values[i]
        mech = (st,dp,rk)
        plt.annotate(f'{mech}', (zd[i], yd[i]), (zd[i]+0.001, yd[i]+0.001), size=4)
    plt.scatter(zd,yd,c=mfd)
    plt.fill_between(hor, e1SVSH, e2SVSH, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])

    plt.subplot(337)
    plt.ylabel('P'); plt.xlabel('SV')
    slope = n[0]/n[1]   # P/SV
    ver = hor*slope
    plt.plot(hor,ver)
    plt.scatter(yd,xd,c=mfd)
    plt.fill_between(hor, e1PSV, e2PSV, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])

    plt.subplot(338)
    plt.ylabel('P'); plt.xlabel('SH')
    slope = n[0]/n[2]   # P/SH
    ver = hor*slope
    plt.plot(hor,ver)
    plt.scatter(zd,xd,c=mfd)
    plt.fill_between(hor, e1PSH, e2PSH, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])
    plt.title('shortest distance')

    plt.subplot(339)
    plt.ylabel('SV'); plt.xlabel('SH')
    slope = n[1]/n[2]   # SV/SH
    ver = hor*slope
    plt.plot(hor,ver)
    plt.scatter(zd, yd,c=mf3D)
    plt.fill_between(hor, e1SVSH, e2SVSH, color='b', alpha=0.1)
    plt.xlim([-3,3]);plt.ylim([-3,3])

    plt.colorbar()
    plt.show()

    #new possible faulting mechanims based on 3D misfit cutoff values
    st_ls = []; dp_ls = []; rk_ls = []; mf_ls = []
    for i in select:
        st = event['Strike'].values[i] ; dp = event['Dip'].values[i] ; rk = event['Rake'].values[i]
        mf = mf3D[i]
        st_ls.append(st); dp_ls.append(dp); rk_ls.append(rk); mf_ls.append(mf)
    data = {'Strike': st_ls, 'Dip': dp_ls, 'Rake': rk_ls, 'Misfit Val': mf_ls}
    misfitfile = pd.DataFrame(data, columns = ['Strike', 'Dip', 'Rake', 'Misfit Val'])
    misfitfile.to_csv(path + 'mf3d_' + str(eventfile) + '.csv', index=False)

print('173a')
misfit3D('S0173a', 195, -176, 144, 18, 47, 47)

print('173ab')
misfit3D('S0173ab', -286, 407, 561, 18, 47, 47)

print('235b')
misfit3D('S0235b', -80, -234, -318, 10, 26, 26)

print('325a')
misfit3D('S0325a', 132, 257, 1, 26, 93, 93)

print('325ab')
misfit3D('S0325ab', 229, 364, -360, 26, 93, 93)
