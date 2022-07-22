import numpy as np

eps=1.e-6
rad=180./np.pi; halfpi = 0.5*np.pi; twopi = 2.0*np.pi

def azdp(v):
       vr=v[0]; vt=v[1]; vp=v[2]
       dparg = np.sqrt(vt*vt + vp*vp)
       if dparg>1.:
          dparg = 1.
          print('argument error for np.arccos: ',dparg,'  Set to 1.')
       if dparg<-1.:
          dparg = -1.
          print('argument error for np.arccos: ',dparg,'  Set to -1.')
       vdp = np.arccos(dparg)
       dp = halfpi - vdp
       vaz = halfpi + np.arctan2(vt,vp)
       if vr>0.: vaz = vaz + np.pi
       st = vaz + halfpi
       if st>=twopi: st = st - twopi
       if vaz>=twopi: vaz = vaz - twopi
       vaz = vaz*rad; st = st*rad; vdp = vdp*rad; dp = dp*rad
       return vaz,vdp,st,dp

def getmt(fault):
    #converts DC to MT
    rad = 180./np.pi; m0 = 1
    st,dp,rk = fault

    st = st/rad; dp = dp/rad; rk = rk/rad
    sd = np.sin(dp); cd = np.cos(dp)
    sd2 = np.sin(2*dp); cd2 = np.cos(2*dp)
    ss = np.sin(st); cs = np.cos(st)
    ss2 = np.sin(2*st); cs2 = np.cos(2*st)
    sr = np.sin(rk); cr = np.cos(rk)

    # formulas from Aki & Richards Box 4.4
    # mt(1-6): Mrr, Mtt, Mff, Mrt, Mrf, Mtf
    # mt(1-6): Mzz, Mxx, Myy, Mzx, -Mzy, -Mxy
    mt = [ sr*sd2, -1.*sd*cr*ss2 - sd2*sr*ss*ss, sd*cr*ss2 - sd2*sr*cs*cs, -1.*cd*cr*cs - cd2*sr*ss, cd*cr*ss - cd2*sr*cs, -1.*sd*cr*cs2 - 0.5*sd2*sr*ss2]
    return mt

def getplanes(xm):
    """
    needs function azdp. - converts MT to DC
    IN: xm = list of moment tensor elements in Harvard order (GCMT)
    OUT: strike dip rake (twice). Also: P and T vectors
    """
    xmatrix = np.array([[xm[0],xm[3],xm[4]],[xm[3],xm[1],xm[5]],[xm[4],xm[5],xm[2]]])
    tr=(xm[0]+xm[1]+xm[2])/3.
    if np.abs(tr)>eps:
       xmatrix[0,0] = xmatrix[0,0] - tr
       xmatrix[1,1] = xmatrix[1,1] - tr
       xmatrix[2,2] = xmatrix[2,2] - tr
    #print('removed isotropic component from Moment Tensor:')
    d,pt = np.linalg.eigh(xmatrix)
    jt = np.argmax(d) ; dmax = d[jt]
    jp = np.argmin(d) ; dmin = d[jp]
    for j in [0,1,2]:
        if j!=jp and j!=jt: jn=j
    if (jn+jp+jt)!=3:
        print('ERROR in axis determination')
        return 0

    p = pt[:,jp]
    t = pt[:,jt]
    n = pt[:,jn]
    if p[0] < 0.: p = -1.*p
    if t[0] < 0.: t = -1.*t
    pole1 = (t+p)/np.sqrt(2.)
    pole2 = (t-p)/np.sqrt(2.)
    if p[0] > t[0]: pole1 = -1.*pole1
    # planes' poles not part of function output, but they could be in future
    azt,dpt,st,dp = azdp(t)
    azn,dpn,st,dp = azdp(n)
    azp,dpp,st,dp = azdp(p)

    az1,dp1,st1,dip1 = azdp(pole1)
    az2,dp2,st2,dip2 = azdp(pole2)

    if -1.*d[jp]>d[jt]:
       djpt = d[jp]
    else:
       djpt = d[jt]
    clvd = d[jn]/djpt

    m0 = 0.5*(np.abs(d[jp])+np.abs(d[jt]))

    x = np.array([0.,-1*np.cos(st1/rad),np.sin(st1/rad)])
    vfin = np.dot(pole2,x)
    if vfin>1.:
       vfin = 1.
    if vfin<-1.:
       vfin = -1.
    rake1 = rad*np.arccos(vfin)
    if pole2[0]<0.: rake1 = -1.*rake1


    x = np.array([0.,-1*np.cos(st2/rad),np.sin(st2/rad)])
    vfin = np.dot(pole1,x)
    if vfin>1.:
       vfin = 1.
    if vfin<-1.:
       vfin = -1.
    rake2 = rad*np.arccos(vfin)
    if pole1[0]<0.: rake2 = -1.*rake2
    return 3*tr,clvd, m0,(azt,dpt),(azn,dpn),(azp, dpp), (st1,dip1,rake1), (st2,dip2,rake2)

def sumtensor(scale1, mt1, scale2, mt2):
    scale_mt1 = [i*scale1 for i in mt1]
    scale_mt2 = [i*scale2 for i in mt2]

    sum_mt = [a+b for a,b in zip(scale_mt1,scale_mt2)]
    # print(sum_mt)
    return sum_mt


# #--173a THRUST - AUG
# mt173a = getmt([2,24,86])
# #strikeslip ish
# mt173a_alt = getmt([170,68,68])
#
# #---173ab STRIKESLIP - OCT
# mt173ab = getmt([30,48,-12])
#
#
# summed = sumtensor(0.4,mt173a,0.6,mt173ab)
# print(getplanes(summed))

#325aa
#---Insight---
mt_in = getmt([166,80,-40])
#----Ours----
mt_our = getmt([348,34,-122])

summed = sumtensor(0.5, mt_in, 0.5, mt_our)
print(getplanes(summed))
