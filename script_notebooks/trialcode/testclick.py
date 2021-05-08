import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

x, y, c, s = np.random.rand(4, 100)
def onpick3(event):
    ind = event.ind
    print('onpick3 scatter:', ind, np.take(x, ind), np.take(y, ind))

fig = plt.figure()
ax1 = fig.add_subplot(111)
col = ax1.scatter(x, y, 100*s, c, picker=True)
#fig.savefig('pscoll.eps')
fig.canvas.mpl_connect('pick_event', onpick3)

plt.show()




# hypothetically observed amplitudes P, SV, SH
vobserved = np.array([195, -176, 144])
vobslength = np.linalg.norm(vobserved)
# and errors (always positive):
eobs = np.array([18, 47, 47])
# normalized:
n = vobserved/vobslength

eca = np.arctan(np.max([eobs[0]*(1-n[0]),eobs[1]*(1-n[1]),eobs[2]*(1-n[2])])/vobslength)

#read in modeled csv file
file = pd.read_csv('S0173a.csv')
xd = file['P']
yd = file['SV']
zd = file['SH']

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


selsiz = 144
selcol = 'red'

plt.subplot(331)
plt.ylabel('P'); plt.xlabel('SV')
plt.plot([0,vobserved[1]],[0,vobserved[0]]) # P/SV
for i in select:
    plt.scatter(yd[i],xd[i],s=selsiz,c=selcol,marker='s')
    st = file['Strike'].values[i] ; dp = file['Dip'].values[i] ; rk = file['Rake'].values[i]
    mech = (st,dp,rk)
    plt.annotate(f'{mech}', (yd[i], xd[i]), (yd[i]+0.001, xd[i]+0.001), size=4)
plt.scatter(yd,xd,c=mf3D)
plt.xlim([-3,3]);plt.ylim([-3,3])


plt.show()

st_ls = []; dp_ls = []; rk_ls = []; mf_ls = []
for i in select:
    st = file['Strike'].values[i] ; dp = file['Dip'].values[i] ; rk = file['Rake'].values[i]
    mf = mf3D[i]
    st_ls.append(st); dp_ls.append(dp); rk_ls.append(rk); mf_ls.append(mf)
data = {'Strike': st_ls, 'Dip': dp_ls, 'Rake': rk_ls, 'Misfit Val': mf_ls}
df = pd.DataFrame(data, columns = ['Strike', 'Dip', 'Rake', 'Misfit Val'])
print(df)
