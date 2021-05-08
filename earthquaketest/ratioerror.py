import numpy as np

def error(r_o, d, n, sig_d, sig_n):
    a_o = np.arctan(r_o)
    t1 = a_o **2

    t2 = (d**2)/(d**2 + n**2)

    e_n = (sig_n**2)/(n**2)
    e_d = (sig_d**2)/(d**2)
    t3 = e_n + e_d

    return t1*t2*t3

def error2(n, d, sign, sigd):
    a = (np.arctan2(n, d))**2
    b = (d**2/(d**2 + n**2))**2
    ab = a*b

    e_n = (sign**2)/(n**2)
    e_d = (sigd**2)/(d**2)
    c = e_n + e_d

    return ab*c

def all_error(P, SH, SV, sigP, sigS):
    PSH = error2(P, SH, sigP, sigS)
    PSV = error2(P, SV, sigP, sigS)
    SHSV = error2(SH, SV, sigS, sigS)

    print('PSH: ', PSH)
    print('PSV: ', PSV)
    print('SHSV: ', SHSV)

print('Uganda')
Puga = 2108; SHuga = -3721; SVuga = -2415
all_error(Puga, SHuga, SVuga, 90, 319)

# PSH:  0.0012766208291823426
# PSV:  0.0049391396619344686
# SHSV:  0.0011917374519192156

# PSH:  0.03625428788888845
# PSV:  0.03647407005519762
# SHSV:  0.010036175551812293
