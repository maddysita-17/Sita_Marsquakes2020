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
    if n==1 or d==1:
        a = (np.arctan2(n, d))**2
        b = (d**2/(d**2 + n**2))**2
        ab = a*b

        t1 = ((n**2)/(d**2))*((sign**2)/(n**2))
        t2 = ((sigd**2)/(d**2))*((n**2)/(d**2))
        m = sqrt(t1+t2)

        return ab*m

    else:
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

print('173a')
P173a = 195; SH173a = 144; SV173a = -176
all_error(P173a, SH173a, SV173a, 18, 47)

print('235b')
P235b = -80; SH235b = -318; SV235b = -234
all_error(P235b, SH235b, SV235b, 10, 26)

print('235bi')
iP235b = -80; iS235b = 1
all_error(iP235b, iS235b, iS235b, 10, 26)

print('325a')
P325a = 132; SH325a = 1; SV325a = 257
all_error(P325a, SH325a, SV325a, 26, 93)

print('325ab')
P325ab = 228; SH325ab = -360; SV325ab = 364
all_error(P325ab, SH325ab, SV325ab, 26, 93)

print('173ab')
P173ab = -286; SH173ab = 561; SV173ab = 407
all_error(P173ab, SH173ab, SV173ab, 18, 47)

# 173a
# PSH:  0.01251767423775118
# PSV:  0.0854829422800005
# SHSV:  0.38486910966388066

# 235b
# PSH:  0.16539870015659908
# PSV:  0.17732555103975575
# SHSV:  0.011418827906242521

# 235bi
# PSH:  4.006466690869815e-05
# PSV:  4.006466690869815e-05
# SHSV:  208.49539297301268

# 325a
# PSH:  6.960859599569113e-05
# PSV:  0.023926528580082356
# SHSV:  0.1309449311459457

# 325ab
# PSH:  0.2697542935333619
# PSV:  0.012644457747317964
# SHSV:  0.020518775419130133

# 173ab
# PSH:  0.0015375149263553327
# PSV:  0.0029083317634612115
# SHSV:  0.0021530793596216697
