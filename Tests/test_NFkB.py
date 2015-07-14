__author__ = 'geena'

from pysb import *
from test_NFkB import model
from pysb.integrate import odesolve
import numpy as np
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *

def test_simulation():
    time = np.linspace(0, 18000, 1801)
    x = odesolve(model, time, verbose=True)
    expectedx_final = {'IkBan_obs': 3145.8443053128276, 'NFkBn_obs': 178.06326194655537, 'IkBapc_NFkBpc_obs': 1.3899039187480429e-07, 'IkBa_on_obs': 0.045201392272254615, 'IKKKa_obs': 0.00054704675357157236, 'A20_off_obs': 1.9547986077277453, 'A20t_obs': 6.2762670477040663, 'IkBan_NFkBn_obs': 28.012691613691331, 'IKKn_obs': 199999.9951281352, 'IKKa_obs': 2.3395704210536575e-08, 'TNFR1a_obs': 5.3362282855621739e-06, 'IkBa_off_obs': 1.9547986077277459, 'TNF_ext_obs': 9.5649122910157076e-08, 'A20_on_obs': 0.045201392272254573, 'IkBa_p_obs': 2.3942600622297971e-09, 'IKKii_obs': 0.0045466692017435365, 'IKKKn_obs': 99999.999452953241, 'IKKi_obs': 0.00032517206132429724, 'TNFR1i_obs': 1999.9999946637713, 'IkBat_obs': 6.2762670477040681, 'NFkB_c_obs': 139.40965738805579, 'IkBa_obs': 8607.2535257899253, 'A20_obs': 5714.7212868294837, 'IkBa_NFkB_obs': 99655.514388905125}
    for i in expectedx_final.keys():
        final_obs = x[i][-1]
        assert abs(final_obs - expectedx_final[i]) < 1e-3

def test_ODEs():
    generate_equations(model, verbose = True)
    one = [(0, ':', -__s0*ki + __s15*__s19*keff), (1, ':', -0.5*__s0**2*__s1*k1 + __s20*k4), (2, ':', 0.5*__s0**2*__s1*k1 - __s2*__s8*k3_div_k2 - __s2*k3), (3, ':', __s2*__s8*k3_div_k2 + __s2*k3 - __s3*k4), (4, ':', __s10*__s2*a2 - __s4*tp), (5, ':', __s13*__s2*a3 - __s5*tp), (6, ':', -__s10*__s6*a1 + __s13*c6a + __s5*tp - __s6*i1), (7, ':', -__s11*__s7*a1_mult_kv + __s6*i1), (8, ':', -__s8*c5 + __s9*c4), (9, ':', __s16*c1 - __s9*c3), (10, ':', -__s10*__s2*a2 - __s10*__s6*a1 - __s10*c5a - __s10*i1a + __s11*e1a + __s12*c4), (11, ':', __s10*i1a - __s11*__s7*a1_mult_kv - __s11*e1a), (12, ':', -__s12*c3 + __s17*c1a), (13, ':', __s10*__s6*a1 - __s13*__s2*a3 - __s13*c6a + __s14*e2a), (14, ':', __s11*__s7*a1_mult_kv - __s14*e2a), (15, ':', -__s15*kf + __s18*__s21*kb), (16, ':', -__s11*__s16*q2 + __s22*__s7*q1), (17, ':', -__s11*__s17*q2 + __s23*__s7*q1), (18, ':', -Tdeg*__s18), (19, ':', __s0*ki - __s15*__s19*keff), (20, ':', -__s20*k4 + __s3*k4), (21, ':', __s15*kf - __s18*__s21*kb), (22, ':', __s11*__s16*q2 - __s22*__s7*q1), (23, ':', __s11*__s17*q2 - __s23*__s7*q1), (24, ':', 0), (25, ':', Tdeg*__s18 + __s10*c5a + __s12*c3 + __s4*tp + __s8*c5 + __s9*c3)]
    for i, odes in enumerate(model.odes):
        if (i,":",odes) == one:
            print(True)
        else:
            print(False)

        # print i,":",odes
        # assert(odes.equal_to(ode))