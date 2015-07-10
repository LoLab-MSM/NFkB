__author__ = 'geena'

from pysb import *
from NFkB_Messy_Model import model
from pysb.integrate import odesolve
import numpy as np
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *

@with_model
def test_NFkB():
    time = np.linspace(0, 18000, 1801)
    x = odesolve(model, time, verbose=True)
    expectedx_final = {'IkBan_obs': 3145.8443053128276, 'NFkBn_obs': 178.06326194655537, 'IkBapc_NFkBpc_obs': 1.3899039187480429e-07, 'IkBa_on_obs': 0.045201392272254615, 'IKKKa_obs': 0.00054704675357157236, 'A20_off_obs': 1.9547986077277453, 'A20t_obs': 6.2762670477040663, 'IkBan_NFkBn_obs': 28.012691613691331, 'IKKn_obs': 199999.9951281352, 'IKKa_obs': 2.3395704210536575e-08, 'TNFR1a_obs': 5.3362282855621739e-06, 'IkBa_off_obs': 1.9547986077277459, 'TNF_ext_obs': 9.5649122910157076e-08, 'A20_on_obs': 0.045201392272254573, 'IkBa_p_obs': 2.3942600622297971e-09, 'IKKii_obs': 0.0045466692017435365, 'IKKKn_obs': 99999.999452953241, 'IKKi_obs': 0.00032517206132429724, 'TNFR1i_obs': 1999.9999946637713, 'IkBat_obs': 6.2762670477040681, 'NFkB_c_obs': 139.40965738805579, 'IkBa_obs': 8607.2535257899253, 'A20_obs': 5714.7212868294837, 'IkBa_NFkB_obs': 99655.514388905125}
    for i in expectedx_final.keys():
        final_obs = x[i][-1]
        assert abs(final_obs - expectedx_final[i]) < 1e-3
