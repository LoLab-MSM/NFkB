
"""Demonstration of the move_connected keyword for Rules
"""

from __future__ import print_function
from pysb import *
import numpy as np
from pysb.simulator import ScipyOdeSimulator
# from pysb.logging import setup_logger
from pysb.core import *
import matplotlib.pyplot as plt


Model()

Monomer('A', ['b'])
Monomer('B', ['a'])

# One parent 3d compartment, a 2D membrane, and a 3D volume
Parameter('Vmain', 1)
Parameter('V2D', 1)
Parameter('V3D', 1)
Compartment('P', None, 3, Vmain)
Compartment('M', P, 2, V2D)
Compartment('V', M, 3, V3D)

# Compartment('cytosol', dimension=3, size=cyto_vol, parent=ec_membrane)

# A and B, both embedded in Parent compartment, bind reversibly and release A
Parameter('Kab_f', 1)
Parameter('Kab_r', 1)
Rule('AP_bind_BP', A(b=None) ** P + B(a=None) ** V <> A(b=1) **M  % B(a=1) ** M,
     Kab_f, Kab_r)
Rule('ABP_unbind_AP',A(b=1) ** P % B(a=1) ** P >> A(b=None) ** P,
     Kab_f)

# The A:B complex is transported back and forth from X to Y
Parameter('Ktrans_f', 1)
Parameter('Ktrans_r', 1)
Parameter('Ktrans_v', 1)
# move_connected is required or B will be "left behind" and BNG will complain
# (change move_connected to False and run pysb.tools.export_bng_net on this file
# and watch for the WARNING line in the output log)
Rule('A1_A2', A(b=ANY) ** P <> A(b=ANY) ** M,
     Ktrans_f, Ktrans_r, move_connected=True)
Rule('A2_A3', A(b=ANY) ** M >> A(b=ANY) ** V,
     Ktrans_v,move_connected=True)
Rule('A_deg', A(b = None)**P >> None, Ktrans_v)

# A and B bind and unbind, and release A in the 3D volume
Parameter('Vab_f', 1)
Parameter('Vab_r', 1)
Rule('AV_bind_BV', A(b=None) ** V + B(a=None) ** V <> A(b=1) ** V % B(a=1) ** V,
     Vab_f, Vab_r)
Rule('ABV_unbind_AV',A(b=1) ** V % B(a=1) ** V >> A(b=None) ** V,
     Vab_f)

Parameter('AP_0', 1)
Initial(A(b=None)**P, AP_0)
Parameter('BP_0', 1)
Initial(B(a=None)**P, BP_0)

Observable('AP_obs', A(b=None) ** P)
Observable('AM_obs', A(b=ANY) ** M)
Observable('AV_obs', A(b=ANY) ** V)
# Observable('Ax_any_obs', A(b=ANY) ** X)
Observable('BP_obs', B(a=None) ** P)
Observable('BV_obs', B(a=None) ** V)
Observable('APBP_obs', A(b=1) ** P % B(a=1) ** P)
Observable('AVBV_obs', A(b=1) ** V % B(a=1) ** V)


tspan = np.linspace(0, 60, 61)
# x = odesolve(model,tspan, verbose=False)
sim = ScipyOdeSimulator(model, tspan=tspan)
sim_result = sim.run()

plt.figure(figsize = (15,5))
# plt.figure()
plt.subplot(131)
plt.plot(tspan, sim_result.observables['AP_obs'], marker = '*',label = 'Acyto', lw = 4)
plt.plot(tspan, sim_result.observables['AV_obs'],marker = '*', label = 'Anuc', lw = 4)
# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in sec)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(132)
plt.plot(tspan, sim_result.observables['BP_obs'],marker = '*',label = 'Bcyto', lw = 4)
plt.plot(tspan, sim_result.observables['BV_obs'],marker = '*', label = 'Bnuc', lw = 4)
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
plt.xlabel("Time (in sec)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(133)
plt.plot(tspan, sim_result.observables['APBP_obs'],marker = '*', label = 'ABcyto', lw = 4)
plt.plot(tspan, sim_result.observables['AVBV_obs'],marker = '*', label = 'ABnuc', lw = 4)
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in sec)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.subplot(234)
# plt.plot(tspan, sim_result.observables['AV_obs'],marker = '*', label = 'AV', lw = 4)
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in sec)", fontsize=10)
# plt.ylabel("Concentration uM", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(235)
# plt.plot(tspan, sim_result.observables['BV_obs'],marker = '*', label = 'BV', lw = 4)
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in sec)", fontsize=10)
# plt.ylabel("Concentration uM", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(236)
# plt.plot(tspan, sim_result.observables['AVBV_obs'],marker = '*', label = 'AVBV', lw = 4)
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in sec)", fontsize=10)
# plt.ylabel("Concentration uM", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)


plt.tight_layout()
plt.show()

