from pysb import *
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import matplotlib.pyplot as plt
import numpy as np
from sympy import sympify
from pysb.bng import generate_equations
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *

Model ()

#Deterministic approximation
#Number of active TNFR1 receptors B:

#Declaring the monomers
Monomer('TNF_ext') #s19
Monomer('TNFR1i') #22
Monomer('TNFR1a') #s16
Monomer('IKKKa') #s1
Monomer('IKKn') #s2
Monomer('IKKa') #s3
Monomer('IKKi') #s4
Monomer('IkBa_p') #5
Monomer('IkBapc_NFkBpc') #s6
Monomer('NFkB_c') #s7
Monomer('NFkBn') #s8
Monomer('A20') #s9
Monomer('A20t') #s10
Monomer('IkBac') #s11
Monomer('IkBan') #s12
Monomer('IkBat') #s13
Monomer('IkBac_NFkBc') #s14
Monomer('IkBan_NFkBn') #s15

Monomer('A20_on') #s17
Monomer('IkBa_on') #s18
Monomer('IKKKn') #s20
Monomer('IKKii') #21

Monomer('A20_off') #23
Monomer('IkBa_off') #24

Parameter('KNN', 2.0e5)
Parameter('KN', 1.0e5)
Parameter('M', 2000)
Parameter('AN', 2)
Parameter('ANa', 2)

#Declaring initial conditions
Initial(IKKn(), Parameter('IKKn_0', 2e5))
Initial(NFkBn(), Parameter('NFkBn_0', 1))
Initial(A20(), Parameter('A20_0', 10000))
Initial(A20t(), Parameter('A20t_0', 10))
Initial(IkBac(), Parameter('IkBa_0', 0.14*100000))
Initial(IkBan(), Parameter('IkBan_0', 0.06*100000))
Initial(IkBat(), Parameter('IkBat_0', 10))
Initial(IkBac_NFkBc(), Parameter('IkBa_NFkB_0', 100000))
Initial(TNF_ext(), Parameter('TNF_ext_0', 0.1)) #1e-1
Initial(IKKKn(), Parameter('IKKKn_0', KN.value))
Initial(IKKii(), Parameter('IKKii_0', KNN.value-IKKn_0.value))
Initial(TNFR1i(), Parameter('TNFR1i_0', M.value))
Initial(A20_off(), Parameter('A20_off_0', AN.value))
Initial(IkBa_off(), Parameter('ANa_0', ANa.value))

#Cell Parameters
Parameter('kv', 5.0) #Nucealr to cytoplasm volume

#Declaring Parameters
Parameter('kb', 1.2e-5) #receptor activation rate
Parameter('kf', 1.2e-3) #receptor inactivation rate
Parameter('Tdeg', 7.7e-4) #TNF loss same as cdeg
Parameter('ka',1e-5) #IKKK kinase activation rate
Parameter('ka20', 1e5) #A20 TNFR1 block
Parameter('ki', 0.01) #IKKK kinase inactivation rate
Parameter('k1', 2*6e-10) #IKKn activation by IKKK
Parameter('k3', 0.002) #IKKa inactivation by A20
Parameter('k2', 10000) #IKKa inactivation
Parameter('k4', 0.001) #IKKii transfer rate
Parameter('c1a', 0.1) #inducible IkBa mRNA synthesis

#IkBa Parameters
Parameter('a2', 1e-7) #IkBa phosphorylation b/c IKKa
Parameter('tp', 0.01) #degradation of phosph-IkBa complex with NFkB
Parameter('a3', 5e-7) #IkBa_NFkB phosphorylation b/c IKKa
Parameter('a1', 5e-7) #IkBa*NFkB association
Parameter('c6a', 0.00002) #spontaneous IkBa_NFkB defg of IkBa complex to NFkB
Parameter('c5a', 0.0001) #IkBa deg rate

#Transport Parameters
Parameter('i1a', 0.002) #IkBa nuclear import
Parameter('e1a', 0.005) #IkBa nuclear export
Parameter('e2a', 0.05) #IkBa_NFkB nuclear export
Parameter('i1', 0.01) #NFkB nuclear import

#A20 and IkBa synthesis Parameters
Parameter('c4', 0.5) #A20, IkBa transformation rate
Parameter('c5', 0.0005) #A20 degredation rate
Parameter('c1', 0.1) #inducible A20 mRNA synthesis
Parameter('c3', 0.00075) #A20 and IkBa mRNA deg rate
Parameter('q1', 4e-7) #NFkB attaching @ A20 and IkBa site
Parameter('q2', 1e-6) #IkBa inducible detaching from A20, IkBa site

Parameter('k3_div_k2', k3.value/k2.value) #for IKKa and A20
Parameter('a1_mult_kv', a1.value*kv.value) #for volume IkBa association NFkB

#Declaring expression
Observable('A20_obs', A20())
Expression('keff', ka*ka20/(ka20+A20_obs)) #10000 #michaelis menten

#Declaring rules
Rule('TNF_ext_and_TNFR1i', TNF_ext() + TNFR1i() >> TNF_ext() + TNFR1a(), kb)
Rule('TNFR1a_and_IKKKn', TNFR1a() + IKKKn() >> TNFR1a() + IKKKa(), keff)
Rule('TNFR1a_to_TNFR1i', TNFR1a() >> TNFR1i(), kf)
Rule('IKKKa_to_IKKKn', IKKKa() >> IKKKn(), ki)
Rule('IKKKa_and_IKKn', IKKKa() + IKKKa() + IKKn() >> IKKKa() + IKKKa() + IKKa(), k1) #proportion of IKKKa changes IKKn to IKKa
Rule('IKKa_to_IKKi', IKKa() >> IKKi(), k3) #IKKa active to IKKi inactive
Rule('IKKa_and_A20', IKKa() + A20() >> A20() + IKKi(), k3_div_k2)
Rule('IKKi_to_IKKii', IKKi() >> IKKii(), k4)
Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)
Rule('IkBapc_NFkBpc_to_NFkB_c', IkBapc_NFkBpc() >> NFkB_c(), tp)
Rule('IkBa_p_deg', IkBa_p() >> None, tp)
Rule('A20_deg', A20() >> None, c5)
Rule('A20t_deg', A20t() >> None, c3)
Rule('IkBa_cyt_deg', IkBac() >> None, c5a)
Rule('IkBac_NFkBc_to_NFkB_c', IkBac_NFkBc() >> NFkB_c(), c6a)
Rule('IkBat_deg', IkBat() >> None, c3)
Rule('NFkB_c_to_NFkBn', NFkB_c() >> NFkBn(), i1)
Rule('NFkB_c_and_IkBa', NFkB_c() + IkBac() >> IkBac_NFkBc(), a1)
Rule('NFkBn_and_IkBan', NFkBn() + IkBan() >> IkBan_NFkBn(), a1_mult_kv)
Rule('IKKa_and_IkBac_NFkBc', IKKa() + IkBac_NFkBc() >> IKKa() + IkBapc_NFkBpc(), a3)
Rule('IKKa_and_IkBa', IKKa() + IkBac() >> IKKa() + IkBa_p(), a2)
Rule('A20t_to_A20', A20t() >> A20t() + A20(), c4)
Rule('A20_on_to_A20t', A20_on() >> A20_on() + A20t(), c1)
Rule('IkBat_to_IkBan', IkBat() >> IkBat() + IkBac(), c4)
Rule('IkBac_to_IkBan', IkBac() >> IkBan(), i1a)
Rule('IkBan_to_IkBac', IkBan() >> IkBac(), e1a)
Rule('IkBa_on_to_IkBat', IkBa_on() >> IkBa_on() + IkBat(), c1a)
Rule('IkBan_NFkBn_to_IkBac_NFkBc', IkBan_NFkBn() >> IkBac_NFkBc(), e2a)
Rule('IkBan_and_A20_on', IkBan() + A20_on() >> IkBan() + A20_off(), q2)
Rule('IkBan_and_IkBa_on', IkBan() + IkBa_on() >> IkBan() + IkBa_off(), q2)
Rule('IKKii_to_IKKn', IKKii() >> IKKn(), k4)
Rule('NFkBn_and_A20_off', NFkBn() + A20_off() >> NFkBn() + A20_on(), q1)
Rule('NFkBn_and_IkBa_off', NFkBn() + IkBa_off() >> NFkBn() + IkBa_on(), q1)


generate_equations(model, verbose=True)
# odes_unstruct = [(i,":",odes)for i,odes in enumerate(model.odes)]

#
# for i,sp in enumerate(model.species):
#     print i,":",sp
# # print
# # for i,rxn in enumerate(model.reactions):
# #     print i,":",rxn
# # print
# myodes = []
# for i,ode in enumerate(model.odes):
#     myodes.append(i, odes)
# #     print myodes # i,":",ode

# Simulate the model
# time = np.linspace(0, 18000, 1801)
# #
# # # ODE simulation
# # #@TODO fix array problem
# # # param_values = np.array([p.value for p in model.parameters])
# # # param_values[18]=0.1
# plt.figure('Unstruct')
# x = odesolve(model, time, verbose=True) #integrator='lsoda',
# for obs in ["IKKKa_obs", "IKKii_obs", "IKKa_obs", "IKKi_obs", "IkBa_p_obs", "IkBan_NFkBn_obs"]:
#     plt.plot(time, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     # re.match(r"(\w+)_obs", obs).group()
#     #plt.plot(time/60, x[obs.name], label=obs.name)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)

# plt.figure()
# for i in range(len(model.species)):
#     plt.plot(time/60, x["__s%d" %i])
# plt.figure(2)
# x = odesolve(model, time, verbose=True)
# for obs in ["TNFR1a_obs", "A20_on_obs"]:
#     plt.plot(time, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# plt.figure(3)
# x = odesolve(model, time, verbose=True)
# for obs in ["A20_obs", "A20_off_obs", "IkBa_obs", "IkBan_obs"]:
#     plt.plot(time, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# x = odesolve(model, time, verbose=True)
# plt.figure(4)
# plt.plot(time/60, x["NFkBn_obs"], label=NFkBn_obs)
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend()
#
# plt.figure(5)
# plt.plot(time/60, x["IkBa_off_obs"], label=IkBa_off_obs)
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend()
#
# plt.figure(6)
# plt.plot(time/60, x["IkBa_on_obs"], label=IkBa_on_obs)
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend()

# plt.show()
# # SSA simulation
# x = run_ssa(model, time, verbose=True)
# plt.figure('SSA')
# plt.plot(time, x['B_obs'], label=B_obs, lw=2)
# plt.xlabel('Time (seconds)')
# plt.ylabel('Amount of B_obs')
# plt.legend(loc=0)

# plt.figure('ODE')
# plt.plot(time, x['IKKKa_obs'], label=IKKKa_obs, lw=2)
# plt.xlabel('Time (seconds)')
# plt.ylabel('Amount of IKKKa_obs')
# plt.legend(loc=0)
#plt.show()
