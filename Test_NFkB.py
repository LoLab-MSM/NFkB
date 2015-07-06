from pysb import *
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import matplotlib.pyplot as plt
import numpy as np
from sympy import sympify
from pysb.bng import generate_equations

Model ()

#Deterministic approximation
#Number of active TNFR1 receptors B:

#Declaring the monomers
Monomer('IKKKa')
Monomer('IKKn')
Monomer('IKKa')
Monomer('IKKi')
Monomer('IkBa_p')
Monomer('IkBapc_NFkBpc')
Monomer('NFkB_c')
Monomer('NFkBn')
Monomer('A20')
Monomer('A20t')
Monomer('IkBa')
Monomer('IkBan')
Monomer('IkBat')
Monomer('IkBa_NFkB')
Monomer('IkBan_NFkBn')
Monomer('B')
Monomer('A20_gs')
Monomer('IkBa_gs')
Monomer('TNF_ext')
Monomer('KN')
Monomer('KNN')
Monomer('M')
Monomer('AN')
Monomer('ANa')

#Declaring initial conditions
Initial(IKKKa(), Parameter('IKKKa_0'))
Initial(IKKn(), Parameter('IKKn_0', 2e5))
Initial(IKKa(), Parameter('IKKa_0'))
Initial(IKKi(), Parameter('IKKi_0'))
Initial(IkBa_p(), Parameter('IkBa_p_0'))
Initial(IkBapc_NFkBpc(), Parameter('IkBapc_NFkBpc_0'))
Initial(NFkB_c(), Parameter('NFkB_c_0'))
Initial(NFkBn(), Parameter('NFkBn_0', 1))
Initial(A20(), Parameter('A20_0', 10000))
Initial(A20t(), Parameter('A20t_0', 10))
Initial(IkBa(), Parameter('IkBa_0', 0.14*2e-5))
Initial(IkBan(), Parameter('IkBan_0', 0.06*2e-5))
Initial(IkBat(), Parameter('IkBat_0', 10))
Initial(IkBa_NFkB(), Parameter('IkBa_NFkB_0', 2e5))
Initial(IkBan_NFkBn(), Parameter('IkBan_NFkBn_0'))
Initial(B(), Parameter('B_0', 2000))
Initial(A20_gs(), Parameter('A20_gs_0'))
Initial(IkBa_gs(), Parameter('IkBa_gs_0'))
Initial(TNF_ext(), Parameter('TNF_ext_0', 1e-1))
Initial(KN(), Parameter('KN_0', 1e5))
Initial(KNN(), Parameter('KNN_0', 1e-9))
Initial(M(), Parameter('M_0', 2000))
Initial(AN(), Parameter('AN_0', 2))
Initial(ANa(), Parameter('ANa_0', 2))

#Parameters for TNFa and TNFR1 activation
# Parameter('kon', 1.83e-7) #value in M^-1*s^-1? @TODO check this value
# Parameter('koff', 3.48e-4)
# Parameter('kint', 7.70e-4)
# Parameter('cdeg', 2e-4) #extracellular TNFa degredation

#Reporter gene parameters
# Parameter('q1r', 10e-7) #NFkB binding at reporter gene promoter
# Parameter('q2r', 10e-7) #IkBa inducible NFkB detaching from reporter gene
# Parameter('q2rr', 10e-3) #spontaneous NFkB detaching from reporter gene
# Parameter('c1r', 5e-2) #inducible reportner mRNA synthesis
# Parameter('c1rr', 10e-3) #reporter mRNA constitutive synthesis
Parameter('c1a', 1*0.1) #inducible IkBa mRNA synthesis
#Parameter('c3r', varies) parameter varies in paper

#Cell Parameters
Parameter('kv', 5) #Nucealr to cytoplasm volume
Parameter('Kn', 1e5) #total number of IKKK kinase molecules
#Parameter('K_nn', 2e-05) #total number of IKK kinase molecules

#Declaring Parameters
Parameter('kb', .000012) #receptor activation rate
Parameter('kf', 1.2e-3) #receptor inactivation rate
Parameter('Tdeg', 2e-4) #TNF loss same as cdeg
Parameter('ka',2e-05) #IKKK kinase activation rate
Parameter('ka20', 1e5) #A20 TNFR1 block
Parameter('ki', 0.01) #IKKK kinase inactivation rate
Parameter('AR', 0) #active receptors
Parameter('k1', 2*6e-10) #IKKn activation by IKKK
Parameter('k3', 0.002) #IKKa inactivation by A20
Parameter('k2', 10000) #IKKa inactivation
Parameter('k4', 1e-3) #IKKii transfer rate

#IkBa Parameters
Parameter('a2', 1e-7) #IkBa phosphorylation b/c IKKa
Parameter('tp', 0.01) #degradation of phosph-IkBa complex with NFkB
Parameter('a3', 5e-7) #IkBa_NFkB phosphorylation b/c IKKa
Parameter('a1', 5e-7) #IkBa*NFkB association
Parameter('c6a', 0.00002) #spontaneous IkBa_NFkB defg of IkBa complex to NFkB
Parameter('c5a', 0.0001) #IkBa deg rate

#Transport Parameters
Parameter('i1a', 2e-3) #IkBa nuclear import
Parameter('e1a', 5e-3) #IkBa nuclear export
Parameter('e2a', 5e-2) #IkBa_NFkB nuclear export
Parameter('i1', 0.01) #NFkB nuclear import

#A20 and IkBa synthesis Parameters
Parameter('c4', 0.5) #A20, IkBa transformation rate
Parameter('c5', 0.0005) #A20 degredation rate
Parameter('c1', 0.01) #inducible A20 mRNA synthesis
Parameter('G', 0) #initial status of A20 promoter
Parameter('c3', 0.00075) #A20 and IkBa mRNA deg rate
Parameter('q1', 4e-7) #NFkB attaching @ A20 and IkBa site
Parameter('q2', 1e-6) #IkBa inducible detaching from A20, IkBa site

Parameter('k3_div_k2', k3.value/k2.value) #for IKKa and A20
Parameter('a1_mult_kv', a1.value*kv.value) #for volume IkBa association NFkB

#Declaring expression
Observable('A20_obs', A20())
Expression('keff', sympify("ka*ka20/(ka20+A20_obs)")) #10000 #michaelis menten

#Rule for TNFa and TNFR1 activation
# Rule('IKKn5_IKKn6_create_IKKn7', IKKn5() + IKKn6() >> IKKn7(), kon) #creating TNFR1|TNFa complexes
# Rule('degrade_IKKn5', IKKn5() >> None, cdeg) #degredation of free TNFa
# Rule('IKKn7_create_IKKn5_IKKn6', IKKn7() >> IKKn5() + IKKn6(), koff) #(TNFR1|TNFa) creating free TNFa and TNFR1
# Rule('IKKn7_create_IKKn6', IKKn7() >> IKKn6(), cdeg) #TNFR1|TNFa complex creates free TNFR1 receptors
# Rule('IKKn7_create_IKKn8', IKKn7() >> IKKn8(), kint) #TNFR1|TNFa complex internalized

#Declaring rules
Rule('IKKKa_to_KN', IKKKa() >> KN(), ki)
Rule('IKKKa_and_IKKn', IKKKa() + IKKKa() + IKKn() >> IKKKa() + IKKKa() + IKKa(), k1)
Rule('IKKa_to_IKKi', IKKa() >> IKKi(), k3)
Rule('IKKa_and_A20', IKKa() + A20() >> A20() + IKKi(), k3_div_k2) #Exp1)
Rule('IKKi_to_KNN', IKKi() >> KNN(), k4)
Rule('IkBa_p_deg', IkBa_p() >> None, tp)
Rule('IkBapc_NFkBpc_to_NFkB_c', IkBapc_NFkBpc() >> NFkB_c(), tp)
Rule('IkBa_NFkB_to_NFkB_c', IkBa_NFkB() >> NFkB_c(), c6a)
Rule('NFkB_c_to_NFkBn', NFkB_c() >> NFkBn(), i1)
Rule('NFkB_c_and_IkBa', NFkB_c() + IkBa() >> IkBa_NFkB(), a1)
Rule('NFkBn_and_IkBan', NFkBn() + IkBan() >> IkBan_NFkBn(), a1_mult_kv)
Rule('A20t_to_A20', A20t() >> A20t() + A20(), c4)
Rule('A20_deg', A20() >> None, c5)
Rule('B_and_KN', B() + KN() >> B() + IKKKa(), keff)
Rule('A20_gs_to_A20t', A20_gs() >> A20_gs() + A20t(), c1)
Rule('A20t_deg', A20t() >> None, c3)
Rule('IKKa_and_IkBa', IKKa() + IkBa() >> IKKa() + IkBa_p(), a2)
Rule('IkBat_to_IkBa', IkBat() >> IkBat() + IkBa(), c4)
Rule('IkBa_deg', IkBa() >> None, c5a)
Rule('IkBa_to_IkBan', IkBa() >> IkBan(), i1a)
Rule('IkBan_to_IkBa', IkBan() >> IkBa(), e1a)
Rule('IkBa_gs_to_IkBat', IkBa_gs() >> IkBa_gs() + IkBat(), c1a)
Rule('IkBat_deg', IkBat() >> None, c3)
Rule('IKKa_and_IkBa_NFkB', IKKa() + IkBa_NFkB() >> IKKa() + IkBapc_NFkBpc(), a3)
Rule('IkBan_NFkBn_to_IkBa_NFkB', IkBan_NFkBn() >> IkBa_NFkB(), e2a)
Rule('B_to_M', B() >> M(), kf)
Rule('IkBan_and_A20_gs', IkBan() + A20_gs() >> IkBan() + AN(), q2)
Rule('IkBan_and_IkBa_gs', IkBan() + IkBa_gs() >> IkBan() + ANa(), q2)
Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)
Rule('KNN_to_IKKn', KNN() >> IKKn(), k4)
Rule('TNF_ext_and_M', TNF_ext() + M() >> TNF_ext() + B(), kb)
Rule('NFkBn_and_AN', NFkBn() + AN() >> NFkBn() + A20_gs(), q1)
Rule('NFkBn_and_ANa', NFkBn() + ANa() >> NFkBn() + IkBa_gs(), q1)

#Declaring observables
Observable('IKKKa_obs', IKKKa())
Observable('IKKn_obs', IKKn())
Observable('IKKa_obs', IKKa())
Observable('IKKi_obs', IKKi())
Observable('IkBa_p_obs', IkBa_p())
Observable('IkBapc_NFkBpc_obs', IkBapc_NFkBpc())
Observable('NFkB_c_obs', NFkB_c())
Observable('NFkBn_obs', NFkBn())
Observable('A20t_obs', A20t())
Observable('IkBa_obs', IkBa())
Observable('IkBan_obs', IkBan())
Observable('IkBat_obs', IkBat())
Observable('IkBa_NFkB_obs', IkBa_NFkB())
Observable('IkBan_NFkBn_obs', IkBan_NFkBn())
Observable('B_obs', B())
Observable('A20_gs_obs', A20_gs())
Observable('IkBa_gs_obs', IkBa_gs())
Observable('TNF_ext_obs', TNF_ext())
Observable('KN_obs', KN())
Observable('KNN_obs', KNN())
Observable('M_obs', M())
Observable('AN_obs', AN())
Observable('ANa_obs', ANa())



generate_equations(model, verbose=True)

# for i,sp in enumerate(model.species):
#     print i,":",sp
# print
# for i,rxn in enumerate(model.reactions):
#     print i,":",rxn
# print
# for i,ode in enumerate(model.odes):
#     print i,":",ode


# Simulate the model 
time = np.linspace(0, 18000, 15)

# ODE simulation
x = odesolve(model, time, verbose=True) #integrator='lsoda',
#obs = []
# for i in model.observables:
#     plt.figure(1)
#     plt.plot(time,x[i.name], label=i.name)
#     plt.xlabel('time')
#     plt.ylabel('# of molecules')
# plt.legend()
# plt.show()


plt.figure('ODE')
plt.plot(time, x['B_obs'], label=B_obs, lw=2)
plt.xlabel('Time (seconds)')
plt.ylabel('Amount of B_obs')
plt.legend(loc=0)

# SSA simulation
x = run_ssa(model, time, verbose=True)
plt.figure('SSA')
plt.plot(time, x['B_obs'], label=B_obs, lw=2)
plt.xlabel('Time (seconds)')
plt.ylabel('Amount of B_obs')
plt.legend(loc=0)

plt.show()


