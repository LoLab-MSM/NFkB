__author__ = 'geena'

from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from sympy import sympify
from pysb.bng import *
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import matplotlib.pyplot as plt
import numpy as np

Model ()

#added monomers
Monomer('TNF_ext')
Monomer('TNFR1', ['state'], {'state': ['i', 'a']})
Monomer('IKKK', ['state'], {'state': ['a', 'n']})
Monomer('IKK', ['state'], {'state': ['n', 'a', 'i', 'ii']})
Monomer('NFkB', ['b','loc'], {'loc': ['n', 'c']})
Monomer('IkBa_gene', ['state'], {'state': ['on', 'off']})
Monomer('IkBa_mRNA')
Monomer('IkBa', ['b', 'phos', 'loc'], {'phos': ['u', 'p'], 'loc': ['n', 'c']})
Monomer('A20_gene', ['state'], {'state': ['on', 'off']})
Monomer('A20_mRNA')
Monomer('A20')

Parameter('KNN', 2.0e5)
Parameter('KN', 1.0e5)
Parameter('M', 2000)
Parameter('AN', 2)
Parameter('ANa', 2)

#Declaring initial conditions
Initial(IKKK(state = 'a'), Parameter('IKKKa_0'))
Initial(IKK(state = 'n'), Parameter('IKKn_0', 2e5))
Initial(IKK(state = 'a'), Parameter('IKKa_0'))
Initial(IKK(state = 'i'), Parameter('IKKi_0'))
Initial(IkBa(b = None, phos = 'p', loc = 'c'), Parameter('IkBap_0'))
Initial(NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c'), Parameter('IkBap_NFkB'))
Initial(NFkB(b = None, loc = 'c'), Parameter('NFkBc'))
Initial(NFkB(b = None, loc = 'n'), Parameter('NFkBn_0', 1))
Initial(A20(), Parameter('A20_0', 10000))
Initial(A20_mRNA(), Parameter('A20t_0', 10))
Initial(IkBa(b = None, phos = 'u', loc = 'c'), Parameter('IkBa_0', 0.14*100000))
Initial(IkBa(b = None, phos = 'u', loc = 'n'), Parameter('IkBan_0', 0.06*100000))
Initial(IkBa_mRNA(), Parameter('IkBat_0', 10))
Initial(NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c'), Parameter('IkBa_NFkB_0', 100000))
Initial(NFkB(b = 1, loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n'), Parameter('IkBan_NFkBn'))
Initial(TNFR1(state = 'a'), Parameter('TNFR1a'))
Initial(A20_gene(state = 'on'), Parameter('A20_on_0'))
Initial(IkBa_gene(state = 'on'), Parameter('IkBa_on_0'))
Initial(TNF_ext(), Parameter('TNF_ext_0', 0.1))
Initial(IKKK(state = 'n'), Parameter('IKKKn_0', KN.value))
Initial(IKK(state = 'ii'), Parameter('IKKii'))
Initial(TNFR1(state = 'i'), Parameter('TNFR1i_0', M.value))
Initial(A20_gene(state = 'off'), Parameter('A20_off_0', AN.value))
Initial(IkBa_gene(state = 'off'), Parameter('ANa_0', ANa.value))


#Cell Parameters
Parameter('kv', 5.0) #Nucealr to cytoplasm volume

#Declaring Parameters
Parameter('kb', 1.2e-5) #receptor activation rate
Parameter('kf', 1.2e-3) #receptor inactivation rate
Parameter('Tdeg', 7.7e-4) #TNF loss same as cdeg
Parameter('ka',1e-5) #IKKK kinase activation rate
Parameter('ka20', 1e5) #A20 TNFR1 block
Parameter('ki', 0.01) #IKKK kinase inactivation rate
Parameter('k1', 2.0*6e-10) #IKKn activation by IKKK
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

#Declaring observables
Observable('IKKKa_obs', IKKK(state = 'a'))
Observable('IKKKn_obs', IKK(state = 'n'))
Observable('IKKa_obs', IKK(state = 'a'))
Observable('IKKn_obs', IKK(state = 'n'))
Observable('IKKi_obs', IKK(state = 'i'))
Observable('IKKii_obs', IKK(state = 'ii'))
Observable('TNF_ext_obs', TNF_ext())
Observable('IkBap_obs', IkBa(b = None, phos = 'p', loc = 'c'))
Observable('IkBan_NFkBn_obs', NFkB(b = 1,loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n'))
Observable('IkBap_NFkBc_obs', NFkB(b = 1,loc = 'n') % IkBa(b = 1, phos = 'p', loc = 'c'))
Observable('NFkB_c_obs', NFkB(b = None, loc = 'c'))
Observable('NFkBn_obs', NFkB(b = None, loc = 'n'))

Observable('Nuclear_NFkBn', NFkB(b = 1,loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n') + NFkB(b = None, loc = 'n'))
#
Observable('A20_obs', A20())
# Observable('IKKKa_obs', IKKKa()) #
# Observable('IKKKn_obs', IKKKn())
# Observable('IKKa_obs', IKKa()) #
# Observable('IKKi_obs', IKKi()) #
# Observable('IKKii_obs', IKKii()) #
# Observable('TNF_ext_obs', TNF_ext()) #
# Observable('IkBa_p_obs', IkBa_p()) #
# Observable('IkBapc_NFkBpc_obs', IkBapc_NFkBpc()) #
# Observable('NFkB_c_obs', NFkB_c())
# Observable('NFkBn_obs', NFkBn())
# Observable('TNFR1a_obs', TNFR1a())
# Observable('TNFR1i_obs', TNFR1i())
# Observable('A20_off_obs', A20_off())
# Observable('IkBa_off_obs', IkBa_off())
# Observable('A20_on_obs', A20_on())
# Observable('IkBa_on_obs', IkBa_on())
# Observable('A20t_obs', A20t())
# Observable('IkBa_obs', IkBa())
# Observable('IkBan_obs', IkBan())
# Observable('IkBat_obs', IkBat())
# Observable('IkBa_NFkB_obs', IkBa_NFkB())
# Observable('IkBan_NFkBn_obs', IkBan_NFkBn())
# Observable('IKKn_obs', IKKn()) #

#Declaring expression
Expression('keff', sympify("ka*ka20/(ka20+A20_obs)")) #10000 #michaelis menten

#TNFa binding to TNFR1
# Rule('TNF_ext_and_TNFR1i', TNF_ext() + TNFR1i() >> TNF_ext() + TNFR1a(), kb)
Rule('TNF_activate_TNFR1', TNF_ext() + TNFR1(state = 'i') >> TNF_ext() + TNFR1(state = 'a'), kb)
# Rule('TNF_ext_and_TNFR1i', TNF_ext(b = None) + TNFR1(b = None, state = 'i') >> TNF_ext(b = 1) % TNFR1(b = 1, state = 'i'), kb)
# Rule('activate_TNFR1', TNF_ext(b = 1) % TNFR1(b = 1, state = 'i') >> TNF_ext(b = 1) % TNFR1(b = None, state = 'a'), kb)
# catalyze_state(TNF_ext(), 'b', TNFR1(state = 'i'), 'b', TNFR1(state = 'a'), kb) #using macros

#TNFR1 activating IKKK
# Rule('TNFR1a_and_IKKKn', TNFR1a() + IKKKn() >> TNFR1a() + IKKKa(), keff)
Rule('TNFR1a_creates_IKKKa', TNFR1(state = 'a') + IKKK(state = 'n') >> TNFR1(state = 'a') + IKKK(state = 'a'), keff)
# Rule('TNFR1a_and_IKKKn', TNFR1(b = None, state = 'a') + IKKK(b = None, state = 'n') >> TNFR1(b = 1, state = 'a') % IKKK(b = 1, state = 'n'), keff)
# Rule('TNFR1a_and_IKKKn', TNFR1(b = 1, state = 'a') % IKKK(b = 1, state = 'n') >> TNFR1(b = 1, state = 'a') + IKKK(b = 1, state = 'a'), keff)
# catalyze_state(TNFR1(), 'b', IKKK(state = 'n'), 'b', IKKK(state = 'a'), keff)

#Inactivation of TNFR1
# Rule('TNFR1a_to_TNFR1i', TNFR1a() >> TNFR1i(), kf)
Rule('inactivate_TNFR1', TNFR1(state = 'a') >> TNFR1(state = 'i'), kf)

#IKKKa creates KN (total amount of IKKK kinase molecules)
# Rule('IKKKa_to_IKKKn', IKKKa() >> IKKKn(), ki)
Rule('IKKKa_to_IKKKn', IKKK(state = 'a') >> IKKK(state = 'n'), ki)

#IKKK activates IKK
# Rule('IKKKa_and_IKKn', IKKKa() + IKKKa() + IKKn() >> IKKKa() + IKKKa() + IKKa(), k1) #proportion of IKKKa changes IKKn to IKKa
Rule('IKKKa_activates_IKK', IKKK(state = 'a') + IKKK(state = 'a') + IKK(state = 'n') >> IKKK(state = 'a') + IKKK(state = 'a') + IKK(state = 'a'), k1)

#IKK inactivation
# Rule('IKKa_to_IKKi', IKKa() >> IKKi(), k3) #IKKa active to IKKi inactive
Rule('IKKa_to_IKKi', IKK(state = 'a') >> IKK(state = 'i'), k3)

#A20 inactivating IKK
# Rule('IKKa_and_A20', IKKa() + A20() >> A20() + IKKi(), k3_div_k2)
Rule('IKKa_and_A20', IKK(state = 'a') + A20() >> IKK(state = 'i') + A20(), k3_div_k2) #Exp1) #A20 mediated IKKa to IKKi

#inactive IKK inactivates immediate IKK
# Rule('IKKi_to_IKKii', IKKi() >> IKKii(), k4)
Rule('IKKi_to_IKKii', IKK(state = 'i') >> IKK(state = 'ii'), k4)

#degredation of TNF
# Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)
Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)

#IkBa|NFkB phos complex cyto to free NFkB in cyto
# Rule('IkBapc_NFkBpc_to_NFkB_c', IkBapc_NFkBpc() >> NFkB_c(), tp)
Rule('IkBapc_NFkBpc_to_NFkB_c', IkBa(b = 1, phos = 'p', loc = 'c') % NFkB(b = 1,  loc = 'c') >> NFkB(b = None, loc = 'c'), tp)

#degredation of phos IkBa
# Rule('IkBa_p_deg', IkBa_p() >> None, tp)
Rule('IkBa_p_deg', IkBa(b = None, phos = 'p', loc = 'c') >> None, tp)

#degrade A20
# Rule('A20_deg', A20() >> None, c5)
Rule('A20_deg', A20() >> None, c5)

#degrade A20t
# Rule('A20t_deg', A20t() >> None, c3)
Rule('A20_mRNA_deg', A20_mRNA() >> None, c3)


#IkBa cyto degredation
# Rule('IkBa_cyt_deg', IkBac() >> None, c5a)
Rule('IkBa_cyt_deg', IkBa(b = None, phos = 'u', loc = 'c') >> None, c5a)

#release of NFkB in cyto
# Rule('IkBac_NFkBc_to_NFkB_c', IkBac_NFkBc() >> NFkB_c(), c6a)
Rule('IkBa_NFkB_to_NFkB_c', IkBa(b = 1, phos = 'u', loc = 'c') % NFkB(b = 1, loc = 'c') >> NFkB(b = None, loc = 'c'), c6a)

#degredation of IkBa_mRNA
# Rule('IkBat_deg', IkBat() >> None, c3)
Rule('IkBa_mRNA_deg', IkBa_mRNA() >> None, c3)

#Transport of NFkB from cyto to nuc
# Rule('NFkB_c_to_NFkBn', NFkB_c() >> NFkBn(), i1)
Rule('NFkB_c_to_NFkBn', NFkB(b = None, loc = 'c') >> NFkB(b = None, loc = 'n'), i1)

#NFkB in cyto and IkBa in cyto create bind in cyto
# Rule('NFkB_c_and_IkBa', NFkB_c() + IkBa() >> IkBa_NFkB(), a1)
Rule('NFkB_c_and_IkBa', NFkB(b = None, loc = 'c') + IkBa(b = None, phos = 'u', loc = 'c') >> NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c') , a1)

#NFkBn and IkBan create IkBan|NFkBn complex in nucleus
# Rule('NFkBn_and_IkBan', NFkBn() + IkBan() >> IkBan_NFkBn(), a1_mult_kv)
Rule('NFkBn_and_IkBan', NFkB(b = None, loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n') >> NFkB(b = 1, loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n'), a1_mult_kv)

#IKKa creates IkBa|NFkB phos complex in cyto
# Rule('IKKa_and_IkBa_NFkB', IKKa() + IkBac_NFkBc() >> IKKa() + IkBapc_NFkBpc(), a3)
Rule('IKKa_and_IkBa_NFkB', IKK(state = 'a') + NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c') >> IKK(state = 'a') + NFkB(b = 1,  loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c'), a3)

#IKKa phos IkBa
# Rule('IKKa_and_IkBa', IKKa() + IkBa() >> IKKa() + IkBa_p(), a2)
Rule('IKKa_and_IkBap', IKK(state = 'a') + IkBa(b = None, phos = 'u', loc = 'c') >> IKK(state  = 'a') + IkBa(b = None, phos = 'p', loc = 'c'), a2)

#A20 translation
# Rule('A20_mRNA_to_A20', A20_mRNA() >> A20_mRNA() + A20(), c4)
Rule('A20_mRNA_to_A20', A20_mRNA() >> A20_mRNA() + A20(), c4)

#A20 transcription
# Rule('A20_on_to_A20_mRNA', A20_on() >> A20_on() + A20_mRNA(), c1)
Rule('A20_gene_to_A20_mRNA', A20_gene(state = 'on') >> A20_gene(state ='on') + A20_mRNA(), c1)


#IkBa translation
# Rule('IkBa_mRNA_to_IkBan', IkBa_mRNA() >> IkBa_mRNA() + IkBan(), c4)
Rule('IkBa_mRNA_to_IkBa', IkBa_mRNA() >> IkBa_mRNA() + IkBa(b = None, phos = 'u', loc = 'c'), c4)

#IkBa cyto to nuc
# Rule('IkBac_to_IkBan', IkBac() >> IkBan(), i1a)
Rule('IkBac_to_IkBan', IkBa(b = None, phos = 'u', loc = 'c') >> IkBa(b = None, phos = 'u', loc = 'n'), i1a)

#IkBa nuc to cyto
# Rule('IkBan_to_IkBac', IkBan() >> IkBac(), e1a)
Rule('IkBan_to_IkBac', IkBa(b = None, phos = 'u', loc = 'n') >> IkBa(b = None, phos = 'u', loc = 'c'), e1a)

#IkBa Transcription
# Rule('IkBa_on_to_IkBat', IkBa_on() >> IkBa_on() + IkBat(), c1a)
Rule('IkBa_on_to_IkBat', IkBa_gene(state = 'on') >> IkBa_gene(state = 'on') + IkBa_mRNA(), c1a)

#IkBa\NFkB from nuc to cyto
# Rule('IkBan_NFkBn_to_IkBac_NFkBc', IkBan_NFkBn() >> IkBac_NFkBc(), e2a)
Rule('IkBan_NFkBn_to_IkBac_NFkBc', NFkB(b = 1, loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n') >> NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c'), e2a)

#turning off A20 gene state
# Rule('IkBan_and_A20_on', IkBan() + A20_on() >> IkBan() + A20_off(), q2)
Rule('IkBan_and_A20_on', IkBa(b = None, phos = 'u', loc = 'n') + A20_gene(state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n') + A20_gene(state = 'off'), q2)

#turning off IkBa gene state
# Rule('IkBan_and_IkBa_on', IkBan() + IkBa_on() >> IkBan() + IkBa_off(), q2)
Rule('IkBan_and_IkBa_on', IkBa(b = None, phos = 'u', loc = 'n') + IkBa_gene(state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n') + IkBa_gene(state = 'off'), q2)

#IKKii creates IKKn
# Rule('IKKii_to_IKKn', IKKii() >> IKKn(), k4)
Rule('IKKii_to_IKKn', IKK(state = 'ii') >> IKK(state = 'n'), k4)

#NFkB turning off A20 gene state
# Rule('NFkBn_and_A20_off', NFkBn() + A20_off() >> NFkBn() + A20_on(), q1)
Rule('NFkBn_and_A20_off', NFkB(b = None, loc = 'n') + A20_gene(state = 'off') >> NFkB(b = None, loc = 'n') + A20_gene(state = 'on'), q1)

#NFkB turning off IkBa gene state
# Rule('NFkBn_and_IkBa_off', NFkBn() + IkBa_off() >> NFkBn() + IkBa_on(), q1)
Rule('NFkBn_and_IkBa_off', NFkB(b = None, loc = 'n') + IkBa_gene(state = 'off') >> NFkB(b = None, loc = 'n') + IkBa_gene(state = 'on'), q1)


generate_equations(model, verbose = True)

time = np.linspace(0, 18000, 1801)
x = odesolve(model, time, verbose=True) #integrator='lsoda',


for obs in ["IKKKa_obs", "IKKii_obs", "IKKa_obs", "IKKi_obs","TNF_ext_obs", "IkBap_obs", "IkBap_NFkBc"]:
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.plot(time/60, x["TNF_ext_obs"], label = 'TNF')
    plt.subplot(2,1,2)
    plt.plot(time/60, x["IKKKa_obs"], label = 'IKKKa')
    # plt.figure(2)
    # plt.subplot(2,1,1)
    # plt.plot(time/60, x["IKKii_obs"], label = 'IKKii')
    # plt.subplot(2,1,2)
    # plt.plot(time/60, x["IKKn_obs"], label = 'IKKn')
    # plt.figure(3)
    # plt.subplot(2,1,1)
    # plt.plot(time/60, x["IKKa_obs"], label = 'IKKa')
    # plt.subplot(2,1,2)
    # plt.plot(time/60, x["IKKi_obs"], label = 'IKKi')
    # plt.figure(4)
    # plt.subplot(2,1,1)
    # plt.plot(time/60, x["IkBap_obs"], label = 'IkBap')
    # plt.subplot(2,1,2)
    # plt.plot(time/60, x["IkBap_NFkBc_obs"], label = 'IkBap_NFkBc')
    plt.legend(loc=0, prop={'size': 16})
    plt.xlabel("Time (in minutes)", fontsize=16)
    plt.ylabel("Concentrations", fontsize=16)

plt.show()

# plt.figure(1)
# for obs in ["IKKKa_obs", "IKKii_obs", "IKKa_obs", "IKKi_obs", "IkBa_p_obs", "IkBan_NFkBn_obs"]:
#     plt.subplot(3,2,1)
#     plt.plot(time, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#
#
#     # re.match(r"(\w+)_obs", obs).group()
#     #plt.plot(time/60, x[obs.name], label=obs.name)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)

# # plt.figure()
# # for i in range(len(model.species)):
# #     plt.plot(time/60, x["__s%d" %i])
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
# plt.plot(time/60, x["NFkBn_obs"], label='NFkBn')
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

# generate_equations(model, verbose = True)
# for i,ode in enumerate(model.odes):
#     print i,":",ode
