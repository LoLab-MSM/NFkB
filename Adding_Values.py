__author__ = 'geena'

from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from sympy import sympify

Model ()

#added monomers
Monomer('TNF_ext')
Monomer('TNFR1', ['state'], {'state': ['i', 'a']})
Monomer('IKKK', ['state'], {'state': ['i', 'a', 'n']})
Monomer('IKK', ['state'], {'state': ['n', 'a', 'i', 'ii']})
Monomer('IkBa', ['b', 'phos', 'loc', 'state'], {'phos': ['u', 'p'], 'loc': ['nuc', 'cyt'], 'state': ['on', 'off']})
Monomer('NFkB', ['b', 'phos', 'loc'], {'phos': ['u', 'p'], 'loc': ['nuc', 'cyt']})
Monomer('A20', ['state'], {'state': ['on', 'off', 'a']})
Monomer('IkBa_mRNA')
Monomer('A20_mRNA')

Parameter('KNN', 2.0e5)
Parameter('KN', 1.0e5)
Parameter('M', 2000)
Parameter('AN', 2)
Parameter('ANa', 2)

#Declaring initial conditions
# @TODO FIGURE OUT INITIAL CONDITION PROBLEM
Initial(IKKK(state = 'a'), Parameter('IKKKa_0'))
Initial(IKK(state = 'n'), Parameter('IKKn_0', 2e5))
Initial(IKK(state = 'a'), Parameter('IKKa_0'))
Initial(IKK(state = 'i'), Parameter('IKKi_0'))
Initial(NFkB(b = 1, phos = 'p', loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c', state = 'off'), Parameter('IkBapc_NFkBpc_0'))
Initial(NFkB(b = None, phos = 'u', loc = 'c'), Parameter('NFkB_c_0'))
Initial(NFkB(b = None, phos = 'u', loc = 'n'), Parameter('NFkBn_0', 1))
Initial(A20(state = 'a'), Parameter('A20_0', 10000))
Initial(A20_mRNA(), Parameter('A20t_0', 10))
Initial(IkBa(b = None, phos = 'u', loc = 'c', state = 'off'), Parameter('IkBa_0', 0.14*100000))
Initial(IkBa(b = None, phos = 'u', loc = 'n', state = 'off'), Parameter('IkBan_0', 0.06*100000))
Initial(IkBa_mRNA(), Parameter('IkBat_0', 10))
Initial(NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off'), Parameter('IkBa_NFkB_0', 100000))
Initial(NFkB(b = 1, phos = 'u', loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n', state = 'off'), Parameter('IkBan_NFkBn_0'))
Initial(TNFR1(state = 'a'), Parameter('TNFR1a_0'))
Initial(A20(state = 'on'), Parameter('A20_on_0'))
Initial(IkBa(b = None, phos = 'u', loc = 'n', state = 'on'), Parameter('IkBa_on_0'))
Initial(TNF_ext(), Parameter('TNF_ext_0', 0.1)) #1e-1
Initial(IKKK(state = 'n'), Parameter('IKKKn_0', KN.value))
Initial(IKK(state = 'ii'), Parameter('IKKii_0', KNN.value-IKKn_0.value))
Initial(TNFR1(state = 'i'), Parameter('TNFR1i_0', M.value))
Initial(A20(state = 'off'), Parameter('A20_off_0', AN.value))
Initial(IkBa(b = None, phos = 'u', loc = 'n', state = 'off'), Parameter('ANa_0', ANa.value))

#Cell Parameters
Parameter('kv', 5.0) #Nucealr to cytoplasm volume
#Parameter('Kn', 1e5) #total number of IKKK kinase molecules
#Parameter('K_nn', 2e-05) #total number of IKK kinase molecules

#Declaring Parameters
Parameter('kb', 1.2e-5) #receptor activation rate
Parameter('kf', 1.2e-3) #receptor inactivation rate
Parameter('Tdeg', 7.7e-4) #TNF loss same as cdeg
Parameter('ka',1e-5) #IKKK kinase activation rate
Parameter('ka20', 1e5) #A20 TNFR1 block
Parameter('ki', 0.01) #IKKK kinase inactivation rate
Parameter('AR', 0) #active receptors
Parameter('k1', 2*6e-10) #IKKn activation by IKKK
Parameter('k3', 0.002) #IKKa inactivation by A20
Parameter('k2', 10000) #IKKa inactivation
Parameter('k4', 0.001) #IKKii transfer rate
Parameter('AA', 1.0) #IkBa on (or off)
Parameter('c0', 0.1) # inducible A20 and IkBa mRNA synthesis
Parameter('c1a', 0.1) #inducible IkBa mRNA synthesis
Parameter('AB', 1.0) #A20 (on or off)

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
Parameter('G', 0) #initial status of A20 promoter
Parameter('c3', 0.00075) #A20 and IkBa mRNA deg rate
Parameter('q1', 4e-7) #NFkB attaching @ A20 and IkBa site
Parameter('q2', 1e-6) #IkBa inducible detaching from A20, IkBa site

Parameter('k3_div_k2', k3.value/k2.value) #for IKKa and A20
Parameter('a1_mult_kv', a1.value*kv.value) #for volume IkBa association NFkB


#Declaring expression
Observable('A20_obs', A20())
Expression('keff', sympify("ka*ka20/(ka20+A20_obs)")) #10000 #michaelis menten

#TNFa binding to TNFR1
# Rule('TNF_ext_and_TNFR1i', TNF_ext() + TNFR1i() >> TNF_ext() + TNFR1a(), kb)
Rule('test', TNF_ext() + TNFR1(state = 'i') >> TNF_ext() + TNFR1(state = 'a'), kb)
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
Rule('IKKa_and_A20', IKK(state = 'a') + A20(state = 'a') >> A20(state = 'a') + IKK(state = 'i'), k3_div_k2) #Exp1) #A20 mediated IKKa to IKKi

#inactive IKK inactivates immediate IKK
# Rule('IKKi_to_IKKii', IKKi() >> IKKii(), k4)
Rule('IKKi_to_IKKii', IKK(state = 'i') >> IKK(state = 'ii'), k4)

#degredation of TNF
# Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)
Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)

#degredation of phos IkBa
# Rule('IkBa_p_deg', IkBa_p() >> None, tp)
Rule('IkBa_p_deg', IkBa(b = None, phos = 'p', loc = 'c', state = 'off') >> None, tp)

#degrade A20 @TODO question about A20!
# Rule('A20_deg', A20() >> None, c5)
Rule('A20_deg', A20(state = 'a') >> None, c5)

#degrade A20t
# Rule('A20t_deg', A20t() >> None, c3)
Rule('A20_mRNA_deg', A20_mRNA() >> None, c3)

#IkBa cyto degredation
# Rule('IkBa_cyt_deg', IkBac() >> None, c5a)
Rule('IkBa_cyt_deg', IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> None, c5a)

#degredation of IkBa_mRNA
# Rule('IkBat_deg', IkBat() >> None, c3)
Rule('IkBa_mRNA_deg', IkBa_mRNA() >> None, c3)

#IkBa|NFkB phos complex cyto to free NFkB in cyto
# Rule('IkBapc_NFkBpc_to_NFkB_c', IkBapc_NFkBpc() >> NFkB_c(), tp)
Rule('IkBapc_NFkBpc_to_NFkB_c', NFkB(b = 1, phos = 'p', loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c', state = 'off') >> NFkB(b = None, phos = 'u', loc = 'c'), tp)

#release of NFkB in cyto
# Rule('IkBac_NFkBc_to_NFkB_c', IkBac_NFkBc() >> NFkB_c(), c6a)
Rule('IkBa_NFkB_to_NFkB_c', NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off') >> NFkB(b = None, phos = 'u', loc = 'c'), c6a)

#Transport of NFkB from cyto to nuc
# Rule('NFkB_c_to_NFkBn', NFkB_c() >> NFkBn(), i1)
Rule('NFkB_c_to_NFkBn', NFkB(b = None, phos = 'u', loc = 'c') >> NFkB(b = None, phos = 'u', loc = 'n'), i1)

#NFkB in cyto and IkBa in cyto create bind in cyto
# Rule('NFkB_c_and_IkBa', NFkB_c() + IkBa() >> IkBa_NFkB(), a1)
Rule('NFkB_c_and_IkBa', NFkB(b = None, phos = 'u', loc = 'c') + IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off') , a1)

#NFkBn and IkBan create IkBan|NFkBn complex in nucleus
# Rule('NFkBn_and_IkBan', NFkBn() + IkBan() >> IkBan_NFkBn(), a1_mult_kv)
Rule('NFkBn_and_IkBan', NFkB(b = None, phos = 'u', loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n', state = 'off') >> NFkB(b = 1, phos = 'u', loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n', state = 'off'), a1_mult_kv)

#IKKa creates IkBa|NFkB phos complex in cyto
# Rule('IKKa_and_IkBa_NFkB', IKKa() + IkBac_NFkBc() >> IKKa() + IkBapc_NFkBpc(), a3)
Rule('IKKa_and_IkBa_NFkB', IKK(state = 'a') + NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off') >> IKK(state = 'a') + NFkB(b = 1, phos = 'p', loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c', state = 'off'), a3)

#A20 translation
# Rule('A20_mRNA_to_A20', A20_mRNA() >> A20_mRNA() + A20(), c4)
Rule('A20_mRNA_to_A20', A20_mRNA() >> A20_mRNA() + A20(state = 'a'), c4)

#A20 transcription
# Rule('A20_on_to_A20_mRNA', A20_on() >> A20_on() + A20_mRNA(), c1)
Rule('A20_on_to_A20_mRNA', A20(state = 'on') >> A20(state ='on') + A20_mRNA(), c1)

#IKKa phos IkBa
# Rule('IKKa_and_IkBa', IKKa() + IkBa() >> IKKa() + IkBa_p(), a2)
Rule('IKKa_and_IkBap', IKK(state = 'a') + IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> IKK(state  = 'a') + IkBa(b = None, phos = 'p', loc = 'c', state = 'off'), a2)

#IkBa translation
# Rule('IkBa_mRNA_to_IkBan', IkBa_mRNA() >> IkBa_mRNA() + IkBan(), c4)
Rule('IkBa_mRNA_to_IkBan', IkBa_mRNA() >> IkBa_mRNA() + IkBa(b = None, phos = 'u', loc = 'n', state = 'off'), c4)

#IkBa cyto to nuc
# Rule('IkBac_to_IkBan', IkBac() >> IkBan(), i1a)
Rule('IkBac_to_IkBan', IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'off'), i1a)

#IkBa nuc to cyto
# Rule('IkBan_to_IkBac', IkBan() >> IkBac(), e1a)
Rule('IkBan_to_IkBac', IkBa(b = None, phos = 'u', loc = 'n', state = 'off') >> IkBa(b = None, phos = 'u', loc = 'c', state = 'off'), e1a)

#IkBa Transcription
# Rule('IkBa_on_to_IkBat', IkBa_on() >> IkBa_on() + IkBat(), c1a)
Rule('IkBa_on_to_IkBat', IkBa(b = None, phos = 'u', loc = 'n', state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'on') + IkBa_mRNA(), c1a)

#IkBa\NFkB from nuc to cyto
# Rule('IkBan_NFkBn_to_IkBac_NFkBc', IkBan_NFkBn() >> IkBac_NFkBc(), e2a)
Rule('IkBan_NFkBn_to_IkBac_NFkBc', NFkB(b = 1, phos = 'u', loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n', state = 'off') >> NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off'), e2a)

#turning off A20 gene state
# Rule('IkBan_and_A20_on', IkBan() + A20_on() >> IkBan() + A20_off(), q2)
Rule('IkBan_and_A20_on', IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + A20(state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + A20(state = 'off'), q2)

#turning off IkBa gene state
# Rule('IkBan_and_IkBa_on', IkBan() + IkBa_on() >> IkBan() + IkBa_off(), q2)
Rule('IkBan_and_IkBa_on', IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + IkBa(b = None, phos = 'u', loc = 'n', state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + IkBa(b = None, phos = 'u', loc = 'n', state = 'off'), q2)

#IKKii creates IKKn
# Rule('IKKii_to_IKKn', IKKii() >> IKKn(), k4)
Rule('IKKii_to_IKKn', IKK(state = 'ii') >> IKK(state = 'n'), k4)

#NFkB turning off A20 gene state
# Rule('NFkBn_and_A20_off', NFkBn() + A20_off() >> NFkBn() + A20_on(), q1)
Rule('NFkBn_and_A20_off', NFkB(b = None, phos = 'u', loc = 'n') + A20(state = 'off') >> NFkB(b = None, phos = 'u', loc = 'n') + A20(state = 'on'), q1)

#NFkB turning off IkBa gene state
# Rule('NFkBn_and_IkBa_off', NFkBn() + IkBa_off() >> NFkBn() + IkBa_on(), q1)
Rule('NFkBn_and_IkBa_off', NFkB(b = None, phos = 'u', loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n', state = 'off') >> NFkB(b = None, phos = 'u', loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n', state = 'on'), q1)


for i,ode in enumerate(model.odes):
    print i,":",ode