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


#TNFa binding to TNFR1
Rule('TNF_ext_and_TNFR1i', TNF_ext() + TNFR1i() >> TNF_ext() + TNFR1a(), kb)
Rule('test', TNF_ext() + TNFR1(state = 'i') >> TNF_ext() + TNFR1(state = 'a'), kb)
# Rule('TNF_ext_and_TNFR1i', TNF_ext(b = None) + TNFR1(b = None, state = 'i') >> TNF_ext(b = 1) % TNFR1(b = 1, state = 'i'), kb)
# Rule('activate_TNFR1', TNF_ext(b = 1) % TNFR1(b = 1, state = 'i') >> TNF_ext(b = 1) % TNFR1(b = None, state = 'a'), kb)
# catalyze_state(TNF_ext(), 'b', TNFR1(state = 'i'), 'b', TNFR1(state = 'a'), kb) #using macros

#TNFR1 activating IKKK
Rule('TNFR1a_and_IKKKn', TNFR1a() + IKKKn() >> TNFR1a() + IKKKa(), keff)
Rule('TNFR1a_creates_IKKKa', TNFR1(state = 'a') + IKKK(b = None, state = 'n') >> TNFR1(state = 'a') + IKKK(b = None, state = 'a'), keff)
# Rule('TNFR1a_and_IKKKn', TNFR1(b = None, state = 'a') + IKKK(b = None, state = 'n') >> TNFR1(b = 1, state = 'a') % IKKK(b = 1, state = 'n'), keff)
# Rule('TNFR1a_and_IKKKn', TNFR1(b = 1, state = 'a') % IKKK(b = 1, state = 'n') >> TNFR1(b = 1, state = 'a') + IKKK(b = 1, state = 'a'), keff)
# catalyze_state(TNFR1(), 'b', IKKK(state = 'n'), 'b', IKKK(state = 'a'), keff)

#Inactivation of TNFR1
Rule('TNFR1a_to_TNFR1i', TNFR1a() >> TNFR1i(), kf)
Rule('deactivate_TNFR1', TNFR1(state = 'a') >> TNFR1(state = 'i'), kf)

#IKKKa creates KN (total amount of IKKK kinase molecules)
Rule('IKKKa_to_IKKKn', IKKKa() >> IKKKn(), ki)
Rule('IKKKa_to_IKKKn', IKKK(b = None, state = 'a') >> IKKK(b = None, state = 'n'), ki)

#IKKK activates IKK
Rule('IKKKa_and_IKKn', IKKKa() + IKKKa() + IKKn() >> IKKKa() + IKKKa() + IKKa(), k1) #proportion of IKKKa changes IKKn to IKKa
Rule('IKKKa_activates_IKK', IKKK(state = 'a') + IKKK(state = 'a') + IKK(state = 'n') >> IKKK(state = 'a') + IKKK(state = 'a') + IKK(state = 'a'), k1)

#IKK inactivation
Rule('IKKa_to_IKKi', IKKa() >> IKKi(), k3) #IKKa active to IKKi inactive
Rule('IKKa_to_IKKi', IKK(state = 'a') >> IKK(state = 'i'), k3)

#A20 inactivating IKK
Rule('IKKa_and_A20', IKKa() + A20() >> A20() + IKKi(), k3_div_k2)
Rule('IKKa_and_A20', IKK(state = 'a') + A20(state = 'a') >> A20(state = 'a') + IKK(state = 'i'), k3_div_k2) #Exp1) #A20 mediated IKKa to IKKi

#inactive IKK inactivates immediate IKK
Rule('IKKi_to_IKKii', IKKi() >> IKKii(), k4) #
Rule('IKKi_to_IKKii', IKK(state = 'i') >> IKK(state = 'ii'), k4)

#degredation of TNF
Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)
Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)

#degredation of phos IkBa
Rule('IkBa_p_deg', IkBa_p() >> None, tp)
Rule('IkBa_p_deg', IkBa_p() >> None, tp)

#degrade A20 @TODO question about A20!
Rule('A20_deg', A20() >> None, c5)
Rule('A20_deg', A20() >> None, c5)

#degrade A20t
Rule('A20t_deg', A20t() >> None, c3)
Rule('A20_mRNA_deg', A20_mRNA() >> None, c3)

#IkBa cyto degredation
Rule('IkBa_cyt_deg', IkBac() >> None, c5a)
Rule('IkBa_cyt_deg', IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> None, c5a)

#degredation of IkBa_mRNA
Rule('IkBat_deg', IkBat() >> None, c3)
Rule('IkBa_mRNA_deg', IkBa_mRNA() >> None, c3)

#IkBa|NFkB phos complex cyto to free NFkB in cyto
Rule('IkBapc_NFkBpc_to_NFkB_c', IkBapc_NFkBpc() >> NFkB_c(), tp)
Rule('IkBapc_NFkBpc_to_NFkB_c', IkBa(b = ) >> NFkB_c(), tp)

#release of NFkB in cyto
Rule('IkBa_NFkB_to_NFkB_c', IkBa_NFkB() >> NFkB_c(), c6a)
Rule('IkBa_NFkB_to_NFkB_c', NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off') >> NFkB(b = None, phos = 'u', loc = 'c'), c6a)

#Transport of NFkB from cyto to nuc
Rule('NFkB_c_to_NFkBn', NFkB_c() >> NFkBn(), i1)
Rule('NFkB_c_to_NFkBn', NFkB(b = None, phos = 'u', loc = 'c') >> NFkB(b = None, phos = 'u', loc = 'c'), i1)

#NFkB in cyto and IkBa in cyto create bind in cyto
Rule('NFkB_c_and_IkBa', NFkB_c() + IkBa() >> IkBa_NFkB(), a1)
Rule('NFkB_c_and_IkBa', NFkB(b = None, phos = 'u', loc = 'c') + IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off') , a1)

#NFkBn and IkBan create IkBan|NFkBn complex in nucleus
Rule('NFkBn_and_IkBan', NFkBn() + IkBan() >> IkBan_NFkBn(), a1_mult_kv)
Rule('NFkBn_and_IkBan', NFkB(b = None, phos = 'u', loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n', state = 'off') >> NFkB(b = 1, phos = 'u', loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n', state = 'off'), a1_mult_kv)

#IKKa creates IkBa|NFkB phos complex in cyto
Rule('IKKa_and_IkBa_NFkB', IKKa() + IkBa_NFkB() >> IKKa() + IkBapc_NFkBpc(), a3)
Rule('IKKa_and_IkBa_NFkB', IKK(state = 'a') + NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off') >> IKK(state = 'a') + NFkB(b = 1, phos = 'p', loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c', state = 'off'), a3)

#A20 translation
Rule('A20_mRNA_to_A20', A20_mRNA() >> A20_mRNA() + A20(), c4)
Rule('A20_mRNA_to_A20', A20_mRNA() >> A20_mRNA() + A20(state = 'a'), c4)

#A20 transcription
Rule('A20_on_to_A20_mRNA', A20_on() >> A20_on() + A20_mRNA(), c1)
Rule('A20_on_to_A20_mRNA', A20(state = 'on') >> A20(state ='on') + A20_mRNA(), c1)

#IKKa phos IkBa
Rule('IKKa_and_IkBa', IKKa() + IkBa() >> IKKa() + IkBa_p(), a2)
Rule('IKKa_and_IkBa', IKK(state = 'a') + IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> IKK(state  = 'a') + IkBa(b = None, phos = 'p', loc = 'c', state = 'off') ), a2)

#IkBa translation
Rule('IkBa_mRNA_to_IkBa', IkBa_mRNA() >> IkBa_mRNA() + IkBa(), c4)
Rule('IkBa_mRNA_to_IkBa', IkBa_mRNA() >> IkBa_mRNA() + IkBa(), c4)

#IkBa cyto to nuc
Rule('IkBa_to_IkBan', IkBa() >> IkBan(), i1a)
Rule('IkBa_to_IkBan', IkBa(b = None, phos = 'u', loc = 'c', state = 'off') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'off'), i1a)

#IkBa nuc to cyto
Rule('IkBan_to_IkBa', IkBan() >> IkBa(), e1a)
Rule('IkBan_to_IkBa', IkBa(b = None, phos = 'u', loc = 'n', state = 'off') >> IkBa(b = None, phos = 'u', loc = 'c', state = 'off'), e1a)

#IkBa Transcription
Rule('IkBa_on_to_IkBat', IkBa_on() >> IkBa_on() + IkBat(), c1a)
Rule('IkBa_on_to_IkBat', IkBa(b = None, phos = 'u', loc = 'n', state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'on') + IkBa_mRNA(), c1a)

#IkBa\NFkB from nuc to cyto
Rule('IkBan_NFkBn_to_IkBa_NFkB', IkBan_NFkBn() >> IkBa_NFkB(), e2a)
Rule('IkBan_NFkBn_to_IkBa_NFkB', NFkB(b = 1, phos = 'u', loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n', state = 'off') >> NFkB(b = 1, phos = 'u', loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c', state = 'off'), e2a)

#turning off A20 gene state
Rule('IkBan_and_A20_on', IkBan() + A20_on() >> IkBan() + A20_off(), q2)
Rule('IkBan_and_A20_on', IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + A20(state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + A20(state = 'off'), q2)

#turning off IkBa gene state
Rule('IkBan_and_IkBa_on', IkBan() + IkBa_on() >> IkBan() + IkBa_off(), q2)
Rule('IkBan_and_IkBa_on', IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + IkBa(b = None, phos = 'u', loc = 'n', state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n', state = 'off') + IkBa(b = None, phos = 'u', loc = 'n', state = 'off'), q2)

#IKKii creates IKKn
Rule('IKKii_to_IKKn', IKKii() >> IKKn(), k4)
Rule('IKKii_to_IKKn', IKK(state = 'ii') >> IKK(state = 'n'), k4)

#NFkB turning off A20 gene state
Rule('NFkBn_and_A20_off', NFkBn() + A20_off() >> NFkBn() + A20_on(), q1)
Rule('NFkBn_and_A20_off', NFkBn(b = None, phos = 'u', loc = 'n') + A20(state = 'off') >> NFkB(b = None, phos = 'u', loc = 'n') + A20(state = 'on'), q1)

#NFkB turning off IkBa gene state
Rule('NFkBn_and_IkBa_off', NFkBn() + IkBa_off() >> NFkBn() + IkBa_on(), q1)
Rule('NFkBn_and_IkBa_off', NFkBn(b = None, phos = 'u', loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n', state = 'off') >> NFkBn(b = None, phos = 'u', loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n', state = 'on'), q1)


# for i,ode in enumerate(model.odes):
#     print i,":",ode