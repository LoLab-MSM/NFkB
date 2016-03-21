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

#Declaration of monomers
Monomer('IkBa', ['b', 'c', 'S'], {'S': ['C', 'N']})
Monomer('IkBb', ['b', 'c', 'S'], {'S': ['C', 'N']})
Monomer('IkBe', ['b', 'c', 'S'], {'S': ['C', 'N']})
Monomer('IkBd', ['b', 'c', 'S'], {'S': ['C', 'N']})
Monomer('IkBa_mRNA')
Monomer('IkBb_mRNA')
Monomer('IkBe_mRNA')
Monomer('IkBd_mRNA')
Monomer('NFkB', ['b', 'S'], {'S': ['C', 'N']})
Monomer('IKK1', ['b', 'c', 'S'], {'S': ['C', 'N']})
Monomer('IKK2', ['b', 'c', 'S'], {'S': ['C', 'N']})

#Declaration of initial conditions

Parameter('IKK1_0', 0.1) #Inhibitor Kinase Kinase a
Parameter('IKK2_0', 0.1) #Inhibitor Kinase Kinase b
Parameter('NFkB_0', 0.125) #Nuclear Factor-kappaB

Initial(IKK1(b=None, c=None, S='C'), IKK1_0)
Initial(IKK2(b=None, c=None, S='C'), IKK2_0)
Initial(NFkB(b=None, S='N'), NFkB_0)

# IkB mRNA and protein synthesis reactions

# Parameter('mRNA_a', 0.035)
# Parameter('mRNA_b', 3e-3)
# Parameter('mRNA_e', 4e-3)
# Parameter('mRNA_d', 2e-3)
# Parameter('synth', 0.25)
# Parameter('synth_a', 2e-4)
# Parameter('synth_b', 1e-5)
# Parameter('synth_e', 3e-6)
# Parameter('synth_d', 1e-7)
# Parameter('psynth_a', 2e-4)
# Parameter('psynth_b', 1e-5)
# Parameter('psynth_e', 3e-6)
# Parameter('psynth_d', 1e-7)
#
# Rule('a_mRNA', IkBa_mRNA() >> None, mRNA_a)
# Rule('b_mRNA', IkBb_mRNA() >> None, mRNA_b)
# Rule('e_mRNA', IkBe_mRNA() >> None, mRNA_e)
# Rule('d_mRNA', IkBd_mRNA() >> None, mRNA_d)
#
# Rule('a_psynth', None >> IkBa(b=None, c=None, S='C'), synth)
# Rule('b_psynth', None >> IkBb(b=None, c=None, S='C'), synth)
# Rule('e_psynth', None >> IkBe(b=None, c=None, S='C'), synth)
# Rule('d_psynth', None >> IkBd(b=None, c=None, S='C'), synth)
#
# Rule('a_synth', None >> IkBa_mRNA(), psynth_a)
# Rule('b_synth', None >> IkBb_mRNA(), psynth_b)
# Rule('e_synth', None >> IkBe_mRNA(), psynth_e)
# Rule('d_synth', None >> IkBd_mRNA(), psynth_d)

#IkB(a,b,e) association and dissociation from IKK2 and IkBd association and dissociation from IKK2

Parameter('a_IKKf', 1.35)
Parameter('a_IKKr', 0.75)
Parameter('b_IKKf', 0.36)
Parameter('IkB_IKK_r', 0.105)
Parameter('ed_IKKf', 0.54)
Parameter('IkB_IKKf', 30)
Parameter('IkB_IKKr', 6e-5)
Parameter('IkB_IKK_NFkBaf', 11.1)
Parameter('IkB_IKK_NFkBbf', 2.88)
Parameter('IkB_IKK_NFkBedf', 4.2)

Rule('a2_adc', IkBa(b=None, c=None, S='C') + IKK2(b=None, c=None, S='C') <> IkBa(b=1, c=None, S='C')%IKK2(b=1, c=None, S='C'), a_IKKf, a_IKKr)
Rule('b2_adc', IkBb(b=None, c=None, S='C') + IKK2(b=None, c=None, S='C') <> IkBb(b=1, c=None, S='C')%IKK2(b=1, c=None, S='C'), b_IKKf, IkB_IKK_r)
Rule('e2_adc', IkBe(b=None, c=None, S='C') + IKK2(b=None, c=None, S='C') <> IkBe(b=1, c=None, S='C')%IKK2(b=1, c=None, S='C'), ed_IKKf, IkB_IKK_r)
Rule('d1_adc', IkBd(b=None, c=None, S='C') + IKK1(b=None, c=None, S='C') <> IkBd(b=1, c=None, S='C')%IKK1(b=1, c=None, S='C'), ed_IKKf, IkB_IKK_r)
#
# Rule('an_adc', IkBa(b=None, c=None, S='C') + NFkB(b=None, S='C') <> IkBa(b=2, c=None, S='C')%NFkB(b=2, S='C'), IkB_IKKf, IkB_IKKr)
# Rule('bn_adc', IkBb(b=None, c=None, S='C') + NFkB(b=None, S='C') <> IkBb(b=2, c=None, S='C')%NFkB(b=2, S='C'), IkB_IKKf, IkB_IKKr)
# Rule('en_adc', IkBe(b=None, c=None, S='C') + NFkB(b=None, S='C') <> IkBe(b=2, c=None, S='C')%NFkB(b=2, S='C'), IkB_IKKf, IkB_IKKr)
# Rule('dn_adc', IkBd(b=None, c=None, S='C') + NFkB(b=None, S='C') <> IkBd(b=2, c=None, S='C')%NFkB(b=2, S='C'), IkB_IKKf, IkB_IKKr)
#
# Rule('an_adn', IkBa(b=None, c=None, S='N') + NFkB(b=None, S='N') <> IkBa(b=3, c=None, S='N')%NFkB(b=3, S='N'), IkB_IKKf, IkB_IKKr)
# Rule('bn_adn', IkBb(b=None, c=None, S='N') + NFkB(b=None, S='N') <> IkBb(b=3, c=None, S='N')%NFkB(b=3, S='N'), IkB_IKKf, IkB_IKKr)
# Rule('en_adn', IkBe(b=None, c=None, S='N') + NFkB(b=None, S='N') <> IkBe(b=3, c=None, S='N')%NFkB(b=3, S='N'), IkB_IKKf, IkB_IKKr)
# Rule('dn_adn', IkBd(b=None, c=None, S='N') + NFkB(b=None, S='N') <> IkBd(b=3, c=None, S='N')%NFkB(b=3, S='N'), IkB_IKKf, IkB_IKKr)

# Rule('a2n_c1', IkBa(b=None, c=1, S='C')%IKK2(b=None, c=1, S='C') + NFkB(b=None, S='C') <> IkBa(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=1, c=1, S='C'), IkB_IKKf, IkB_IKKr)
# Rule('b2n_c1', IkBb(b=None, c=1, S='C')%IKK2(b=None, c=1, S='C') + NFkB(b=None, S='C') <> IkBb(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=1, c=1, S='C'), IkB_IKKf, IkB_IKKr)
# Rule('e2n_c1', IkBe(b=None, c=1, S='C')%IKK2(b=None, c=1, S='C') + NFkB(b=None, S='C') <> IkBe(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=1, c=1, S='C'), IkB_IKKf, IkB_IKKr)
# Rule('d1n_c1', IkBd(b=None, c=1, S='C')%IKK1(b=None, c=1, S='C') + NFkB(b=None, S='C') <> IkBd(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK1(b=1, c=1, S='C'), IkB_IKKf, IkB_IKKr)
#
# Rule('a2n_c', IkBa(b=1, c=None, S='C')%NFkB(b=1, S='C') + IKK2(b=None, c=None, S='C') <> IkBa(b=1, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=None, c=1, S='C'), IkB_IKK_NFkBaf, a_IKKr)
# Rule('b2n_c', IkBb(b=1, c=None, S='C')%NFkB(b=1, S='C') + IKK2(b=None, c=None, S='C') <> IkBb(b=1, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=None, c=1, S='C'), IkB_IKK_NFkBbf, IkB_IKKr)
# Rule('e2n_c', IkBe(b=1, c=None, S='C')%NFkB(b=1, S='C') + IKK2(b=None, c=None, S='C') <> IkBe(b=1, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=None, c=1, S='C'), IkB_IKK_NFkBedf, IkB_IKKr)
# Rule('d1n_c', IkBd(b=1, c=None, S='C')%NFkB(b=1, S='C') + IKK1(b=None, c=None, S='C') <> IkBd(b=1, c=1, S='C')%NFkB(b=1, S='C')%IKK1(b=None, c=1, S='C'), IkB_IKK_NFkBedf, IkB_IKKr)


# #IkB adn NFkB cellular localization reactions
#
# Parameter('af', 0.09)
# Parameter('ncf', 0.012)
# Parameter('bf', 0.009)
# Parameter('ef', 0.045)
# Parameter('anf', 0.276)
# Parameter('anr', 0.828)
# Parameter('bnf', 0.0276)
# Parameter('bnr', 0.414)
# Parameter('enf', 0.138)
# Parameter('dnf', 0.276)
# Parameter('nf', 5.4)
# Parameter('nr', 0.0048)
#
# Rule('a_nc', IkBa(b=None, c=None, S='C') <> IkBa(b=None, c=None, S='N'), af, ncf)
# Rule('b_nc', IkBb(b=None, c=None, S='C') <> IkBb(b=None, c=None, S='N'), bf, ncf)
# Rule('e_nc', IkBe(b=None, c=None, S='C') <> IkBe(b=None, c=None, S='N'), ef, ncf)
# Rule('d_nc', IkBd(b=None, c=None, S='C') <> IkBd(b=None, c=None, S='N'), ef, ncf)
#
# Rule('an_nc', IkBa(b=1, c=None, S='C')%NFkB(b=1, S='C') <> IkBa(b=1, c=None, S='N')%NFkB(b=1, S='N'), anf, anr)
# Rule('bn_nc', IkBb(b=1, c=None, S='C')%NFkB(b=1, S='C') <> IkBb(b=1, c=None, S='N')%NFkB(b=1, S='N'), bnf, bnr)
# Rule('en_nc', IkBe(b=1, c=None, S='C')%NFkB(b=1, S='C') <> IkBe(b=1, c=None, S='N')%NFkB(b=1, S='N'), enf, bnr)
# Rule('dn_nc', IkBd(b=1, c=None, S='C')%NFkB(b=1, S='C') <> IkBd(b=1, c=None, S='N')%NFkB(b=1, S='N'), dnf, bnr)
#
# Rule('n_nc', NFkB(b=None, S='C') <> NFkB(b=None, S='N'), nf, nr)

#IkB Protein Degradation Reactions
# Rule('ad_c', IkBa(b=None, S='C') >> None, 0.12)
# Rule('bd_c', IkBb(b=None, S='C') >> None, 0.18)
# Rule('ed_c', IkBe(b=None, S='C') >> None, 0.18)
# Rule('dd_c', IkBd(b=None, S='C') >> None, 0.0014)
#
# Rule('ad_n', IkBa(b=None, S='N') >> None, 0.12)
# Rule('bd_n', IkBb(b=None, S='N') >> None, 0.18)
# Rule('ed_n', IkBe(b=None, S='N') >> None, 0.18)
# Rule('dd_n', IkBd(b=None, S='N') >> None, 0.0014)
#
# Rule('an_n', IkBa(b=1, S='C')%NFkB(b=1, S='C') >> NFkB(b=None, S='C'), 0.00006)
# Rule('bn_n', IkBb(b=1, S='C')%NFkB(b=1, S='C') >> NFkB(b=None, S='C'), 0.00006)
# Rule('en_n', IkBe(b=1, S='C')%NFkB(b=1, S='C') >> NFkB(b=None, S='C'), 0.00006)
# Rule('dn_n', IkBd(b=1, S='C')%NFkB(b=1, S='C') >> NFkB(b=None, S='C'), 0.00006)
#
# Rule('an_n', IkBa(b=1, S='N')%NFkB(b=1, S='N') >> NFkB(b=None, S='N'), 0.00006)
# Rule('bn_n', IkBb(b=1, S='N')%NFkB(b=1, S='N') >> NFkB(b=None, S='N'), 0.00006)
# Rule('en_n', IkBe(b=1, S='N')%NFkB(b=1, S='N') >> NFkB(b=None, S='N'), 0.00006)
# Rule('dn_n', IkBd(b=1, S='N')%NFkB(b=1, S='N') >> NFkB(b=None, S='N'), 0.00006)
#
# Rule('a2_c', IkBa(b=1, S='C')%IKK2(b=1, S='C') >> IKK2(b=None, S='C'), 0.0018)
# Rule('b2_c', IkBb(b=1, S='C')%IKK2(b=1, S='C') >> IKK2(b=None, S='C'), 0.0006)
# Rule('e2_c', IkBe(b=1, S='C')%IKK2(b=1, S='C') >> IKK2(b=None, S='C'), 0.0012)
# Rule('d1_c', IkBd(b=1, S='C')%IKK1(b=1, S='C') >> IKK1(b=None, S='C'), 0.0012)
#
# Rule('a2n_c', IkBa(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=1, c=1, S='C') >> IKK2(b=None, c=None, S='C') + NFkB(b=None, S='C'), (0.36))
# Rule('b2n_c', IkBb(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=1, c=1, S='C') >> IKK2(b=None, c=None, S='C') + NFkB(b=None, S='C'), (0.12))
# Rule('e2n_c', IkBe(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK2(b=1, c=1, S='C') >> IKK2(b=None, c=None, S='C') + NFkB(b=None, S='C'), (0.18))
# Rule('d1n_c', IkBd(b=None, c=1, S='C')%NFkB(b=1, S='C')%IKK1(b=1, c=1, S='C') >> IKK1(b=None, c=None, S='C') + NFkB(b=None, S='C'), (0.18))

generate_network(model)
generate_equations(model)

for i,ode in enumerate(model.odes):
     print i,":",ode

for i,reactions in enumerate(model.reactions):
     print i,":",reactions

