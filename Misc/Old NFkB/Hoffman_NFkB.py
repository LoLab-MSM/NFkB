__author__ = 'geena'

from pysb import *
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import scipy.interpolate

Model ()

result = []
with open('Old NFkB/_diswitch.txt','r') as inputFile:
    for line in inputFile:
        result.append(float(line))
ikk2 = np.asarray(result)
# print(result)
# quit()
#Declaration of monomers

""" Delcares IkBa, IkBb, IkBe, and IkBd monomers. All four of these monomers have binding sites to IKK2
    (for IkBa, IkBb, IkBe), and IKK1 (for IkBd). The monomers import and export between the Cytoplasm
    and the Nucleus.  Declares the IkBa, IkBb, IkBe, and IkBd mRNA monomers. The monomers go through constitutive RNA synthesis
    without the presence of NFkB in the nucleus,  and inducible RNA synthesis in the presence of NFkB in the nucleus. """

# def ikba_and_mRNA_monomers():
Monomer('IkBa', ['ikk', 'nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBa_mRNA')

# def ikbb_and_mRNA_monomers():
Monomer('IkBb', ['ikk', 'nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBb_mRNA')

# def ikbe_and_mRNA_monomers():
Monomer('IkBe', ['ikk', 'nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBe_mRNA')

# def ikbd_and_mRNA_monomers():
Monomer('IkBd', ['ikk', 'nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBd_mRNA')


# Declares NFkB, IKK1, and IKK2 which having a binding site to IkB and NFkB can transport between the Cytoplasm
# and the Nucleus. IKK1 and IKK2 only exist in the Cytoplasm.

# def nfkb_and_ikk_monomers():
Monomer('NFkB', ['ikb', 'S'], {'S': ['C', 'N']})

# def ikk1_monomer():
Monomer('IKK1', ['ikb', 'S'], {'S': ['C', 'N']})

# def ikk2_monomer():
Monomer('IKK2', ['ikb', 'S'], {'S': ['C', 'N']})

#Declaration of initial conditions
# def initial_ikk1_conditions():
Parameter('IKK1_0', 0.1) #Inhibitor Kinase Kinase a
Initial(IKK1(ikb=None, S='C'), IKK1_0)

# def initial_ikk2_conditions():
Parameter('IKK2_0', 0.1) #Inhibitor Kinase Kinase b
Initial(IKK2(ikb=None, S='C'), IKK2_0)

# def initial nfkb_conditions():
Parameter('NFkB_0', 0.125) #Nuclear Factor-kappaB
Initial(NFkB(ikb=None, S='N'), NFkB_0)

#RNA synthesis by NFkBn and Hill Coefficient

# Observable('IkBa_mRNA_obs', IkBa_mRNA())
# Observable('IkBa_obs', IkBa(b=ANY, c=ANY, S=ANY))
#Observables for NFkB
Observable('NFkBn_free', NFkB(ikb=None, S='N'))
Observable('IKK2_obs', IKK2(ikb = None, S='C'))
Observable('NFkBc_free', NFkB(ikb=None, S='C'))
Observable('NFkBn_obs', NFkB(ikb=WILD, S='N'))
Observable('NFkBn_bound', NFkB(ikb=ANY, S='N'))

Observable('IkBa_mRNA_obs', IkBa_mRNA())
Observable('IkBb_mRNA_obs', IkBb_mRNA())
Observable('IkBe_mRNA_obs', IkBe_mRNA())
Observable('IkBd_mRNA_obs', IkBd_mRNA())

Observable('IkBa_obs', IkBa(ikk = None, nfkb = None, S='C'))
Observable('IkBb_obs', IkBb(ikk = WILD, nfkb = WILD, S=WILD))
Observable('IkBe_obs', IkBe(ikk = WILD, nfkb = WILD, S=WILD))
Observable('IkBd_obs', IkBd(ikk = WILD, nfkb = WILD, S=WILD))


Parameter('hill', 3)
Parameter('a', 6)
Parameter('b', 0.25)
Parameter('e', 0.5)
Parameter('d', 0.025)

Expression('a_NFkBn', a*(NFkBn_free)**(hill))
Expression('b_NFkBn', b*(NFkBn_free)**(hill))
Expression('e_NFkBn', e*(NFkBn_free)**(hill))
Expression('d_NFkBn', d*(NFkBn_free)**(hill))

Rule('an_mRNA', None >> IkBa_mRNA(), a_NFkBn)
Rule('bn_mRNA', None >> IkBb_mRNA(), b_NFkBn)
Rule('en_mRNA', None >> IkBe_mRNA(), e_NFkBn)
Rule('dn_mRNA', None >> IkBd_mRNA(), d_NFkBn)

# IkB mRNA and protein synthesis reactions

Parameter('mRNA_a', 0.035)
Rule('a_mRNA', IkBa_mRNA() >> None, mRNA_a)


Parameter('mRNA_b', 3e-3)
Parameter('mRNA_e', 4e-3)
Parameter('mRNA_d', 2e-3)

Rule('b_mRNA', IkBb_mRNA() >> None, mRNA_b)
Rule('e_mRNA', IkBe_mRNA() >> None, mRNA_e)
Rule('d_mRNA', IkBd_mRNA() >> None, mRNA_d)


# Parameter('synth', 0.2448)
# Rule('a_psynth', None >> IkBa(ikk=None, nfkb=None, S='C'), synth)
# Rule('b_psynth', None >> IkBb(ikk=None, nfkb=None, S='C'), synth)
# Rule('e_psynth', None >> IkBe(ikk=None, nfkb=None, S='C'), synth)
# Rule('d_psynth', None >> IkBd(ikk=None, nfkb=None, S='C'), synth)

Parameter('syntha', 0.2448)
Rule('a_psynth', IkBa_mRNA() >> IkBa(ikk=None, nfkb=None, S='C') + IkBa_mRNA(), syntha)

Parameter('synthb', 0.2448)
Rule('b_psynth', IkBb_mRNA() >> IkBb(ikk=None, nfkb=None, S='C') + IkBb_mRNA(), synthb)

Parameter('synthe', 0.2448)
Rule('e_psynth', IkBe_mRNA() >> IkBe(ikk=None, nfkb=None, S='C') + IkBe_mRNA(), synthe)

Parameter('synthd', 0.2448)
Rule('d_psynth', IkBd_mRNA() >> IkBd(ikk=None, nfkb=None, S='C') + IkBd_mRNA(), synthd)


Parameter('psynth_a', 2e-4)
Parameter('psynth_b', 1e-5)
Parameter('psynth_e', 3e-6)
Parameter('psynth_d', 1e-7)
Rule('a_synth', None >> IkBa_mRNA(), psynth_a)
Rule('b_synth', None >> IkBb_mRNA(), psynth_b)
Rule('e_synth', None >> IkBe_mRNA(), psynth_e)
Rule('d_synth', None >> IkBd_mRNA(), psynth_d)

#IkB(a,b,e) association and dissociation from IKK2 and IkBd association and dissociation from IKK2

Parameter('a_2f', 1.35)
Parameter('a_2r', 0.075)
Parameter('b_kf', 0.36)
Parameter('b_2r', 0.105)
Parameter('e_2r', 0.105)
Parameter('d_1r', 0.105)
Parameter('e_2f', 0.54)
Parameter('d_1f', 0.54)
Rule('a2_adc', IkBa(ikk=None, nfkb=None, S='C') + IKK2(ikb=None, S='C') <> IkBa(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C'), a_2f, a_2r)
Rule('b2_adc', IkBb(ikk=None, nfkb=None, S='C') + IKK2(ikb=None, S='C') <> IkBb(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C'), b_kf, b_2r)
Rule('e2_adc', IkBe(ikk=None, nfkb=None, S='C') + IKK2(ikb=None, S='C') <> IkBe(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C'), e_2f, e_2r)
Rule('d1_adc', IkBd(ikk=None, nfkb=None, S='C') + IKK1(ikb=None, S='C') <> IkBd(ikk=1, nfkb=None, S='C')%IKK1(ikb=1, S='C'), d_1f, d_1r)



Parameter('IkB_IKKf', 30)
Parameter('IkB_IKKr', 6e-5)
Rule('a2n_c1', IkBa(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C') + NFkB(ikb=None, S='C') <> IkBa(ikk=1, nfkb=2, S='C')%IKK2(ikb=1, S='C')%NFkB(ikb=2, S='C'), IkB_IKKf, IkB_IKKr)
Rule('b2n_c1', IkBb(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C') + NFkB(ikb=None, S='C') <> IkBb(ikk=1, nfkb=2, S='C')%IKK2(ikb=1, S='C')%NFkB(ikb=2, S='C'), IkB_IKKf, IkB_IKKr)
Rule('e2n_c1', IkBe(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C') + NFkB(ikb=None, S='C') <> IkBe(ikk=1, nfkb=2, S='C')%IKK2(ikb=1, S='C')%NFkB(ikb=2, S='C'), IkB_IKKf, IkB_IKKr)
Rule('d1n_c1', IkBd(ikk=1, nfkb=None, S='C')%IKK1(ikb=1, S='C') + NFkB(ikb=None, S='C') <> IkBd(ikk=1, nfkb=2, S='C')%IKK1(ikb=1, S='C')%NFkB(ikb=2, S='C'), IkB_IKKf, IkB_IKKr)


Parameter('an_2f', 11.1)
Parameter('bn_2f', 2.88)
Parameter('en_2f', 4.2)
Parameter('dn_1f', 4.2)
Parameter('an_2r', 0.075)
Parameter('bn_2r', 0.105)
Parameter('en_2r', 0.105)
Parameter('dn_1r', 0.105)
Rule('a2n_c', IkBa(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK2(ikb=None, S='C') <> IkBa(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK2(ikb=2, S='C'), an_2f, an_2r)
Rule('b2n_c', IkBb(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK2(ikb=None, S='C') <> IkBb(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK2(ikb=2, S='C'), bn_2f, bn_2r)
Rule('e2n_c', IkBe(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK2(ikb=None, S='C') <> IkBe(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK2(ikb=2, S='C'), en_2f, en_2r)
Rule('d1n_c', IkBd(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK1(ikb=None, S='C') <> IkBd(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK1(ikb=2, S='C'), dn_1f, dn_1r)


# Parameter('IkB_IKKf', 30)
# Parameter('IkB_IKKr', 6e-5)
Rule('an_adc', IkBa(ikk=None, nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBa(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('bn_adc', IkBb(ikk=None, nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBb(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('en_adc', IkBe(ikk=None, nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBe(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('dn_adc', IkBd(ikk=None, nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBd(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)

Rule('an_adn', IkBa(ikk=None, nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBa(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('bn_adn', IkBb(ikk=None, nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBb(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('en_adn', IkBe(ikk=None, nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBe(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('dn_adn', IkBd(ikk=None, nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBd(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)


# #IkB and NFkB cellular localization reactions

Parameter('af', 0.09)
Parameter('bf', 0.009)
Parameter('ef', 0.045)
Parameter('df', 0.045)
Parameter('ancf', 0.012)
Parameter('bncf', 0.012)
Parameter('encf', 0.012)
Parameter('dncf', 0.012)
Rule('a_nc', IkBa(ikk=None, nfkb=None, S='C') <> IkBa(ikk=None, nfkb=None, S='N'), af, ancf)
Rule('b_nc', IkBb(ikk=None, nfkb=None, S='C') <> IkBb(ikk=None, nfkb=None, S='N'), bf, bncf)
Rule('e_nc', IkBe(ikk=None, nfkb=None, S='C') <> IkBe(ikk=None, nfkb=None, S='N'), ef, encf)
Rule('d_nc', IkBd(ikk=None, nfkb=None, S='C') <> IkBd(ikk=None, nfkb=None, S='N'), ef, dncf)

Parameter('anf', 0.276)
Parameter('bnf', 0.0276)
Parameter('enf', 0.138)
Parameter('dnf', 0.276)
Parameter('anr', 0.828)
Parameter('bnr', 0.414)
Parameter('enr', 0.414)
Parameter('dnr', 0.414)
Rule('an_nc', IkBa(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBa(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), anf, anr)
Rule('bn_nc', IkBb(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBb(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), bnf, bnr)
Rule('en_nc', IkBe(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBe(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), enf, enr)
Rule('dn_nc', IkBd(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBd(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N'), dnf, dnr)

Parameter('nf', 5.4)
Parameter('nr', 0.0048)
Rule('n_nc', NFkB(ikb=None, S='C') <> NFkB(ikb=None, S='N'), nf, nr)

# IkB Protein Degradation Reactions

Parameter('a_dc', 0.12)
Parameter('b_dc', 0.18)
Parameter('e_dc', 0.18)
Parameter('d_dc', 0.0014)
Parameter('a_dn', 0.12)
Parameter('b_dn', 0.18)
Parameter('e_dn', 0.18)
Parameter('d_dn', 0.0014)

Rule('ad_c', IkBa(ikk=None, nfkb=None, S='C') >> None, a_dc)
Rule('bd_c', IkBb(ikk=None, nfkb=None, S='C') >> None, b_dc)
Rule('ed_c', IkBe(ikk=None, nfkb=None, S='C') >> None, e_dc)
Rule('dd_c', IkBd(ikk=None, nfkb=None, S='C') >> None, d_dc)

Rule('ad_n', IkBa(ikk=None, nfkb=None, S='N') >> None, a_dn)
Rule('bd_n', IkBb(ikk=None, nfkb=None, S='N') >> None, b_dn)
Rule('ed_n', IkBe(ikk=None, nfkb=None, S='N') >> None, e_dn)
Rule('dd_n', IkBd(ikk=None, nfkb=None, S='N') >> None, d_dn)

Parameter('c_bn', 0.00006)
Parameter('n_bn', 0.00006)

Rule('an_c', IkBa(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('bn_c', IkBb(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('en_c', IkBe(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('dn_c', IkBd(ikk=None, nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)

Rule('an_n', IkBa(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('bn_n', IkBb(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('en_n', IkBe(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('dn_n', IkBd(ikk=None, nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)

Parameter('a2_2', 0.0018)
Parameter('b2_2', 0.0006)
Parameter('ed21_21', 0.0012)

Rule('a2_c', IkBa(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C') >> IKK2(ikb=None, S='C'), a2_2)
Rule('b2_c', IkBb(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C') >> IKK2(ikb=None, S='C'), b2_2)
Rule('e2_c', IkBe(ikk=1, nfkb=None, S='C')%IKK2(ikb=1, S='C') >> IKK2(ikb=None, S='C'), ed21_21)
Rule('d1_c', IkBd(ikk=1, nfkb=None, S='C')%IKK1(ikb=1, S='C') >> IKK1(ikb=None, S='C'), ed21_21)

Parameter('an2_2n', 0.36)
Parameter('bn2_2n', 0.12)
Parameter('edn21_21n', 0.18)

Rule('a2n_2nc', IkBa(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK2(ikb=2, S='C') >> IKK2(ikb=None, S='C') + NFkB(ikb=None, S='C'), an2_2n)
Rule('b2n_2nc', IkBb(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK2(ikb=2, S='C') >> IKK2(ikb=None, S='C') + NFkB(ikb=None, S='C'), bn2_2n)
Rule('e2n_2nc', IkBe(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK2(ikb=2, S='C') >> IKK2(ikb=None, S='C') + NFkB(ikb=None, S='C'), edn21_21n)
Rule('d1n_1nc', IkBd(ikk=2, nfkb=1, S='C')%NFkB(ikb=1, S='C')%IKK1(ikb=2, S='C') >> IKK1(ikb=None, S='C') + NFkB(ikb=None, S='C'), edn21_21n)



#Dictionary to substitute in species names to match matlab files
species_dict = {
    0: 'IKK1',
    1: 'IKK',
    2: 'NFkBn',
    3: 'SOURCE',
    4: 'IkBat',
    5: 'IkBbt',
    6: 'IkBet',
    7: 'IkBdt',
    8: 'NFkB',
    9: 'SINK',
    10: 'IkBa',
    11: 'IkBb',
    12: 'IkBe',
    13: 'IkBd',
    14: 'IkBaIKK',
    15: 'IkBbIKK',
    16: 'IkBeIKK',
    17: 'IkBdIKK1',
    18: 'IkBaNFkB',
    19: 'IkBbNFkB',
    20: 'IkBeNFkB',
    21: 'IkBdNFkB',
    22: 'IkBan',
    23: 'IkBbn',
    24: 'IkBen',
    25: 'IkBdn',
    26: 'IkBaIKKNFkB',
    27: 'IkBbIKKNFkB',
    28: 'IkBeIKKNFkB',
    29: 'IkBdIKK1NFkB',
    30: 'IkBaNFkBn',
    31: 'IkBbNFkBn',
    32: 'IkBeNFkBn',
    33: 'IkBdNFkBn'
}

#updating the species names in the odes
for  ode in model.odes:
    for i in range(len(model.species)):
       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
    print ode

#equilibrium phase 1 of model using odesolve
time_equil = np.linspace(-8000, 0, 8001)
equil = odesolve(model, time_equil, verbose=True)

plt.figure(1)
plt.plot(time_equil/60., equil['NFkBn_free']*1000, label= NFkBn_free.name)
plt.xlabel("Time (in hours)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
plt.ylim(ymin = 0, ymax =100)
plt.legend(loc=0)

last_conc = [equil[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]]
print(last_conc)#setting previous conc species equal to none

#phase 2 of model using Solver
tspan = np.linspace(0, 720, 721)
solver = Solver(model,tspan,verbose=True)

# ikk2_vals = []
nsims = len(tspan) - 1 #simulation time is 720
# yobs = np.empty((nsims, len(model.observables)))
rows = len(tspan) # nSims
cols = len(model.observables) #length of observables
yobs = np.empty((rows, cols)) #creating empty array of R and C

assert model.observables[1].name == 'IKK2_obs'
assert str(model.species[1]) == "IKK2(ikb=None, S='C')"

# print(model.parameters['IKK2_0'].value)
#using imported IKK2 values during phase 2 simulation
for i in range(nsims):
    print(i) #printing each simulation iteration
    last_conc[1] = ikk2[i]
    solver.tspan = [tspan[i], tspan[i+1]] #simulating in 2 step intervals
    solver.run(y0 = last_conc) #setting initial concentrations to previous simulation value (equil phase 1)
    # print(last_conc)
    if i == 0:
        yobs[0, :] = solver.yobs_view[0]
    yobs[i + 1 , :] = solver.yobs_view[1] #taking each simulation index and all obs and keeping the last obs in last time point of simulation
    last_conc = solver.y[1,:] # running solver taking first row and all simulation species


#plotting NFkBn free
plt.figure(1)
plt.plot(tspan/60., yobs[:,0]*1000, label= NFkBn_free.name)
plt.xlabel("Time (in hours)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
plt.ylim(ymin = 0, ymax =100)
plt.legend(loc=0)

# #plotting NFkBc free
# plt.figure(2)
# plt.plot(tspan/60., yobs[:,2], label= NFkBc_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)
#
# #plotting IkBa unbound cytosol
# plt.figure(2)
# plt.plot(tspan/60., yobs[:,10], label= IkBa_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)


plt.show()



