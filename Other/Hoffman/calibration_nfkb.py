# __author__ = 'geena'

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
#from pysb.simulator import ScipyOdeSimulator
import numpy as np


def run_sim_with_perturbations(sim, tspan, perturbations):
    perturbation_times = sorted(perturbations.keys())

    for t in perturbation_times:
        if not any(np.isclose(t, tspan)):
            raise ValueError('Perturbation time %f is not defined in tspan'
                             % t)

    # Create the tspans for each simulation stage
    tspans = np.split(tspan, [np.abs(tspan-value).argmin() for value in
                              perturbation_times])
    for i in range(len(tspans) - 1):
        tspans[i] = np.append(tspans[i], tspans[i + 1][0])

    # Do an initial simulation with the above setup until the point at which
    # we want to change a parameter (perturbation_time)
    res = sim.run(tspan=tspans[0])

    df_out = res.dataframe

    for i, t in enumerate(perturbation_times, 1):
        res = sim.run(initials=res.species[-1],
                      param_values=perturbation_params[t],
                      tspan=tspans[i])

        df_out = pd.concat([df_out, res.dataframe.iloc[1:]])

    return df_out

data = pd.read_csv('/Users/geenaildefonso/data_ikk.csv')
IKK_flux = tuple(data.loc[:, "IKK_flux"])
print(np.transpose(IKK_flux))
# print(IKK.shape)
# print(IKK)
# my_dict = {}
#
# # you have ten values, so loop using range
# for i in range(721):
# 	my_dict[i] = {'some_key': IKK_flux[i]}
#
# # print("MY DICT IN SORTED TUPLE FORM\n")
# print(sorted(my_dict.items(), key=lambda x: x[0]))

# -----------------------------
# convert outer key values into inner key values
# {int: {'some_key': float} } --> {float: {'some_key': int} }


# make new dict
# my_new_dict = {}
#
# for keys in my_dict:
# 	my_new_key = my_dict[keys]['some_key']
# 	my_new_dict[my_new_key] = {'some_key': keys}
#
# print("\n\nMY NEW DICT IN SORTED TUPLE FORM\n")
# print(sorted(my_new_dict.items(), key=lambda x: x[0]))





Model ()

#Declaration of monomers

""" Delcares IkBa, IkBb, IkBe, and IkBd monomers. All four of these monomers have binding sites to IKK2
    (for IkBa, IkBb, IkBe), and IKK1 (for IkBd). The monomers import and export between the Cytoplasm
    and the Nucleus.  Declares the IkBa, IkBb, IkBe, and IkBd mRNA monomers. The monomers go through constitutive RNA synthesis
    without the presence of NFkB in the nucleus,  and inducible RNA synthesis in the presence of NFkB in the nucleus. """

# def ikba_and_mRNA_monomers():
Monomer('IkBa', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBa_mRNA')

# def ikbb_and_mRNA_monomers():
Monomer('IkBb', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBb_mRNA')

# def ikbe_and_mRNA_monomers():
Monomer('IkBe', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBe_mRNA')

# def ikbd_and_mRNA_monomers():
Monomer('IkBd', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBd_mRNA')


# Declares NFkB, IKK1, and IKK2 which having a binding site to IkB and NFkB can transport between the Cytoplasm
# and the Nucleus. IKK1 and IKK2 only exist in the Cytoplasm.

# def nfkb_and_ikk_monomers():
Monomer('NFkB', ['ikb', 'S'], {'S': ['C', 'N']})

# def ikk1_monomer():
# Monomer('IKK1', ['S'], {'S': ['C', 'N']})
#
# # def ikk2_monomer():
# Monomer('IKK2', ['S'], {'S': ['C', 'N']})
# def a20_monomers():
Monomer('A20')
Monomer('A20t')

# def ligand_to_receptor_monomers():
Monomer('TNF', ['c', 'tnfr'])
Monomer('TNFRM')
Monomer('TNFR', ['tnf'])

# def complex_monomers():
Monomer('C1', ['tnf', 'state'], {'state': ['a', 'i']})
Monomer('TTR')

# def kinase_monomers():
Monomer('IKKK', ['state'], {'state': ['a', 'i']})
Monomer('IKK', ['state'], {'state': ['a', 'i', 'ai']})


#Declaration of initial conditions
# def initial_ikk1_conditions():
# def initial_conditions():
# Parameter('IKKK_0', 0.1) #Inhibitor Kinase Kinase Kinase
# Initial(IKKK(state = 'i'), IKKK_0)

# def initial_ikk2_conditions():
# Parameter('IKK_0', 0.1) #Inhibitor Kinase Kinase
# Initial(IKK(state = 'i'), IKK_0)

# def initial nfkb_conditions():
# Parameter('NFkB_0', 0.125) #Nuclear Factor-kappaB
# Initial(NFkB(ikb=None, S='N'), NFkB_0)

# def initial TTR_conditions():
# Parameter('TTR_0', 8.3e-4) #TRADD-TRAF-RIP
# Initial(TTR(), TTR_0)

# #IkBa initial conditions
Parameter('IkBa_0', 0.0025) #IkBa
Initial(IkBa(nfkb = None, S = 'C'), IkBa_0)

Parameter('IkBan_0', 0.0013) #IkBa
Initial(IkBa(nfkb = None, S = 'N'), IkBan_0)

Parameter('IkBa_NFkB_0', 0.0621) #IkBa
Initial(IkBa(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBa_NFkB_0)

Parameter('IkBa_NFkBn_0', 0.0208) #IkBa
Initial(IkBa(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBa_NFkBn_0)

Parameter('IkBat_0', 0.0020) #IkBa
Initial(IkBa_mRNA(), IkBat_0)
#
# #IkBb initial conditions
Parameter('IkBb_0', 0.0044) #IkBb
Initial(IkBb(nfkb = None, S = 'C'), IkBb_0)

Parameter('IkBbn_0', 0.0002) #IkBa
Initial(IkBb(nfkb = None, S = 'N'), IkBbn_0)

Parameter('IkBb_NFkB_0', 0.0237) #IkBa
Initial(IkBb(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBb_NFkB_0)

Parameter('IkBb_NFkBn_0', 0.0016) #IkBb
Initial(IkBb(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBb_NFkBn_0)

Parameter('IkBbt_0', 0.0033) #IkBa
Initial(IkBb_mRNA(), IkBbt_0)
#
# #IkBe initial conditions
Parameter('IkBe_0', 0.0003) #IkBa
Initial(IkBe(nfkb = None, S = 'C'), IkBe_0)

Parameter('IkBen_0', 0.0001) #IkBa
Initial(IkBe(nfkb = None, S = 'N'), IkBen_0)

Parameter('IkBeNFkB_0', 0.0045) #IkBa
Initial(IkBe(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBeNFkB_0)

Parameter('IkBeNFkBn_0', 0.0015) #IkBa
Initial(IkBe(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBeNFkBn_0)

Parameter('IkBet_0', 0.0003) #IkBa
Initial(IkBe_mRNA(), IkBet_0)

#IkBd initial conditions
Parameter('IkBd_0', 0.0003) #IkBa
Initial(IkBd(nfkb = None, S = 'C'), IkBd_0)

Parameter('IkBdn_0', 0.0003) #IkBa
Initial(IkBd(nfkb = None, S = 'N'), IkBdn_0)

Parameter('IkBdNFkB_0', 0.0058) #IkBa
Initial(IkBd(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBdNFkB_0)

Parameter('IkBdNFkBn_0', 0.0039) #IkBa
Initial(IkBd(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBdNFkBn_0)

Parameter('IkBdt_0', 0.0001) #IkBa
Initial(IkBd_mRNA(), IkBdt_0)

Parameter('NFkB_0', 0.000) #Nuclear Factor-kappaB
Initial(NFkB(ikb=None, S='C'), NFkB_0)

Parameter('NFkBn_0', 0.0012) #Nuclear Factor-kappaB in nucleus
Initial(NFkB(ikb=None, S='N'), NFkBn_0)

Parameter('IKKK_off_0', 0.1000) #Inhibitor Kinase Kinase Kinase
Initial(IKKK(state = 'i'), IKKK_off_0)

Parameter('IKKK_0', 0.000) #Inhibitor Kinase Kinase Kinase
Initial(IKKK(state = 'a'), IKKK_0)

Parameter('IKK_off_0', 0.0985) #TRADD-TRAF-RIP
Initial(IKK(state = 'i'), IKK_off_0)

Parameter('IKK_0', 0.0002) #TRADD-TRAF-RIP
Initial(IKK(state = 'a'), IKK_0)

Parameter('IKKi_0', 0.0013) #TRADD-TRAF-RIP
Initial(IKK(state = 'ai'), IKKi_0)


Parameter('TNF_0', .2) #TNF
Initial(TNF(c = None, tnfr = None), TNF_0)

Parameter('TNFRM_0', 0.0) #TRADD-TRAF-RIP
Initial(TNFRM(), TNFRM_0)

Parameter('TNFR_0', 0.0) #TNFR Tumor Necrosis Factor Receptor 1
Initial(TNFR(tnf = None), TNFR_0)

Parameter('TNFRtnf_0', 0.0) #TNF binding TNFR
Initial(TNF(c = None, tnfr = 1)%TNFR(tnf = 1), TNFRtnf_0)

Parameter('C1_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = None, state = 'a'), C1_0)

Parameter('C1_off_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = None, state = 'i'), C1_off_0)

Parameter('C1tnf_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = 1, state = 'a')%TNF(c =  1, tnfr = None), C1tnf_0)

Parameter('C1tnf_off_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = 1, state = 'i')%TNF(c =  1, tnfr = None), C1tnf_off_0)

Parameter('TTR_0', 8.3e-4) #TRADD-TRAF-RIP
Initial(TTR(), TTR_0)

Parameter('a20_0', 0.0049) #TRADD-TRAF-RIP
Initial(A20(), a20_0)

Parameter('a20t_0', 0.0001) #TRADD-TRAF-RIP
Initial(A20t(), a20t_0)

# Parameter('TNF_0', 0.2)
# Initial(TNF(c = None, tnfr = None), TNF_0)


# def  ikb_mrna_to_ikb():
Parameter('psynth_a', 7e-5)
Parameter('psynth_b', 1e-5)
Parameter('psynth_e', 1e-6)
Parameter('psynth_d', 1e-7)
Rule('a_synth', None >> IkBa_mRNA(), psynth_a)
Rule('b_synth', None >> IkBb_mRNA(), psynth_b)
Rule('e_synth', None >> IkBe_mRNA(), psynth_e)
Rule('d_synth', None >> IkBd_mRNA(), psynth_d)

Parameter('hill', 3)
# Parameter('a', 0)
# Parameter('b', 0)
# Parameter('e', 0)
# Parameter('d', 0)
# Parameter('a', 8)
# Parameter('b', 0.02)
# Parameter('e', 0.3)
# Parameter('d', 0.025)

Observable('NFkBn_free', NFkB(ikb=None, S='N'))

# Parameter('zero', 0)
# Parameter('a_NFkBn', 0)
# Parameter('b_NFkBn', 0)
# Parameter('e_NFkBn', 0)
# Parameter('d_NFkBn', 0)

Parameter('a', 0)
Parameter('b', 0)
Parameter('e', 0)
Parameter('d', 0)

Expression('a_NFkBn', a*model.observables['NFkBn_free']**(hill))
Expression('b_NFkBn', b*model.observables['NFkBn_free']**(hill))
Expression('e_NFkBn', e*model.observables['NFkBn_free']**(hill))
Expression('d_NFkBn', d*model.observables['NFkBn_free']**(hill))

Rule('an_mRNA', None >> IkBa_mRNA(), a_NFkBn)
Rule('bn_mRNA', None >> IkBb_mRNA(), b_NFkBn)
Rule('en_mRNA', None >> IkBe_mRNA(), e_NFkBn)
Rule('dn_mRNA', None >> IkBd_mRNA(), d_NFkBn)
#
# my_dict = {}
# #
# # # you have ten values, so loop using range
# for i in range(721):
# 	my_dict[i] = {'': IKK_flux[i]}
#
# #
# # # print("MY DICT IN SORTED TUPLE FORM\n")
# # print(sorted(my_dict.items(), key=lambda x: x[0]))
# #
# perturbation_params = {sorted(my_dict.items(), key=lambda x: x[0])}



# perturbation_params = {
#     10: {'a_NFkBn': 8*model.observables['NFkBn_free']**3},
#     37: {'b_NFkBn': .02*model.observables['NFkBn_free']**3, 'e_NFkBn': .3*model.observables['NFkBn_free']**3},
#     90: {'d_NFkBn': .025*model.observables['NFkBn_free']**3}
# }


#A20 mRNA and Protein Synthesis and Degradation Reactions
# def a20_mrna_to_a20():
Parameter('A20_mRNA', 2e-6)
Parameter('A20n', 0.4)
Parameter('A20_mRNA_c_deg', 0.035)
Parameter('a1d_c_deg', 0.36)
Parameter('A20_synth', 0.25)
Parameter('A20_deg', 0.0029)

Observable('obs_A20t', A20t())

Expression('A20t_NFkBn', A20n*(NFkBn_free)**(hill))
Expression('A20_synthesis', A20_synth*model.observables['obs_A20t'])

Rule('A20t_synth', None >> A20t(), A20_mRNA)
Rule('A20t_mediated_nfkbn', None >> A20t(), A20t_NFkBn)
Rule('A20t_deg', A20t() >> None, A20_mRNA_c_deg)

Rule('synth_A20', None >> A20(), A20_synthesis)
Rule('deg_A20', A20() >> None, A20_deg)



Parameter('nfkb_synthb', 0)
Parameter('nfkb_synthe', 0)
Parameter('nfkb_synthd', 0)



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
Rule('a_psynth', IkBa_mRNA() >> IkBa(nfkb=None, S='C') + IkBa_mRNA(), syntha)

Parameter('synthb', 0.2448)
Rule('b_psynth', IkBb_mRNA() >> IkBb(nfkb=None, S='C') + IkBb_mRNA(), synthb)

Parameter('synthe', 0.2448)
Rule('e_psynth', IkBe_mRNA() >> IkBe(nfkb=None, S='C') + IkBe_mRNA(), synthe)

Parameter('synthd', 0.2448)
Rule('d_psynth', IkBd_mRNA() >> IkBd(nfkb=None, S='C') + IkBd_mRNA(), synthd)

#IkB(a,b,e) association and dissociation from IKK2 and IkBd association and dissociation from IKK2
# def ikb_assoc_diss_nfkb():

Parameter('IkB_IKKf', 30)
Parameter('IkB_IKKr', 6e-5)
Rule('an_adc', IkBa(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('bn_adc', IkBb(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('en_adc', IkBe(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('dn_adc', IkBd(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)

Rule('an_adn', IkBa(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('bn_adn', IkBb(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('en_adn', IkBe(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('dn_adn', IkBd(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)


# #IkB and NFkB cellular localization reactions
# def ikb_nfkb_localization():
Parameter('af', 0.09)
Parameter('bf', 0.009)
Parameter('ef', 0.045)
Parameter('df', 0.045)
Parameter('ancf', 0.012)
Parameter('bncf', 0.012)
Parameter('encf', 0.012)
Parameter('dncf', 0.012)
Rule('a_nc', IkBa(nfkb=None, S='C') <> IkBa(nfkb=None, S='N'), af, ancf)
Rule('b_nc', IkBb(nfkb=None, S='C') <> IkBb(nfkb=None, S='N'), bf, bncf)
Rule('e_nc', IkBe(nfkb=None, S='C') <> IkBe(nfkb=None, S='N'), ef, encf)
Rule('d_nc', IkBd(nfkb=None, S='C') <> IkBd(nfkb=None, S='N'), ef, dncf)

Parameter('anf', 0.276)
Parameter('bnf', 0.0276)
Parameter('enf', 0.138)
Parameter('dnf', 0.276)
Parameter('anr', 0.828)
Parameter('bnr', 0.414)
Parameter('enr', 0.414)
Parameter('dnr', 0.414)
Rule('an_nc', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'), anf, anr)
Rule('bn_nc', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'), bnf, bnr)
Rule('en_nc', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'), enf, enr)
Rule('dn_nc', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'), dnf, dnr)

Parameter('nf', 5.4)
Parameter('nr', 0.0048)
Rule('n_nc', NFkB(ikb=None, S='C') <> NFkB(ikb=None, S='N'), nf, nr)

# IkB Protein Degradation Reactions
# def ikb_deg_reactions():
Parameter('a_dc', 0.12)
Parameter('b_dc', 0.18)
Parameter('e_dc', 0.18)
Parameter('d_dc', 0.0014)
Parameter('a_dn', 0.12)
Parameter('b_dn', 0.18)
Parameter('e_dn', 0.18)
Parameter('d_dn', 0.0014)

Rule('ad_c', IkBa(nfkb=None, S='C') >> None, a_dc)
Rule('bd_c', IkBb(nfkb=None, S='C') >> None, b_dc)
Rule('ed_c', IkBe(nfkb=None, S='C') >> None, e_dc)
Rule('dd_c', IkBd(nfkb=None, S='C') >> None, d_dc)

Rule('ad_n', IkBa(nfkb=None, S='N') >> None, a_dn)
Rule('bd_n', IkBb(nfkb=None, S='N') >> None, b_dn)
Rule('ed_n', IkBe(nfkb=None, S='N') >> None, e_dn)
Rule('dd_n', IkBd(nfkb=None, S='N') >> None, d_dn)

Parameter('c_bn', 0.00006)
Parameter('n_bn', 0.00006)

Rule('an_c', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C') , c_bn)
Rule('bn_c', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('en_c', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('dn_c', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)

Rule('an_n', IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('bn_n', IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('en_n', IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('dn_n', IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)

#Added Reactions

#IKK-mediated IkB degradation reactions
# def ikb_ikk_mediated_deg():
Observable('IKKa_obs', IKK(state = 'a'))

Parameter('a_f_deg', 0.36/.1)
Parameter('b_f_deg', 0.12/.1)
Parameter('e_f_deg', 0.18/.1)
Parameter('d_f_deg', 0.18*.01)

Parameter('and_c_n', 0.36)
Parameter('bnd_c_n', 0.12)
Parameter('end_c_n', 0.18)
Parameter('dnd_c_n', 0.18*.01)
# Parameter('IKK_boom', 0)

# Expression('IKK_ikba_flux', IKK_boom*a_f_deg)
# Expression('IKK_ikbb_flux', IKK_boom*b_f_deg)
# Expression('IKK_ikbe_flux', IKK_boom*e_f_deg)

Expression('IKK_ikba_flux', model.observables['IKKa_obs']*a_f_deg)
Expression('IKK_ikbb_flux', model.observables['IKKa_obs']*b_f_deg)
Expression('IKK_ikbe_flux', model.observables['IKKa_obs']*e_f_deg)
#
# my_dict = {}
# #
# # # you have ten values, so loop using range
# for i in range(721):
#     my_dict[i] = {'IKK_boom': IKK_flux[i]}
#
# #
# # # print("MY DICT IN SORTED TUPLE FORM\n")
# # print(sorted(my_dict.items(), key=lambda x: x[0]))
# perturbation_params = sorted(my_dict.items(), key=lambda x: x[0])
# for i in tuple(range(721)):
#     perturbation_params = {
#         i: {'a_f_deg': [w*.36 for w in IKK_flux]}
#     }

# for i in tuple(range(720)):
# perturbation_params = {
#     1: {'a_f_deg': .36}
# }

# IkBa >> None
# IkBb >> None
# IkBe >> None
# IkBd >> None

# IkBa : NFkB >> NFkB
# IkBb : NFkB >> NFkB
# IkBe : NFkB >> NFkB
# IkBd : NFkB >> NFkB

Rule('a_c_deg', IkBa(nfkb=None, S='C') >> None, IKK_ikba_flux)
Rule('b_c_deg', IkBb(nfkb=None, S='C') >> None, IKK_ikbb_flux)
Rule('e_c_deg', IkBe(nfkb=None, S='C') >> None, IKK_ikbe_flux)
Rule('d_c_deg', IkBd(nfkb=None, S='C') >> None, d_f_deg)

Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikba_flux)
Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikbb_flux)
Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikbe_flux)
Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), dnd_c_n)


# Rule('a_c_deg', IkBa(nfkb=None, S='C') >> None, IKK_ikba_flux)
# Rule('b_c_deg', IkBb(nfkb=None, S='C') >> None, b_f_deg)
# Rule('e_c_deg', IkBe(nfkb=None, S='C') >> None, e_f_deg)
# Rule('d_c_deg', IkBd(nfkb=None, S='C') >> None, d_f_deg)
#
# Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), and_c_n)
# Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), bnd_c_n)
# Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), end_c_n)
# Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), dnd_c_n)

#
# Rule('a_c_deg', IkBa(nfkb=None, S='C') + IKK(state = 'a') >> IKK(state = 'a'), a_f_deg)
# Rule('b_c_deg', IkBb(nfkb=None, S='C') + IKK(state = 'a') >> IKK(state = 'a'), b_f_deg)
# Rule('e_c_deg', IkBe(nfkb=None, S='C') + IKK(state = 'a') >> IKK(state = 'a'), e_f_deg)
# Rule('d_c_deg', IkBd(nfkb=None, S='C') + IKK(state = 'a') >> IKK(state = 'a'), d_f_deg)
#
# Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK(state = 'a') >> NFkB(ikb=None, S='C') + IKK(state = 'a'), and_c_n)
# Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK(state = 'a')>> NFkB(ikb=None, S='C') + IKK(state = 'a'), bnd_c_n)
# Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK(state = 'a') >> NFkB(ikb=None, S='C') + IKK(state = 'a'), end_c_n)
# Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') + IKK(state = 'a')>> NFkB(ikb=None, S='C') + IKK(state = 'a'), dnd_c_n)



#IKK Activation Module
#TNF-Independent Complex 1 Activity Reactions
# def tnf_independent_to_c1():
Parameter('synth_tnfrm', 2e-7)
Parameter('deg_tnfrm', 0.0058)
Parameter('tnfr_f_tnfrm', 1e-5*3) # times 6
Parameter('tnfr_r_tnfrm', 0.1*3)
Parameter('deg_TNFR', 0.023)
Parameter('TNFR_TTR_f_C1', 100.0)
Parameter('TNFR_TTR_r_C1', 0.1)
Parameter('C1_f_a', 30.0)
Parameter('C1_r_a', 2.0)
Parameter('C1_f_A20', 1000.0)
Parameter('C1_f_TNFR_TTR', 0.1)
Parameter('C1_i_deg', 0.023)
Parameter('C1_a_deg', 0.023)
# Observable('TNFRM_obs', TNFRM())

# Expression('tnfrm_rule', tnfr_f_tnfrm/(TNFRM_obs)**2)
# None >> tnfrm
# tnfrm >> None
# 3*tnfrm <> tnfr
# tnfr >> None

# Rule('tnfrm_synth', None >> TNFRM(), synth_tnfrm)
# Rule('tnfrm_deg', TNFRM() >> None, deg_tnfrm)
# Rule('TNFR_3tnfrm', TNFRM() + TNFRM() + TNFRM() <> TNFR(tnf = None), tnfr_f_tnfrm, tnfr_r_tnfrm)

# Parameter('tnfrm_none', 3*1e-05)
# Parameter('none_tnfr', 1e-05)
# Parameter('none_tnfrm', 3*.1)
# Parameter('tnfr_none', .1)
# Parameter('tnfrm_tnf_none', 3*1100)
# Parameter('tnf_tnf_none', 1100)
# Parameter('none_tnfrtnf', 1100)
# Rule('one', TNFRM() >> None, tnfrm_none)
# Rule('two', None >> TNFR(tnf = None), none_tnfr)
# Rule('three', None >> TNFRM(), none_tnfrm)
# Rule('four', TNFR(tnf = None) >> None, tnfr_none)
# Rule('five', TNFRM() >> None, tnfrm_tnf_none)
# Rule('six', TNF(c = None, tnfr = None) >> None, tnf_tnf_none)
# Rule('seven', None >> TNFR(tnf = 1)%TNF(c=None, tnfr = 1), none_tnfrtnf)

# Rule('tnfrm3_to_TNFR', TNFRM() + TNFRM() + TNFRM() >> TNFR(tnf = None), tnfr_f_tnfrm)
# Rule('TNFR_to_tnfrm3', TNFR(tnf = None) >> TNFRM() + TNFRM() + TNFRM(), tnfr_r_tnfrm)

Rule('tnfrm3_to_TNFR', TNFRM()  >> TNFR(tnf = None), tnfr_f_tnfrm)
Rule('TNFR_to_tnfrm3', TNFR(tnf = None) >> TNFRM() , tnfr_r_tnfrm)



Rule('TNFR_deg', TNFR(tnf = None) >> None, deg_TNFR)

# tnfr + ttr <> C1_off[state = i]
# C1_off[state = i] <> C1_off[state = a]
# C1 + A20 >> C1_off[state = i]
# C1 >> tnfr + ttr
# C1_off[state = i] >> None
# C1 >> None

Rule('TNFR_TTR_C1i', TNFR(tnf = None) + TTR() <> C1(tnf = None, state = 'i'), TNFR_TTR_f_C1, TNFR_TTR_r_C1)
Rule('a_C1_i', C1(tnf = None, state = 'i') <> C1(tnf = None, state = 'a'), C1_f_a, C1_r_a)
Rule('C1_i_A20', C1(tnf = None, state = 'a') + A20() >> C1(tnf = None, state = 'i') + A20(), C1_f_A20)
Rule('C1_TNFR_TTR', C1(tnf = None, state = 'a') >> TNFR(tnf = None) + TTR(), C1_f_TNFR_TTR)
Rule('C1i_deg', C1(tnf = None, state = 'i') >> None, C1_i_deg)
Rule('C1a_deg', C1(tnf = None, state = 'a') >> None, C1_a_deg)

#TNF-Dependent Complex 1 Activity Reactions
# def tnf_dependent_to_c1():
Parameter('tnf_deg', 0.0154)
Parameter('tnf_tnfrm_f_TNFRtnf', 1100.0*3) # times 2
Parameter('tnf_TNFR_f_TNFRtnf', 1100.0)
Parameter('tnf_TNFR_r_TNFRtnf', 0.021)
Parameter('deg_TNFRtnf', 0.023)
Parameter('TNFRtnf_TTR_f_C1itnf', 100)
Parameter('TNFRtnf_TTR_r_C1itnf', 0.1)
Parameter('C1itnf_f_C1atnf', 30.0)
Parameter('C1itnf_r_C1atnf', 2.0)

# tnf >> None
# tnf + 3*tnfrm >> tnfr : tnf
# tnfr + tnf <> tnfr : tnf
# tnfr : tnf >> None
# tnfr : tnf + ttr <> C1_off[state = i] : tnf
# C1_off[state = i] : tnf <> C1 : tnf

Rule('deg_tnf', TNF(c = None, tnfr = None) >> None, tnf_deg)

# Expression('tnfrmtnf_rule', tnf_tnfrm_f_TNFRtnf/(TNFRM_obs)**2)
# Rule('tnf_tnfrm_TNFRtnf', TNF(c = None, tnfr = None) + TNFRM() + TNFRM() + TNFRM() >> TNFR(tnf = 1)%TNF(c=None, tnfr = 1), tnf_tnfrm_f_TNFRtnf)
Rule('tnf_tnfrm_TNFRtnf', TNF(c = None, tnfr = None) + TNFRM()  >> TNFR(tnf = 1)%TNF(c=None, tnfr = 1), tnf_tnfrm_f_TNFRtnf)



Rule('tnf_TNFR_TNFRtnf', TNFR(tnf = None) + TNF(c = None, tnfr = None)  <> TNFR(tnf = 1)%TNF(c=None, tnfr = 1), tnf_TNFR_f_TNFRtnf, tnf_TNFR_r_TNFRtnf)

Rule('TNFRtnf_deg', TNFR(tnf = 1)%TNF(c =None, tnfr = 1) >> None, deg_TNFRtnf)

Rule ('TNFRtnf_TTR_C1itnf', TNFR(tnf = 1) % TNF(c = None, tnfr = 1) + TTR() <> C1(tnf = 1, state = 'i') % TNF(c = 1, tnfr =None), TNFRtnf_TTR_f_C1itnf, TNFRtnf_TTR_r_C1itnf)
Rule('C1itnf_C1atnf', C1(tnf = 1, state = 'i') % TNF(c = 1, tnfr = None) <> C1(tnf = 1, state = 'a') % TNF(c = 1, tnfr = None), C1itnf_f_C1atnf, C1itnf_r_C1atnf)



Parameter('C1atnf_f_C1itnf', 1000.0)
Parameter('C1tnf_f_TNFRtnf_TTR', 0.1)
Parameter('deg_C1itnf', 0.023)
Parameter('deg_C1atnf', 0.023)

# C1 : tnf + A20 >> C1_off[state = i] : tnf
# C1 : tnf >> tnfr : tnf + ttr
# C1_off[state = i] : tnf >> None
# C1 : tnf >> None

Rule('C1atnf_C1itnf',C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) + A20() >> C1(tnf = 1, state = 'i')%TNF(c = 1, tnfr = None) + A20(), C1atnf_f_C1itnf)
Rule('C1tnf_TNFRtnf_TTR', C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) >> TNFR(tnf = 1)%TNF(c = None, tnfr = 1) + TTR(), C1tnf_f_TNFRtnf_TTR)
Rule('C1itnf_deg', C1(tnf = 1, state = 'i')%TNF(c = 1, tnfr = None) >> None, deg_C1itnf)
Rule('C1atnf_deg', C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) >> None, deg_C1atnf)


Parameter('C1i_f_tnf', 0.021)
Parameter('C1i_r_tnf', 1100.0)
Parameter('C1a_f_tnf', 0.021)
Parameter('C1a_r_tnf', 1100.0)

# C1_off[state = i] : tnf <> C1_off[state =i] + tnf
# C1 : tnf <> C1 + tnf

Rule('C1itnf', C1(tnf = 1, state = 'i')%TNF(c = 1, tnfr = None) <> C1(tnf = None, state = 'i') + TNF(c = None, tnfr = None), C1i_f_tnf, C1i_r_tnf)
Rule('C1atnf', C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) <> C1(tnf = None, state = 'a') + TNF(c = None, tnfr = None), C1a_f_tnf, C1a_r_tnf)

#IKKK (TAB1/2-TAK1 complex) Activity Reactions
# def ikkk_to_ikk_complex():
Parameter('IKKKi_f_IKKK', 5e-7)
Parameter('IKKKa_r_IKKK', 0.25)
Parameter('IKKKi_f_IKKKaC1', 500.0)
Parameter('IKKKi_f_IKKKaC1tnf', 500.0)

# IKKKi <> IKKKa
# IKKK + C1 >> IKKKa
# IKKK + C1 : tnf >> IKKKa

Rule('IKKKi_IKKKa', IKKK(state = 'i') <> IKKK(state = 'a'), IKKKi_f_IKKK, IKKKa_r_IKKK)
Rule('IKKKi_IKKKaC1', IKKK(state = 'i') + C1(tnf = None, state = 'a') >> IKKK(state = 'a') + C1(tnf = None, state = 'a'), IKKKi_f_IKKKaC1)
Rule('IKKKi_IKKKaC1tnf', IKKK(state = 'i') + C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) >> IKKK(state = 'a') + C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None), IKKKi_f_IKKKaC1tnf)


#IKK Activity Reactions
Parameter('IKKi_f_IKK', 5e-5)
Parameter('IKKa_r_IKK', 0.02)
Parameter('IKKi_f_IKKaIKKK', 520.0)
Parameter('IKKa_f_IKKai', 0.15)
Parameter('IKKai_f_IKKi', 0.02)

# IKKi <> IKKa
# IKKi + IKKKa >> IKKa
# IKKa >> IKKai
# IKKai >> IKKi

Rule('IKKi_IKKa', IKK(state = 'i') <> IKK(state = 'a'), IKKi_f_IKK, IKKa_r_IKK)
Rule('IKKi_IKKaIKKK', IKK(state = 'i') + IKKK(state = 'a') >> IKK(state = 'a') + IKKK(state = 'a'), IKKi_f_IKKaIKKK)
Rule('IKKi_IKKai', IKK(state = 'a') >> IKK(state = 'ai'), IKKa_f_IKKai)
Rule('IKKai_IKKi', IKK(state = 'ai') >> IKK(state = 'i'), IKKai_f_IKKi)



#Dictionary to substitute in species names to match matlab files
# def species_dict():
species_dict = {
    0: 'IkBa',
    1: 'IkBan',
    2: 'IkBaNFkB',
    3: 'IkBaNFkBn',
    4: 'IkBat',
    5: 'IkBb',
    6: 'IkBbn',
    7: 'IkBbNFkB',
    8: 'IkBbNFkBn',
    9: 'IkBbt',
    10: 'IkBe',
    11: 'IkBen',
    12: 'IkBeNFkB',
    13: 'IkBeNFkBn',
    14: 'IkBet',
    15: 'IkBd',
    16: 'IkBdn',
    17: 'IkBdNFkB',
    18: 'IkBdNFkBn',
    19: 'IkBdt',
    20: 'NFkB',
    21: 'NFkBn',
    22: 'IKKK_off',
    23: 'IKKK',
    24: 'IKK_off',
    25: 'IKK',
    26: 'IKK_i',
    27: 'TNF',
    28: 'tnfrm',
    29: 'TNFR',
    30: 'TNFRtnf',
    31: 'C1',
    32: 'C1_off',
    33: 'C1tnf',
    34: 'C1tnf_off',
    35: 'TTR',
    36: 'a20',
    37: 'a20t',
    38: 'SOURCE',
    39: 'SINK'
}

# Observable('TNF_obs', TNF(c = None, tnfr = None))
# Observable('IKKK_obs', IKKK(state = 'a'))
# Observable('IKKKoff_obs', IKKK(state = 'i'))
# # Observable('IkBd_obs', IkBd(nfkb=None, S='C'))
# Observable('IKK_obs', IKK(state = 'a'))
Observable('NFkBn_obs', NFkB(ikb=None, S='N'))
# Observable('NFkB_obs', NFkB(ikb=None, S='C'))
# Observable('A20t_obs', A20t())
# Observable('A20_obs', A20())
Observable('TNFRM_obs', TNFRM())
# Observable('TNFRTNF_obs', TNF(c = None, tnfr = 1)%TNFR(tnf = 1))
Observable('TNFR_obs', TNFR(tnf = None))
# Observable('TTR_obs', TTR())
# Observable('C1_obs', C1(tnf = None, state = 'a'))
# Observable('C1off_obs', C1(tnf = None, state = 'i'))
# Observable('C1tnf_obs', C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None))
# Observable('C1tnfoff_obs', C1(tnf = 1, state = 'i')%TNF(c = 1, tnfr = None))
#
#
#
# Observable('IKKa_obs', IKK(state = 'a'))
# Observable('IKKi_obs', IKK(state = 'i'))
# Observable('IKKai_obs', IKK(state = 'ai'))
# Observable('IKKKa_obs', IKKK(state = 'a'))
# #RNA synthesis by NFkBn and Hill Coefficient
#
# # Observable('IkBa_mRNA_obs', IkBa_mRNA())
# # Observable('IkBa_obs', IkBa(b=ANY, c=ANY, S=ANY))
# #Observables for NFkB
# # def observables():
#
# # Observable('IKK2_obs', IKK2(ikb = None, S='C'))
#
# # Observable('NFkBn_obs', NFkB(ikb=WILD, S='N'))
# Observable('NFkBn_bound', NFkB(ikb=ANY, S='N'))

Observable('IkBa_mRNA_obs', IkBa_mRNA())
Observable('IkBb_mRNA_obs', IkBb_mRNA())
Observable('IkBe_mRNA_obs', IkBe_mRNA())
Observable('IkBd_mRNA_obs', IkBd_mRNA())

# Observable('TNF_obs', TNF(c = None, tnfr = None))
# Observable('A20_obs', A20())

Observable('IkBa_obs', IkBa(nfkb = None, S='C'))
# Observable('IkBan_obs', IkBa(nfkb = None, S='N'))
# Observable('IkBaNFkB_obs', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
# Observable('IkBaNFkBn_obs', IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'))

Observable('IkBb_obs', IkBb(nfkb = None, S='C'))
# Observable('IkBbn_obs', IkBb(nfkb = None, S='N'))
# Observable('IkBbNFkB_obs', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
# Observable('IkBbNFkBn_obs', IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'))

Observable('IkBe_obs', IkBe(nfkb = None, S='C'))
# Observable('IkBen_obs', IkBe(nfkb = None, S='N'))
# Observable('IkBeNFkB_obs', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
# Observable('IkBeNFkBn_obs', IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'))


Observable('IkBd_obs', IkBd(nfkb = None, S='C'))
# Observable('IkBdn_obs', IkBd(nfkb = None, S='N'))
# Observable('IkBdNFkB_obs', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
# Observable('IkBdNFkBn_obs', IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'))

tspan = np.linspace(0, 720, 721)
# print(len(tspan))
# print(tspan)
sim = ScipyOdeSimulator(model, tspan = tspan)
simulation_result = sim.run()
# df = run_sim_with_perturbations(sim, tspan, perturbation_params)
# print(df.shape)
# print(len(df['__s21']))
# print(df['__s21'])

for  j,ode in enumerate(model.odes):
    for i in range(len(model.species)):
       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
    print j,":",ode


print(model.parameters)
print(len(model.reactions))
#
# plt.figure()
# plt.plot(tspan/60, df['__s21'], label = 'NFkBn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# plt.show()


# tspan = np.linspace(0, 720, 721)
# x = odesolve(model,tspan,verbose=True)
# print("Initial Conditions")
# y = model.initial_conditions
# print(np.transpose(y))
# #
# print('TTR conc')
# print(x['TTR_obs'])

plt.figure()
plt.plot(tspan, simulation_result.observables['NFkBn_obs'], label = 'NFkBn_obs')
plt.xlabel("Time (in min)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)


plt.figure()
plt.plot(tspan, simulation_result.observables['TNFRM_obs'], label = 'TNFRM_obs')
plt.xlabel("Time (in min)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.figure()
plt.plot(tspan, simulation_result.observables['TNFR_obs'], label = 'TNFR_obs')
plt.xlabel("Time (in min)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan, simulation_result.observables['IKKa_obs'], label = 'IKK_obs')
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBaNFkB_obs'], label = 'IkBaNFkB_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBaNFkBn_obs'], label = 'IkBaNFkBn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBan_obs'], label = 'IkBan_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan/60, x['IkBb_mRNA_obs'], label = 'IkBb_mRNA_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBb_obs'], label = 'IkBb_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBbNFkB_obs'], label = 'IkBbNFkB_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBbNFkBn_obs'], label = 'IkBbNFkBn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBbn_obs'], label = 'IkBbn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan/60, x['IkBe_mRNA_obs'], label = 'IkBe_mRNA_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBe_obs'], label = 'IkBe_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBeNFkB_obs'], label = 'IkBeNFkB_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBeNFkBn_obs'], label = 'IkBeNFkBn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBen_obs'], label = 'IkBen_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan/60, x['IkBd_mRNA_obs'], label = 'IkBd_mRNA_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBd_obs'], label = 'IkBd_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBdNFkB_obs'], label = 'IkBdNFkB_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBdNFkBn_obs'], label = 'IkBdNFkBn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBdn_obs'], label = 'IkBdn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)


# plt.figure()
# plt.plot(tspan/60, x['NFkB_obs'], label = 'NFkB_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()

# plt.figure()
# plt.plot(tspan/60, x['NFkBn_obs'], label = 'NFkBn_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # # plt.show()
# #
# # plt.figure()
# # plt.plot(tspan/60, x['IKKKoff_obs'], label = 'IKKKoff_obs')
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# # # plt.show()
# #
# plt.figure()
# plt.plot(tspan/60, x['IKKK_obs'], label = 'IKKK_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # # plt.show()
# #
# # plt.figure()
# # plt.plot(tspan/60, x['IKKi_obs'], label = 'IKKi_obs')
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# #
# plt.figure()
# plt.plot(tspan/60, x['IKKa_obs'], label = 'IKKa_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['IKKai_obs'], label = 'IKKai_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['TNF_obs'], label = 'TNF_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['TNFRM_obs'], label = 'TNFRM_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['TNFR_obs'], label = 'TNFR_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# # plt.plot(tspan/60, x['A20t_obs'], label = 'A20t_obs')
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# # # plt.show()
# #
# # plt.figure()
# # plt.plot(tspan/60, x['A20_obs'], label = 'A20_obs')
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# # plt.show()
# #
# plt.figure()
# plt.plot(tspan/60, x['TTR_obs'], label = 'TTR_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['TNFRM_obs'], label = 'TNFRM_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan/60, x['TNFR_obs'], label = 'TNFR_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan/60, x['TNFRTNF_obs'], label = 'TNFRTNF_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan/60, x['C1_obs'], label = 'C1_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['C1off_obs'], label = 'C1off_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['C1tnf_obs'], label = 'C1tnf_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['C1tnfoff_obs'], label = 'C1tnfoff_obs')
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)



# plt.show()



# for obs in ["TNF_obs", "IKKK_obs","IkBd_obs", "IKK_obs","NFkBn_obs", "a20t_obs","a20_obs"]:
# plt.figure(1)
# # plt.plot(time/60, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
# plt.subplot(2,1,1)
# plt.plot(tspan/60, x["TNF_obs"], label = 'TNF')
# plt.subplot(2,1,2)
# plt.plot(tspan/60, x["IKKK_obs"], label = 'IKKKa')
# plt.figure(2)
# plt.subplot(2,1,1)
# plt.plot(tspan/60, x["IkBd_obs"], label = 'IkBd')
# plt.subplot(2,1,2)
# plt.plot(tspan/60, x["IKK_obs"], label = 'IKKa')
# plt.figure(3)
# plt.subplot(2,1,1)
# plt.plot(tspan/60, x["NFkBn_obs"], label = 'NFkBn')
# plt.subplot(2,1,2)
# plt.plot(tspan/60, x["a20t_obs"], label = 'A20t')
# plt.figure(4)
# # plt.subplot(2,2,3)
# plt.plot(tspan/60, x["a20_obs"], label = 'A20')
# # plt.subplot(2,1,2)
# # plt.plot(time/60, x["IkBap_NFkBc_obs"], label = 'IkBap_NFkBc')
# plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xlim(xmin = 0)
# plt.ylim(ymin = 0)
# plt.show()

# #equilibrium phase 1 of model using odesolve
# # time_equil = np.linspace(-8000, 0, 8001)
# # equil = odesolve(model, time_equil, verbose=True)
#
# # plt.figure(1)
# # plt.plot(time_equil/60., equil['NFkBn_free']*1000, label= NFkBn_free.name)
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax =100)
# # plt.legend(loc=0)
#
# # last_conc = [equil[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]]
# # print(last_conc)#setting previous conc species equal to none
#
#phase 2 of model using Solver
# tspan = np.linspace(0, 720, 721)
# x = odesolve(model,tspan,verbose=True)
# #
#updating the species names in the odes

# time_equil = np.linspace(-20004, 0, 2001)
# equil = odesolve(model, time_equil, verbose=True)
#
# # plt.figure(1)
# # plt.plot(time_equil/60., equil['NFkBn_free']*1000, label= NFkBn_free.name)
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax =100)
# # plt.legend(loc=0)
#
# last_conc = [equil[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]]
# print(last_conc)


# phase 2 of model using Solver
# tspan = np.linspace(0, 720, 721)
# x = odesolve(model,tspan,verbose=True)


# for  j,ode in enumerate(model.odes):
#     for i in range(len(model.species)):
#        ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
#     print j,":",ode
# #
# # # for i,ode in enumerate(model.odes):
# # #     print i,":", ode
# # #
# for i,sp in enumerate(model.species):
#     print i,":", sp
# #
# #
# plt.figure()
# plt.plot(tspan/60, x['IKKa_obs'], label = IKKa_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax =0.04)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['IKKi_obs'], label = IKKi_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax =0.04)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['IKKai_obs'], label = IKKai_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax =0.04)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['IKKKa_obs'], label = IKKKa_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax =0.04)
# plt.legend(loc=0)

# plt.show()
# # #
# plt.figure()
# plt.plot(tspan/60, x['NFkBn_free'], label = NFkBn_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # plt.show()
#
# # plt.figure(1)
# # plt.plot(tspan/60, x['NFkBc_free'], label = NFkBc_free.name)
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# # # plt.show()
# #
# #
# # plt.figure(1)
# # plt.plot(tspan/60, x['NFkBn_free'], label = NFkBn_free.name)
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# #
# # plt.figure(1)
# # plt.plot(tspan/60, x['NFkBn_free'], label = NFkBn_free.name)
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# #
# # plt.figure(1)
# # plt.plot(tspan/60, x['NFkBn_bound'], label = NFkBn_bound.name)
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['IKK_obs'], label = IKK.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['A20_obs'], label = A20.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# plt.show()
# #
# # ikk2_vals = []
# nsims = len(tspan) - 1 #simulation time is 720
# # yobs = np.empty((nsims, len(model.observables)))
# rows = len(tspan) # nSims
# cols = len(model.observables) #length of observables
# yobs = np.empty((rows, cols)) #creating empty array of R and C
#
# assert model.observables[1].name == 'IKK2_obs'
# assert str(model.species[1]) == "IKK2(ikb=None, S='C')"
#
# # print(model.parameters['IKK2_0'].value)
# #using imported IKK2 values during phase 2 simulation
# for i in range(nsims):
#     print(i) #printing each simulation iteration
#     last_conc[1] = ikk2[i]
#     solver.tspan = [tspan[i], tspan[i+1]] #simulating in 2 step intervals
#     solver.run(y0 = last_conc) #setting initial concentrations to previous simulation value (equil phase 1)
#     # print(last_conc)
#     if i == 0:
#         yobs[0, :] = solver.yobs_view[0]
#     yobs[i + 1 , :] = solver.yobs_view[1] #taking each simulation index and all obs and keeping the last obs in last time point of simulation
#     last_conc = solver.y[1,:] # running solver taking first row and all simulation species


#plotting NFkBn free
# plt.figure(1)
# plt.plot(tspan/60., yobs[:,0]*1000, label= NFkBn_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0, ymax =100)
# plt.legend(loc=0)

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







# # -*- coding: utf-8 -*-
#
#
# import matplotlib, os, sys
# import numpy as np
# import pandas as pd
# # sys.path.append('/Users/geenaildefonso/Projects/ParticleSwarmOptimization')
# from simplepso.pso import PSO
#
# r = os.system('python -c "import matplotlib.pyplot as plt;plt.figure()"')
# if r != 0:
#     matplotlib.use('Agg')
#     import matplotlib.pyplot as plt
#
#     show = False
# else:
#     import matplotlib.pyplot as plt
#
#     show = True
#
# # from exported_nnrm import model
# from old_hoffman_nfkb import model
# from pysb.integrate import Solver
# from pysb.simulator.scipyode import ScipyOdeSimulator
# show = False
# # you can do this in pycharm by right clicking on folder and
# # do "Mark directory as" and do "Source"
#
#
# # defining observable names from model to be used for calibration
# # data_names and var_names from csv data file for all time points
# obs_names = ['IkBa_obs', 'IkBb_obs', 'IkBe_obs', 'IkBd_obs']
# data_names = ['ikba_avg', 'ikbb_avg', 'ikbe_avg', 'ikbd_avg']
# var_names = ['ikba_var', 'ikbb_var', 'ikbe_var', 'ikbd_var']
#
# # Total starting amounts of proteins in obs_names, for normalizing simulations
# obs_totals = [model.parameters['IkBa_0'].value,
#               model.parameters['IkBb_0'].value,
#               model.parameters['IkBe_0'].value,
#               model.parameters['IkBd_0'].value]
#
# # handle experimental data
# exp_data = pd.read_csv('nfkbn_data_jp.csv', delimiter='\t', )
# # all_points = [0.0, 10., 30., 60., 90., 120., 180., 300., 480., 960., 1440.]
#
# # this is if you want to train to subsets of the points
# # exclude_last_points = all_points[:-4]
# # exp_data = exp_data[exp_data['Time'].isin(exclude_last_points)]
#
# exp_obs_norm = exp_data[data_names].as_matrix().T
# var_norm = exp_data[var_names].as_matrix().T
# std_norm = var_norm ** 0.5
#
# tspan = np.array(exp_data['Time'])
#
# # following 3 lines can be deleted once you see what it was doing below
# tmul = 10
# ntimes = len(tspan)
# old_tspan = np.linspace(tspan[0], tspan[-1], (ntimes - 1) * tmul + 1)[::tmul]
#
# t_to_plot = tspan
# t_exp = tspan
# # before they would not match up! So when you compare two vectors, they were
# # not aligned, thus give crazy errors! Only first and last time point aligned
# print("original sim time {}".format(old_tspan))
# print("algined sim time {}".format(t_to_plot))
# print("exp time {}".format(t_exp))
#
# # deprecated solver
# solver = Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6)
#
# # this is the solver you want to move towards
# # the integrate.solver is being deprecated
#
# # solver = ScipyOdeSimulator(model=model, tspan=tspan, integrator='vode',
# #                            integrator_options={'rtol':1e-6,
# #                                                'atol':1e-6})
#
# rate_params = model.parameters_rules()
# param_values = np.array([p.value for p in model.parameters])
#
# rate_mask = np.array(
#     [i for i, p in enumerate(model.parameters) if p in rate_params])
#
# k_ids = [p.value for p in model.parameters_rules()]
# nominal_values = np.array([p.value for p in model.parameters])
# xnominal = np.log10(nominal_values[rate_mask])
# bounds_radius = 2
#
#
# def display_jp(position, save_name='nfkb_untrained'):
#     Y = np.copy(position)
#     param_values[rate_mask] = 10 ** Y
#     solver.run(param_values)
#
#     obs1 = solver.yobs['IkBa_obs']
#     obs1 = obs1 / model.parameters['IkBa_0'].value
#
#     obs2 = solver.yobs['IkBb_obs']
#     obs2 = obs2 / model.parameters['IkBb_0'].value
#
#     obs3 = solver.yobs['IkBe_obs']
#     obs3 /= model.parameters['IkBe_0'].value
#
#     obs4 = solver.yobs['IkBd_obs']
#     obs4 /= model.parameters['IkBd_0'].value
#
#     # create subplots of each to see what is going on
#
#     plt.figure(figsize=(10, 8))
#     plt.subplot(411)
#     plt.plot(t_to_plot, obs1, 'b-o', label='IkBa_obs')
#     plt.errorbar(t_exp, exp_obs_norm[0], yerr=std_norm[0], ecolor='r',
#                  color='r', elinewidth=0.5, capsize=0, label=data_names[0])
#     plt.legend(loc=0)
#
#     plt.subplot(412)
#     plt.plot(t_to_plot, obs2, 'b-o', label='IkBb_obs')
#     plt.errorbar(t_exp, exp_obs_norm[1], yerr=std_norm[1], ecolor='r',
#                  color='r', elinewidth=0.5, capsize=0, label=data_names[1])
#     plt.legend(loc=0)
#
#     plt.subplot(413)
#     plt.plot(t_to_plot, obs3, 'b-o', label='IkBe_obs')
#     plt.errorbar(t_exp, exp_obs_norm[2], yerr=std_norm[2], ecolor='r',
#                  color='r', elinewidth=0.5, capsize=0, label=data_names[2])
#     plt.legend(loc=0)
#
#     plt.subplot(414)
#     plt.plot(t_to_plot, obs4, 'b-o', label='IkBd_obs')
#     plt.errorbar(t_exp, exp_obs_norm[3], yerr=std_norm[3], ecolor='r',
#                  color='r', elinewidth=0.5, capsize=0, label=data_names[3])
#     plt.legend(loc=0)
#
#     plt.tight_layout()
#     plt.savefig('{}.png'.format(save_name), bbox_inches='tight')
#     if show:
#         plt.show()
#     plt.close()
#     observable = [obs1, obs2, obs3, obs4]
#
#     error = calculate_error(observable)
#
#     print('')
#     for n, name in enumerate(obs_names):
#         print("{} error is {}".format(name, error[n]))
#
#
# def calculate_error(obs_list):
#     error = []
#     for ysim_norm, data_name, var_name in zip(obs_list, data_names, var_names):
#         ydata = exp_data[data_name]
#         yvar = exp_data[var_name]
#         error.append(
#             np.sum((ydata[1:] - ysim_norm[1:]) ** 2 / (2 * yvar[1:])) / len(
#                 ydata[1:]))
#
#     return error
#
#
# def likelihood(position):
#     Y = np.copy(position)
#     param_values[rate_mask] = 10 ** Y
#     solver.run(param_values=param_values)
#     obs1 = solver.yobs['IkBa_obs']
#     obs1 /= model.parameters['IkBa_0'].value
#
#     obs2 = solver.yobs['IkBb_obs']
#     obs2 /= model.parameters['IkBb_0'].value
#
#     obs3 = solver.yobs['IkBe_obs']
#     obs3 /= model.parameters['IkBe_0'].value
#
#     obs4 = solver.yobs['IkBd_obs']
#     obs4 /= model.parameters['IkBd_0'].value
#
#     observable = [obs1, obs2, obs3, obs4]
#     error = calculate_error(observable)
#     error = np.array(error).sum()
#     return error,
#
#
# if __name__ == '__main__':
#
#     display_jp(position=xnominal)
#     # quit()
#
#     out_dir = 'PSO_results'
#     if not os.path.exists(out_dir):
#         os.mkdir(out_dir)
#
#     # create 100 best fits parameter sets
#     for i in range(100):
#         pso_fn = PSO(save_sampled=False, verbose=False, num_proc=16)
#         print("running {}".format(i))
#         pso_fn.set_cost_function(likelihood)
#         pso_fn.set_start_position(xnominal)
#         pso_fn.set_bounds(3)
#         pso_fn.set_speed(-.25, .25)
#         pso_fn.run(100, 500)  # particles in swarm and iterations
#         print("Iteration = {} : best value = {}".format(i, pso_fn.best.fitness.values[0]))
#         csv_out = os.path.join(out_dir, 'best_fit_{}.csv'.format(i))
#         fig_out = os.path.join(out_dir, 'best_fit_{}'.format(i))
#         np.savetxt(csv_out, pso_fn.best)
#         display_jp(pso_fn.best, fig_out)
