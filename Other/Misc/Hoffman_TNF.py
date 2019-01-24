__author__ = 'geena'

from pysb import *
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.kappa import contact_map, set_kappa_path, influence_map
from pysb.tools.render_reactions import run
import pygraphviz as pyg
from matplotlib.backends.backend_pdf import PdfPages

Model ()


xlsx_file = pd.ExcelFile('data.xlsx')
pandas_df = xlsx_file.parse('data')

# print('IKKK_off values')
# print(df['IKKK_off values'])

# from xlrd import open_workbook
#
# book = open_workbook('matlabdata.xlsx')
# sheet = book.sheet_by_index(101)
#
# # read header values into the list
# keys = [sheet.cell(0, col_index).value for col_index in xrange(sheet.ncols)]
#
# dict_list = []
# for row_index in xrange(1, sheet.nrows):
#     d = {keys[col_index]: sheet.cell(row_index, col_index).value
#          for col_index in xrange(sheet.ncols)}
#     dict_list.append(d)
#
# print dict_list
# workbook = xlrd.open_workbook('matlabdata.xlsx')
#
# if sheet.cell()
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

Monomer('NEMO',['brip','btak', 'bikk', 'state'], {'state':['unmod', 'K63ub']})

# Parameter('_0', 9000)
# Parameter('IKK_0', 9000)  # molecules per cell

#
# Initial(TNF(brec=None), TNF_0)
# Initial(TNFR(blig=None, bscaf=None, brip=None, state='unmod'), TNFR_0)
# Initial(TRADD(brec=None, btraf=None, brip=None, state='unmod'), TRADD_0)
# Initial(RIP1(bDD=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'), RIP1_0)
# Initial(TRAF(brip=None, bciap=None, state='unmod'), TRAF_0)
# Initial(cIAP(btraf=None), cIAP_0)
# Initial(A20(brip=None), A20_0)
# Initial(CYLD(brip=None, state='inactive'), CYLD_0)
# Initial(TAK1(brip=None, bmapk=None), TAK1_0)
# Initial(NEMO(brip=None, btak=None, bikk=None, state='unmod'), NEMO_0)
# Initial(LUBAC(brip=None), LUBAC_0)
# Initial(IKK(bnemo=None, bnfkb=None, state='inactive'), IKK_0)


#Declaration of initial conditions
# def initial_ikk1_conditions():
# def initial_conditions():
Parameter('IKKK_0', 0.1) #Inhibitor Kinase Kinase Kinase
Initial(IKKK(state = 'i'), IKKK_0)

# def initial_ikk2_conditions():
Parameter('IKK_0', 0.1) #Inhibitor Kinase Kinase
Initial(IKK(state = 'i'), IKK_0)

# def initial nfkb_conditions():
Parameter('NFkB_0', 0.125) #Nuclear Factor-kappaB
Initial(NFkB(ikb=None, S='N'), NFkB_0)

# def initial TTR_conditions():
Parameter('TTR_0', 8.3e-4) #TRADD-TRAF-RIP
Initial(TTR(), TTR_0)

Parameter('TNF_0', 0.2)
Initial(TNF(c = None, tnfr = None), TNF_0)

Observable('NFkBn_free', NFkB(ikb=None, S='N'))
Observable('IKKK_off_obs', IKKK(state = 'i'))
Observable('IKKK_obs', IKKK(state = 'a'))
Observable('IKK_off_obs', IKK(state = 'i'))
Observable('IKK_obs', IKK(state = 'a'))
Observable('IKK_i_obs', IKK(state = 'ai'))

Observable('NFkBn_obs', NFkB(ikb=None, S='N'))
Observable('NFkB_obs', NFkB(ikb=None, S='C'))

Observable('IkBa_obs', IkBa(nfkb = None, S = 'C'))
Observable('IkBb_obs', IkBb(nfkb = None, S = 'C'))
Observable('IkBe_obs', IkBe(nfkb = None, S = 'C'))
Observable('IkBd_obs', IkBd(nfkb = None, S = 'C'))

Observable('IkBan_obs', IkBa(nfkb = None, S = 'N'))
Observable('IkBbn_obs', IkBb(nfkb = None, S = 'N'))
Observable('IkBen_obs', IkBe(nfkb = None, S = 'N'))
Observable('IkBdn_obs', IkBd(nfkb = None, S = 'N'))

Observable('IkBaNFkB_obs', IkBa(nfkb = 1, S = 'C')% NFkB(ikb=1, S='C'))
Observable('IkBbNFkB_obs', IkBb(nfkb = 1, S = 'C')% NFkB(ikb=1, S='C'))
Observable('IkBeNFkB_obs', IkBe(nfkb = 1, S = 'C')% NFkB(ikb=1, S='C'))
Observable('IkBdNFkB_obs', IkBd(nfkb = 1, S = 'C')% NFkB(ikb=1, S='C'))

Observable('IkBaNFkBn_obs', IkBa(nfkb = 1, S = 'N')% NFkB(ikb=1, S='N'))
Observable('IkBbNFkBn_obs', IkBb(nfkb = 1, S = 'N')% NFkB(ikb=1, S='N'))
Observable('IkBeNFkBn_obs', IkBe(nfkb = 1, S = 'N')% NFkB(ikb=1, S='N'))
Observable('IkBdNFkBn_obs', IkBd(nfkb = 1, S = 'N')% NFkB(ikb=1, S='N'))

Observable('TNFR_obs', TNFR(tnf = None))
Observable('TNFRtnf_obs', TNFR(tnf = 1)%TNF(c = None, tnfr = 1))
# Observable('TNF_obs', IkBd_mRNA())
Observable('tnfrm_obs', TNFRM())

Observable('C1_obs', C1(tnf = None, state ='a'))
Observable('C1_off_obs', C1(tnf = None, state ='i'))
Observable('C1tnf_obs', C1(tnf =1, state ='a')%TNF(c = 1, tnfr = None))
Observable('C1tnf_off_obs', C1(tnf =1, state ='i')%TNF(c = 1, tnfr = None))

Observable('TNF_obs', TNF(c = None, tnfr = None))
Observable('A20_obs', A20())
Observable('A20t_obs', A20t())

Observable('IkBat_obs', IkBa_mRNA())
Observable('IkBbt_obs', IkBb_mRNA())
Observable('IkBet_obs', IkBe_mRNA())
Observable('IkBdt_obs', IkBd_mRNA())

#RNA synthesis by NFkBn and Hill Coefficient

# Observable('IkBa_mRNA_obs', IkBa_mRNA())
# Observable('IkBa_obs', IkBa(b=ANY, c=ANY, S=ANY))
#Observables for NFkB
# def observables():

# Observable('IkBa_obs', IkBa(ikk = None, nfkb = None, S='C'))
# Observable('IkBb_obs', IkBb(ikk = WILD, nfkb = WILD, S=WILD))
# Observable('IkBe_obs', IkBe(ikk = WILD, nfkb = WILD, S=WILD))
# Observable('IkBd_obs', IkBd(ikk = WILD, nfkb = WILD, S=WILD))

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
Parameter('a', 8)
Parameter('b', 0.02)
Parameter('e1', 0.3)
Parameter('d', 0.025)

Expression('a_NFkBn', a*(NFkBn_free)**(hill))
Expression('b_NFkBn', b*(NFkBn_free)**(hill))
Expression('e_NFkBn', e1*(NFkBn_free)**(hill))
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

Rule('an_c', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
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
Parameter('a_f_deg', 0.36)
Parameter('b_f_deg', 0.12)
Parameter('e_f_deg', 0.18)
Parameter('d_f_deg', 0.18)

Parameter('and_c_n', 0.36)
Parameter('bnd_c_n', 0.12)
Parameter('end_c_n', 0.18)
Parameter('dnd_c_n', 0.18)

# IkBa >> None
# IkBb >> None
# IkBe >> None
# IkBd >> None

# IkBa : NFkB >> NFkB
# IkBb : NFkB >> NFkB
# IkBe : NFkB >> NFkB
# IkBd : NFkB >> NFkB

Rule('a_c_deg', IkBa(nfkb=None, S='C') >> None, a_f_deg)
Rule('b_c_deg', IkBb(nfkb=None, S='C') >> None, b_f_deg)
Rule('e_c_deg', IkBe(nfkb=None, S='C') >> None, e_f_deg)
Rule('d_c_deg', IkBd(nfkb=None, S='C') >> None, d_f_deg)

Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), and_c_n)
Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), bnd_c_n)
Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), end_c_n)
Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), dnd_c_n)

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
Expression('A20_synthesis', A20_synth*obs_A20t)

Rule('A20t_synth', None >> A20t(), A20_mRNA)
Rule('A20t_mediated_nfkbn', None >> A20t(), A20t_NFkBn)
Rule('A20t_deg', A20t() >> None, A20_mRNA_c_deg)

Rule('synth_A20', None >> A20(), A20_synthesis)
Rule('deg_A20', A20() >> None, A20_deg)

#IKK Activation Module
#TNF-Independent Complex 1 Activity Reactions
# def tnf_independent_to_c1():
Parameter('synth_tnfrm', 2e-7)
Parameter('deg_tnfrm', 0.0058)
Parameter('tnfr_f_tnfrm', 1e-5)
# Parameter('tnfr_f_tnfrm', 6*1e-5)
Parameter('tnfr_r_tnfrm', 0.1)
Parameter('deg_TNFR', 0.023)
Parameter('TNFR_TTR_f_C1', 100.0)
Parameter('TNFR_TTR_r_C1', 0.1)
Parameter('C1_f_a', 30.0)
Parameter('C1_r_a', 2.0)
Parameter('C1_f_A20', 1000.0)
Parameter('C1_f_TNFR_TTR', 0.1)
Parameter('C1_i_deg', 0.023)
Parameter('C1_a_deg', 0.023)

# None >> tnfrm
# tnfrm >> None
# 3*tnfrm <> tnfr
# tnfr >> None

Rule('tnfrm_synth', None >> TNFRM(), synth_tnfrm)
Rule('tnfrm_deg', TNFRM() >> None, deg_tnfrm)
Rule('TNFR_3tnfrm', TNFRM() <> TNFR(tnf = None), tnfr_f_tnfrm, tnfr_r_tnfrm)
# Rule('TNFR_3tnfrm', TNFRM() + TNFRM() + TNFRM() <> TNFR(tnf = None), tnfr_f_tnfrm, tnfr_r_tnfrm)

# Rule('tnfrm3_to_TNFR', TNFRM() + TNFRM() + TNFRM() >> TNFR(tnf = None), tnfr_f_tnfrm)
# Rule('TNFR_to_tnfrm3', TNFR(tnf = None) >> TNFRM() + TNFRM() + TNFRM(), tnfr_r_tnfrm)
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
Parameter('tnf_tnfrm_f_TNFRtnf', 1100.0)
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
Rule('tnf_tnfrm_TNFRtnf', TNF(c = None, tnfr = None) + TNFRM() >> TNFR(tnf = 1)%TNF(c=None, tnfr = 1), tnf_tnfrm_f_TNFRtnf)
# Rule('tnf_tnfrm_TNFRtnf', TNF(c = None, tnfr = None) + TNFRM() + TNFRM() + TNFRM() >> TNFR(tnf = 1)%TNF(c=None, tnfr = 1), tnf_tnfrm_f_TNFRtnf)
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

bind(RIP1(bDD=ANY, bub2=None, state='K63ub'), 'bub2', NEMO(brip=None), 'brip', [1e-6, 1e-3])


"""TAK1 and IKK recruitment to NEMO-K63ub
1. NEMO(brip=ANY)%cIAP(btraf=ANY) >> NEMO(K63ub)
2. NEMO(K63ub) + TAK1 <> NEMO(K63ub):TAK1
3. NEMO(K63ub) + IKK <> NEMO(K63ub):IKK"""

Rule('NEMO_ubiquitylation', NEMO(brip=ANY, state='unmod') % cIAP(btraf=ANY) >> NEMO(brip=ANY, state='K63ub') % cIAP(btraf=ANY), Parameter('NEMO_K63ub', 1e-2))


bind(NEMO(bikk=None, state='K63ub'), 'bikk', IKK(bnemo=None), 'bnemo', [1e-6, 1e-3])

#Dictionary to substitute in species names to match matlab files
# def species_dict():
species_dict = {
    0: 'IKKK_off',
    1: 'IKK_off',
    2: 'NFkBn',
    3: 'TTR',
    4: 'TNF',
    5: 'SOURCE',
    6: 'IkBat',
    7: 'IkBbt',
    8: 'IkBet',
    9: 'IkBdt',
    10: 'NFkB',
    11: 'A20t',
    12: 'A20',
    13: 'tnfrm',
    14: 'SINK',
    15: 'IKKK',
    16: 'IKK',
    17: 'IkBa',
    18: 'IkBb',
    19: 'IkBe',
    20: 'IkBd',
    21: 'TNFR',
    22: 'TNFRtnf',
    23: 'IKK_i',
    24: 'IkBaNFkB',
    25: 'IkBbNFkB',
    26: 'IkBeNFkB',
    27: 'IkBdNFkB',
    28: 'IkBan',
    29: 'IkBbn',
    30: 'IkBen',
    31: 'IkBdn',
    32: 'C1_off',
    33: 'C1tnf_off',
    34: 'IkBaNFkBn',
    35: 'IkBbNFkBn',
    36: 'IkBeNFkBn',
    37: 'IkBdNFkBn',
    38: 'C1',
    39: 'C1tnf'
}

#phase 2 of model using Solver
tspan = np.linspace(0, 720, 721)
x = odesolve(model,tspan,verbose=True)






