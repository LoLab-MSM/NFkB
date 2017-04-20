from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

Model()

Monomer('TNF', ['brec'])
Monomer('TNFR', ['blig', 'brip', 'bDD'])
Monomer('TRADD', ['brec', 'brip', 'state','bDD1', 'bDD2'], {'state': ['unmod', 'K63ub']})
Monomer('RIP1', ['bscf', 'bub1', 'bub2', 'bub3','bDD', 'btraf', 'bRHIM', 'bMLKL', 'state'], {'state': ['unmod', 'K63ub', 'deub','ub', 'po4', 'trunc']})
Monomer('TRAF', ['brip', 'bciap', 'bcyld', 'state'], {'state': ['unmod', 'K63ub']})
Monomer('cIAP', ['btraf'])
Monomer('A20', ['brip'])
Monomer('A20t')
Monomer('CYLD', ['brip','btraf'])
Monomer('FADD', ['bDD', 'bDED1', 'bDED2'])
Monomer('proC8', ['bDED'])
Monomer('C8', ['bf', 'state'], {'state': ['A', 'I']})
Monomer('flip_L', ['bDED'])
# Monomer('flip_S', ['bDED'])
Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
Monomer('TAK1', ['brip', 'bmapk'])
Monomer('NEMO', ['brip', 'btak', 'bikk'])
Monomer('LUBAC', ['brip'])
Monomer('IKK', ['bind', 'bnemo', 'state'], {'state': ['I', 'A', 'AI']})

Parameter('TNF_0', .869) #698 is 30ng/ml of TNF
Parameter('TNFR_0', 5.8515)
Parameter('TRADD_0', 8.715)
Parameter('RIP1_0', 5.8515) #47000
Parameter('TRAF_0', 5.8515)
Parameter('cIAP_0', 1.245) #10000
Parameter('A20_0', 6.225) #2256
Parameter('CYLD_0', 6.225) #50000
Parameter('FADD_0', 1)
Parameter('flip_L_0', 4.858)
Parameter('Lubac_0', 4.858)
# Parameter('proC8_0', 16057.0)
Parameter('C8_0', 1.245) #10000
Parameter('RIP3_0', 2.49) #20000
# Parameter('MLKL_0', 10000.0) # 1000000
Parameter('MLKLa_0', 6.225) # 100000

Initial(TNF(brec=None), TNF_0)
Initial(TNFR(blig=None, brip=None, bDD = None), TNFR_0)
Initial(TRADD(brec=None, brip=None, state='unmod', bDD1 = None, bDD2 = None), TRADD_0)
Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'), RIP1_0)
Initial(TRAF(brip=None, bciap=None, bcyld = None, state='unmod'), TRAF_0)
Initial(cIAP(btraf=None), cIAP_0)
Initial(A20(brip = None), A20_0)
Initial(CYLD(brip=None, btraf = None), CYLD_0)
Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
Initial(RIP3(bRHIM=None, bDD = None, state='unmod'), RIP3_0)
Initial(flip_L(bDED=None), flip_L_0)
Initial(LUBAC(brip=None), Lubac_0)
# Initial(proC8(bDED=None), proC8_0)
Initial(C8(bf=None, state='A'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)
# Initial(MLKL(bRHIM=None, state='active'), MLKLa_0)


Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf', 1e-06)
Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr', 0.001)
Parameter('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc', 0.1)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kf', 1e-06)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kr', 0.001)
Parameter('catalyze_C8ACYLDU_to_C8A_CYLDT_kc', 0.1)
# Parameter('TNF_0', 23500.0)

Parameter('bind_RIP1K63ubANY_TAK1_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_TAK1_kr', 0.001)
Parameter('bind_RIP1K63ubANY_NEMO_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_NEMO_kr', 0.001)
Parameter('bind_RIP1K63ubANYNEMO_LUBAC_kf', 1e-06)
Parameter('bind_RIP1K63ubANYNEMO_LUBAC_kr', 0.001)
Parameter('k_A20_2', 10.0)
Parameter('k_A20_3', 10.0)
Parameter('k_A20_4', 10.0)
Parameter('k_CYLD_1', 10.0)
Parameter('k_CYLD_2', 10.0)
Parameter('k_CYLD_3', 10.0)
Parameter('k_CYLD_4', 10.0)
Parameter('bind_NEMORIP1TAK1_IKKinactive_kf', 1e-06)
Parameter('bind_NEMORIP1TAK1_IKKinactive_kr', 0.001)
Parameter('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc', 0.1)
Parameter('bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kf', 1e-06)
Parameter('bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kr', 0.001)
Parameter('catalyze_TAK1ANYMAPKinactive_to_TAK1ANY_MAPKactive_kc', 0.1)
Parameter('k_CYLD_deac', 0.01)
Parameter('k_ComplexII_1', 10.0)
Parameter('k_ComplexII_2', 10.0)
Parameter('TNFR1_FADD_kc_2', 0.1)
Parameter('TRADD_FADD_kc_2', 0.1)
Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 0.1)
Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 0.01)
Parameter('bind_RIP1ANYunmod_RIP3unmod_kf', 1e-06)
Parameter('bind_RIP1ANYunmod_RIP3unmod_kr', 0.001)
Parameter('k19', 0.01)
Parameter('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kr', 0.001)
Parameter('catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc_kc', 0.1)
Parameter('bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kf', 1e-06)
Parameter('bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kr', 0.001)
Parameter('catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc_kc', 0.1)
Parameter('bind_FADDproC8proC8_RIP1unmod_kf', 1e-06)
Parameter('bind_FADDproC8proC8_RIP1unmod_kr', 0.001)
Parameter('RIP1_trunc_kc3', 0.1)
Parameter('bind_TRADDFADDproC8proC8_RIP1unmod_kf', 1e-06)
Parameter('bind_TRADDFADDproC8proC8_RIP1unmod_kr', 0.001)
Parameter('RIP1_trunc_kc4', 0.1)
Parameter('bind_C8A_RIP3unmod_to_C8ARIP3unmod_kf', 1e-06)
Parameter('bind_C8A_RIP3unmod_to_C8ARIP3unmod_kr', 0.001)
Parameter('catalyze_C8ARIP3unmod_to_C8A_RIP3trunc_kc', 0.1)
Parameter('bind_FADD_flip_L_2_kf', 7.27e-06)
Parameter('bind_FADD_flip_L_2_kr', 0.018)
Parameter('kc_c8_1', 1.0)
Parameter('bind_TRADDunmodANY_RIP1unmod_kf', 1e-06)
Parameter('bind_TRADDunmodANY_RIP1unmod_kr', 0.001)

Parameter('catalyze_cIAPTRAFunmod_to_cIAP_TRAFK63ub_kc', 0.1)

Parameter('catalyze_CYLDTRAFK63ub_to_CYLD_TRAFunmod_kc', 1.0)


#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)
Parameter('bind_TNF_TNFR_kf', 1.245e-10)
Parameter('bind_TNF_TNFR_kr', 1.245e-7)
Parameter('bind_TNFRANY_TRADD_kf', 1.245e-10)
Parameter('bind_TNFRANY_TRADD_kr', 1.245e-7)
Parameter('bind_TNFRANY_RIP1unmod_kf', 1.245e-10)
Parameter('bind_TNFRANY_RIP1unmod_kr', 1.245e-7)
Parameter('bind_RIP1ANY_TRAFunmod_kf', 1.245e-10)
Parameter('bind_RIP1ANY_TRAFunmod_kr', 1.245e-7)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf', 1.245e-10)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr', 1.245e-7)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf', 1.245e-10)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr', 1.245e-7)

Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)

Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)

Rule('bind_TNFRANY_RIP1unmod', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)

# Rule('Complex_II_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod') <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod'), k_ComplexII_1)
Rule('Complex_I_ubiquitylation1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)

# Rule('Complex_II_2', TNFR(brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TNFR(brip=None), k_ComplexII_2)
Rule('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod', cIAP(btraf=None) + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') <> cIAP(btraf=1) % TRAF(brip=None, bciap=1, bcyld = None, state='unmod'), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)

Rule('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub', CYLD(brip=None) + TRAF(brip=None, bciap=None, bcyld=None, state='unmod') <> CYLD(brip=1) % TRAF(brip=None, bciap=None, bcyld=1, state='unmod'), bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf, bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr)

Rule('Complex_I_ubiquitylation', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)



Parameter('CompI_UB1', 0.001245)
Parameter('CompI_UB2', 0.001245)
Parameter('bind_RIP1K63ubANY_A20_kf', 1.245e-10)
Parameter('bind_RIP1K63ubANY_A20_kr', 1.245e-7)
Parameter('k_A20_1', 0.001245)
Parameter('bind_RIP1K63ubANY_CYLD_kf', 1.245e-10)
Parameter('bind_RIP1K63ubANY_CYLD_kr', 1.245e-7)

Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
     >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), CompI_UB2)


Rule('ComplexI_Lubac', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
     bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)


#RIP1 K63ub to be deub by A20 or CYLD
#
Rule('bind_RIP1K63ubANY_A20', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
     bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)

# Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
#      >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None),
#      k_A20_1)


Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
     >>  TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub'), k_A20_1)



Rule('bind_RIP1K63ubANY_CYLD', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
     bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)

Rule('A20_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
     >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
     k_A20_1)



Parameter('bind_FADD_proC8_2_kf', 9.05115e-10)
Parameter('bind_FADD_proC8_2_kr', 2.241e-6)
Parameter('bind_FADDANY_flip_L_kf', 9.05115e-10)
Parameter('bind_FADDANY_flip_L_kr', 2.241e-6)
Parameter('bind_FADDANY_proC8_kf', 7.27e-06)
Parameter('bind_FADDANY_proC8_kr', 2.241e-6)
Parameter('kc_c8_2', 0.0001245)
Parameter('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf', 1.245e-6)
Parameter('k20', 1.245e-7)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 1.245e-10)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 1.245e-7)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 1.245e-5)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf', 1e-07)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr', 4.98e-5)
Parameter('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc', 1.245e-6)

#RIP1 deub and necrosome formation

Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)

#DOES THIS HAPPEN?!?!??! @TODO
# Rule('bind_FADD_flip_L_2', C8(bf=None, state='I') + flip_L(bDED=None) <> C8(bf=1, state='I') % flip_L(bDED=1), bind_FADD_flip_L_2_kf, bind_FADD_flip_L_2_kr)
# Rule('C8_activation1', C8(bf=1, state='I') % flip_L(bDED=1) >> flip_L(bDED=None) + C8(bf=None, state='A'), kc_c8_1)


Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, state='A')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='A'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)


Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='A') + flip_L(bDED=None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)


TNF
Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) + RIP3(bRHIM=None, bDD = None, state='unmod')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_proC8_kf, bind_FADDANY_proC8_kr)


# Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') +  <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)

Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)

Parameter('TNF_deg1', 0.01245)
Rule('TNF_deg', TNF(brec = None) >> None, TNF_deg1)

Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf)

Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), k20)

Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
     <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)

Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod')
     >>  MLKL(bRHIM=None, state='active') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)

Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)

Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)



Observable('RIP1_obs',RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'))
Observable('RIP3_obs', RIP3(bRHIM=None, bDD = None, state='unmod'))
Observable('MLKL_obs', MLKL(bRHIM=None, state='unmod'))
Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('RIP13_obs', RIP1(bDD=ANY, bRHIM=1, state='unmod') % RIP3(bRHIM=1, bDD = None, state='unmod'))
Observable('RIP13po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('TNF_obs', TNF(brec = ANY))
Observable('RIPk63_obs',  RIP1(bscf=None, bub1=None, state='K63ub'))
Observable('TNFR_TRADD', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None))
Observable('CI_k63_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) %
           RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))
Observable('CI', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None, state='K63ub') % TRAF(brip=4, bciap=5, state='unmod') % cIAP(btraf = 5))



tspan = np.linspace(0, 1440, 1441)
# print(len(tspan))
# print(tspan)
sim = ScipyOdeSimulator(model, tspan = tspan)
simulation_result = sim.run()


for i,sp in enumerate(model.species):
    print i,":",sp

print(len(model.species))

plt.figure(figsize = (15,5))
# plt.figure()
plt.subplot(131)
plt.plot(tspan/60, simulation_result.observables['TNF_obs'], marker = '*',label = 'TNF')
# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(132)
plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],marker = '*',label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(133)
plt.plot(tspan/60, simulation_result.observables['MLKLa_obs'],marker = '*', label = 'MLKLa_obs')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.tight_layout()
plt.show()