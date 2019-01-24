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
Monomer('flip_L', ['bDED', 'state'], {'state': ['A', 'I']})
Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
Monomer('TAK1', ['brip', 'bmapk'])
Monomer('NEMO', ['brip', 'btak', 'bikk', 'state'], {'state': ['I', 'A']})
Monomer('LUBAC', ['brip'])

Parameter('TNF_0', 1.96e-3) #698 is 30ng/ml of TNF
Parameter('TNFR_0', 0.01)
Parameter('TRADD_0', 8.3e-4)
Parameter('RIP1_0', 0.04) #47000
Parameter('TRAF_0', 8.3e-4) # 8.3e-4
Parameter('cIAP_0', 8.3e-4) #10000 8.3e-4
Parameter('A20_0', 0.0049) #2256
Parameter('CYLD_0', 0.004) #50000
Parameter('FADD_0', 0.0033)
Parameter('flip_L_0', 0.004) # 0.004 # 0.09
Parameter('Lubac_0', 0.003)
Parameter('C8_0', 0.093) #10000 # 0.033 # 0.093
Parameter('RIP3_0', 0.004) #20000
Parameter('NEMO_0', 0.1) # 1000000
Parameter('MLKLa_0', 0.0034) # 100000 #0.0034

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
Initial(flip_L(bDED=None, state = 'I'), flip_L_0)
Initial(LUBAC(brip=None), Lubac_0)
Initial(NEMO(brip = None, btak = None, bikk = None, state = 'I'), NEMO_0)
Initial(C8(bf=None, state='A'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)


#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)
Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 100)
Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 0.75)
Parameter('bind_TNF_TNFR_kf', 800)
Parameter('bind_TNF_TNFR_kr', 0.021)
Parameter('TNF_deg1', 0.09)
Parameter('bind_TNFRANY_TRADD_kf', 100)
Parameter('bind_TNFRANY_TRADD_kr', 0.75)
Parameter('bind_TNFRANY_RIP1unmod_kf', 100)
Parameter('bind_TNFRANY_RIP1unmod_kr', 0.75)
Parameter('bind_RIP1ANY_TRAFunmod_kf', 100)
Parameter('bind_RIP1ANY_TRAFunmod_kr', 0.75)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf', 100)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr', 0.75)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf', 100)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr', 0.75)

Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)

Rule('TNF_deg', TNF(brec = None) >> None, TNF_deg1)

Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <>
     TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)

Rule('bind_TNFRANY_RIP1unmod', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)

# Rule('Complex_II_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod') <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod'), k_ComplexII_1)
Rule('Complex_I_ubiquitylation1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)

# Rule('Complex_II_2', TNFR(brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TNFR(brip=None), k_ComplexII_2)
# Rule('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod', cIAP(btraf=None) + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') <>
#      cIAP(btraf=1) % TRAF(brip=None, bciap=1, bcyld = None, state='unmod'), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)
#
# Rule('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub', CYLD(brip=None) + TRAF(brip=None, bciap=None, bcyld=None, state='unmod') <>
#      CYLD(brip=1) % TRAF(brip=None, bciap=None, bcyld=1, state='unmod'), bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf, bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr)

Rule('Complex_I_ubiquitylation', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)



Parameter('CompI_UB1', 30)
Parameter('CompI_UB2', 30)
Parameter('bind_RIP1K63ubANY_A20_kf', 500)
Parameter('bind_RIP1K63ubANY_A20_kr', 0.75)
Parameter('k_A20_1', 100)
Parameter('bind_RIP1K63ubANY_CYLD_kf', 500)
Parameter('bind_RIP1K63ubANY_CYLD_kr', 0.75)

Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
     >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), CompI_UB2)


Rule('ComplexI_Lubac', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
     bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)


#RIP1 K63ub to be deub by A20 or CYLD

Rule('bind_RIP1K63ubANY_A20', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
     bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)

# Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
#      >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None),
#      k_A20_1)


Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
     >>  TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None), k_A20_1)



Rule('bind_RIP1K63ubANY_CYLD', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
     bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)

Rule('A20_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
     >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
     k_A20_1)

#Initiating Necroptosis

Parameter('bind_FADD_proC8_2_kf', 500)
Parameter('bind_FADD_proC8_2_kr', 0.018)
Parameter('bind_FADDANY_flip_L_kf', 500)
Parameter('bind_FADDANY_flip_L_kr', 0.018)
Parameter('bind_FADDANY_proC8_kf', 100)
Parameter('bind_FADDANY_proC8_kr', 0.018)
Parameter('kc_c8_2', 100)
Parameter('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf', 30)
Parameter('k20', 0.001)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 30)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 0.75)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 30)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf', 30)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr', 0.75)
Parameter('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc', 0.75)

#RIP1 deub and necrosome formation

Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)

#DOES THIS HAPPEN?!?!??! @TODO
# Rule('bind_FADD_flip_L_2', C8(bf=None, state='I') + flip_L(bDED=None) <> C8(bf=1, state='I') % flip_L(bDED=1), bind_FADD_flip_L_2_kf, bind_FADD_flip_L_2_kr)
# Rule('C8_activation1', C8(bf=1, state='I') % flip_L(bDED=1) >> flip_L(bDED=None) + C8(bf=None, state='A'), kc_c8_1)


Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, state='I')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='I'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)


Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='I') + flip_L(bDED=None, state = 'I')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='I') % flip_L(bDED=4, state = 'I'), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)



#@TODO ADDED THESE REACTION
Rule('catalyze_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='I') % flip_L(bDED=4, state = 'I') >>
     TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None, state='A') + flip_L(bDED=None, state = 'A') , bind_FADDANY_flip_L_kf)




Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + RIP3(bRHIM=None, bDD = None, state='unmod')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_proC8_kf, bind_FADDANY_proC8_kr)



# Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') +  <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)

Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod')
     >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + FADD(bDD=None,bDED1 = None, bDED2 = None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)

#@TODO ORIGINAL RULE
# Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)

Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf)

Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), k20)

Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
     <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)

Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod')
     >>  MLKL(bRHIM=None, state='active') + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)

# Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)
#
# Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)
# Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)


Observable('RIP1_obs',RIP1(state='unmod'))
Observable('RIP1k63_obs',RIP1(state='K63ub'))
Observable('RIP1deub_obs',RIP1(state='deub'))
Observable('RIP1po4_obs',RIP1(state='po4'))
Observable('RIP3_obs', RIP3(state='unmod'))
Observable('MLKL_obs', MLKL(bRHIM=None, state='unmod'))
Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('RIP13_obs', RIP1(bDD=ANY, bRHIM=1, state='deub') % RIP3(bRHIM=1, bDD = None, state='unmod'))
Observable('RIP13po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('RIP1deub3po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('TNF_obs', TNF(brec = ANY))
Observable('RIPk63_obs',  RIP1(bscf=None, bub1=None, state='K63ub'))
Observable('TNFR_TRADD', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None))
Observable('CI_k63_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) %
           RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))
Observable('CI', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None, state='K63ub') % TRAF(brip=4, bciap=5, state='unmod') % cIAP(btraf = 5))
Observable('A20_obs', A20(brip = None))
Observable('C8_obs', C8())
Observable('flip_obs', flip_L())


# mlklnum = list(np.linspace(0.0009, 0.09, num = 10))
c8num = list(np.linspace(0.033, 0.1, 10))
tspan = np.linspace(0, 1440, 1441)
sim = ScipyOdeSimulator(model, tspan)
sim_result = sim.run(param_values= {'C8_0': c8num})
sim_df = sim_result.dataframe

# print('tnf')
# print(sim_df)
# # quit()
print(sim_df.loc[0, ['C8_obs', 'MLKLa_obs']])
print(sim_df.loc[9, ['C8_obs', 'MLKLa_obs']])
quit()
print(sim_df.iloc[0]['C8_obs'].iloc[:])
print('bhasfadkj')
print(sim_df.iloc[0]['MLKLa_obs'].iloc[:])
# print('mlkl2')
# print(sim_df[0]['MLKLa_obs'].iloc[:])
# print('mlkl3')
# print(sim_df[6]['MLKLa_obs'].iloc[:])
quit()


for n in range(0,10):
    plt.plot(tspan/60, sim_df.loc[n]['MLKLa_obs'].iloc[:], lw =2)
plt.xlabel('Time in Hr', fontsize=16)
plt.ylabel('MLKL uM Concentration', fontsize=16)
plt.title('TNF:10 ng/ml')
plt.legend(c8num, loc='best')
# plt.legend(loc = 'best')
plt.show()
quit()

sim1 = ScipyOdeSimulator(model, tspan)
sim2 =ScipyOdeSimulator(model, tspan)
sim3 = ScipyOdeSimulator(model, tspan)
sim4 =ScipyOdeSimulator(model, tspan)

L1 = sim1.run(param_values={'TNF_0':1.96e-5})
L2 = sim2.run(param_values={'TNF_0':1.96e-4})
L3 = sim3.run(param_values={'TNF_0':1.96e-3})
L4 = sim4.run(param_values={'TNF_0':1.96e-2})


# print('mlkl')
# print(L1.observables['MLKLa_obs'])
# print(L2.observables['MLKLa_obs'])
# print(L3.observables['MLKLa_obs'])
# print(L4.observables['MLKLa_obs'])
#
#
# # plt.figure()
# # plt.plot(tspan/60, L1.observables['C8_obs'],label = 'C8.1')
# # plt.plot(tspan/60, L2.observables['C8_obs'],label = 'C81')
# # plt.plot(tspan/60, L3.observables['C8_obs'],label = 'C810')
# # plt.plot(tspan/60, L4.observables['C8_obs'],label = 'C8100')
# # plt.xlabel("Time (in hr)", fontsize=15)
# # plt.ylabel("Concentration uM", fontsize=15)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# #
# # plt.figure()
# # plt.plot(tspan/60, L1.observables['flip_obs'],label = 'flip.1')
# # plt.plot(tspan/60, L2.observables['flip_obs'],label = 'flip1')
# # plt.plot(tspan/60, L3.observables['flip_obs'],label = 'flip10')
# # plt.plot(tspan/60, L4.observables['flip_obs'],label = 'flip100')
# # plt.xlabel("Time (in hr)", fontsize=15)
# # plt.ylabel("Concentration uM", fontsize=15)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
#
plt.figure(figsize = (20,7))
# plt.figure()
plt.subplot(221)
plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
plt.plot(tspan/60, L3.observables['TNF_obs'],label = 'TNF10')
plt.plot(tspan/60, L4.observables['TNF_obs'],label = 'TNF100')

# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Concentration uM", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(222)
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
plt.plot(tspan/60, L1.observables['CI_k63_obs'],label = 'CI_k63.1')
plt.plot(tspan/60, L2.observables['CI_k63_obs'],label = 'CI_k631')
plt.plot(tspan/60, L3.observables['CI_k63_obs'],label = 'CI_k6310')
plt.plot(tspan/60, L4.observables['CI_k63_obs'],label = 'CI_k63100')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Concentration uM", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
#
# plt.figure()
plt.subplot(223)
plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
plt.plot(tspan/60, L2.observables['RIP13po4_obs'],label = 'RIP13po41')
plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po410')
plt.plot(tspan/60, L4.observables['RIP13po4_obs'],label = 'RIP13po4100')

# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Concentration uM", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(224)
plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
plt.plot(tspan/60, L3.observables['MLKLa_obs'],label = 'MLKLa10')
plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLa100')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Concentration uM", fontsize=15)
plt.legend(loc = 0)
#
# plt.figure(figsize = (30, 10))
# # plt.figure()
# plt.subplot(141)
# plt.plot(tspan/60, L1.observables['RIP1_obs'],label = 'RIP1.1')
# plt.plot(tspan/60, L2.observables['RIP1_obs'],label = 'RIP11')
# plt.plot(tspan/60, L3.observables['RIP1_obs'],label = 'RIP110')
# plt.plot(tspan/60, L4.observables['RIP1_obs'],label = 'RIP1100')
# # plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymax =9e-4)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(142)
# plt.plot(tspan/60, L1.observables['RIP1k63_obs'],label = 'RIP1k63.1')
# plt.plot(tspan/60, L2.observables['RIP1k63_obs'],label = 'RIP1k631')
# plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k6310')
# plt.plot(tspan/60, L4.observables['RIP1k63_obs'],label = 'RIP1k63100')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(143)
# plt.plot(tspan/60, L1.observables['RIP1deub_obs'],label = 'RIP1deub.1')
# plt.plot(tspan/60, L2.observables['RIP1deub_obs'],label = 'RIP1deub1')
# plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub10')
# plt.plot(tspan/60, L4.observables['RIP1deub_obs'],label = 'RIP1deub100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(144)
# plt.plot(tspan/60, L1.observables['RIP1po4_obs'],label = 'RIP1po4.1')
# plt.plot(tspan/60, L2.observables['RIP1po4_obs'],label = 'RIP1po41')
# plt.plot(tspan/60, L3.observables['RIP1po4_obs'],label = 'RIP1po410')
# plt.plot(tspan/60, L4.observables['RIP1po4_obs'],label = 'RIP1po4100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure(figsize = (23,5))
# # plt.figure()
# plt.subplot(141)
# plt.plot(tspan/60, L1.observables['RIP13_obs'],label = 'RIP13.1')
# plt.plot(tspan/60, L2.observables['RIP13_obs'],label = 'RIP131')
# plt.plot(tspan/60, L3.observables['RIP13_obs'],label = 'RIP1310')
# plt.plot(tspan/60, L4.observables['RIP13_obs'],label = 'RIP13100')
# #
# # plt.plot(tspan/60, L1.observables['RIP1k63_obs'],label = 'RIP1k63.1')
# # plt.plot(tspan/60, L2.observables['RIP1k63_obs'],label = 'RIP1k631')
# # plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k6310')
# # plt.plot(tspan/60, L4.observables['RIP1k63_obs'],label = 'RIP1k63100')
# # plt.plot(tspan/60, L1.observables['RIP1deub_obs'],label = 'RIP1deub.1')
# # plt.plot(tspan/60, L2.observables['RIP1deub_obs'],label = 'RIP1deub1')
# # plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub10')
# # plt.plot(tspan/60, L4.observables['RIP1deub_obs'],label = 'RIP1deub100')
# # plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymax =9e-4)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(142)
# plt.plot(tspan/60, L1.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po4.1')
# plt.plot(tspan/60, L2.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po41')
# plt.plot(tspan/60, L3.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po410')
# plt.plot(tspan/60, L4.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po4100')
# # plt.plot(tspan/60, L1.observables['RIP1_obs'],label = 'RIP1.1')
# # plt.plot(tspan/60, L2.observables['RIP1_obs'],label = 'RIP11')
# # plt.plot(tspan/60, L3.observables['RIP1_obs'],label = 'RIP110')
# # plt.plot(tspan/60, L4.observables['RIP1_obs'],label = 'RIP1100')
# # plt.plot(tspan/60, L1.observables['RIP3_obs'],label = 'RIP3.1')
# # plt.plot(tspan/60, L2.observables['RIP3_obs'],label = 'RIP31')
# # plt.plot(tspan/60, L3.observables['RIP3_obs'],label = 'RIP310')
# # plt.plot(tspan/60, L4.observables['RIP3_obs'],label = 'RIP3100')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(143)
# plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
# plt.plot(tspan/60, L2.observables['RIP13po4_obs'],label = 'RIP13po41')
# plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po410')
# plt.plot(tspan/60, L4.observables['RIP13po4_obs'],label = 'RIP13po4100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(144)
# plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
# plt.plot(tspan/60, L3.observables['MLKLa_obs'],label = 'MLKLa10')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLa100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

plt.tight_layout()
plt.show()