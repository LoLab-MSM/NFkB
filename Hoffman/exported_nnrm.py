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
Monomer('IKK', ['bind'])
Monomer('MAPK', ['btak', 'state'], {'state': ['inactive', 'active']})
Monomer('NFkB', ['ikb', 'loc'], {'loc': ['C', 'N']})
Monomer('IkBa', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
Monomer('IkBb', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
Monomer('IkBe', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
Monomer('IkBd', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
Monomer('IkBa_mRNA')
Monomer('IkBb_mRNA')
Monomer('IkBe_mRNA')
Monomer('IkBd_mRNA')


Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf', 1e-06)
Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr', 0.001)
Parameter('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc', 0.1)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kf', 1e-06)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kr', 0.001)
Parameter('catalyze_C8ACYLDU_to_C8A_CYLDT_kc', 0.1)
# Parameter('TNF_0', 23500.0)
Parameter('TNF_0', 200) #698 is 30ng/ml of TNF
Parameter('TNFR_0', 47000.0)
Parameter('TRADD_0', 70000.0)
Parameter('RIP1_0', 47000.0) #47000
Parameter('TRAF_0', 47000.0)
Parameter('cIAP_0', 10000.0)
Parameter('A20_0', 2256.0)
Parameter('CYLD_0', 50000.0) #50000
Parameter('TAK1_0', 9000.0)
Parameter('NEMO_0', 10000.0)
Parameter('LUBAC_0', 10000.0)
Parameter('IKK_0', 75000.0)
Parameter('NFkBn_0', 1000.0) #100,000
Parameter('IkBa_0', 1000.0)
Parameter('IkBb_0', 1000.0)
Parameter('IkBe_0', 1000.0)
Parameter('IkBd_0', 1000.0)
Parameter('MAPK_0', 10000.0)
Parameter('bind_TNF_TNFR_kf', 1e-06)
Parameter('bind_TNF_TNFR_kr', 0.001)
Parameter('bind_TNFRANY_TRADD_kf', 1e-06)
Parameter('bind_TNFRANY_TRADD_kr', 0.001)
Parameter('bind_TNFRANY_RIP1unmod_kf', 1e-06)
Parameter('bind_TNFRANY_RIP1unmod_kr', 0.001)
Parameter('bind_TRADDunmodANY_RIP1unmod_kf', 1e-06)
Parameter('bind_TRADDunmodANY_RIP1unmod_kr', 0.001)
Parameter('bind_RIP1ANY_TRAFunmod_kf', 1e-06)
Parameter('bind_RIP1ANY_TRAFunmod_kr', 0.001)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf', 1e-06)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr', 0.001)
Parameter('catalyze_cIAPTRAFunmod_to_cIAP_TRAFK63ub_kc', 0.1)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf', 1e-06)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr', 0.001)
Parameter('catalyze_CYLDTRAFK63ub_to_CYLD_TRAFunmod_kc', 1.0)
Parameter('CompI_UB1', 10.0)
Parameter('CompI_UB2', 10.0)
Parameter('bind_RIP1K63ubANY_A20_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_A20_kr', 0.001)
Parameter('bind_RIP1K63ubANY_CYLD_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_CYLD_kr', 0.001)
Parameter('bind_RIP1K63ubANY_TAK1_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_TAK1_kr', 0.001)
Parameter('bind_RIP1K63ubANY_NEMO_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_NEMO_kr', 0.001)
Parameter('bind_RIP1K63ubANYNEMO_LUBAC_kf', 1e-06)
Parameter('bind_RIP1K63ubANYNEMO_LUBAC_kr', 0.001)
Parameter('k_A20_1', 10.0)
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
Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 0.1)
Parameter('bind_RIP1ANYunmod_RIP3unmod_kf', 1e-06)
Parameter('bind_RIP1ANYunmod_RIP3unmod_kr', 0.001)
Parameter('k19', 0.01)
Parameter('k20', 0.001)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 1e-06)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 0.001)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 0.1)
Parameter('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf', 1e-06)
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
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf', 1e-07)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr', 0.4)
Parameter('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc', 0.01)
Parameter('bind_FADDANY_proC8_kf', 7.27e-06)
Parameter('bind_FADDANY_proC8_kr', 0.018)
Parameter('bind_FADD_proC8_2_kf', 7.27e-06)
Parameter('bind_FADD_proC8_2_kr', 0.018)
Parameter('bind_FADDANY_flip_L_kf', 7.27e-05)
Parameter('bind_FADDANY_flip_L_kr', 0.018)
Parameter('bind_FADD_flip_L_2_kf', 7.27e-06)
Parameter('bind_FADD_flip_L_2_kr', 0.018)
# Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf', 1e-07)
# Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr', 0.4)
# Parameter('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc', 0.01)
Parameter('kc_c8_1', 1.0)
Parameter('kc_c8_2', 1.0)


Parameter('IkBat_0', 20) #IkBa .002
Initial(IkBa_mRNA(), IkBat_0)
Parameter('IkBbt_0', 33) #IkBa .0033
Initial(IkBb_mRNA(), IkBbt_0)
Parameter('IkBet_0', 3) #IkBa .0003
Initial(IkBe_mRNA(), IkBet_0)
Parameter('IkBdt_0', 1) #IkBa .0001
Initial(IkBd_mRNA(), IkBdt_0)

Initial(TNF(brec=None), TNF_0)
Initial(TNFR(blig=None, brip=None, bDD = None), TNFR_0)
Initial(TRADD(brec=None, brip=None, state='unmod', bDD1 = None, bDD2 = None), TRADD_0)
Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'), RIP1_0)
Initial(TRAF(brip=None, bciap=None, bcyld = None, state='unmod'), TRAF_0)
Initial(cIAP(btraf=None), cIAP_0)
Initial(A20(brip=None), A20_0)
Initial(CYLD(brip=None, btraf = None), CYLD_0)
Initial(TAK1(brip=None, bmapk=None), TAK1_0)
Initial(NEMO(brip=None, btak=None, bikk=None), NEMO_0)
Initial(LUBAC(brip=None), LUBAC_0)
Initial(IKK(bind=None), IKK_0)
Initial(NFkB(ikb=None, loc='N'), NFkBn_0)
Initial(IkBa(ikk=None, nfkb=None, loc='C', state='active'), IkBa_0)
Initial(IkBb(ikk=None, nfkb=None, loc='C', state='active'), IkBb_0)
Initial(IkBe(ikk=None, nfkb=None, loc='C', state='active'), IkBe_0)
Initial(IkBd(ikk=None, nfkb=None, loc='C', state='active'), IkBd_0)
# Initial(MAPK(btak=None, state='inactive'), MAPK_0)

# Parameter('TNFa_0', 600.0)
# Parameter('TNFR1_0', 4800.0)
# Parameter('TRADD_0', 9000.0)
# Parameter('RIP1_0', 12044.0)
# Parameter('TRAF_0', 9000.0)
# Parameter('cIAP_0', 9000.0)
# Parameter('NSC_0', 0.0)
# Parameter('NFkB_0', 50000.0)
# Parameter('CYLD_0', 9000.0)
Parameter('FADD_0', 8030.0)
Parameter('flip_L_0', 39023.0)
# Parameter('flip_S_0', 39023.0)
# Parameter('proC8_0', 16057.0)
Parameter('C8_0', 100) #10000
Parameter('RIP3_0', 20000.0) #20000
# Parameter('MLKL_0', 10000.0) # 1000000
# Parameter('MLKLa_0', 10000.0) # 100000
Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
Initial(RIP3(bRHIM=None, bDD = None, state='unmod'), RIP3_0)
Initial(flip_L(bDED=None), flip_L_0)
# Initial(flip_S(bDED=None), flip_S_0)
# Initial(proC8(bDED=None), proC8_0)
Initial(C8(bf=None, state='I'), C8_0)
# Initial(MLKL(bRHIM=None, state='unmod'), MLKL_0)
# Initial(MLKL(bRHIM=None, state='active'), MLKLa_0)


# TNF + TNFR <> TNF:TNFR
# TNFR + TRADD <> TNFR:TRADD
# TNFR + RIP1 <> TNFR:RIP1
# TNFR:TRADD:RIP1(deub) >> TRADD:RIP1(deub) + TNFR
# TNFR:RIP1(deub) >> TNFR + RIP1(deub)


# Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)
# Rule('bind_TNFRANY_TRADD', TNFR(blig=None, brip=None) + TRADD(brec=None, brip=None) <> TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)
# Rule('bind_TNFRANY_RIP1unmod', TNFR(blig=None, brip=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') <> TNFR(blig=None, brip=1) % RIP1(bscf=1, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)
# Rule('Complex_II_1', TNFR(brip=3) % TRADD(brec=3, brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') >> TRADD(brec=None, brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TNFR(brip=None), k_ComplexII_1)
# Rule('Complex_II_2', TNFR(brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TNFR(brip=None), k_ComplexII_2)


# Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)
# Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec=None, brip=None) <> TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)
# Rule('bind_TNFRANY_RIP1unmod', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') <> TNFR(blig=None, brip=1) % TRADD(brec=1, brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)
# Rule('Complex_II_1', TNFR(brip=3) % TRADD(brec=3, brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') >> TRADD(brec=None, brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'), k_ComplexII_1)
# Rule('Complex_II_2', TNFR(brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TNFR(brip=None), k_ComplexII_2)

#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)

Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)
Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)
Rule('bind_TNFRANY_RIP1unmod', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod') <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)
# Rule('Complex_II_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod') <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod'), k_ComplexII_1)
Rule('Complex_I_ubiquitylation1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod') <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)
# Rule('Complex_II_2', TNFR(brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TNFR(brip=None), k_ComplexII_2)
Rule('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod', cIAP(btraf=None) + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') <> cIAP(btraf=1) % TRAF(brip=None, bciap=1, bcyld = None, state='unmod'), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)
Rule('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub', CYLD(brip=None) + TRAF(brip=None, bciap=None, bcyld=None, state='unmod') <> CYLD(brip=1) % TRAF(brip=None, bciap=None, bcyld=1, state='unmod'), bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf, bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr)
Rule('Complex_I_ubiquitylation', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)
Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None), CompI_UB2)



Rule('bind_RIP1K63ubANY_A20', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='K63ub') + A20(brip=None) <> RIP1(bscf=1, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='K63ub') % A20(brip=1), bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)
Rule('A20_2', RIP1(bscf=1, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='K63ub') % A20(brip=1)  >> RIP1(bscf=None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='deub') + A20(brip=None), k_A20_1)
Rule('bind_RIP1K63ubANY_CYLD', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM= None, state='K63ub') + CYLD(brip=None) <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM= 1, state='K63ub') % CYLD(brip=1), bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)
Rule('A20_1', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM = 1, state='K63ub') % CYLD(brip=1) >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM = None, state='deub') + CYLD(brip=None), k_A20_1)
Rule('bind_RIP1K63ubANY_NEMO', RIP1(bub1 = None, bub2=None, bub3 = None, bDD = None, state='K63ub') + NEMO(brip=None) <> RIP1(bub1 = None, bub2=1, bub3 = None, bDD = None, state='K63ub') % NEMO(brip=1), bind_RIP1K63ubANY_NEMO_kf, bind_RIP1K63ubANY_NEMO_kr)
Rule('bind_RIP1K63ubANY_TAK1', RIP1(bub1 = None, bub2=1, bub3 = None, bDD = None, state='K63ub') % NEMO(brip=1) + TAK1(brip=None) <> RIP1(bub1 = None, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2), bind_RIP1K63ubANY_TAK1_kf, bind_RIP1K63ubANY_TAK1_kr)
Rule('bind_RIP1K63ubANYNEMO_LUBAC', RIP1(bub1 = None, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) + LUBAC(brip=None) <> RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3), bind_RIP1K63ubANYNEMO_LUBAC_kf, bind_RIP1K63ubANYNEMO_LUBAC_kr)
Rule('bind_NEMORIP1TAK1_IKKinactive', RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) + IKK(bind=None) <> RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = 4, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) % IKK(bind=4), bind_NEMORIP1TAK1_IKKinactive_kf, bind_NEMORIP1TAK1_IKKinactive_kr)
# Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = 4, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) % IKK(bind=4) >> RIP1(bub1 = None, bub2=None, bub3 = None, bDD = None, state='K63ub') + NEMO(brip=None) + TAK1(brip=None) + LUBAC(brip=None) + IKK(bind=None), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)
Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = 4, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) % IKK(bind=4) >> IKK(bind=None), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)


Rule('bind_TRADDANYRIP1ANY_FADD', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') + FADD(bDD=None, bDED1 = None, bDED2 = None) <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)
Rule('bind_FADDANY_proC8', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1, bDED1 = None, bDED2 = None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1, bDED1 = 2, bDED2 = None) % TRADD(brec = None, brip = None, bDD1=2, bDD2=None), bind_FADDANY_proC8_kf, bind_FADDANY_proC8_kr)
# Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') +  <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)
Rule('bind_FADD_flip_L_2', C8(bf=None, state='I') + flip_L(bDED=None) <> C8(bf=1, state='I') % flip_L(bDED=1), bind_FADD_flip_L_2_kf, bind_FADD_flip_L_2_kr)
Rule('C8_activation1', C8(bf=1, state='I') % flip_L(bDED=1) >> flip_L(bDED=None) + C8(bf=None, state='A'), kc_c8_1)
Rule('bind_FADD_proC8_2', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1, bDED1 = 2) % TRADD(brec = None, brip = None, bDD1=2, bDD2=None) + C8(bf=None, state='A') <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1, bDED1 = 2) % TRADD(brec = None, brip = None, bDD1=2, bDD2=3) % C8(bf=3, state='A'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)
Rule('bind_FADDANY_flip_L', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1, bDED1 = 2) % TRADD(brec = None, brip = None, bDD1=2, bDD2=3) % C8(bf=3, state='A') + RIP3(bRHIM=None, bDD = None, state='unmod') <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=4, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1, bDED1 = 2) % TRADD(brec = None, brip = None, bDD1=2, bDD2=3) % C8(bf=3, state='A') % RIP3(bRHIM=4, bDD = None, state='unmod'), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)
Rule('C8_activation2', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=1, btraf=4, bMLKL = None, bRHIM = None, state = 'deub') % FADD(bDD=1, bDED1 = 2) % TRADD(brec = None, brip = None, bDD1=2, bDD2=3) % C8(bf=3, state='A') % RIP3(bRHIM=4, bDD = None, state='unmod') >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=4, bDD = None, state='unmod') + FADD(bDD=None, bDED1 = None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + C8(bf=None, state='A') , kc_c8_2)



Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=4, bDD = None, state='unmod') <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=4, bDD = None, state='po4'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=4, bDD = None, state='po4') >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = None, bRHIM = None, state = 'po4')% RIP3(bRHIM=4, bDD = None, state='po4'), k20)
Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = None, bRHIM = None, state = 'po4')% RIP3(bRHIM=4, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod') <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = 1, bRHIM = None, state = 'po4')% RIP3(bRHIM=4, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf, bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)
Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=4, bMLKL = 1, bRHIM = None, state = 'po4')% RIP3(bRHIM=4, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod') >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state='po4') + RIP3(bRHIM=None,bDD = None, state='po4') +  MLKL(bRHIM=None, state='active')+  MLKL(bRHIM=None, state='active'), catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)
# Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)
# Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)


Parameter('kr', 0.01)
Parameter('kf', 0.1)
# Rule('IKK_activate_IkBa', IKK(bind=None) + IkBa(ikk=None, nfkb=None, loc='C', state='active') <> IKK(bind=1) % IkBa(ikk=1, nfkb=None, loc='C', state='active'), kf, kr)
# Rule('IKK_activate_IkBb', IKK(bind=None) + IkBb(ikk=None, nfkb=None, loc='C', state='active') <> IKK(bind=1) % IkBb(ikk=1, nfkb=None, loc='C', state='active'), kf, kr)
# Rule('IKK_activate_IkBe', IKK(bind=None) + IkBe(ikk=None, nfkb=None, loc='C', state='active') <> IKK(bind=1) % IkBe(ikk=1, nfkb=None, loc='C', state='active'), kf, kr)
# Rule('IKK_activate_IkBd', IKK(bind=None) + IkBd(ikk=None, nfkb=None, loc='C', state='active') <> IKK(bind=1) % IkBd(ikk=1, nfkb=None, loc='C', state='active'), kf, kr)
# Rule('IKKIkBa_active', IKK(bind=1) % IkBa(ikk=1, nfkb=None, loc='C', state='active') >> IKK(bind=None, state='active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'), kf)
# Rule('IKKIkBb_active', IKK(bind=1) % IkBa(ikk=1, nfkb=None, loc='C', state='active') >> IKK(bind=None, state='active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'), kf)
# Rule('IKKIkBe_active', IKK(bind=1) % IkBa(ikk=1, nfkb=None, loc='C', state='active') >> IKK(bind=None, state='active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'), kf)
# Rule('IKKIkBd_active', IKK(bind=1) % IkBa(ikk=1, nfkb=None, loc='C', state='active') >> IKK(bind=None, state='active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'), kf)

Rule('IKK_activate_IkBa', IKK(bind=None) + IkBa(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IKK(bind=2) % IkBa(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), kf, kr)
Rule('IKK_activate_IkBb', IKK(bind=None) + IkBb(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IKK(bind=2) % IkBb(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), kf, kr)
Rule('IKK_activate_IkBe', IKK(bind=None) + IkBe(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IKK(bind=2) % IkBe(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), kf, kr)
Rule('IKK_activate_IkBd', IKK(bind=None) + IkBd(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IKK(bind=2) % IkBd(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), kf, kr)
Rule('IKKIkBa_active', IKK(bind=2) % IkBa(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> IKK(bind=None) + IkBa(ikk=None, nfkb=None, loc='C', state='phos') + NFkB(ikb=None, loc='C'), kf)
Rule('IKKIkBb_active', IKK(bind=2) % IkBa(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> IKK(bind=None) + IkBa(ikk=None, nfkb=None, loc='C', state='phos') + NFkB(ikb=None, loc='C'), kf)
Rule('IKKIkBe_active', IKK(bind=2) % IkBa(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> IKK(bind=None) + IkBa(ikk=None, nfkb=None, loc='C', state='phos') + NFkB(ikb=None, loc='C'), kf)
Rule('IKKIkBd_active', IKK(bind=2) % IkBa(ikk=2, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> IKK(bind=None) + IkBa(ikk=None, nfkb=None, loc='C', state='phos') + NFkB(ikb=None, loc='C'), kf)


Parameter('psynth_a', 7e-5)
Parameter('psynth_b', 1e-5)
Parameter('psynth_e', 1e-6)
Parameter('psynth_d', 1e-7)
Rule('a_synth', None >> IkBa_mRNA(), psynth_a)
Rule('b_synth', None >> IkBb_mRNA(), psynth_b)
Rule('e_synth', None >> IkBe_mRNA(), psynth_e)
Rule('d_synth', None >> IkBd_mRNA(), psynth_d)

# Parameter('hill', 3)
# Parameter('a', 8)
# Parameter('b', 0.02)
# Parameter('e', 0.3)
# Parameter('d', 0.025)
#
#
# Observable('NFkBn_free', NFkB(ikb=None, loc='N'))
# Expression('a_NFkBn', a*(NFkBn_free)**(hill))
# Expression('b_NFkBn', b*(NFkBn_free)**(hill))
# Expression('e_NFkBn', e*(NFkBn_free)**(hill))
# Expression('d_NFkBn', d*(NFkBn_free)**(hill))
#
# Rule('an_mRNA', None >> IkBa_mRNA(), a_NFkBn)
# Rule('bn_mRNA', None >> IkBb_mRNA(), b_NFkBn)
# Rule('en_mRNA', None >> IkBe_mRNA(), e_NFkBn)
# Rule('dn_mRNA', None >> IkBd_mRNA(), d_NFkBn)

# IkB mRNA and protein synthesis reactions

Parameter('mRNA_a', 0.035)
Parameter('mRNA_b', 3e-3)
Parameter('mRNA_e', 4e-3)
Parameter('mRNA_d', 2e-3)

Rule('a_mRNA', IkBa_mRNA() >> None, mRNA_a)
Rule('b_mRNA', IkBb_mRNA() >> None, mRNA_b)
Rule('e_mRNA', IkBe_mRNA() >> None, mRNA_e)
Rule('d_mRNA', IkBd_mRNA() >> None, mRNA_d)


Parameter('synth', 0.2448)
# Rule('a_psynth', None >> IkBa(ikk=None, nfkb=None, S='C'), synth)
# Rule('b_psynth', None >> IkBb(ikk=None, nfkb=None, S='C'), synth)
# Rule('e_psynth', None >> IkBe(ikk=None, nfkb=None, S='C'), synth)
# Rule('d_psynth', None >> IkBd(ikk=None, nfkb=None, S='C'), synth)

Rule('a_psynth', IkBa_mRNA() >> IkBa(ikk=None, nfkb=None, loc='C', state='active') + IkBa_mRNA(), synth)
Rule('b_psynth', IkBb_mRNA() >> IkBb(ikk=None, nfkb=None, loc='C', state='active') + IkBb_mRNA(), synth)
Rule('e_psynth', IkBe_mRNA() >> IkBe(ikk=None, nfkb=None, loc='C', state='active') + IkBe_mRNA(), synth)
Rule('d_psynth', IkBd_mRNA() >> IkBd(ikk=None, nfkb=None, loc='C', state='active') + IkBd_mRNA(), synth)


# IkBa + NFkB <> IkBa:NFkB
# IkBb + NFkB <> IkBa:NFkB
# IkBe + NFkB <> IkBa:NFkB
# IkBd + NFkB <> IkBa:NFkB

# IkBan + NFkBn <> IkBan:NFkBn
# IkBbn + NFkBn <> IkBan:NFkBn
# IkBen + NFkBn <> IkBan:NFkBn
# IkBdn + NFkBn <> IkBan:NFkBn

Parameter('IkB_IKKf', 30.0)
Parameter('IkB_IKKr', 6e-05)
Rule('an_adc', IkBa(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBa(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)
Rule('bn_adc', IkBb(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBb(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)
Rule('en_adc', IkBe(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBe(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)
Rule('dn_adc', IkBd(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBd(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)
Rule('an_adn', IkBa(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBa(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)
Rule('bn_adn', IkBb(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBb(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)
Rule('en_adn', IkBe(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBe(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)
Rule('dn_adn', IkBd(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBd(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)


#IkB shuttling
# IkBa <> IkBan
# IkBb <> IkBbn
# IkBe <> IkBen
# IkBd <> IkBdn

Parameter('af', 0.09)
Parameter('bf', 0.009)
Parameter('ef', 0.045)
Parameter('df', 0.045)
Parameter('ancf', 0.012)
Parameter('bncf', 0.012)
Parameter('encf', 0.012)
Parameter('dncf', 0.012)
Rule('a_nc', IkBa(ikk=None, nfkb=None, loc='C', state='active') <> IkBa(ikk=None, nfkb=None, loc='N', state='active'), af, ancf)
Rule('b_nc', IkBb(ikk=None, nfkb=None, loc='C', state='active') <> IkBb(ikk=None, nfkb=None, loc='N', state='active'), bf, bncf)
Rule('e_nc', IkBe(ikk=None, nfkb=None, loc='C', state='active') <> IkBe(ikk=None, nfkb=None, loc='N', state='active'), ef, encf)
Rule('d_nc', IkBd(ikk=None, nfkb=None, loc='C', state='active') <> IkBd(ikk=None, nfkb=None, loc='N', state='active'), ef, dncf)



# IkBa:NFkB <> IkBan:NFkBn
# IkBb:NFkB <> IkBbn:NFkBn
# IkBe:NFkB <> IkBen:NFkBn
# IkBd:NFkB <> IkBdn:NFkBn

# NFkB <> NFkBn

Parameter('anf', 0.276)
Parameter('bnf', 0.0276)
Parameter('enf', 0.138)
Parameter('dnf', 0.276)
Parameter('anr', 0.828)
Parameter('bnr', 0.414)
Parameter('enr', 0.414)
Parameter('dnr', 0.414)
Parameter('nf', 5.4)
Parameter('nr', 0.0048)
Rule('an_nc', IkBa(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IkBa(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), anf, anr)
Rule('bn_nc', IkBb(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IkBb(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), bnf, bnr)
Rule('en_nc', IkBe(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IkBe(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), enf, enr)
Rule('dn_nc', IkBd(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') <> IkBd(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'), dnf, dnr)
Rule('n_nc', NFkB(ikb=None, loc='C') <> NFkB(ikb=None, loc='N'), nf, nr)



# IkBa >> None
# IkBb >> None
# IkBe >> None
# IkBd >> None
# IkBa >> None
# IkBb >> None
# IkBe >> None
# IkBd >> None

Parameter('a_dc', 0.12)
Parameter('b_dc', 0.18)
Parameter('e_dc', 0.18)
Parameter('d_dc', 0.0014)
Parameter('a_dn', 0.12)
Parameter('b_dn', 0.18)
Parameter('e_dn', 0.18)
Parameter('d_dn', 0.0014)
Rule('ad_c', IkBa(ikk=None, nfkb=None, loc='C', state='phos') >> None, a_dc)
Rule('bd_c', IkBb(ikk=None, nfkb=None, loc='C', state='phos') >> None, b_dc)
Rule('ed_c', IkBe(ikk=None, nfkb=None, loc='C', state='phos') >> None, e_dc)
Rule('dd_c', IkBd(ikk=None, nfkb=None, loc='C', state='phos') >> None, d_dc)
Rule('ad_n', IkBa(ikk=None, nfkb=None, loc='N', state='active') >> None, a_dn)
Rule('bd_n', IkBb(ikk=None, nfkb=None, loc='N', state='active') >> None, b_dn)
Rule('ed_n', IkBe(ikk=None, nfkb=None, loc='N', state='active') >> None, e_dn)
Rule('dd_n', IkBd(ikk=None, nfkb=None, loc='N', state='active') >> None, d_dn)


# IkBa : NFkB >> NFkB
# IkBb : NFkB >> NFkB
# IkBe : NFkB >> NFkB
# IkBd : NFkB >> NFkB
# IkBa : NFkB >> NFkB
# IkBb : NFkB >> NFkB
# IkBe : NFkB >> NFkB
# IkBd : NFkB >> NFkB

Parameter('c_bn', 6e-05)
Parameter('n_bn', 6e-05)
Rule('an_c', IkBa(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> NFkB(ikb=None, loc='C'), c_bn)
Rule('bn_c', IkBb(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> NFkB(ikb=None, loc='C'), c_bn)
Rule('en_c', IkBe(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> NFkB(ikb=None, loc='C'), c_bn)
Rule('dn_c', IkBd(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C') >> NFkB(ikb=None, loc='C'), c_bn)
Rule('an_n', IkBa(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N') >> NFkB(ikb=None, loc='N'), n_bn)
Rule('bn_n', IkBb(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N') >> NFkB(ikb=None, loc='N'), n_bn)
Rule('en_n', IkBe(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N') >> NFkB(ikb=None, loc='N'), n_bn)
Rule('dn_n', IkBd(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N') >> NFkB(ikb=None, loc='N'), n_bn)


# TRADD + RIP1 <> TRADD:RIP1
# RIP1 + TRAF <> RIP1:TRAF
# cIAP + TRAF <> cIAP:TRAF
# cIAP:TRAF >> cIAP + TRAF(K63ub)
# CYLD + TRAF(K63ub) <> CYLD:TRAF(K63ub)
# CYLD:TRAF(K63ub) >> CYLD + TRAF(unmod)
# TNFR:RIP1:TRAF:cIAP >> TNFR:RIP1(K63ub):TRAF(K63ub):cIAP
# TRADD:RIP1:TRAF:cIAP >> TRADD(K63ub):RIP1(K63ub):TRAF(K63ub):cIAP


# Rule('bind_TRADDunmodANY_RIP1unmod', TRADD(brec=None, brip=None, state='unmod') + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') <> TRADD(brec=None, brip=1, state='unmod') % RIP1(bscf=1, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'), bind_TRADDunmodANY_RIP1unmod_kf, bind_TRADDunmodANY_RIP1unmod_kr)
# Rule('bind_RIP1ANY_TRAFunmod', RIP1(bscf=None, btraf=None) + TRAF(brip=None, state='unmod') <> RIP1(bscf=None, btraf=1) % TRAF(brip=1, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)
# Rule('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod', cIAP(btraf=None) + TRAF(brip=None, bciap=None, state='unmod') <> cIAP(btraf=1) % TRAF(brip=None, bciap=1, state='unmod'), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)
# Rule('catalyze_cIAPTRAFunmod_to_cIAP_TRAFK63ub', cIAP(btraf=1) % TRAF(brip=None, bciap=1, state='unmod') >> cIAP(btraf=None) + TRAF(brip=None, bciap=None, state='K63ub'), catalyze_cIAPTRAFunmod_to_cIAP_TRAFK63ub_kc)
# Rule('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub', CYLD(brip=None) + TRAF(brip=None, state='K63ub') <> CYLD(brip=1) % TRAF(brip=1, state='K63ub'), bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf, bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr)
# Rule('catalyze_CYLDTRAFK63ub_to_CYLD_TRAFunmod', CYLD(brip=1) % TRAF(brip=1, state='K63ub') >> CYLD(brip=None) + TRAF(brip=None, state='unmod'), catalyze_CYLDTRAFK63ub_to_CYLD_TRAFunmod_kc)
# Rule('Complex_I_ubiquitylation1', TNFR(blig=None, brip=2) % RIP1(bscf=2, btraf=3, state='unmod') % TRAF(brip=3, bciap=4, state='unmod') % cIAP(btraf=4) >> TNFR(blig=None, brip=2) % RIP1(bscf=2, btraf=3, state='K63ub') % TRAF(brip=3, bciap=4, state='K63ub') % cIAP(btraf=4), CompI_UB1)
# Rule('Complex_I_ubiquitylation2', TRADD(brec=None, brip=1, state='unmod') % RIP1(bscf=1, btraf=2, state='unmod') % TRAF(brip=2, bciap=3, state='unmod') % cIAP(btraf=3) >> TRADD(brec=None, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=3, state='K63ub') % cIAP(btraf=3), CompI_UB2)


# Rule('bind_TRADDunmodANY_RIP1unmod', TRADD(brec=None, brip=None, state='unmod') + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') <> TRADD(brec=None, brip=1, state='unmod') % RIP1(bscf=1, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'), bind_TRADDunmodANY_RIP1unmod_kf, bind_TRADDunmodANY_RIP1unmod_kr)
# Rule('bind_RIP1ANY_TRAFunmod', RIP1(bscf=None, btraf=None) + TRAF(brip=None, state='unmod') <> RIP1(bscf=None, btraf=1) % TRAF(brip=1, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)

# Rule('catalyze_cIAPTRAFunmod_to_cIAP_TRAFK63ub', cIAP(btraf=1) % TRAF(brip=None, bciap=1, state='unmod') >> cIAP(btraf=None) + TRAF(brip=None, bciap=None, state='K63ub'), catalyze_cIAPTRAFunmod_to_cIAP_TRAFK63ub_kc)

# Rule('catalyze_CYLDTRAFK63ub_to_CYLD_TRAFunmod', CYLD(brip=1) % TRAF(brip=1, state='K63ub') >> CYLD(brip=None) + TRAF(brip=None, state='unmod'), catalyze_CYLDTRAFK63ub_to_CYLD_TRAFunmod_kc)








# RIP1(K63ub) + A20 <> RIP1(K63ub):A20
# RIP1(K63ub) + CYLD <> RIP1(K63ub):CYLD
# RIP1(K63ub) + TAK1 <> RIP1(K63ub):TAK1
# RIP1(K63ub) + NEMO <> RIP1(K63ub):NEMO
# RIP1(K63ub) + LUBAC <> RIP1(K63ub):LUBAC




# RIP1(K63ub):A20:CYLD:TRAF(K63ub) >> A20 + CYLD + TRAF(unmod)
# RIP1(K63ub):A20:NEMO:TRAF(K63ub) >> A20 + NEMO + TRAF(unmod)
# RIP1(K63ub):A20:TRAF(K63ub) >> A20 + TRAF(unmod)
# RIP1(K63ub):TAK1:CYLD:TRAF(K63ub) >> RIP(deub) + CYLD + TRAF(unmod)
# RIP1(K63ub):A20:CYLD:TRAF(K63ub) >> RIP(deub) + A20 + TRAF(unmod)
# RIP1(K63ub):TAK1:CYLD:TRAF(K63ub) >> RIP(deub) + CYLD + TRAF(unmod)
# RIP1(K63ub):CYLD:TRAF(K63ub) >> RIP(deub) + CYLD + TRAF(unmod)



#
# Rule('A20_2', TRADD(brec=4, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=3, state='K63ub')% A20(brip=3) % NEMO(brip=4) >> A20(brip=None) + NEMO(brip=None) + TRAF(brip=None, bciap = None, state='unmod') + TRADD(brec=None, brip=None, state='unmod') + RIP1(bscf=None, btraf=None, state='deub'), k_A20_2)
# Rule('A20_3', TRADD(brec=4, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=None, state='K63ub')% A20(brip=4)  >> A20(brip=None) + TRAF(brip=None, bciap = None, state='unmod') + TRADD(brec=None, brip=None, state='unmod') + RIP1(bscf=None, btraf=None, state='deub'), k_A20_3)
# # Rule('A20_4', RIP1(bscf=None, btraf=2, bub1=1, bub2=None, bub3=None, state='K63ub') % A20(brip=1) % TRAF(brip=2, state='K63ub') >> A20(brip=None) + TRAF(brip=None, state='unmod'), k_A20_4)
# Rule('CYLD_1', TRADD(brec=4, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=3, state='K63ub') % TAK1(brip=3) % CYLD(brip=4)  >> TAK1(brip=None) + CYLD(brip=None) + TRAF(brip=None, bciap = None, state='unmod') + TRADD(brec=None, brip=None, state='unmod') + RIP1(bscf=None, btraf=None, state='deub'), k_CYLD_1)
# # Rule('CYLD_2', TRADD(brec=4, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=3, state='K63ub') % A20(brip=3) % CYLD(brip=4) >> A20(brip=None) + TRAF(brip=None, state='unmod'), k_CYLD_2)
# Rule('CYLD_3', TRADD(brec=4, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=None, state='K63ub') % CYLD(brip=4) >> TRAF(brip=None, bciap = None, state='unmod') + TRADD(brec=None, brip=None, state='unmod') + RIP1(bscf=None, btraf=None, state='deub') + CYLD(brip=None), k_CYLD_3)
# Rule('CYLD_4', RIP1(bscf=None, btraf=3, bub1=None, bub2=1, bub3=None, state='K63ub') % TRAF(brip=3, state='K63ub') % CYLD(brip=1) >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TRAF(brip=None, state='unmod') + CYLD(brip=None), k_CYLD_4)


# Rule('A20_1', RIP1(bscf=None, btraf=3, bub1=1, bub2=2, bub3=None, state='K63ub') % A20(brip=1) % CYLD(brip=2) % TRAF(brip=3, state='K63ub') >> A20(brip=None) + CYLD(brip=None) + TRAF(brip=None, state='unmod'), k_A20_1)
# Rule('A20_2', RIP1(bscf=None, btraf=3, bub1=1, bub2=2, bub3=None, state='K63ub') % A20(brip=1) % NEMO(brip=2) % TRAF(brip=3, state='K63ub') >> A20(brip=None) + NEMO(brip=None) + TRAF(brip=None, state='unmod'), k_A20_2)
# Rule('A20_3', RIP1(bscf=None, btraf=2, bub1=1, bub2=None, bub3=None, state='K63ub') % A20(brip=1) % TRAF(brip=2, state='K63ub') >> A20(brip=None) + TRAF(brip=None, state='unmod'), k_A20_3)
# # Rule('A20_4', RIP1(bscf=None, btraf=2, bub1=1, bub2=None, bub3=None, state='K63ub') % A20(brip=1) % TRAF(brip=2, state='K63ub') >> A20(brip=None) + TRAF(brip=None, state='unmod'), k_A20_4)
# Rule('CYLD_1', RIP1(bscf=None, btraf=3, bub1=1, bub2=2, bub3=None, state='K63ub') % TAK1(brip=1) % CYLD(brip=2) % TRAF(brip=3, state='K63ub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TAK1(brip=None) + CYLD(brip=None) + TRAF(brip=None, state='unmod'), k_CYLD_1)
# Rule('CYLD_2', RIP1(bscf=None, btraf=3, bub1=1, bub2=2, bub3=None, state='K63ub') % A20(brip=1) % CYLD(brip=2) % TRAF(brip=3, state='K63ub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + A20(brip=None) + TRAF(brip=None, state='unmod'), k_CYLD_2)
# Rule('CYLD_3', RIP1(bscf=None, btraf=3, bub1=None, bub2=1, bub3=None, state='K63ub') % TRAF(brip=3, state='K63ub') % CYLD(brip=1) >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TRAF(brip=None, state='unmod') + CYLD(brip=None), k_CYLD_3)
# # Rule('CYLD_4', RIP1(bscf=None, btraf=3, bub1=None, bub2=1, bub3=None, state='K63ub') % TRAF(brip=3, state='K63ub') % CYLD(brip=1) >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TRAF(brip=None, state='unmod') + CYLD(brip=None), k_CYLD_4)



# NEMO:RIP1:TAK1 + IKK(inactive) <> NEMO:RIP1:TAK1:IKK(inactive)
# NEMO:RIP1:TAK1:IKK(inactive) >> NEMO:RIP1:TAK1 + IKK(active)
# TAK1 + MAPK(inactive) <> TAK1:MAPK(inactive)
# TAK1:MAPK(inactive) >> TAK1 + MAPK(active)
# RIP1(K63ub):TAK1:CYLD >> RIP1(K63ub):TAK1



# Rule('bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive', TAK1(brip=None, bmapk=None) + MAPK(btak=None, state='inactive') <> TAK1(brip=None, bmapk=1) % MAPK(btak=1, state='inactive'), bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kf, bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kr)
# Rule('catalyze_TAK1ANYMAPKinactive_to_TAK1ANY_MAPKactive', TAK1(brip=None, bmapk=1) % MAPK(btak=1, state='inactive') >> TAK1(brip=None, bmapk=None) + MAPK(btak=None, state='active'), catalyze_TAK1ANYMAPKinactive_to_TAK1ANY_MAPKactive_kc)
# Rule('CYLD_deactivation', RIP1(bscf=None, bub1=1, bub2=2, state='K63ub') % TAK1(brip=1) % CYLD(brip=2) >> RIP1(bscf=None, bub1=1, bub2=None, state='K63ub') % TAK1(brip=1), k_CYLD_deac)








# TNFR:TRADD:RIP1 + FADD >> TNFR:TRADD + RIP1:FADD
# TRADD:RIP1 + FADD >> TRADD + RIP1:FADD
# TRADD:RIP1 + FADD <> TRADD:RIP1:FADD
# RIP1(unmod) + RIP3(unmod) <> RIP1(unmod):RIP3(unmod)
# RIP1(unmod):RIP3(unmod) >> RIP1(unmod):RIP3(po4)
# RIP1(unmod):RIP3(po4) >> RIP1(po4):RIP3(po4)

# Rule('TNFR1_to_FADD_2', TNFR(bDD=2) % TRADD(bDD1=2, bDD2=3) % RIP1(bDD=3, btraf=None) + FADD(bDD=None) >> TNFR(bDD=2) % TRADD(bDD1=2, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1), TNFR1_FADD_kc_2)
# Rule('TRADD_to_FADD_2', TRADD(bDD1=None, bDD2=2) % RIP1(bDD=2, btraf=None) + FADD(bDD=None) >> TRADD(bDD1=None, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1), TRADD_FADD_kc_2)






# Rule('bind_RIP1ANYunmod_RIP3unmod', RIP1(bDD=None, bRHIM=None, state='unmod') + RIP3(bRHIM=None, bDD = None, state='unmod') <> RIP1(bDD=None, bRHIM=1, state='unmod') % RIP3(bRHIM=1, bDD = None, state='unmod'), bind_RIP1ANYunmod_RIP3unmod_kf, bind_RIP1ANYunmod_RIP3unmod_kr)
# Rule('Rip3_PO4lation', RIP1(bRHIM=4, state='deub') % RIP3(bRHIM=4, bDD = None,state='unmod') >> RIP1(bRHIM=4, state='unmod') % RIP3(bRHIM=4,bDD = None, state='po4'), k19)
# Rule('bind_RIP1ANYunmod_RIP3unmod', RIP1(bDD=None, bRHIM=None, state='deub') + RIP3(bRHIM=None, bDD = None, state='unmod') <> RIP1(bDD=None, bRHIM=1, state='deub') % RIP3(bRHIM=1, bDD = None, state='unmod'), bind_RIP1ANYunmod_RIP3unmod_kf, bind_RIP1ANYunmod_RIP3unmod_kr)
# Rule('Rip3_PO4lation', RIP1(bRHIM=1, state='deub') % RIP3(bRHIM=3, bDD = None,state='unmod') >> RIP1(bRHIM=1, state='unmod') % RIP3(bRHIM=3,bDD = None, state='po4'), k19)


# Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='unmod') <> RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf, bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)

# Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod') >> RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='active'), catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)






# MLKL(active) + MLKL(unmod) <> MLKL(active):MLKL(unmod)
# MLKL(active):MLKL(unmod) >> MLKL(active) + MLKL(active)
# RIP1(po4):RIP3(po4) + MLKL(unmod) <> RIP1(po4):RIP3(po4):MLKL(unmod)
# RIP1(po4):RIP3(po4):MLKL(unmod) >> RIP1(po4) + RIP3(po4) + MLKL(active)
# FADD:FLIP:pC8 + RIP1 <> FADD:FLIP:pC8:RIP1
# FADD:FLIP:pC8:RIP1 >> FADD:FLIP:pC8 + RIP1(trunc)



# Rule('catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc', FADD(bDD=50, bDED1=1, bDED2=2) % flip_L(bDED=1) % proC8(bDED=2) % RIP1(bDD=50, state='unmod') >> FADD(bDD=None, bDED1=1, bDED2=2) % flip_L(bDED=1) % proC8(bDED=2) + RIP1(bDD=None, state='trunc'), catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)



# TRADD:FADD:FLIP:pC8 + RIP1(unmod) <> TRADD:FADD:FLIP:pC8:RIP1(unmod)
# TRADD:FADD:FLIP:pC8:RIP1(unmod) >> TRADD:FADD:FLIP:pC8 + RIP1(trunc)
# FADD:pC8:pC8 + RIP1(unmod) <> FADD:pC8:pC8:RIP1(unmod)
# FADD:pC8:pC8:RIP1(unmod) >> FADD + C8(active) + RIP1(trunc)
# TRADD:FADD:pC8:pC8 + RIP1(unmod) <>  TRADD:FADD:pC8:pC8:RIP1(unmod)
# TRADD:FADD:pC8:pC8:RIP1(unmod) >> TRADD:FADD + C8(active) + RIP1(trunc)
# C8 + RIP1 <> C8:RIP1
# C8:RIP1 >> C8 + RIP1(trunc)


# Rule('bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=1, bDED2=3) % flip_L(bDED=1) % proC8(bDED=3) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=1, bDED2=3) % flip_L(bDED=1) % proC8(bDED=3) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
# Rule('catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc', TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=1, bDED2=3) % flip_L(bDED=1) % proC8(bDED=3) % RIP1(bDD=50, state='unmod') >> TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=1, bDED2=3) % flip_L(bDED=1) % proC8(bDED=3) + RIP1(bDD=None, state='trunc'), catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)
# Rule('bind_FADDproC8proC8_RIP1unmod', FADD(bDD=None, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) + RIP1(bDD=None, state='unmod') <> FADD(bDD=50, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) % RIP1(bDD=50, state='unmod'), bind_FADDproC8proC8_RIP1unmod_kf, bind_FADDproC8proC8_RIP1unmod_kr)
# Rule('RIP1_trunc_3b', RIP1(bDD=1, state='unmod') % FADD(bDD=1, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) >> RIP1(bDD=None, state='trunc') + FADD(bDD=None, bDED1=None, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc3)
# Rule('bind_TRADDFADDproC8proC8_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDproC8proC8_RIP1unmod_kf, bind_TRADDFADDproC8proC8_RIP1unmod_kr)
# Rule('RIP1_trunc_4b', RIP1(bDD=1, state='unmod') % TRADD(bDD1=2, bDD2=3) % FADD(bDD=1, bDED1=2, bDED2=4) % proC8(bDED=4) % proC8(bDED=3) >> RIP1(bDD=None, state='trunc') + TRADD(bDD1=2, bDD2=None) % FADD(bDD=None, bDED1=2, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc4)
# Rule('bind_C8A_RIP1unmod_to_C8ARIP1unmod', C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='unmod') <> C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod'), bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf, bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr)
# Rule('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc', C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod') >> C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='trunc'), catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc)




# Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf', 1e-06)
# Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr', 0.001)
# Parameter('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc', 0.1)
# Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kf', 1e-06)
# Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kr', 0.001)
# Parameter('catalyze_C8ACYLDU_to_C8A_CYLDT_kc', 0.1)


# FADD + pC8 <> FADD:pC8
# FADD + pC8 <> FADD:pC8
# FADD + FLIP <> FADD:FLIP
# FADD + FLIP <> FADD:FLIP ( need?!)
# FADD:pC8:pC8 >> FADD + C8(active)
# TRADD:FADD:pC8:pC8 >> TRADD:FADD + C8(active)
# C8(active) + RIP3(unmod) <> C8(active):RIP3(unmod)
# C8(active):RIP3(unmod) >> C8(active) + RIP3(trunc)
# C8(active) + CYLD(U) <> C8(active):CYLD(U)
# C8(active):CYLD(U) >> C8(active) + CYLD(T)




# Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)
# Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)
# Rule('C8_activation1', FADD(bDD=None, bDED1=1, bDED2=2) % proC8(bDED=1) % proC8(bDED=2) >> FADD(bDD=None, bDED1=None, bDED2=None) + C8(bf=None, state='A'), kc_c8_1)
# Rule('C8_activation2', TRADD(bDD1=None, bDD2=1) % FADD(bDD=1, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) >> TRADD(bDD1=None, bDD2=1) % FADD(bDD=1, bDED1=None, bDED2=None) + C8(bf=None, state='A'), kc_c8_2)
# Rule('bind_C8A_RIP3unmod_to_C8ARIP3unmod', C8(bf=None, state='A') + RIP3(bRHIM=None, bDD = None,state='unmod') <> C8(bf=1, state='A') % RIP3(bRHIM=1,bDD = None, state='unmod'), bind_C8A_RIP3unmod_to_C8ARIP3unmod_kf, bind_C8A_RIP3unmod_to_C8ARIP3unmod_kr)
# Rule('catalyze_C8ARIP3unmod_to_C8A_RIP3trunc', C8(bf=1, state='A') % RIP3(bRHIM=1, bDD = None,state='unmod') >> C8(bf=None, state='A') + RIP3(bRHIM=None,bDD = None, state='trunc'), catalyze_C8ARIP3unmod_to_C8A_RIP3trunc_kc)
# # Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') + CYLD(btraf=None, state='U') <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)
# Rule('catalyze_C8ACYLDU_to_C8A_CYLDT', C8(bf=1, state='A') % CYLD(btraf=1, state='U') >> C8(bf=None, state='A') + CYLD(btraf=None, state='T'), catalyze_C8ACYLDU_to_C8A_CYLDT_kc)





# Defining Observables
Observable('TNFR_TNF', TNFR(blig=ANY))
Observable('TRADD_CompI', TRADD(brec=ANY))
Observable('RIP1_CompI', TNFR() % RIP1())
Observable('RIP1ub', RIP1(state='K63ub'))
# Observable('TRAF_CompI', TRAF(brip=ANY) % RIP1(btraf=ANY))
# Observable('Obs_TRAF_ComplexI', TRAF(brip=ANY) % RIP1(btraf=ANY))
Observable('Ubiquitylated_ComplexI', RIP1(state='K63ub'))
Observable('MAPK_activity', MAPK(state='active'))
Observable('RIP1deub', RIP1(state='deub'))
Observable('Obs_ComplexII', TRADD(brec=None, brip=1) % RIP1(bscf=1, btraf=None, bub1=None, bub2=None, bub3=None, state='deub'))
Observable('TNFR_endocytosis', TNF(brec=3) % TNFR(blig=3) % RIP1(btraf=None, bub1=None, bub2=None, bub3=None, state='deub'))
Observable('CYLD_CompI_1', RIP1(bscf=ANY, bub2=2, bub3=None, state='K63ub') % CYLD(brip=2))
Observable('CYLD_CompI_2', RIP1(bscf=ANY, btraf=3, bub2=2, bub3=None, state='K63ub') % CYLD(brip=2) % TRAF(brip=3, state='K63ub'))
Observable('NFkBn_obs', NFkB(ikb=ANY, loc='N'))
Observable('NFkB_obs', NFkB(ikb=ANY, loc='C'))
Observable('IkBa_obs', IkBa(ikk=None, nfkb=None, loc='C', state='active'))
Observable('IkBb_obs', IkBb(ikk=None, nfkb=None, loc='C', state='active'))
Observable('IkBe_obs', IkBe(ikk=None, nfkb=None, loc='C', state='active'))
Observable('IkBd_obs', IkBd(ikk=None, nfkb=None, loc='C', state='active'))
Observable('IkBan_obs', IkBa(ikk=None, nfkb=None, loc='N', state='active'))
Observable('IkBbn_obs', IkBb(ikk=None, nfkb=None, loc='N', state='active'))
Observable('IkBen_obs', IkBe(ikk=None, nfkb=None, loc='N', state='active'))
Observable('IkBdn_obs', IkBd(ikk=None, nfkb=None, loc='N', state='active'))
Observable('IkBaNFkB_obs', IkBa(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'))
Observable('IkBbNFkB_obs', IkBb(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'))
Observable('IkBeNFkB_obs', IkBe(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'))
Observable('IkBdNFkB_obs', IkBd(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C'))
Observable('IkBaNFkBn_obs', IkBa(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'))
Observable('IkBbNFkBn_obs', IkBb(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'))
Observable('IkBeNFkBn_obs', IkBe(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'))
Observable('IkBdNFkBn_obs', IkBd(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N'))
Observable('RIP1_obs',RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'))
Observable('RIP3_obs', RIP3(bRHIM=None, bDD = None, state='unmod'))
Observable('MLKL_obs', MLKL(bRHIM=None, state='unmod'))
Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('RIP13_obs', RIP1(bDD=ANY, bRHIM=1, state='unmod') % RIP3(bRHIM=1, bDD = None, state='unmod'))
Observable('RIP13po4_obs', RIP1(bDD=None, btraf=4, state = 'po4') % RIP3(bRHIM=4, bDD = None, state='po4'))
Observable('TNF_obs', TNF(brec = ANY))
Observable('RIPk63_obs',  RIP1(bscf=None, bub1=None, state='K63ub'))
Observable('TNFR_TRADD', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None))
Observable('CI_k63_obs', TRADD(brec=None, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=3, state='K63ub') % cIAP(btraf=3))
Observable('IkBamrna', IkBa_mRNA())
Observable('IKK_obs', IKK(bind=None))
# Observable('IKKany', IKK(bind = ANY))
Observable('CI', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None, state='K63ub') % TRAF(brip=4, bciap=5, state='unmod') % cIAP(btraf = 5))


# generate_network(model)
# generate_equations(model)

# print(model.initials)

# tspan = np.linspace(0, 1440, 1441)
# x = odesolve(model,tspan,verbose=True)
#
# print("printing species length")
# print(len(model.species))
#
# print("printing reactions length")
# print(len(model.reactions))
#
# last_conc = [x[y][-1] for y in ['__s%d' %i for i in np.arange(len(model.species))]]
# print(last_conc)
#
# first_conc = [x[y][0] for y in ['__s%d' %i for i in np.arange(len(model.species))]]
# print(first_conc)
#
# for i,sp in enumerate(model.species):
#     print i,":", sp

# print("IkBmrna conc")
# print(x['IkBamrna'])
#
# print('IKK conc')
# print(x['IKK_obs'])
# print(x['IKKany'])

# # with PdfPages('multipage_pdf.pdf') as pdf:
# plt.figure()
# plt.plot(tspan/60, x['CI'], label="CI")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('CI')
#
# plt.figure()
# plt.plot(tspan/60, x['CI_k63_obs'], label="CI_k63")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('CI_k63')
#
#
# plt.figure()
# plt.plot(tspan/60, x['TNF_obs'], label="TNF")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('TNF')
# # pdf.savefig()
# # plt.close()
#
# plt.figure()
# plt.plot(tspan/60, x['RIP1ub'], label="RIP1ub")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('RIP1')
# # pdf.savefig()
# # plt.close()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBa_obs'], label="IkBa")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('IkBa')
# # pdf.savefig()
# # plt.close()
#
# plt.figure()
# plt.plot(tspan/60, x['RIP13po4_obs'], label="RIP13po4")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('RIP13po4')
# # pdf.savefig()
# # plt.close()
#
# plt.figure()
# plt.plot(tspan/60, x['NFkB_obs'], label="NFkB")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = .0002)
# plt.legend(loc=0)
# plt.title('NFkB')
# # pdf.savefig()
# # plt.close()
#
# plt.figure()
# plt.plot(tspan/60, x['IkBaNFkB_obs'], label="IkBaNFkB")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = .0002)
# plt.legend(loc=0)
# plt.title('IkBaNFkB')
# # pdf.savefig()
# # plt.close()
#
# plt.figure()
# plt.plot(tspan/60, x['NFkBn_obs'], label="NFkBn")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('NFkBn')
# # pdf.savefig()
# # plt.close()
#
# #
# plt.figure()
# plt.plot(tspan/60, x['RIPk63_obs'], label="RIPk63")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('RIPk63')
# #
# plt.figure()
# plt.plot(tspan/60, x['IKK'], label="IKK")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('IKK')
#
# plt.figure()
# plt.plot(tspan/60, x['IkBamrna'], label="IkBamrna")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('mrna')
# # pdf.savefig()
# # plt.close()
#
# plt.figure()
# plt.plot(tspan/60, x['MLKLa_obs'], label="MLKLa")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('MLKLa')
# # pdf.savefig()
# # plt.close()
#
# plt.show()





