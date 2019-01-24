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
# Monomer('MAPK', ['btak', 'state'], {'state': ['inactive', 'active']})
# Monomer('NFkB', ['ikb', 'loc'], {'loc': ['C', 'N']})
# Monomer('IkBa', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
# Monomer('IkBb', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
# Monomer('IkBe', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
# Monomer('IkBd', ['ikk', 'nfkb', 'loc', 'state'], {'loc': ['C', 'N'], 'state': ['inactive', 'active', 'phos']})
# Monomer('IkBa_mRNA')
# Monomer('IkBb_mRNA')
# Monomer('IkBe_mRNA')
# Monomer('IkBd_mRNA')


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
# Monomer('gene', ['type','state'], {'type': ['a', 'b', 'e', 'd'],'state':['off', 'on']})

# Declares NFkB, IKK1, and IKK2 which having a binding site to IkB and NFkB can transport between the Cytoplasm
# and the Nucleus. IKK1 and IKK2 only exist in the Cytoplasm.

# def nfkb_and_ikk_monomers():
Monomer('NFkB', ['ikb', 'S'], {'S': ['C', 'N']})



Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf', 1e-06)
Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr', 0.001)
Parameter('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc', 0.1)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kf', 1e-06)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kr', 0.001)
Parameter('catalyze_C8ACYLDU_to_C8A_CYLDT_kc', 0.1)
# Parameter('TNF_0', 23500.0)
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
Parameter('bind_RIP1K63ubANY_NEMO_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_NEMO_kr', 0.001)
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
Parameter('k20', 0.001)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 1e-06)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 1e-03)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 0.1)
Parameter('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf', 1e-02) #1e-06
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

Parameter('IKK_off_0', 0.0985) #TRADD-TRAF-RIP
Initial(IKK(bind=None, bnemo = None, state = 'I'), IKK_off_0)

Parameter('IKK_0', 0.0002) #TRADD-TRAF-RIP
Initial(IKK(bind=None, bnemo = None, state = 'A'), IKK_0)

Parameter('IKKi_0', 0.0013) #TRADD-TRAF-RIP
Initial(IKK(bind=None, bnemo = None, state = 'AI'), IKKi_0)

Parameter('TNF_0', .2) #698 is 30ng/ml of TNF
Initial(TNF(brec=None), TNF_0)

Parameter('a20_0', 0.0049) #TRADD-TRAF-RIP
Initial(A20(brip = None), a20_0)

Parameter('a20t_0', 0.0001) #TRADD-TRAF-RIP
Initial(A20t(), a20t_0)
Parameter('TNFR_0', .02)
Parameter('TRADD_0', .01)
Parameter('RIP1_0', .1) #47000
Parameter('TRAF_0', .02)
Parameter('cIAP_0', .02) #10000
# # Parameter('A20_0', 50000.0) #2256
Parameter('CYLD_0',.02) #50000
Parameter('TAK1_0', .002)
Parameter('NEMO_0', .1)
Parameter('LUBAC_0', .02)
# Parameter('IKK_0', 75000.0)
# Parameter('NFkBn_0', 100000.0) #100,000
# Parameter('IkBa_0', 10000)
# Parameter('IkBb_0', 10000)
# Parameter('IkBe_0', 10000)
# Parameter('IkBd_0', 10000)

# Parameter('TNF_0', ) #698 is 30ng/ml of TNF
# Parameter('TNFR_0', 47000.0)
# Parameter('TRADD_0', 70000.0)
# Parameter('RIP1_0', 47000.0) #47000
# Parameter('TRAF_0', 47000.0)
# # Parameter('cIAP_0', 10000.0) #10000
# # Parameter('A20_0', 50000.0) #2256
# # Parameter('CYLD_0', 50000.0) #50000
# Parameter('TAK1_0', 9000.0)
# Parameter('NEMO_0', 10000.0)
# Parameter('LUBAC_0', 10000.0)
# Parameter('IKK_0', 75000.0)
# Parameter('NFkBn_0', 100000.0) #100,000
# Parameter('IkBa_0', 10000)
# Parameter('IkBb_0', 10000)
# Parameter('IkBe_0', 10000)
# Parameter('IkBd_0', 10000)
# Parameter('MAPK_0', 10000.0)
#
Initial(TNFR(blig=None, brip=None, bDD = None), TNFR_0)
Initial(TRADD(brec=None, brip=None, state='unmod', bDD1 = None, bDD2 = None), TRADD_0)
Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'), RIP1_0)
Initial(TRAF(brip=None, bciap=None, bcyld = None, state='unmod'), TRAF_0)
Initial(cIAP(btraf=None), cIAP_0)
Initial(LUBAC(brip=None), LUBAC_0)

# Initial(A20(brip=None), A20_0)
Initial(CYLD(brip=None, btraf = None), CYLD_0)


#
Initial(TAK1(brip=None, bmapk=None), TAK1_0)
Initial(NEMO(brip=None, btak=None, bikk=None), NEMO_0)
# Initial(IKK(bind=None, bnemo = None, state = 'I'), IKK_0)
# Initial(NFkB(ikb=None, loc='N'), NFkBn_0)
# Initial(IkBa(ikk=None, nfkb=None, loc='C', state='active'), IkBa_0)
# Initial(IkBb(ikk=None, nfkb=None, loc='C', state='active'), IkBb_0)
# Initial(IkBe(ikk=None, nfkb=None, loc='C', state='active'), IkBe_0)
# Initial(IkBd(ikk=None, nfkb=None, loc='C', state='active'), IkBd_0)
# # Initial(MAPK(btak=None, state='inactive'), MAPK_0)
#
# # Parameter('TNFa_0', 600.0)
# # Parameter('TNFR1_0', 4800.0)
# # Parameter('TRADD_0', 9000.0)
# # Parameter('RIP1_0', 12044.0)
# # Parameter('TRAF_0', 9000.0)
# # Parameter('cIAP_0', 9000.0)
# # Parameter('NSC_0', 0.0)
# # Parameter('NFkB_0', 50000.0)
# # Parameter('CYLD_0', 9000.0)
# Parameter('FADD_0', 8030.0)
# # Parameter('flip_L_0', 39023.0)
# # Parameter('flip_S_0', 39023.0)
# # Parameter('proC8_0', 16057.0)
# Parameter('C8_0', 10000) #10000
# Parameter('RIP3_0', 20000.0) #20000
# # Parameter('MLKL_0', 10000.0) # 1000000
# Parameter('MLKLa_0', 50000.0) # 100000
# Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
# Initial(RIP3(bRHIM=None, bDD = None, state='unmod'), RIP3_0)
# # Initial(flip_L(bDED=None), flip_L_0)
# # Initial(flip_S(bDED=None), flip_S_0)
# # Initial(proC8(bDED=None), proC8_0)
# Initial(C8(bf=None, state='A'), C8_0)
# Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)
# Initial(MLKL(bRHIM=None, state='active'), MLKLa_0)


# TNF + TNFR <> TNF:TNFR
# TNFR + TRADD <> TNFR:TRADD
# TNFR + RIP1 <> TNFR:RIP1
# TNFR:TRADD:RIP1(deub) >> TRADD:RIP1(deub) + TNFR
# TNFR:RIP1(deub) >> TNFR + RIP1(deub)

#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)

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


Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
     >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), CompI_UB2)



Rule('ComplexI_Lubac', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
     bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)


#RIP1 K63ub to be deub by A20 or CYLD
#
# Rule('bind_RIP1K63ubANY_A20', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
#      bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)
#
# # Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
# #      >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None),
# #      k_A20_1)
#
#
# Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
#      >>  TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub'), k_A20_1)
#
#
#
# Rule('bind_RIP1K63ubANY_CYLD', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
#      bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)
#
# Rule('A20_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
#      >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
#      k_A20_1)


# #RIP1 deub and necrosome formation
#
# Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)
#
# #DOES THIS HAPPEN?!?!??! @TODO
# # Rule('bind_FADD_flip_L_2', C8(bf=None, state='I') + flip_L(bDED=None) <> C8(bf=1, state='I') % flip_L(bDED=1), bind_FADD_flip_L_2_kf, bind_FADD_flip_L_2_kr)
# # Rule('C8_activation1', C8(bf=1, state='I') % flip_L(bDED=1) >> flip_L(bDED=None) + C8(bf=None, state='A'), kc_c8_1)
#
#
# Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, state='A')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='A'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)
#
#
# Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='A') + flip_L(bDED=None)
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)
#
#
# TNF
# Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) + RIP3(bRHIM=None, bDD = None, state='unmod')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_proC8_kf, bind_FADDANY_proC8_kr)
#
#
# # Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') +  <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)
#
#
#
# Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)
#
#
# Parameter('TNF_deg1', .1)
# Rule('TNF_deg', TNF(brec = None) >> None, TNF_deg1)
#
# Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf)
#
# Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), k20)
#
# Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
#      <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)
#
# Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod')
#      >>  MLKL(bRHIM=None, state='active') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)
#
# Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)
# Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)

# #RIP1K63ub recruits TAK1, NEMO, and IKK for NFkB pathway

Parameter('bind_RIP1K63ubANY_NEMO_kf', 1e-06)
Parameter('bind_RIP1K63ubANY_NEMO_kr', 0.001)
Parameter('bind_NEMORIP1TAK1_IKKinactive_kf', 1e-06)
Parameter('bind_NEMORIP1TAK1_IKKinactive_kr', 0.001)
Parameter('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc', 0.1)
Parameter('bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kf', 1e-06)
Parameter('bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kr', 0.001)
Parameter('catalyze_TAK1ANYMAPKinactive_to_TAK1ANY_MAPKactive_kc', 0.1)


Rule('bind_RIP1K63ubANY_NEMO', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + NEMO(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7), bind_RIP1K63ubANY_NEMO_kf, bind_RIP1K63ubANY_NEMO_kr)



Rule('bind_RIP1K63ubANY_TAK1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) + TAK1(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8), bind_RIP1K63ubANY_TAK1_kf, bind_RIP1K63ubANY_TAK1_kr)



# Rule('bind_RIP1K63ubANYNEMO_LUBAC', RIP1(bub1 = None, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) + LUBAC(brip=None) <> RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3), bind_RIP1K63ubANYNEMO_LUBAC_kf, bind_RIP1K63ubANYNEMO_LUBAC_kr)




Rule('bind_NEMORIP1TAK1_IKKinactive', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) + IKK(bind=None, bnemo = None, state = 'I')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = 9,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) % IKK(bind=9, bnemo = None, state = 'I'), bind_NEMORIP1TAK1_IKKinactive_kf, bind_NEMORIP1TAK1_IKKinactive_kr)



Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = 9,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) % IKK(bind=9, bnemo = None, state = 'I')
     >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) + IKK(bind=None, bnemo = None, state = 'A'), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)

Parameter('IKKa_f_IKKai', 0.15)
Parameter('IKKai_f_IKKi', 0.02)
Rule('IKKi_IKKai', IKK(bind=None, bnemo = None, state = 'A') >> IKK(bind=None, bnemo = None, state = 'AI'), IKKa_f_IKKai)
Rule('IKKai_IKKi', IKK(bind=None, bnemo = None, state = 'AI') >> IKK(bind=None, bnemo = None, state = 'I'), IKKai_f_IKKi)

# # Rule('bind_RIP1K63ubANY_A20', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='K63ub') + A20(brip=None) <> RIP1(bscf=1, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='K63ub') % A20(brip=1), bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)
# # Rule('A20_2', RIP1(bscf=1, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='K63ub') % A20(brip=1)  >> RIP1(bscf=None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None,state='deub') + A20(brip=None), k_A20_1)
# # Rule('bind_RIP1K63ubANY_CYLD', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM= None, state='K63ub') + CYLD(brip=None) <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM= 1, state='K63ub') % CYLD(brip=1), bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)
# # Rule('A20_1', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM = 1, state='K63ub') % CYLD(brip=1) >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None,bRHIM = None, state='deub') + CYLD(brip=None), k_A20_1)
# # Rule('bind_RIP1K63ubANY_NEMO', RIP1(bub1 = None, bub2=None, bub3 = None, bDD = None, state='K63ub') + NEMO(brip=None) <> RIP1(bub1 = None, bub2=1, bub3 = None, bDD = None, state='K63ub') % NEMO(brip=1), bind_RIP1K63ubANY_NEMO_kf, bind_RIP1K63ubANY_NEMO_kr)
# # Rule('bind_RIP1K63ubANY_TAK1', RIP1(bub1 = None, bub2=1, bub3 = None, bDD = None, state='K63ub') % NEMO(brip=1) + TAK1(brip=None) <> RIP1(bub1 = None, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2), bind_RIP1K63ubANY_TAK1_kf, bind_RIP1K63ubANY_TAK1_kr)
# # Rule('bind_RIP1K63ubANYNEMO_LUBAC', RIP1(bub1 = None, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) + LUBAC(brip=None) <> RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3), bind_RIP1K63ubANYNEMO_LUBAC_kf, bind_RIP1K63ubANYNEMO_LUBAC_kr)
# # Rule('bind_NEMORIP1TAK1_IKKinactive', RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = None, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) + IKK(bind=None, bnemo = None) <> RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = 4, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) % IKK(bind=4, bnemo = None), bind_NEMORIP1TAK1_IKKinactive_kf, bind_NEMORIP1TAK1_IKKinactive_kr)
# # # Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = 4, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) % IKK(bind=4) >> RIP1(bub1 = None, bub2=None, bub3 = None, bDD = None, state='K63ub') + NEMO(brip=None) + TAK1(brip=None) + LUBAC(brip=None) + IKK(bind=None), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)
# # # Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = 4, state='K63ub') % NEMO(brip=1) % TAK1(brip=2) % LUBAC(brip=3) % IKK(bind=4) >> IKK(bind=None), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)
# # Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', RIP1(bub1 = 3, bub2=1, bub3 = 2, bDD = 4, state='K63ub') % NEMO(brip=1, bikk = None) % TAK1(brip=2) % LUBAC(brip=3) % IKK(bind=4, bnemo = None) >> NEMO(brip = None, bikk = 5) % IKK(bind=None, bnemo = 5), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)
#


# def  ikb_mrna_to_ikb():
Parameter('psynth_a', 7e-5)
Parameter('psynth_b', 1e-5)
Parameter('psynth_e', 1e-6)
Parameter('psynth_d', 1e-7)
Rule('a_synth', None >> IkBa_mRNA(), psynth_a)
Rule('b_synth', None >> IkBb_mRNA(), psynth_b)
Rule('e_synth', None >> IkBe_mRNA(), psynth_e)
Rule('d_synth', None >> IkBd_mRNA(), psynth_d)
#
# Rule('a_synth', NFkB(ikb=None, S='N') >> IkBa_mRNA() + NFkB(ikb=None, S='N'), psynth_a)
# Rule('b_synth', NFkB(ikb=None, S='N') >> IkBb_mRNA() + NFkB(ikb=None, S='N'), psynth_b)
# Rule('e_synth', NFkB(ikb=None, S='N') >> IkBe_mRNA() + NFkB(ikb=None, S='N'), psynth_e)
# Rule('d_synth', NFkB(ikb=None, S='N') >> IkBd_mRNA() + NFkB(ikb=None, S='N'), psynth_d)

#
# Parameter('a_delay', 0.1)
# Parameter('be_delay', 0.027)
# Parameter('d_delay', 0.011)
#
# Rule('genea_off_on', gene(type = 'a', state = 'off') >> gene(type = 'a',state = 'on'), a_delay)
# Rule('geneb_off_on', gene(type = 'b', state = 'off') >> gene(type = 'b',state = 'on'), be_delay)
# Rule('genee_off_on', gene(type = 'e', state = 'off') >> gene(type = 'e',state = 'on'), be_delay)
# Rule('gened_off_on', gene(type = 'd', state = 'off') >> gene(type = 'd',state = 'on'), d_delay)


Parameter('hill', 3)
Parameter('a', 8)
Parameter('b', 0.02)
Parameter('e', 0.3)
Parameter('d', 0.025)

Observable('NFkBn_free', NFkB(ikb=None, S='N'))

Expression('a_NFkBn', a*(NFkBn_free)**(hill))
Expression('b_NFkBn', b*(NFkBn_free)**(hill))
Expression('e_NFkBn', e*(NFkBn_free)**(hill))
Expression('d_NFkBn', d*(NFkBn_free)**(hill))
#
# Rule('an_mRNA', gene(type = 'a',state = 'on') >> IkBa_mRNA() + gene(type = 'a',state = 'on'), a_NFkBn)
# Rule('bn_mRNA', gene(type = 'b',state = 'on') >> IkBb_mRNA() + gene(type = 'b',state = 'on'), a_NFkBn)
# Rule('en_mRNA', gene(type = 'e',state = 'on') >> IkBe_mRNA() + gene(type = 'e',state = 'on'), a_NFkBn)
# Rule('dn_mRNA', gene(type = 'd',state = 'on') >> IkBd_mRNA() + gene(type = 'd',state = 'on'), a_NFkBn)


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
# Parameter('a_f_deg', 0.36)
# Parameter('b_f_deg', 0.12)
# Parameter('e_f_deg', 0.18)
# Parameter('d_f_deg', 0.18)
#
# Parameter('and_c_n', 0.36)
# Parameter('bnd_c_n', 0.12)
# Parameter('end_c_n', 0.18)
# Parameter('dnd_c_n', 0.18)
#
# # IkBa >> None
# # IkBb >> None
# # IkBe >> None
# # IkBd >> None
#
# # IkBa : NFkB >> NFkB
# # IkBb : NFkB >> NFkB
# # IkBe : NFkB >> NFkB
# # IkBd : NFkB >> NFkB
#
# Rule('a_c_deg', IkBa(nfkb=None, S='C') >> None, a_f_deg)
# Rule('b_c_deg', IkBb(nfkb=None, S='C') >> None, b_f_deg)
# Rule('e_c_deg', IkBe(nfkb=None, S='C') >> None, e_f_deg)
# Rule('d_c_deg', IkBd(nfkb=None, S='C') >> None, d_f_deg)
#
# Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), and_c_n)
# Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), bnd_c_n)
# Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), end_c_n)
# Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), dnd_c_n)


Observable('IKKa_obs', IKK(bind=None, bnemo = None, state = 'A'))

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

Rule('synth_A20', None >> A20(brip = None), A20_synthesis)
Rule('deg_A20', A20(brip = None) >> None, A20_deg)


# Observable('TNF_obs', TNF(brec = None))
Observable('one', TNF(brec=1) % TNFR(blig=1, brip=None))

Observable('two', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None))

Observable('three', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None)
           % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'))

Observable('four', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) %
           RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'))

Observable('five', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) %
           RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))

Observable('six', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub')
           % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6))

Observable('RIP1deub', RIP1(state='deub'))
# Observable('TAK1',TAK1(brip = ANY))
Observable('TAK1b', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = None,bRHIM=None,bMLKL=None, state='K63ub')
           % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8))
# Observable('IKK_obs', IKK(bind=None, bnemo = None, state = 'A'))
Observable('IKKb', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = 9,bRHIM=None,bMLKL=None, state='K63ub') %
           TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) % IKK(bind=9, bnemo = None, state = 'I'))


Observable('RIP1RIP3po4', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')
           % RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('MLKLa', MLKL(bRHIM=None, state='active'))

Observable('RIPMLKL', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')
           % RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod'))



# Defining Observables
Observable('NEMO_obs', NEMO(brip = ANY, btak = ANY, bikk = ANY))
Observable('A20_obs', A20(brip= ANY))
Observable('TNFR_TNF', TNFR(blig=ANY))
Observable('TRADD_CompI', TRADD(brec=ANY))
Observable('RIP1_CompI', TNFR() % RIP1())
Observable('RIP1ub', RIP1(state='K63ub'))
# Observable('TRAF_CompI', TRAF(brip=ANY) % RIP1(btraf=ANY))
# Observable('Obs_TRAF_ComplexI', TRAF(brip=ANY) % RIP1(btraf=ANY))
Observable('Ubiquitylated_ComplexI', RIP1(state='K63ub'))
# Observable('MAPK_activity', MAPK(state='active'))
# Observable('RIP1deub', RIP1(state='deub'))
Observable('Obs_ComplexII', TRADD(brec=None, brip=1) % RIP1(bscf=1, btraf=None, bub1=None, bub2=None, bub3=None, state='deub'))
Observable('TNFR_endocytosis', TNF(brec=3) % TNFR(blig=3) % RIP1(btraf=None, bub1=None, bub2=None, bub3=None, state='deub'))
Observable('CYLD_CompI_1', RIP1(bscf=ANY, bub2=2, bub3=None, state='K63ub') % CYLD(brip=2))
Observable('CYLD_CompI_2', RIP1(bscf=ANY, btraf=3, bub2=2, bub3=None, state='K63ub') % CYLD(brip=2) % TRAF(brip=3, state='K63ub'))
Observable('NFkBn_obs', NFkB(ikb=None, S='N'))
Observable('IkBa_obs', IkBa(nfkb = None, S='C'))
Observable('IkBan_obs', IkBa(nfkb = None, S='N'))
Observable('IkBaNFkB_obs', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBaNFkBn_obs', IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'))

Observable('IkBb_obs', IkBb(nfkb = None, S='C'))
Observable('IkBbn_obs', IkBb(nfkb = None, S='N'))
Observable('IkBbNFkB_obs', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBbNFkBn_obs', IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'))

Observable('IkBe_obs', IkBe(nfkb = None, S='C'))
Observable('IkBen_obs', IkBe(nfkb = None, S='N'))
Observable('IkBeNFkB_obs', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBeNFkBn_obs', IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'))


Observable('IkBd_obs', IkBd(nfkb = None, S='C'))
Observable('IkBdn_obs', IkBd(nfkb = None, S='N'))
Observable('IkBdNFkB_obs', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBdNFkBn_obs', IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'))


Observable('RIP1_obs',RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'))
Observable('RIP3_obs', RIP3(bRHIM=None, bDD = None, state='unmod'))
Observable('MLKL_obs', MLKL(bRHIM=None, state='unmod'))
Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('RIP13_obs', RIP1(bDD=ANY, bRHIM=1, state='unmod') % RIP3(bRHIM=1, bDD = None, state='unmod'))
Observable('RIP13po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('TNF_obs', TNF(brec = ANY))
Observable('RIPk63_obs',  RIP1(bscf=None, bub1=None, state='K63ub'))
Observable('TNFR_TRADD', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None))
Observable('CI_k63_obs', TRADD(brec=None, brip=1, state='K63ub') % RIP1(bscf=1, btraf=2, state='K63ub') % TRAF(brip=2, bciap=3, state='K63ub') % cIAP(btraf=3))
# Observable('IkBamrna', IkBa_mRNA())
# Observable('IKK_obs', IKK(bind=None))
# Observable('IKKany', IKK(bind = ANY))
Observable('CI', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None, state='K63ub') % TRAF(brip=4, bciap=5, state='unmod') % cIAP(btraf = 5))
#
# generate_network(model)
# generate_equations(model)
# #
# # # # print(model.initials)
#
# for i,sp in enumerate(model.species):
#     print i,":", sp
#
#

tspan = np.linspace(0, 720, 721)
# print(len(tspan))
# print(tspan)
sim = ScipyOdeSimulator(model, tspan = tspan)
simulation_result = sim.run()


for i,sp in enumerate(model.species):
    print i,":",sp

print(len(model.species))

plt.figure(figsize = (15,10))
# plt.figure()
plt.subplot(231)
plt.plot(tspan/60, simulation_result.observables['TNF_obs'], marker = '*',label = 'TNF')
# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(232)
plt.plot(tspan/60, simulation_result.observables['NEMO_obs'],marker = '*',label = 'NEMO')
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(233)
plt.plot(tspan/60, simulation_result.observables['IKKa_obs'],marker = '*', label = 'IKKa')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(234)
plt.plot(tspan/60, simulation_result.observables['NFkBn_obs'], marker = '*', label = 'NFkBn')
# plt.plot(tspan/60, simulation_result.observables['NFkBn_obs'],  color = 'r',label = 'NFkBn_mat')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(235)
plt.plot(tspan/60, simulation_result.observables['A20_obs'],marker = '*', label = 'A20')
# plt.plot(tspan/60, simulation_result.observables['obs_A20t'], color = 'r',label = 'A20_mat')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(236)
plt.plot(tspan/60, simulation_result.observables['IkBd_obs'], marker = '*',label = 'IkBd')
# plt.plot(tspan/60, simulation_result.observables['IkBa_obs'], color = 'r', label = 'IkBa_mat')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.tight_layout()
plt.show()



tspan = np.linspace(0, 1440, 1441)
x = odesolve(model,tspan,verbose=True)
#
# plt.figure()
# plt.plot(tspan/60, x['RIP1deub'], label="RIP1deub")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('RIP13po4')
#
plt.figure()
plt.plot(tspan/60, x['IKK_obs'], label="IKK")
# plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
plt.xlabel("Time (in hr)", fontsize=16)
plt.ylabel("Molecules per Cell", fontsize=16)
# plt.ylim(ymin = 0, ymax = 0.00000008)
plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['RIPMLKL'], label="RIPMLKL")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['RIP1RIP3po4'], label="RIP1RIP3po4")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)

plt.figure()
plt.plot(tspan/60, x['MLKLa'], label="MLKLa")
# plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
plt.xlabel("Time (in hr)", fontsize=16)
plt.ylabel("Molecules per Cell", fontsize=16)
# plt.ylim(ymin = 0, ymax = 0.00000008)
plt.legend(loc=0)
#
# plt.show()
# plt.figure()
# plt.plot(tspan/60, x['IkBa_obs'], label="IkBa")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('IkBa')
#
# plt.figure()
# plt.plot(tspan/60, x['IkBb_obs'], label="IkBb")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('IkBb')
#
# plt.figure()
# plt.plot(tspan/60, x['IkBe_obs'], label="IkBe")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('IkBe')
#
# plt.figure()
# plt.plot(tspan/60, x['IkBd_obs'], label="IkBd")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('IkBd')
# # # plt.show()
# # #
#
# #
# # last_conc = [x[y][-1] for y in ['__s%d' %i for i in np.arange(len(model.species))]]
# # print(last_conc)
# #
# # first_conc = [x[y][0] for y in ['__s%d' %i for i in np.arange(len(model.species))]]
# # print(first_conc)
# #
# #
# # # print("IkBmrna conc")
# # # print(x['IkBamrna'])
# # #
# # # print('IKK conc')
# # # print(x['IKK_obs'])
# # # print(x['IKKany'])
# #
# # # # with PdfPages('multipage_pdf.pdf') as pdf:
# plt.figure()
# plt.plot(tspan/60, x['CI'], label="CI")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('CI')
# #
# # plt.figure()
# # plt.plot(tspan/60, x['CI_k63_obs'], label="CI_k63")
# # # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# # plt.xlabel("Time (in hr)", fontsize=16)
# # plt.ylabel("Molecules per Cell", fontsize=16)
# # # plt.ylim(ymin = 0, ymax = 0.00000008)
# # plt.legend(loc=0)
# # plt.title('CI_k63')
# #
# #
plt.figure()
plt.plot(tspan/60, x['TNF_obs'], label="TNF")
# plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
plt.xlabel("Time (in hr)", fontsize=16)
plt.ylabel("Molecules per Cell", fontsize=16)
# plt.ylim(ymin = 0, ymax = 0.00000008)
plt.legend(loc=0)
plt.title('TNF')
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
#
# plt.figure()
# plt.plot(tspan/60, x['RIP13po4_obs'], label="RIP13po4")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('RIP13po4')
# # # # pdf.savefig()
# # # # plt.close()
# # #
plt.figure()
plt.plot(tspan/60, x['NFkB_obs'], label="NFkB")
# plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
plt.xlabel("Time (in hr)", fontsize=16)
plt.ylabel("Molecules per Cell", fontsize=16)
plt.ylim(ymin = 0, ymax = 1000)
plt.legend(loc=0)
plt.title('NFkB')
#
plt.figure()
plt.plot(tspan/60, x['NFkBn_obs'], label="NFkBn")
# plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
plt.xlabel("Time (in hr)", fontsize=16)
plt.ylabel("Molecules per Cell", fontsize=16)
# plt.ylim(ymin = 0, ymax = .0002)
plt.legend(loc=0)
plt.title('NFkBn')
# # # pdf.savefig()
# # # plt.close()
# #
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
# plt.plot(tspan/60, x['IKK_obs'], label="IKK")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('IKK')
# #
# plt.figure()
# plt.plot(tspan/60, x['IkBamrna'], label="IkBamrna")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Molecules per Cell", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('mrna')
# # # # pdf.savefig()
# # # # plt.close()
# # #
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
plt.show()