# exported from PySB model 'model'

from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, Annotation, ANY, WILD
from pysb.bng import generate_equations, generate_network
from pysb.kappa import contact_map, set_kappa_path, influence_map
from pysb.tools.render_reactions import run
import pygraphviz as pyg

Model()

Monomer('TNFa', ['brec'])
Monomer('TNFR1', ['blig', 'bDD'])
Monomer('TRADD', ['bDD1', 'bDD2'])
Monomer('RIP1', ['bDD', 'btraf', 'bRHIM', 'bMLKL', 'state'], {'state': ['unmod', 'ub', 'po4', 'trunc']})
Monomer('TRAF', ['brip', 'bciap', 'zfinger'])
Monomer('cIAP', ['btraf'])
Monomer('NSC', ['bnfkb'])
Monomer('NFkB', ['bnsc', 'state'], {'state': ['I', 'A']})
Monomer('CYLD', ['btraf', 'state'], {'state': ['U', 'T']})
Monomer('FADD', ['bDD', 'bDED1', 'bDED2'])
Monomer('proC8', ['bDED'])
Monomer('C8', ['bf', 'state'], {'state': ['A', 'I']})
Monomer('flip_L', ['bDED'])
Monomer('flip_S', ['bDED'])
Monomer('RIP3', ['bRHIM', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})

Parameter('TNFa_0', 600.0)
Parameter('TNFR1_0', 4800.0)
Parameter('TRADD_0', 9000.0)
Parameter('RIP1_0', 12044.0)
Parameter('TRAF_0', 9000.0)
Parameter('cIAP_0', 9000.0)
Parameter('NSC_0', 0.0)
Parameter('NFkB_0', 50000.0)
Parameter('CYLD_0', 9000.0)
Parameter('FADD_0', 8030.0)
Parameter('bind_TNFa_TNFR1_kf', 1e-06)
Parameter('bind_TNFa_TNFR1_kr', 0.001)
Parameter('bind_TNFR1ANY_TRADD_kf', 1e-06)
Parameter('bind_TNFR1ANY_TRADD_kr', 0.001)
Parameter('bind_TRADDANYTNFR1ANYANY_RIP1unmod_kf', 1e-06)
Parameter('bind_TRADDANYTNFR1ANYANY_RIP1unmod_kr', 0.001)
Parameter('bind_RIP1ANYunmodTNFR1ANYANYTRADDANYANY_TRAF_kf', 1e-06)
Parameter('bind_RIP1ANYunmodTNFR1ANYANYTRADDANYANY_TRAF_kr', 0.001)
Parameter('bind_TRAF_cIAP_kf', 1e-06)
Parameter('bind_TRAF_cIAP_kr', 0.001)
Parameter('bind_TRAF_CYLD_kf', 1e-06)
Parameter('bind_TRAF_CYLD_kr', 0.001)
Parameter('RIP1_ubiq_kc', 0.1)
Parameter('RIP1_deub_kc', 0.1)
Parameter('RIP1_degr_kc', 0.1)
Parameter('NSC_esta_kc', 1e-08)
Parameter('bind_TNFR1ANYTNFaANY_TRADDANYRIP1ANYunmod_kf', 1e-06)
Parameter('bind_TNFR1ANYTNFaANY_TRADDANYRIP1ANYunmod_kr', 0.001)
Parameter('TNFR1_FADD_kc_1', 0.1)

Parameter('TRADD_FADD_kc_1', 0.1)

Parameter('flip_L_0', 39023.0)
Parameter('flip_S_0', 39023.0)
Parameter('proC8_0', 16057.0)
Parameter('C8_0', 0.0)
Parameter('RIP3_0', 20000.0)
Parameter('MLKL_0', 1000000.0)

Parameter('bind_FADDANY_flip_S_kf', 7.27e-06)
Parameter('bind_FADDANY_flip_S_kr', 0.018)
Parameter('bind_FADD_flip_S_2_kf', 7.27e-06)
Parameter('bind_FADD_flip_S_2_kr', 0.018)






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


Rule('TNFR1_to_FADD_2', TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=3) % RIP1(bDD=3, btraf=None) + FADD(bDD=None) >> TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1), TNFR1_FADD_kc_2)
Rule('TRADD_to_FADD_2', TRADD(bDD1=None, bDD2=2) % RIP1(bDD=2, btraf=None) + FADD(bDD=None) >> TRADD(bDD1=None, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1), TRADD_FADD_kc_2)
Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY, btraf=None) + FADD(bDD=None) <> TRADD(bDD1=50, bDD2=ANY) % RIP1(bDD=ANY, btraf=None) % FADD(bDD=50), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)
Rule('bind_RIP1ANYunmod_RIP3unmod', RIP1(bDD=ANY, bRHIM=None, state='unmod') + RIP3(bRHIM=None, state='unmod') <> RIP1(bDD=ANY, bRHIM=1, state='unmod') % RIP3(bRHIM=1, state='unmod'), bind_RIP1ANYunmod_RIP3unmod_kf, bind_RIP1ANYunmod_RIP3unmod_kr)
Rule('Rip3_PO4lation', RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='po4'), k19)
Rule('Rip1_PO4lation', RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='po4') >> RIP1(bRHIM=ANY, state='po4') % RIP3(bRHIM=ANY, state='po4'), k20)
Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='unmod') <> RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf, bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)
Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod') >> RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='active'), catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)
Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', FADD(bDD=None, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='unmod') <> FADD(bDD=50, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
Rule('catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc', FADD(bDD=50, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod') >> FADD(bDD=None, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='trunc'), catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)
Rule('bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
Rule('catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc', TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod') >> TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='trunc'), catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)
Rule('bind_FADDproC8proC8_RIP1unmod', FADD(bDD=None, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) + RIP1(bDD=None, state='unmod') <> FADD(bDD=50, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) % RIP1(bDD=50, state='unmod'), bind_FADDproC8proC8_RIP1unmod_kf, bind_FADDproC8proC8_RIP1unmod_kr)
Rule('RIP1_trunc_3b', RIP1(bDD=ANY, state='unmod') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> RIP1(bDD=None, state='trunc') + FADD(bDD=None, bDED1=None, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc3)
Rule('bind_TRADDFADDproC8proC8_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDproC8proC8_RIP1unmod_kf, bind_TRADDFADDproC8proC8_RIP1unmod_kr)
Rule('RIP1_trunc_4b', RIP1(bDD=ANY, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY) % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> RIP1(bDD=None, state='trunc') + TRADD(bDD1=ANY, bDD2=None) % FADD(bDD=ANY, bDED1=None, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc4)
Rule('bind_C8A_RIP1unmod_to_C8ARIP1unmod', C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='unmod') <> C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod'), bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf, bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr)
Rule('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc', C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod') >> C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='trunc'), catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc)












Parameter('bind_FADDANY_proC8_kf', 7.27e-06)
Parameter('bind_FADDANY_proC8_kr', 0.018)
Parameter('bind_FADD_proC8_2_kf', 7.27e-06)
Parameter('bind_FADD_proC8_2_kr', 0.018)
Parameter('bind_FADDANY_flip_L_kf', 7.27e-05)
Parameter('bind_FADDANY_flip_L_kr', 0.018)
Parameter('bind_FADD_flip_L_2_kf', 7.27e-06)
Parameter('bind_FADD_flip_L_2_kr', 0.018)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf', 1e-07)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr', 0.4)
Parameter('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc', 0.01)
# Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)
# Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)
Parameter('kc_c8_1', 1.0)
Parameter('kc_c8_2', 1.0)
Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf', 1e-06)
Parameter('bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr', 0.001)
Parameter('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc', 0.1)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kf', 1e-06)
Parameter('bind_C8A_CYLDU_to_C8ACYLDU_kr', 0.001)
Parameter('catalyze_C8ACYLDU_to_C8A_CYLDT_kc', 0.1)
Rule('bind_FADDANY_proC8', FADD(bDD=ANY, bDED1=None) + proC8(bDED=None) <> FADD(bDD=ANY, bDED1=1) % proC8(bDED=1), bind_FADDANY_proC8_kf, bind_FADDANY_proC8_kr)
Rule('bind_FADD_proC8_2', FADD(bDD=ANY, bDED2=None) + proC8(bDED=None) <> FADD(bDD=ANY, bDED2=1) % proC8(bDED=1), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)
Rule('bind_FADDANY_flip_L', FADD(bDD=ANY, bDED1=None) + flip_L(bDED=None) <> FADD(bDD=ANY, bDED1=1) % flip_L(bDED=1), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)
Rule('bind_FADD_flip_L_2', FADD(bDD=ANY, bDED2=None) + flip_L(bDED=None) <> FADD(bDD=ANY, bDED2=1) % flip_L(bDED=1), bind_FADD_flip_L_2_kf, bind_FADD_flip_L_2_kr)
Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)
Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)
Rule('C8_activation1', FADD(bDD=None, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> FADD(bDD=None, bDED1=None, bDED2=None) + C8(bf=None, state='A'), kc_c8_1)
Rule('C8_activation2', TRADD(bDD1=None, bDD2=1) % FADD(bDD=1, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> TRADD(bDD1=None, bDD2=1) % FADD(bDD=1, bDED1=None, bDED2=None) + C8(bf=None, state='A'), kc_c8_2)
Rule('bind_C8A_RIP3unmod_to_C8ARIP3unmod', C8(bf=None, state='A') + RIP3(bRHIM=None, state='unmod') <> C8(bf=1, state='A') % RIP3(bRHIM=1, state='unmod'), bind_C8A_RIP3unmod_to_C8ARIP3unmod_kf, bind_C8A_RIP3unmod_to_C8ARIP3unmod_kr)
Rule('catalyze_C8ARIP3unmod_to_C8A_RIP3trunc', C8(bf=1, state='A') % RIP3(bRHIM=1, state='unmod') >> C8(bf=None, state='A') + RIP3(bRHIM=None, state='trunc'), catalyze_C8ARIP3unmod_to_C8A_RIP3trunc_kc)
Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') + CYLD(btraf=None, state='U') <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)
Rule('catalyze_C8ACYLDU_to_C8A_CYLDT', C8(bf=1, state='A') % CYLD(btraf=1, state='U') >> C8(bf=None, state='A') + CYLD(btraf=None, state='T'), catalyze_C8ACYLDU_to_C8A_CYLDT_kc)


Initial(TNFa(brec=None), TNFa_0)
Initial(TNFR1(blig=None, bDD=None), TNFR1_0)
Initial(TRADD(bDD1=None, bDD2=None), TRADD_0)
Initial(RIP1(bDD=None, btraf=None, bRHIM=None, bMLKL=None, state='unmod'), RIP1_0)
Initial(TRAF(brip=None, bciap=None, zfinger=None), TRAF_0)
Initial(cIAP(btraf=None), cIAP_0)
Initial(NSC(bnfkb=None), NSC_0)
Initial(NFkB(bnsc=None, state='I'), NFkB_0)
Initial(CYLD(btraf=None, state='U'), CYLD_0)
Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
Initial(RIP3(bRHIM=None, state='unmod'), RIP3_0)
Initial(flip_L(bDED=None), flip_L_0)
Initial(flip_S(bDED=None), flip_S_0)
Initial(proC8(bDED=None), proC8_0)
Initial(C8(bf=None, state='I'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKL_0)


Rule('bind_TNFa_TNFR1', TNFa(brec=None) + TNFR1(blig=None) <> TNFa(brec=1) % TNFR1(blig=1), bind_TNFa_TNFR1_kf, bind_TNFa_TNFR1_kr)
Rule('bind_TNFR1ANY_TRADD', TNFR1(blig=ANY, bDD=None) + TRADD(bDD1=None, bDD2=None) <> TNFR1(blig=ANY, bDD=1) % TRADD(bDD1=1, bDD2=None), bind_TNFR1ANY_TRADD_kf, bind_TNFR1ANY_TRADD_kr)
Rule('bind_TRADDANYTNFR1ANYANY_RIP1unmod', TRADD(bDD1=ANY, bDD2=None) % TNFR1(blig=ANY, bDD=ANY) + RIP1(bDD=None, btraf=None, state='unmod') <> TRADD(bDD1=ANY, bDD2=50) % TNFR1(blig=ANY, bDD=ANY) % RIP1(bDD=50, btraf=None, state='unmod'), bind_TRADDANYTNFR1ANYANY_RIP1unmod_kf, bind_TRADDANYTNFR1ANYANY_RIP1unmod_kr)
Rule('bind_RIP1ANYunmodTNFR1ANYANYTRADDANYANY_TRAF', RIP1(bDD=ANY, btraf=None, state='unmod') % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) + TRAF(brip=None) <> RIP1(bDD=ANY, btraf=50, state='unmod') % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % TRAF(brip=50), bind_RIP1ANYunmodTNFR1ANYANYTRADDANYANY_TRAF_kf, bind_RIP1ANYunmodTNFR1ANYANYTRADDANYANY_TRAF_kr)
Rule('bind_TRAF_cIAP', TRAF(bciap=None) + cIAP(btraf=None) <> TRAF(bciap=1) % cIAP(btraf=1), bind_TRAF_cIAP_kf, bind_TRAF_cIAP_kr)
Rule('bind_TRAF_CYLD', TRAF(zfinger=None) + CYLD(btraf=None) <> TRAF(zfinger=1) % CYLD(btraf=1), bind_TRAF_CYLD_kf, bind_TRAF_CYLD_kr)
Rule('RIP1_Ubiquitination', TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=ANY, state='unmod') % TRAF(brip=ANY, bciap=ANY) >> TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=None, state='ub') % TRAF(brip=ANY, bciap=ANY), RIP1_ubiq_kc)
Rule('RIP1_Deubiquitination', TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=ANY, state='ub') % TRAF(brip=ANY, zfinger=ANY) >> TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=None, state='unmod') % TRAF(brip=ANY, zfinger=ANY), RIP1_deub_kc)
Rule('RIP1_Degradation', TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=ANY, state='ub') % TRAF(brip=ANY) >> TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) + TRAF(brip=None), RIP1_degr_kc)
Rule('Establish_NFkB_Signaling_Complex', TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=ANY, state='ub') % TRAF(brip=ANY, bciap=ANY) >> NSC(bnfkb=None), NSC_esta_kc)
Rule('bind_TNFR1ANYTNFaANY_TRADDANYRIP1ANYunmod', TNFR1(blig=ANY, bDD=None) % TNFa(brec=ANY) + TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY, state='unmod') <> TNFR1(blig=ANY, bDD=50) % TNFa(brec=ANY) % TRADD(bDD1=50, bDD2=ANY) % RIP1(bDD=ANY, state='unmod'), bind_TNFR1ANYTNFaANY_TRADDANYRIP1ANYunmod_kf, bind_TNFR1ANYTNFaANY_TRADDANYRIP1ANYunmod_kr)



# Rule('TNFR1_to_FADD_1', TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=3) % RIP1(bDD=3, btraf=4) % TRAF(brip=4) + FADD(bDD=None) >> TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1) + TRAF(brip=None), TNFR1_FADD_kc_1)



#
# Rule('TNFR1_to_FADD_2', TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=3) % RIP1(bDD=3, btraf=None) + FADD(bDD=None) >> TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1), TNFR1_FADD_kc_2)
# Rule('TRADD_to_FADD_2', TRADD(bDD1=None, bDD2=2) % RIP1(bDD=2, btraf=None) + FADD(bDD=None) >> TRADD(bDD1=None, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1), TRADD_FADD_kc_2)
# Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY, btraf=None) + FADD(bDD=None) <> TRADD(bDD1=50, bDD2=ANY) % RIP1(bDD=ANY, btraf=None) % FADD(bDD=50), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)
# Rule('bind_RIP1ANYunmod_RIP3unmod', RIP1(bDD=ANY, bRHIM=None, state='unmod') + RIP3(bRHIM=None, state='unmod') <> RIP1(bDD=ANY, bRHIM=1, state='unmod') % RIP3(bRHIM=1, state='unmod'), bind_RIP1ANYunmod_RIP3unmod_kf, bind_RIP1ANYunmod_RIP3unmod_kr)
# Rule('Rip3_PO4lation', RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='po4'), k19)
# Rule('Rip1_PO4lation', RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='po4') >> RIP1(bRHIM=ANY, state='po4') % RIP3(bRHIM=ANY, state='po4'), k20)
# Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='unmod') <> RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf, bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)
# Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod') >> RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='active'), catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)
# Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', FADD(bDD=None, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='unmod') <> FADD(bDD=50, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
# Rule('catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc', FADD(bDD=50, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod') >> FADD(bDD=None, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='trunc'), catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)
# Rule('bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
# Rule('catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc', TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod') >> TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='trunc'), catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)
# Rule('bind_FADDproC8proC8_RIP1unmod', FADD(bDD=None, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) + RIP1(bDD=None, state='unmod') <> FADD(bDD=50, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) % RIP1(bDD=50, state='unmod'), bind_FADDproC8proC8_RIP1unmod_kf, bind_FADDproC8proC8_RIP1unmod_kr)
# Rule('RIP1_trunc_3b', RIP1(bDD=ANY, state='unmod') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> RIP1(bDD=None, state='trunc') + FADD(bDD=None, bDED1=None, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc3)
# Rule('bind_TRADDFADDproC8proC8_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDproC8proC8_RIP1unmod_kf, bind_TRADDFADDproC8proC8_RIP1unmod_kr)
# Rule('RIP1_trunc_4b', RIP1(bDD=ANY, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY) % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> RIP1(bDD=None, state='trunc') + TRADD(bDD1=ANY, bDD2=None) % FADD(bDD=ANY, bDED1=None, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc4)
# Rule('bind_C8A_RIP1unmod_to_C8ARIP1unmod', C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='unmod') <> C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod'), bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf, bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr)
# Rule('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc', C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod') >> C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='trunc'), catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc)
























# Rule('TRADD_to_FADD_1', TRADD(bDD1=None, bDD2=2) % RIP1(bDD=2) % TRAF(brip=3) + FADD(bDD=None) >> TRADD(bDD1=None, bDD2=None) + RIP1(bDD=1) % FADD(bDD=1) + TRAF(brip=None), TRADD_FADD_kc_1)


# Rule('bind_FADDANY_flip_S', FADD(bDD=ANY, bDED1=None) + flip_S(bDED=None) <> FADD(bDD=ANY, bDED1=1) % flip_S(bDED=1), bind_FADDANY_flip_S_kf, bind_FADDANY_flip_S_kr)
# Rule('bind_FADD_flip_S_2', FADD(bDD=ANY, bDED2=None) + flip_S(bDED=None) <> FADD(bDD=ANY, bDED2=1) % flip_S(bDED=1), bind_FADD_flip_S_2_kf, bind_FADD_flip_S_2_kr)






Observable('Obs_TNFa', TNFa(brec=None))
Observable('ComplexI', TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(state='unmod'))
Observable('Obs_RIP1', RIP1(state='unmod') + RIP1(state='po4') + RIP1(state='ub'))
Observable('Obs_MLKL', MLKL(state='active'))

generate_equations(model)
generate_network(model)

print("DONE")