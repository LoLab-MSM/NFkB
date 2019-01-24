"""
    Overview
    ========

    PySB implementations of the apoptosis-necroptosis reaction model version 1.0
    (ANRM 2.0) originally published in [Irvin,NECRO2016]_.

    This model provides information about the dynamic biomolecular events that
    commit a cell to apoptosis or necroptosis in response to TNFa signalling.
    PARP1 cleavage or activity serves, in this case, as a marker for apoptosis
    or necrosis respectively. An alternate marker for necroptosis

    This file contains functions that implement the parts of the extrinsic apoptosis
    pathway and activation of a necroptosis reporter proteins, PARP and MLKL.

    The modules are organized into three overall sections:
    - ...

    These sections are further devided into...
    For receptor signalling and Bid activation:
    --...

    These modules are to be run in cooperation with Lopez Momp_monomers and Albeck Apaf_to_Parp modules.

    Parameters and Initial Conditions were taken from Jurket Cell et. al. studies:
    1. Micheau O, Tschopp J. Cell. 2003 Jul 25;114(2):181-90.
    2. http://atlasgeneticsoncology.org/Genes/TNFaID319.html (Mass of TNFa is 25.6kDa)
    3. http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=104668&ver=2 (volume of 8e-13L)
    4. Cho KH, Kolch W, Wolkenhauer O, Simulation 2003 vol. 79 no. 12 726-739
    5. Maiti S, Dai W, Alaniz RC, Hahn J, Jayaraman J, Processes 2015, 3, 1-18; doi:10.3390/pr3010001
    6. Rangamani P, Sirovich L, Biotechnology and Bioengineering, Vol.97, No. 5, 2007
    7. Pekalski J, Pawel JK, Kochanczyk N, Junkin M, Kellogg R, Tay S, Lipniacki T PLOS one http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0078887#pone.0078887.s002


    Hua, F., M. Cornejo, M. Cardone, C. Stokes, D. Lauffenburger
    J Immunol 2005; 175:985-995
    --concentrations listed in this were converted to copies per cell using a cell
    --
    """

from pysb import *
from pysb.util import alias_model_components
from pysb.macros import *
from pysb.bng import generate_equations, generate_network
import matplotlib.pyplot as plt
import numpy as np

from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np

# SECTION ONE: TNF Complex I
# ===========================
# This contains FasL binding and assembly of the DISC

# def TNF_ComplexI_Monomers():

""" Declares TNFa, TNFR1, TRADD, TRAF2, RIP1, cIAP, NEMO, LUBAC, A20, TAK1, CYLD, IKKab, IKK
NFkB, MAPK. Binding sites are named after (when known) subdomains that
facilitate binding. Otherwise the binding site will reflect what it binds with.

The 'state' site denotes various localization and/or activity states of a
Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
localization.
"""

Model()

Monomer('TNF'   ,   ['brec'])           #TNFa
Monomer('TNFR'  ,   ['blig','brip'])    #TNFR1
Monomer('TRADD' ,   ['brec','brip', 'state'], {'state':['unmod', 'K63ub']} )   #TRADD
Monomer('RIP1'  ,   ['bscf', 'btraf', 'bub1', 'bub2', 'bub3', 'state'], {'state':['unmod', 'K63ub', 'deub']})   #RIP1
Monomer('TRAF'  ,   ['brip', 'bciap', 'state'], {'state':['unmod', 'K63ub']})   #TRAF2
Monomer('cIAP'  ,   ['btraf'])          #cIAP 1 and 2
Monomer('A20'   ,   ['brip'])
Monomer('CYLD'  ,   ['brip'])
Monomer('TAK1'  ,   ['brip', 'bmapk'])
Monomer('NEMO'  ,   ['brip','btak', 'bikk'])
Monomer('LUBAC' ,   ['brip'])
Monomer('IKK'   ,   ['bind','state'], {'state':['inactive', 'active']})
Monomer('MAPK'  ,   ['btak','state'], {'state':['inactive', 'active']})
Monomer('NFkB', ['ikb', 'loc', 'state'], {'state': ['inactive', 'active'], 'loc': ['C', 'N']})  # NFkB
Monomer('IkBa', ['ikk', 'nfkb', 'loc', 'state'], {'state': ['inactive', 'active', 'phos'], 'loc': ['C', 'N']})
Monomer('IkBb', ['ikk', 'nfkb', 'loc', 'state'], {'state': ['inactive', 'active', 'phos'], 'loc': ['C', 'N']})
Monomer('IkBe', ['ikk', 'nfkb', 'loc', 'state'], {'state': ['inactive', 'active', 'phos'], 'loc': ['C', 'N']})
Monomer('IkBd', ['ikk', 'nfkb', 'loc', 'state'], {'state': ['inactive', 'active', 'phos'], 'loc': ['C', 'N']})
alias_model_components()

# def TNF_ComplexI_Initials():
Parameter('TNF_0'   , 23500) # [1-3] 5e7 cells per mL treated with 50ng/mL TNF corresponds with 23500 copies per cell or 50nM (cellular concentration)
Parameter('TNFR_0'  , 47000) # [5-6] receptors per cell (100nM)
Parameter('TRADD_0' , 70000) # [5,6] molecules per cell
Parameter('RIP1_0'  , 47000) # [4,6] molecules per cell
Parameter('TRAF_0'  , 47000) # [4,6] molecules per cell
Parameter('cIAP_0'  , 10000) # molecules per cell (arbitrarily assigned) 10000
Parameter('A20_0'   ,  2256) # [7] complexes per cell (4.8nM)
Parameter('CYLD_0'  , 50000) # molecules per cell (arbitrarily assigned) 50000
Parameter('TAK1_0'  ,  9000) # molecules per cell
Parameter('NEMO_0'  , 10000) # [7] molecules per cell
Parameter('LUBAC_0' , 10000) # complexes per cell (arbitrarily assigned)
Parameter('IKK_0'   , 75000) # [5-6] molecules per cell
Parameter('NFkB_0'  ,100000) # [6-7] molecules per cell (250nM)
Parameter('MAPK_0'  , 10000) # molecules per cell (arbitrarily assigned)
alias_model_components()

Initial(TNF(brec=None), TNF_0)           #TNFa
Initial(TNFR(blig=None, brip=None), TNFR_0)    #TNFR1
Initial(TRADD(brec=None, brip=None, state='unmod'), TRADD_0)   #TRADD
Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'), RIP1_0)   #RIP1
Initial(TRAF(brip=None, bciap=None, state='K63ub'), TRAF_0)   #TRAF2
Initial(cIAP(btraf=None), cIAP_0)          #cIAP 1 and 2
Initial(A20(brip=None), A20_0)
Initial(CYLD(brip=None), CYLD_0)
Initial(TAK1(brip=None, bmapk=None), TAK1_0)
Initial(NEMO(brip=None,btak=None, bikk=None), NEMO_0)
Initial(LUBAC(brip=None), LUBAC_0)
Initial(IKK(bind=None, state='inactive'), IKK_0)
Initial(NFkB(ikb=None, loc = 'C', state='inactive'), NFkB_0)       #NFkB
Initial(MAPK(btak=None,state='inactive'), MAPK_0)
alias_model_components()

# def TNF_ComplexI_Assembly():
"""Reactions that assemble the TNF complex I are listed in ANRM_Reactions.xlsx. The species involved these reactions are listed in Apoptotic Species (Autosaved).xlsx. All reactions supported by published findings.

    1. TNF + TNFR <> TNF:TNFR
    2. TNF:TNFR + TRADD(unmod) <> TNF:TNFR:TRADD(unmod)
    3. TNF:TNFR + RIP1(unmod) <> TNF:TNFR:RIP1(unmod)
    4. TNFR:TRADD(unmod) + RIP1(unmod) <> TNFR:TRADD(unmod):RIP1(unmod)
    5. (TNFR or TRADD):RIP1 + TRAF <> (TNFR or TRADD):RIP1:TRAF(unmod) #CYLD de-ubiquitinates TRAF2 allowing its dissociation from MLKL. (Does this free TRAF to interact with CPLX1?
    6. TRAF + cIAP <> TRAF:cIAP
    7. TRAF(brip=None, state=K63ub) + CYLD <> TRAF(brip=1, state=K63ub):CYLD >> TRAF(brip=None, state=unmod) + CYLD
"""
bind(TNF(brec = None),'brec', TNFR(blig = None, brip = None), 'blig', [1e-6, 1e-3])
bind(TNFR(blig = ANY, brip = None), 'brip', TRADD(brec = None, brip=None), 'brec', [1e-6, 1e-3])
bind(TNFR(blig = ANY, brip = None), 'brip', RIP1(bscf = None, btraf=None, bub1 = None, bub2 = None, bub3 = None, state = 'unmod'), 'bscf', [1e-6, 1e-3])
bind(TRADD(brec=ANY, brip = None, state = 'unmod'), 'brip', RIP1(bscf = None, btraf=None, bub1 = None, bub2 = None, bub3 = None, state = 'unmod'), 'bscf',[1e-6, 1e-3])
bind(RIP1(bscf=ANY, btraf=None), 'btraf',TRAF(brip = None, state='unmod'), 'brip', [1e-6, 1e-3])

#TRAF availability is regulated by cIAP and CYLD
catalyze_state(cIAP(btraf=None), 'btraf', TRAF(brip = None), 'bciap', 'state', 'unmod', 'K63ub', [1e-6, 1e-3, 1e-1])
catalyze_state(CYLD(brip=None), 'brip', TRAF(), 'brip', 'state', 'K63ub', 'unmod', [1e-6, 1e-3, 1])

# def TNF_ComplexI_K63_Polyubiquitylation():
"""Lys-63-Poly-ubiquitylation reactions that occur in TNFR1
1. TNFR(unmod):...:cIAP >> TNFR(K63ub):...:cIAP
2. RIP1(unmod):...:cIAP >> RIP1(K63ub):...:cIAP
3. TRAF(unmod):cIAP >> TRAF(K63ub):cIAP
"""
Rule('Complex_I_ubiquitylation1', TNFR(blig = ANY, brip = 2)%RIP1(bscf = 2, btraf = 3, state = 'unmod')%TRAF(brip=3, bciap = 4, state = 'unmod')%cIAP(btraf=4) >> TNFR(blig = ANY, brip = 2)%RIP1(bscf = 2, btraf = 3, state = 'K63ub')%TRAF(bciap = 4, brip=3, state = 'K63ub')%cIAP(btraf=4), Parameter('CompI_UB1', 10))
Rule('Complex_I_ubiquitylation2',TRADD(brec = ANY, brip = 1, state='unmod')%RIP1(bscf = 1, btraf = 2, state = 'unmod')%TRAF(brip=2, bciap = 3, state = 'unmod')%cIAP(btraf=3) >> TRADD(brec = ANY, brip = 1, state='K63ub')%RIP1(bscf = 1, btraf = 2, state = 'K63ub')%TRAF(brip=2, bciap = 3, state = 'K63ub')%cIAP(btraf=3), Parameter('CompI_UB2', 10))

# def RIP1ub_complex_assembly():
"""Recruitment of A20, CYLD, NEMO, etc to Ubiquitylated Complex_I
1. RIP1(K63ub) + A20 <> RIP1(K63ub):A20
2. RIP1(K63ub) + CYLD <> RIP1(K63ub):CYLD
3. RIP1(K63ub) + NEMO <> RIP1(K63ub):NEMO
4. RIP1(K63ub) + LUBAC <> RIP1(K63ub):LUBAC
5. RIP1(K63ub) + TAK1 <> RIP1(K63ub):TAK1
    """
bind(RIP1(bscf = ANY, bub1 = None, state='K63ub'), 'bub1', A20(brip=None), 'brip', [1e-6, 1e-3])
bind(RIP1(bscf = ANY, bub2 = None, state='K63ub'), 'bub2', CYLD(brip=None), 'brip', [1e-6, 1e-3])
bind(RIP1(bscf = ANY, bub1 = None, state='K63ub'), 'bub1', TAK1(brip=None), 'brip', [1e-6, 1e-3])
bind(RIP1(bscf = ANY, bub2 = None, bub3=None, state='K63ub'), 'bub2', NEMO(brip=None), 'brip', [1e-6, 1e-3])
bind(RIP1(bscf = ANY, bub2 = 2, bub3 = None, state='K63ub')%NEMO(brip=2), 'bub3', LUBAC(brip=None), 'brip', [1e-6, 1e-3])

# def ComplexI_PCD_Reactions():
"""Pro cell death reactions A20 and CYLD based reactions. The complexes that undergo A20 and CYLD mediated reactions are written explicitly to avoid a clusterfuck of possible binding-complex interactions

    1. RIP1(:A20, :CYLD):TRAF >> A20 + CYLD + TRAF
    2. RIP1(:A20, :NEMO):TRAF >> A20 + NEMO + TRAF
    3. RIP1:A20:TRAF >> A20 + TRAF
    4. RIP1:A20:TRAF >> A20 + TRAF

    5. RIP1(:TAK, :CYLD):TRAF(ub) >> RIP1(deub) + TRAF + TAK + CYLD
    6. RIP1(:A20, :CYLD):TRAF(ub) >> RIP1(deub) + A20 + CYLD + TRAF
    7. RIP1:CYLD:TRAF(ub) >>  RIP1(deub) + CYLD + TRAF
    8. RIP1:CYLD:TRAF(ub) >>  RIP1(deub) + CYLD + TRAF

    """
Rule('A20_1', RIP1(bscf = ANY, bub1=1, bub2=2, bub3=None, btraf=3, state='K63ub')%A20(brip=1)%CYLD(brip=2)%TRAF(brip=3, state='K63ub') >> A20(brip=None) + CYLD(brip=None)+ TRAF(brip=None, state='unmod'), Parameter('k_A20_1', 10))
Rule('A20_2', RIP1(bscf = ANY, bub1=1, bub2=2, bub3=None, btraf=3, state='K63ub')%A20(brip=1)%NEMO(brip=2)%TRAF(brip=3, state='K63ub') >> A20(brip=None) + NEMO(brip=None)+ TRAF(brip=None, state='unmod'), Parameter('k_A20_2', 10))
Rule('A20_3', RIP1(bscf = ANY, bub1=1, bub2=None, bub3=None, btraf=2, state='K63ub')%A20(brip=1)%TRAF(brip=2, state='K63ub') >> A20(brip=None) + TRAF(brip=None, state = 'unmod'), Parameter('k_A20_3', 10))
Rule('A20_4', RIP1(bscf = ANY, bub1=1, bub2=None, bub3=None, btraf=2, state='K63ub')%A20(brip=1)%TRAF(brip=2, state='K63ub') >> A20(brip=None) + TRAF(brip=None, state = 'unmod'), Parameter('k_A20_4', 10))

Rule('CYLD_1', RIP1(bscf = ANY, bub1=1, bub2=2, bub3=None, btraf=3, state='K63ub')%TAK1(brip=1)%CYLD(brip=2)%TRAF(brip=3, state='K63ub') >> RIP1(bscf = ANY, bub1=None, bub2=None, bub3=None, btraf=None, state='deub') + TAK1(brip=None) + CYLD(brip=None)+ TRAF(brip=None, state='unmod'), Parameter('k_CYLD_1', 10))
Rule('CYLD_2', RIP1(bscf = ANY, bub1=1, bub2=2, bub3=None, btraf=3, state='K63ub')%A20(brip=1)%CYLD(brip=2)%TRAF(brip=3, state='K63ub') >> RIP1(bscf = ANY, bub1=None, bub2=None, bub3=None, btraf=None, state='deub')+ A20(brip=None) + TRAF(brip=None, state='unmod'), Parameter('k_CYLD_2', 10))
Rule('CYLD_3', RIP1(bscf = ANY, bub1=None, bub2=1, bub3=None, btraf=3, state='K63ub')%TRAF(brip=3, state='K63ub')%CYLD(brip=1) >> RIP1(bscf = ANY, bub1=None, bub2=None, bub3=None, btraf=None, state='deub') + TRAF(brip=None, state='unmod') + CYLD(brip=None), Parameter('k_CYLD_3', 10))
Rule('CYLD_4', RIP1(bscf = ANY, bub1=None, bub2=1, bub3=None, btraf=3, state='K63ub')%TRAF(brip=3, state='K63ub')%CYLD(brip=1) >> RIP1(bscf = ANY, bub1=None, bub2=None, bub3=None, btraf=None, state='deub') + TRAF(brip=None, state='unmod')+ CYLD(brip=None), Parameter('k_CYLD_4', 10))

# def ComplexI_Amplication_Reactions():
"""Downstream Reactions
    NFkB and MAPK activation
    1. Complex I (containing TAK1 and NEMO) + IKK(inactive) <> Complex I:IKK(inactive) >>  Complex I + IKK(active)
    2. IKK(active) + NFkB(inactive) <> IKK(active):NFkB(inactive) >> IKK(active) + NFkB(active)
    3. Complex I (containing TAK1) + MAPK(inactive) <> Complex I:MAPK(inactive) >> Complex I:MAPK(active)
    4. RIP1(K63ub):TAK1:CYLD >> RIP1(K63ub):TAK #CYLD deactivation is equavent to deletion since there is not mechanism to reactivate CYLD
    5. ComplexI (with RIP1-deub and TAK1, TRAF2 etc missing) >> TNFR + TRADD:RIP1:..Secondary Complex
    """
catalyze_complex(RIP1()%TAK1()%NEMO(bikk=None), 'bikk', IKK(bind=None, state='inactive'), 'bind', IKK(bind=None, state = 'active'), [1e-6, 1e-3, 1e-1])
# catalyze_state(IKK(bind=None, state='active'), 'bind', NFkB(), 'bikk', 'state', 'inactive', 'active', [1e-6, 1e-3, 1e-1])
catalyze_state(TAK1(brip = ANY, bmapk=None), 'bmapk', MAPK(), 'btak', 'state', 'inactive', 'active', [1e-6, 1e-3, 1e-1])

Rule('CYLD_deactivation', RIP1(bscf = ANY, bub1=1, bub2=2, state='K63ub')%TAK1(brip=1)%CYLD(brip=2) >> RIP1(bscf = ANY, bub1=1, bub2=None, state='K63ub')%TAK1(brip=1), Parameter('k_CYLD_deac', 1e-2))
Rule('Complex_II_1', TNFR(brip=3)%TRADD(brec=3, brip=2)%RIP1(bscf=2, bub1=None, bub2=None, bub3=None, btraf=None, state='deub')>> TRADD(brec=None, brip=2)%RIP1(bscf=2, bub1=None, bub2=None, bub3=None, btraf=None, state='deub') + TNFR(brip=None), Parameter('k_ComplexII_1', 10))
Rule('Complex_II_2', TNFR(brip=2)%RIP1(bscf=2, bub1=None, bub2=None, bub3=None, btraf=None, state='deub')>> RIP1(bscf=None, bub1=None, bub2=None, bub3=None, btraf=None, state='deub') + TNFR(brip=None), Parameter('k_ComplexII_2', 10))

# def ComplexI_Amplication_Reactions_2_Monomers():
#     Monomer('TNF_endocytosed')
# def ComplexI_Amplication_Reactions_2_Initials():
#     Parameter('TNF_endocytosed_0', 0)
#     alias_model_components()
#     Initial(TNF_endocytosed(), TNF_endocytosed_0)

# def ComplexI_Amplication_Reactions_2():
#     """Alternate set of Downstream Reactions that includes Receptor Recycle
#         NFkB and MAPK activation
#         1. Complex I (containing TAK1 and NEMO) + IKK(inactive) <> Complex I:IKK(inactive) >>  Complex I + IKK(active)
#         2. IKK(active) + NFkB(inactive) <> IKK(active):NFkB(inactive) >> IKK(active) + NFkB(active)
#         3. Complex I (containing TAK1) + MAPK(inactive) <> Complex I:MAPK(inactive) >> Complex I:MAPK(active)
#         4. RIP1(K63ub):TAK1:CYLD >> RIP1(K63ub):TAK #CYLD deactivation is equavent to deletion since there is not mechanism to reactivate CYLD
#         5. ComplexI (with RIP1-deub and TAK1, TRAF2 etc missing) >> TNFR_end + TRADD:RIP1:..Secondary Complex faciliated by receptor endocytosis.
#         """
#     catalyze_complex(RIP1()%TAK1()%NEMO(bikk=None), 'bikk', IKK(bind=None, state='inactive'), 'bind', IKK(bind=None, state = 'active'), [1e-6, 1e-3, 1e-1])
#     # catalyze_state(IKK(bind=None, state='active'), 'bind', NFkB(), 'bikk', 'state', 'inactive', 'active', [1e-6, 1e-3, 1e-1])
#     catalyze_state(TAK1(brip = ANY, bmapk=None), 'bmapk', MAPK(), 'btak', 'state', 'inactive', 'active', [1e-6, 1e-3, 1e-1])
#
#     Rule('Complex_II_1', TNF(brec=4)%TNFR(blig = 4, brip=3)%TRADD(brec=3, brip=2)%RIP1(bscf=2, bub1=None, bub2=None, bub3=None, btraf=None, state='deub')>> TRADD(brec=None, brip=2)%RIP1(bscf=2, bub1=None, bub2=None, bub3=None, btraf=None, state='deub') + TNF_endocytosed() + TNF(brec=None), Parameter('CPLXII_1', 10))
#     Rule('Complex_II_2', TNF(brec=3)%TNFR(blig = 3, brip=2)%RIP1(bscf=2, bub1=None, bub2=None, bub3=None, btraf=None, state='deub')>> RIP1(bscf=None, bub1=None, bub2=None, bub3=None, btraf=None, state='deub') + TNF_endocytosed()+ TNF(brec=None), Parameter('CPLXII_2', 10))
#
#     Rule('Receptor_Endocytosis', TNF(brec=3)%TNFR(blig = 3, brip=None)>> TNF_endocytosed()+ TNF(brec=None), Parameter('Endocytosis', 1e-2))
#     Rule('Receptor_Recycle', TNF_endocytosed() >> TNFR(blig=None, brip=None), Parameter('k_Receptor_Recycle_1', 1e-2))
#     Rule('CYLD_deactivation', RIP1(bscf = ANY, bub1=1, bub2=2, state='K63ub')%TAK1(brip=1)%CYLD(brip=2) >> RIP1(bscf = ANY, bub1=1, bub2=None, state='K63ub')%TAK1(brip=1), Parameter('k_CYLD_deac', 1e-2))

# def Necrosome_Initials():
Parameter('FADD_0', 8030.0)
Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
Initial(RIP3(bRHIM=None, state='unmod'), RIP3_0)
Initial(flip_L(bDED=None), flip_L_0)
Initial(flip_S(bDED=None), flip_S_0)
Initial(proC8(bDED=None), proC8_0)
Initial(C8(bf=None, state='I'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKL_0)



Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY, btraf=None) + FADD(bDD=None) <> TRADD(bDD1=50, bDD2=ANY) % RIP1(bDD=ANY, btraf=None) % FADD(bDD=50), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)
Rule('bind_FADDANY_proC8', FADD(bDD=ANY, bDED1=None) + proC8(bDED=None) <> FADD(bDD=ANY, bDED1=1) % proC8(bDED=1), bind_FADDANY_proC8_kf, bind_FADDANY_proC8_kr)
Rule('bind_FADD_proC8_2', FADD(bDD=ANY, bDED2=None) + proC8(bDED=None) <> FADD(bDD=ANY, bDED2=1) % proC8(bDED=1), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)
Rule('bind_FADDANY_flip_L', FADD(bDD=ANY, bDED1=None) + flip_L(bDED=None) <> FADD(bDD=ANY, bDED1=1) % flip_L(bDED=1), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)
Rule('bind_FADD_flip_L_2', FADD(bDD=ANY, bDED2=None) + flip_L(bDED=None) <> FADD(bDD=ANY, bDED2=1) % flip_L(bDED=1), bind_FADD_flip_L_2_kf, bind_FADD_flip_L_2_kr)
Rule('bind_FADDANY_flip_S', FADD(bDD=ANY, bDED1=None) + flip_S(bDED=None) <> FADD(bDD=ANY, bDED1=1) % flip_S(bDED=1), bind_FADDANY_flip_S_kf, bind_FADDANY_flip_S_kr)
Rule('bind_FADD_flip_S_2', FADD(bDD=ANY, bDED2=None) + flip_S(bDED=None) <> FADD(bDD=ANY, bDED2=1) % flip_S(bDED=1), bind_FADD_flip_S_2_kf, bind_FADD_flip_S_2_kr)

Rule('bind_RIP1ANYunmod_RIP3unmod', RIP1(bDD=ANY, bRHIM=None, state='unmod') + RIP3(bRHIM=None, state='unmod') <> RIP1(bDD=ANY, bRHIM=1, state='unmod') % RIP3(bRHIM=1, state='unmod'), bind_RIP1ANYunmod_RIP3unmod_kf, bind_RIP1ANYunmod_RIP3unmod_kr)
Rule('Rip3_PO4lation', RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='po4'), k19)
Rule('Rip1_PO4lation', RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='po4') >> RIP1(bRHIM=ANY, state='po4') % RIP3(bRHIM=ANY, state='po4'), k20)
Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='unmod') <> RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf, bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)
Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bMLKL=1, state='po4') % MLKL(bRHIM=1, state='unmod') >> RIP1(bMLKL=None, state='po4') + MLKL(bRHIM=None, state='active'), catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)
Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)
Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)


Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', FADD(bDD=None, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='unmod') <> FADD(bDD=50, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
Rule('catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc', FADD(bDD=50, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod') >> FADD(bDD=None, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='trunc'), catalyze_FADDANYANYflip_LANYproC8ANYRIP1unmod_to_FADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)
Rule('bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kf, bind_TRADDFADDANYANYflip_LANYproC8ANY_RIP1unmod_kr)
Rule('catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc', TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) % RIP1(bDD=50, state='unmod') >> TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=ANY, bDED2=ANY) % flip_L(bDED=ANY) % proC8(bDED=ANY) + RIP1(bDD=None, state='trunc'), catalyze_TRADDFADDANYANYflip_LANYproC8ANYRIP1unmod_to_TRADDFADDANYANYflip_LANYproC8ANY_RIP1trunc_kc)
Rule('bind_FADDproC8proC8_RIP1unmod', FADD(bDD=None, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) + RIP1(bDD=None, state='unmod') <> FADD(bDD=50, bDED1=2, bDED2=3) % proC8(bDED=2) % proC8(bDED=3) % RIP1(bDD=50, state='unmod'), bind_FADDproC8proC8_RIP1unmod_kf, bind_FADDproC8proC8_RIP1unmod_kr)


Rule('RIP1_trunc_3b', RIP1(bDD=ANY, state='unmod') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> RIP1(bDD=None, state='trunc') + FADD(bDD=None, bDED1=None, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc3)

Rule('bind_TRADDFADDproC8proC8_RIP1unmod', TRADD(bDD1=2, bDD2=None) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) + RIP1(bDD=None, state='unmod') <> TRADD(bDD1=2, bDD2=50) % FADD(bDD=2, bDED1=3, bDED2=4) % proC8(bDED=3) % proC8(bDED=4) % RIP1(bDD=50, state='unmod'), bind_TRADDFADDproC8proC8_RIP1unmod_kf, bind_TRADDFADDproC8proC8_RIP1unmod_kr)

Rule('RIP1_trunc_4b', RIP1(bDD=ANY, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY) % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> RIP1(bDD=None, state='trunc') + TRADD(bDD1=ANY, bDD2=None) % FADD(bDD=ANY, bDED1=None, bDED2=None) + C8(bf=None, state='A'), RIP1_trunc_kc4)

Rule('C8_activation1', FADD(bDD=None, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> FADD(bDD=None, bDED1=None, bDED2=None) + C8(bf=None, state='A'), kc_c8_1)

Rule('C8_activation2', TRADD(bDD1=None, bDD2=1) % FADD(bDD=1, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY) % proC8(bDED=ANY) >> TRADD(bDD1=None, bDD2=1) % FADD(bDD=1, bDED1=None, bDED2=None) + C8(bf=None, state='A'), kc_c8_2)

Rule('bind_C8A_RIP1unmod_to_C8ARIP1unmod', C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='unmod') <> C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod'), bind_C8A_RIP1unmod_to_C8ARIP1unmod_kf, bind_C8A_RIP1unmod_to_C8ARIP1unmod_kr)
Rule('catalyze_C8ARIP1unmod_to_C8A_RIP1trunc', C8(bf=1, state='A') % RIP1(bDD=None, bRHIM=1, state='unmod') >> C8(bf=None, state='A') + RIP1(bDD=None, bRHIM=None, state='trunc'), catalyze_C8ARIP1unmod_to_C8A_RIP1trunc_kc)
Rule('bind_C8A_RIP3unmod_to_C8ARIP3unmod', C8(bf=None, state='A') + RIP3(bRHIM=None, state='unmod') <> C8(bf=1, state='A') % RIP3(bRHIM=1, state='unmod'), bind_C8A_RIP3unmod_to_C8ARIP3unmod_kf, bind_C8A_RIP3unmod_to_C8ARIP3unmod_kr)
Rule('catalyze_C8ARIP3unmod_to_C8A_RIP3trunc', C8(bf=1, state='A') % RIP3(bRHIM=1, state='unmod') >> C8(bf=None, state='A') + RIP3(bRHIM=None, state='trunc'), catalyze_C8ARIP3unmod_to_C8A_RIP3trunc_kc)
Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') + CYLD(btraf=None, state='U') <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)
Rule('catalyze_C8ACYLDU_to_C8A_CYLDT', C8(bf=1, state='A') % CYLD(btraf=1, state='U') >> C8(bf=None, state='A') + CYLD(btraf=None, state='T'), catalyze_C8ACYLDU_to_C8A_CYLDT_kc)


# def survival():
Parameter('kf', .10)
Parameter('kr', .10)

# NEMO:IKK(active) + IkBa(inactive) <> IKK(active):IkBa(active)
# NEMO:IKK(active) + IkBb(inactive) <> IKK(active):IkBb(active)
# NEMO:IKK(active) + IkBe(inactive) <> IKK(active):IkBe(active)
# NEMO:IKK(active) + IkBd(inactive) <> IKK(active):IkBd(active)

Rule('IKK_activate_IkBa', IKK(bind= None, state = 'active') + IkBa(ikk = None, nfkb = None, loc = 'C', state = 'inactive') <> IKK(bind= 1,state = 'active')%IkBa(ikk = 1, nfkb = None, loc = 'C', state = 'inactive'), kf,kr)
Rule('IKK_activate_IkBb', IKK(bind= None, state='active') + IkBb(ikk=None, nfkb=None, loc='C', state='inactive') <> IKK(bind= 1,state='active')%IkBb(ikk=1, nfkb=None, loc='C', state='inactive'),kf,kr)
Rule('IKK_activate_IkBe', IKK(bind= None,state='active') + IkBe(ikk=None, nfkb=None, loc='C', state='inactive') <> IKK(bind= 1,state='active')%IkBe(ikk=1, nfkb=None, loc='C', state='inactive'),kf,kr)
Rule('IKK_activate_IkBd', IKK(bind= None,state='active') + IkBd(ikk=None, nfkb=None, loc='C', state='inactive') <> IKK(bind=1,state='active')%IkBd(ikk=1, nfkb=None, loc='C', state='inactive'), kf,kr)


Rule('IKKIkBa_active', IKK(bind= 1, state = 'active')%IkBa(ikk = 1, nfkb = None, loc = 'C', state = 'inactive') >> IKK(bind= None, state = 'active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'), kf,kr)
Rule('IKKIkBb_active', IKK(bind= 1, state = 'active')%IkBa(ikk = 1, nfkb = None, loc = 'C', state = 'inactive') >> IKK(bind= None, state = 'active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'), kf,kr)
Rule('IKKIkBe_active', IKK(bind= 1, state = 'active')%IkBa(ikk = 1, nfkb = None, loc = 'C', state = 'inactive') >> IKK(bind= None, state = 'active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'),kf,kr)
Rule('IKKIkBd_active', IKK(bind= 1, state = 'active')%IkBa(ikk = 1, nfkb = None, loc = 'C', state = 'inactive') >> IKK(bind= None, state = 'active') + IkBa(ikk=None, nfkb=None, loc='C', state='phos'), kf,kr)


Parameter('IkB_IKKf', 30)
Parameter('IkB_IKKr', 6e-5)
Rule('an_adc', IkBa(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBa(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)
Rule('bn_adc', IkBb(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBb(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)
Rule('en_adc', IkBe(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBe(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)
Rule('dn_adc', IkBd(ikk=None, nfkb=None, loc='C', state='active') + NFkB(ikb=None, loc='C') <> IkBd(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C'), IkB_IKKf, IkB_IKKr)

Rule('an_adn', IkBa(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBa(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)
Rule('bn_adn', IkBb(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBb(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)
Rule('en_adn', IkBe(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBe(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)
Rule('dn_adn', IkBd(ikk=None, nfkb=None, loc='N', state='active') + NFkB(ikb=None, loc='N') <> IkBd(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N'), IkB_IKKf, IkB_IKKr)


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
Rule('a_nc', IkBa(ikk=None, nfkb=None, loc='C', state='active') <> IkBa(ikk=None, nfkb=None, loc='N', state='active'), af, ancf)
Rule('b_nc', IkBb(ikk=None, nfkb=None, loc='C', state='active') <> IkBb(ikk=None, nfkb=None, loc='N', state='active'), bf, bncf)
Rule('e_nc', IkBe(ikk=None, nfkb=None, loc='C', state='active') <> IkBe(ikk=None, nfkb=None, loc='N', state='active'), ef, encf)
Rule('d_nc', IkBd(ikk=None, nfkb=None, loc='C', state='active') <> IkBd(ikk=None, nfkb=None, loc='N', state='active'), ef, dncf)

Parameter('anf', 0.276)
Parameter('bnf', 0.0276)
Parameter('enf', 0.138)
Parameter('dnf', 0.276)
Parameter('anr', 0.828)
Parameter('bnr', 0.414)
Parameter('enr', 0.414)
Parameter('dnr', 0.414)
Rule('an_nc', IkBa(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') <> IkBa(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active'), anf, anr)
Rule('bn_nc', IkBb(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') <> IkBb(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active'), bnf, bnr)
Rule('en_nc', IkBe(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') <> IkBe(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active'), enf, enr)
Rule('dn_nc', IkBd(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') <> IkBd(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active'), dnf, dnr)

Parameter('nf', 5.4)
Parameter('nr', 0.0048)
Rule('n_nc', NFkB(ikb=None, loc='C', state = 'active') <> NFkB(ikb=None, loc='N', state = 'active'), nf, nr)

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

Rule('ad_c', IkBa(ikk=None, nfkb=None, loc='C', state='phos') >> None, a_dc)
Rule('bd_c', IkBb(ikk=None, nfkb=None, loc='C', state='phos') >> None, b_dc)
Rule('ed_c', IkBe(ikk=None, nfkb=None, loc='C', state='phos') >> None, e_dc)
Rule('dd_c', IkBd(ikk=None, nfkb=None, loc='C', state='phos') >> None, d_dc)

Rule('ad_n', IkBa(ikk=None, nfkb=None, loc='N', state='active') >> None, a_dn)
Rule('bd_n', IkBb(ikk=None, nfkb=None, loc='N', state='active') >> None, b_dn)
Rule('ed_n', IkBe(ikk=None, nfkb=None, loc='N', state='active') >> None, e_dn)
Rule('dd_n', IkBd(ikk=None, nfkb=None, loc='N', state='active') >> None, d_dn)

Parameter('c_bn', 0.00006)
Parameter('n_bn', 0.00006)

Rule('an_c', IkBa(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') >> NFkB(ikb=None, loc='C', state = 'active'), c_bn)
Rule('bn_c', IkBb(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') >> NFkB(ikb=None, loc='C', state = 'active'), c_bn)
Rule('en_c', IkBe(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') >> NFkB(ikb=None, loc='C', state = 'active'), c_bn)
Rule('dn_c', IkBd(ikk=None, nfkb=1, loc='C', state='active')%NFkB(ikb=1, loc='C', state = 'active') >> NFkB(ikb=None, loc='C', state = 'active'), c_bn)

Rule('an_n', IkBa(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active') >> NFkB(ikb=None, loc='N', state = 'active'), n_bn)
Rule('bn_n', IkBb(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active') >> NFkB(ikb=None, loc='N', state = 'active'), n_bn)
Rule('en_n', IkBe(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active') >> NFkB(ikb=None, loc='N', state = 'active'), n_bn)
Rule('dn_n', IkBd(ikk=None, nfkb=1, loc='N', state='active')%NFkB(ikb=1, loc='N', state = 'active') >> NFkB(ikb=None, loc='N', state = 'active'), n_bn)

# def TNF_ComplexI_Observables():
Observable('TNFR_TNF', TNFR(blig=ANY))
Observable('TRADD_CompI', TRADD(brec = ANY))
Observable('RIP1_CompI', TNFR()%RIP1())
Observable('RIP1ub', RIP1(state='K63ub'))
Observable('TRAF_CompI', TRAF(brip=ANY)%RIP1(btraf=ANY))
Observable('Obs_TRAF_ComplexI', TRAF(brip=ANY)%RIP1(btraf=ANY))
Observable('Ubiquitylated_ComplexI', RIP1(state ='K63ub'))
Observable('MAPK_activity', MAPK(state = 'active'))
Observable('RIP1deub', RIP1(state = 'deub'))
Observable('Obs_ComplexII', TRADD(brec=None, brip=1)%RIP1(bscf=1, bub1=None, bub2=None, bub3=None, btraf=None, state='deub'))
Observable('TNFR_endocytosis', TNF(brec=3)%TNFR(blig = 3)%RIP1(bub1=None, bub2=None, bub3=None, btraf=None, state='deub'))
# Observable('TNF_endocytosed_', TNF_endocytosed())
Observable('CYLD_CompI_1', RIP1(bscf = ANY, bub2=2, bub3=None, state='K63ub')%CYLD(brip=2))
Observable('CYLD_CompI_2', RIP1(bscf = ANY, bub2=2, bub3=None, btraf=3, state='K63ub')%CYLD(brip=2)%TRAF(brip=3, state='K63ub'))
Observable('NFkBn_obs', NFkB(ikb=None, loc='N', state='active'))
Observable('NFkB_obs', NFkB(ikb=None, loc='C', state='active'))

Observable('IkBa_obs', IkBa(ikk=None, nfkb=None, loc='C', state='active'))
Observable('IkBb_obs', IkBb(ikk=None, nfkb=None, loc='C', state='active'))
Observable('IkBe_obs', IkBe(ikk=None, nfkb=None, loc='C', state='active'))
Observable('IkBd_obs', IkBd(ikk=None, nfkb=None, loc='C', state='active'))

Observable('IkBan_obs', IkBa(ikk=None, nfkb=None, loc='N', state='active'))
Observable('IkBbn_obs', IkBb(ikk=None, nfkb=None, loc='N', state='active'))
Observable('IkBen_obs', IkBe(ikk=None, nfkb=None, loc='N', state='active'))
Observable('IkBdn_obs', IkBd(ikk=None, nfkb=None, loc='N', state='active'))

Observable('IkBaNFkB_obs', IkBa(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C', state='active'))
Observable('IkBbNFkB_obs', IkBb(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C', state='active'))
Observable('IkBeNFkB_obs', IkBe(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C', state='active'))
Observable('IkBdNFkB_obs', IkBd(ikk=None, nfkb=1, loc='C', state='active') % NFkB(ikb=1, loc='C', state='active'))

Observable('IkBaNFkBn_obs', IkBa(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N', state='active'))
Observable('IkBbNFkBn_obs', IkBb(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N', state='active'))
Observable('IkBeNFkBn_obs', IkBe(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N', state='active'))
Observable('IkBdNFkBn_obs', IkBd(ikk=None, nfkb=1, loc='N', state='active') % NFkB(ikb=1, loc='N', state='active'))


from pysb.export import pysb_flat
exp = pysb_flat.PysbFlatExporter(model)

print(exp.export())

# generate_equations(model)
# generate_network(model)
#
# print("printing species length")
# print(len(model.species))
#
#
# tspan = np.linspace(0, 720, 721)
# x = odesolve(model,tspan,verbose=True)
#
# plt.figure()
# plt.plot(tspan/60, x['RIP1ub'], label="RIP1ub")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('RIPub')
#
# plt.figure()
# plt.plot(tspan, x['NFkBn_obs'], label="NFkBn")
# # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
# plt.xlabel("Time (in hr)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax = 0.00000008)
# plt.legend(loc=0)
# plt.title('NFkBn')
#
# plt.show()