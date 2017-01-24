import numpy
from pysb import *
from pysb.util import alias_model_components
from pysb.macros import *

def TNFa_to_ComplexI_Monomers():
    """ Declares TNFa, TNFR1, TRADD, RIP1, TRAF2, IAP, NKS (NFkB signaling complex),
    NFkB, CYLD and FADD. Binding sites are named after (when known) subdomains that
    facilitate binding. Otherwise the binding site will reflect what it binds with.

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """
    Monomer('TNFa', ['brec'])  # TNFa
    Monomer('TNFR1', ['blig', 'bDD'])  # TNFR1
    Monomer('TRADD', ['bDD1', 'bDD2'])  # TRADD
    Monomer('RIP1', ['bDD', 'btraf', 'bRHIM', 'bMLKL', 'state'], {'state': ['unmod', 'ub', 'po4', 'trunc']})  # RIP1
    Monomer('TRAF', ['brip', 'bciap', 'zfinger'])  # TRAF2
    Monomer('cIAP', ['btraf'])  # cIAP 1 and 2
    Monomer('NSC', ['bnfkb'])  # Recruitement of LUBAC commits complex I to the NFkB activation pathway. This pathway is approximated by a single activation step which is carried out by a NFkB signaling complex (NSC).
    Monomer('NFkB', ['bnsc', 'state'], {'state': ['I', 'A']})  # NFkB
    Monomer('CYLD', ['btraf', 'state'], {'state': ['U', 'T']})  # CYLD
    Monomer('FADD', ['bDD', 'bDED1', 'bDED2'])  # FADD
    alias_model_components()


def TNFa_to_ComplexI_Initials():
    Parameter('TNFa_0', 600)  # 6000 corresponds to 100ng/ml TNFa
    Parameter('TNFR1_0', 4800)  # 4800 receptors per cell
    Parameter('TRADD_0', 9000)  # molecules per cell (arbitrarily assigned) 9000
    Parameter('RIP1_0', 12044)  # molecules per cell (arbitrarily assigned) 12044
    Parameter('TRAF_0', 9000)  # molecules per cell (arbitrarily assigned) 9000
    Parameter('cIAP_0', 9000)  # molecules per cell (arbitrarily assigned) 9000
    Parameter('NSC_0', 0)  # complexes per cell
    Parameter('NFkB_0', 50000)  # molecules per cell (arbitrarily assigned) 50000
    Parameter('CYLD_0', 9000)  # molecules per cell
    Parameter('FADD_0', 8030)  # 8300 molecules per cell (J Immunol 2005)
    alias_model_components()

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


def TNFa_to_ComplexI():
    """Reaction network that produces Complex I and activates NFkB [1]. Recruitment of LUBAC and NFkB pathway components to complex I is approximated by a one-step conversion to "NFkB Signaling Complex" (reaction 7). A20 ubiquitylates RIP1 and targets it for degradation. The action of A20 is apporixmated at a single reaction (reaction 8).

        1. TNFa + TNFR1 <> TNFa:TNFR1
        2. TNFa:TNFR1 + TRADD <> TNFa:TNFR1:TRADD
        3. TNFa:TNFR1:TRADD + RIP1(unmod) <> TNFa:TNFR1:TRADD:RIP1(unmod)
        4. TNFa:TNFR1:TRADD:RIP1(unmod) + TRAF2 <> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2
        5. TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2 + cIAP <> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:cIAP

        6. TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:cIAP >>TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:cIAP
        7. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 >> NFkB Signaling Complex
        8. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 >> TNFa:TNFR1:TRADD + TRAF2
        9. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 + CYLD <> TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:CYLD
        10.TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:CYLD >> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:CYLD

        1. Olivier Micheau, Jurg Tschopp, Induction of TNF Receptor I-Mediated Apoptosis via Two Sequential Signaling Complexes, Cell, Volume 114, Issue 2, 25 July 2003, Pages 181-190
    """
    bind(TNFa(brec=None), 'brec', TNFR1(blig=None), 'blig', [1e-6, 1e-3])
    bind(TNFR1(blig=ANY, bDD=None), 'bDD', TRADD(bDD1=None, bDD2=None), 'bDD1', [1e-6, 1e-3])
    bind_complex(TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=None), 'bDD2',
                 RIP1(bDD=None, btraf=None, state='unmod'), 'bDD', [1e-6, 1e-3])
    bind_complex(TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=None, state='unmod'),
                 'btraf', TRAF(brip=None), 'brip', [1e-6, 1e-3])
    # we need to get bind_complex working!
    bind(TRAF(bciap=None), 'bciap', cIAP(btraf=None), 'btraf', [1e-6, 1e-3])
    bind(TRAF(zfinger=None), 'zfinger', CYLD(btraf=None), 'btraf', [1e-6, 1e-3])

    ComplexI = TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY)

    # No one step conversions with complexes
    Rule('RIP1_Ubiquitination',
         ComplexI % RIP1(bDD=ANY, btraf=ANY, state='unmod') % TRAF(brip=ANY, bciap=ANY) >> ComplexI % RIP1(bDD=ANY,
                                                                                                           btraf=None,
                                                                                                           state='ub') % TRAF(
             brip=ANY, bciap=ANY), Parameter('RIP1_ubiq_kc', 1e-1))
    Rule('RIP1_Deubiquitination',
         ComplexI % RIP1(bDD=ANY, btraf=ANY, state='ub') % TRAF(brip=ANY, zfinger=ANY) >> ComplexI % RIP1(bDD=ANY,
                                                                                                          btraf=None,
                                                                                                          state='unmod') % TRAF(
             brip=ANY, zfinger=ANY), Parameter('RIP1_deub_kc', 1e-1))
    Rule('RIP1_Degradation',
         ComplexI % RIP1(bDD=ANY, btraf=ANY, state='ub') % TRAF(brip=ANY) >> ComplexI + TRAF(brip=None),
         Parameter('RIP1_degr_kc', 1e-1))
    Rule('Establish_NFkB_Signaling_Complex',
         ComplexI % RIP1(bDD=ANY, btraf=ANY, state='ub') % TRAF(brip=ANY, bciap=ANY) >> NSC(bnfkb=None),
         Parameter('NSC_esta_kc', 1e-8))


def CompI_TRADD_RIP1_Dissociation():
    """Some believe that TRADD, RIP1 and others are released into the cytoplasm after RIP1 deubiquitination or after TNFR endocytosis. This mechanism presents a problem: The RIP1 released into the cytoplasm is in theory the same as RIP1 that existed in the cytoplasm prior to TNF stimulation. Meaning apoptosis whould spontaneously occur. Dissociation of TRADD:RIP1 from complex I is a hypothetical pathway that could distinguish pre- and post- TNF RIP1. It could also serve as a mechanism for FADD independent necrosome formation. This hypothesis is supported by a suggested mechanism, to explain signal transduction from TNFa to Riptosome formation[1].

        TNFa:TNFR1:TRADD:RIP1(unmod) <> TNFa:TNFR1 + TRADD:RIP1(unmod)

        1. Dickens, LS., IR Powley, MA Hughes, M MacFarlane, Exp. Cell. Res. 318 (2012) 1269-1277"""
    bind_complex(TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=None), 'bDD',
                 TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY, state='unmod'), 'bDD1', [1e-6, 1e-3])


def CompII_Hypothesis_1_FADD_CompI_interaction():
    """Michaeu and Tschopp, showed that FADD is absent from Complex I[1]. The presence of TRADD in Complex II is debated [2]. As explained in CompI_TRADD_RIP1_Dissociation(), simply releasing RIP1 from Complex I so that it can bind FADD creates a condition where cell death spontaneously occurs. Here, I hypothesize that FADD transiently associates with Complex I in order to retrieve RIP1 from Complex I. This association, being transient, would have been difficult to detect.

        1. Olivier Micheau, Jurg Tschopp, Induction of TNF Receptor I-Mediated Apoptosis via Two Sequential Signaling Complexes, Cell, Volume 114, Issue 2, 25 July 2003, Pages 181-190

        2. Dickens, LS., IR Powley, MA Hughes, M MacFarlane, "The 'complexities' of life and death: Death receptor signalling platforms" Exp. Cell. Res. 318 (2012) 1269-1277"""

    Rule('TNFR1_to_FADD_1',
         TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=3) % RIP1(bDD=3, btraf=4) % TRAF(brip=4) + FADD(bDD=None) >> TNFR1(
             bDD=2) % TRADD(bDD1=2, bDD2=None) + RIP1(bDD=1, btraf=None) % FADD(bDD=1) + TRAF(brip=None),
         Parameter('TNFR1_FADD_kc_1', 1e-1))

    Rule('TNFR1_to_FADD_2',
         TNFR1(bDD=2) % TRADD(bDD1=2, bDD2=3) % RIP1(bDD=3, btraf=None) + FADD(bDD=None) >> TNFR1(bDD=2) % TRADD(bDD1=2,
                                                                                                                 bDD2=None) + RIP1(
             bDD=1, btraf=None) % FADD(bDD=1), Parameter('TNFR1_FADD_kc_2', 1e-1))


def CompII_Hypothesis_2_FADD_displaces_TRADD():
    """Dickens et al, speculated that TRADD:RIP1 is released into the cytoplasm following RIP1 deubiquitination by CYLD[1]. Furter, Dinckens question presence of TRADD in Complex II[1], but it is established that Complex II is contains FADD. To model these findings, we suggest that FADD replaces the TRADD in the cytoplasmic TRADD:RIP1 complex. """

    Rule('TRADD_to_FADD_1',
         TRADD(bDD1=None, bDD2=2) % RIP1(bDD=2) % TRAF(brip=3) + FADD(bDD=None) >> TRADD(bDD1=None, bDD2=None) + RIP1(
             bDD=1) % FADD(bDD=1) + TRAF(brip=None), Parameter('TRADD_FADD_kc_1', 1e-1))

    Rule('TRADD_to_FADD_2',
         TRADD(bDD1=None, bDD2=2) % RIP1(bDD=2, btraf=None) + FADD(bDD=None) >> TRADD(bDD1=None, bDD2=None) + RIP1(
             bDD=1, btraf=None) % FADD(bDD=1), Parameter('TRADD_FADD_kc_2', 1e-1))


def CompII_Hypothesis_3_FADD_binds_TRADD():
    """Dickens et al, speculated that TRADD:RIP1 is released into the cytoplasm following RIP1 deubiquitination by CYLD[1]. FADD is recruited to RIP1 to form a cytoplasmic complex II. The presence of TRADD in complex II is debated. Here, we hypothesize that it is present. """

    bind_complex(FADD(bDD=None), 'bDD', TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY, btraf=None), 'bDD1', [1e-1, 1e-1])


def ComplexII_to_Bid_Monomers():
    """ Bid is declared in Albeck Modules"""
    Monomer('proC8', ['bDED'])
    Monomer('C8', ['bf', 'state'], {'state': ['A', 'I']})  # active caspase 8
    Monomer('flip_L', ['bDED'])
    Monomer('flip_S', ['bDED'])
    Monomer('RIP3', ['bRHIM', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
    Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})


def ComplexII_to_Bid_Initials():
    Parameter('flip_L_0', 39023)  # 39022 molecules per cell (J Immunol 2005)
    Parameter('flip_S_0', 39023)  # molecules per cell assummed both isoforms w/conc.
    Parameter('proC8_0', 16057)  # procaspase 8 molecules per cell 16057 (J Immunol 2005)
    Parameter('C8_0', 0)  # active caspase 8 dimers per cell.
    Parameter('RIP3_0', 2.0e4)  # molecules per cell
    Parameter('MLKL_0', 1.0e6)  # molecules per cell

    alias_model_components()
    Initial(RIP3(bRHIM=None, state='unmod'), RIP3_0)  # RIP3
    Initial(flip_L(bDED=None), flip_L_0)
    Initial(flip_S(bDED=None), flip_S_0)
    Initial(proC8(bDED=None), proC8_0)
    Initial(C8(bf=None, state='I'), C8_0)
    Initial(MLKL(bRHIM=None, state='unmod'), MLKL_0)


def ComplexIIa_Assembly():
    """Defines the interactions of procaspase 8 and/or cFlip recruitment to FADD as per
    ANRM 1.0.

    Uses FADD, proC8 and flip_L and flip_S monomers and parameters to generate the rules.

    FADD is a scaffold protein that gains the ability to recruit cFlip and procaspase 8 when
    it's death domain DD is occupied. Fas, Rip1 and TRADD are capable of binding FADD.
    """

    # ================================================
    # c-Flip and procaspase 8 recruitment to FADD
    # ------------------------------------------------
    #   X:FADD:proC8 + proC8  <-> X:FADD:proC8:proC8
    #   X:FADD:proC8 + flip_L <-> X:FADD:proC8:flip_L
    #   X:FADD:proC8 + flip_S <-> X:FADD:proC8:flip_S
    #
    #   X = RIP1, TRADD

    # ----------procaspase8 and cFlip recruitment-----
    bind(FADD(bDD=ANY, bDED1=None), 'bDED1', proC8(bDED=None), 'bDED', [7.27e-06, 0.018])  # (J Immunol 2005)
    Rule('bind_FADD_proC8_2', FADD(bDD=ANY, bDED2=None) + proC8(bDED=None) <> FADD(bDD=ANY, bDED2=1) % proC8(bDED=1),
         Parameter('bind_FADD_proC8_2_kf', 7.27e-06), Parameter('bind_FADD_proC8_2_kr',
                                                                0.018))  # (J Immunol 2005) FIX: Bind function fails when one molecule can bind another at more than one binding site. The failure occurs because _rule_name_generic does not distinguish between binding sites. NAMING MACRO SHOULD number repeated names.

    bind(FADD(bDD=ANY, bDED1=None), 'bDED1', flip_L(bDED=None), 'bDED', [7.27e-05, 0.018])  # (J Immunol 2005)
    Rule('bind_FADD_flip_L_2', FADD(bDD=ANY, bDED2=None) + flip_L(bDED=None) <> FADD(bDD=ANY, bDED2=1) % flip_L(bDED=1),
         Parameter('bind_FADD_flip_L_2_kf', 7.27e-06), Parameter('bind_FADD_flip_L_2_kr',
                                                                 0.018))  # (J Immunol 2005)[...source] says cFlip-L has greater affinity for proC8 than proC8 itself.

    bind(FADD(bDD=ANY, bDED1=None), 'bDED1', flip_S(bDED=None), 'bDED', [7.27e-06, 0.018])  # (J Immunol 2005)
    Rule('bind_FADD_flip_S_2', FADD(bDD=ANY, bDED2=None) + flip_S(bDED=None) <> FADD(bDD=ANY, bDED2=1) % flip_S(bDED=1),
         Parameter('bind_FADD_flip_S_2_kf', 7.27e-06), Parameter('bind_FADD_flip_S_2_kr', 0.018))  # (J Immunol 2005)


def ComplexIIb_to_MLKL():
    """Necrosome formation and MLKL activation"""
    bind(RIP1(bDD=ANY, bRHIM=None, state='unmod'), 'bRHIM', RIP3(bRHIM=None, state='unmod'), 'bRHIM', [1e-6, 1e-3])

    Rule('Rip3_PO4lation',
         RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state='unmod') % RIP3(
             bRHIM=ANY, state='po4'), Parameter('k19', 1e-2))
    Rule('Rip1_PO4lation',
         RIP1(bRHIM=ANY, state='unmod') % RIP3(bRHIM=ANY, state='po4') >> RIP1(bRHIM=ANY, state='po4') % RIP3(bRHIM=ANY,
                                                                                                              state='po4'),
         Parameter('k20', 1e-3))

    catalyze_state(RIP1(state='po4'), 'bMLKL', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-6, 1e-3, 1e-1])
    catalyze_state(MLKL(state='active'), 'bRHIM', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-7, 0.4, 0.01])


def observables():
    Observable('Obs_TNFa', TNFa(brec=None))
    Observable('ComplexI', TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(state='unmod'))
    Observable('ComplexI_ub', TNFa(brec=ANY) % TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(state='ub'))
    Observable('ComplexI_TRAF',
               TNFR1(blig=ANY, bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=1, state='unmod') % TRAF(
                   brip=1))
    Observable('TRADD_RIP1', TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY))
    Observable('TRADD_RIP1_2', TRADD(bDD1=None, bDD2=ANY) % RIP1(bDD=ANY, btraf=None))
    Observable('ComplexII', FADD(bDD=ANY) % TRADD(bDD1=ANY, bDD2=ANY) % RIP1(bDD=ANY, btraf=None))
    Observable('ComplexIIa', FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Obs_NFkB', NFkB(state='A'))
    Observable('Obs_FADD_Sole', FADD(bDD=None))
    Observable('Bid_Trunc', Bid(state='T'))
    Observable('Obs_Bid', Bid(bf=None, state='U'))
    Observable('Bid_PO4', Bid(state='po4'))
    Observable('Obs_RIP1', RIP1(state='unmod') + RIP1(state='po4') + RIP1(state='ub'))
    Observable('RIP1_Trunc', RIP1(state='trunc'))
    Observable('RIP3_Trunc', RIP3(state='trunc'))
    Observable('Necrosome', RIP1(bRHIM=ANY, state='po4') % RIP3(bRHIM=ANY, state='po4'))
    Observable('Obs_proC8', proC8())
    Observable('Obs_C8', C8())
    Observable('Obs_C3ub', C3(state='ub'))
    Observable('Obs_C3', C3(state='A'))
    Observable('Obs_pC3', C3(state='pro'))
    Observable('RIP1_FADD', FADD(bDD=ANY) % RIP1(bDD=ANY))
    # Observable('RIP1_nucl', RIP1(bDD = None, bRHIM=ANY, state = 'N')%RIP3(bRHIM=ANY, state = 'N'))
    Observable('Obs_cPARP', PARP(state='C'))
    Observable('Obs_PARP', PARP(state='U'))
    Observable('Obs_MLKL', MLKL(state='active'))
    Observable('Obs_CytoC', CytoC(state='C') + CytoC(state='A'))
    Observable('RIP1_po4', RIP1(state='po4'))
    Observable('PCD_Markers', PARP(state='C') + MLKL(state='active'))
    Observable('Obs_CYLD', CYLD(state='U'))
    Observable('Obs_tCYLD', CYLD(state='T'))