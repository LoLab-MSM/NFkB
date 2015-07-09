from pysb import *
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import matplotlib.pyplot as plt
import numpy as np
from sympy import sympify
from pysb.bng import generate_equations
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.util import alias_model_components

#Model ()

def declare_monomers():
    Monomer('IKKKa') #s1
    Monomer('IKKn') #s2
    Monomer('IKKa') #s3
    Monomer('IKKi') #s4
    Monomer('IkBa_p') #5
    Monomer('IkBapc_NFkBpc') #s6
    Monomer('NFkB_c') #s7
    Monomer('NFkBn') #s8
    Monomer('A20') #s9
    Monomer('A20t') #s10
    Monomer('IkBa') #s11
    Monomer('IkBan') #s12
    Monomer('IkBat') #s13
    Monomer('IkBa_NFkB') #s14
    Monomer('IkBan_NFkBn') #s15
    Monomer('TNFR1a') #s16
    Monomer('A20_on') #s17
    Monomer('IkBa_on') #s18
    Monomer('TNF_ext') #s19
    Monomer('IKKKn') #s20
    Monomer('IKKii') #21
    Monomer('TNFR1i') #22
    Monomer('A20_off') #23
    Monomer('IkBa_off') #24

def declare_parameters():
    Parameter('KNN', 2.0e5)
    Parameter('KN', 1.0e5)
    Parameter('M', 2000)
    Parameter('AN', 2)
    Parameter('ANa', 2)

    #Cell Parameters
    Parameter('kv', 5.0) #Nucealr to cytoplasm volume

    #Declaring Parameters
    Parameter('kb', 1.2e-5) #receptor activation rate
    Parameter('kf', .0012) #receptor inactivation rate
    Parameter('Tdeg', 7.7e-4) #TNF loss same as cdeg
    Parameter('ka',1e-5) #IKKK kinase activation rate
    Parameter('ka20', 1e5) #A20 TNFR1 block
    Parameter('ki', 0.01) #IKKK kinase inactivation rate
    Parameter('AR', 0) #active receptors
    Parameter('k1', 2.0*6e-10) #IKKn activation by IKKK
    Parameter('k3', 0.002) #IKKa inactivation by A20
    Parameter('k2', 10000) #IKKa inactivation
    Parameter('k3_div_k2', k3.value/k2.value) #for IKKa and A20
    Parameter('k4', 0.001) #IKKii transfer rate
    Parameter('AA', 1.0) #IkBa on (or off)
    Parameter('c0', 0.1) # inducible A20 and IkBa mRNA synthesis
    Parameter('c1a', 0.1) #inducible IkBa mRNA synthesis
    Parameter('AB', 1.0) #A20 (on or off)

    #IkBa Parameters
    Parameter('a2', 1e-7) #IkBa phosphorylation b/c IKKa
    Parameter('tp', 0.01) #degradation of phosph-IkBa complex with NFkB
    Parameter('a3', 5e-7) #IkBa_NFkB phosphorylation b/c IKKa
    Parameter('a1', 5e-7) #IkBa*NFkB association
    Parameter('c6a', 0.00002) #spontaneous IkBa_NFkB defg of IkBa complex to NFkB
    Parameter('c5a', 0.0001) #IkBa deg rate

    #Transport Parameters
    Parameter('i1a', 0.002) #IkBa nuclear import
    Parameter('e1a', 0.005) #IkBa nuclear export
    Parameter('e2a', 0.05) #IkBa_NFkB nuclear export
    Parameter('i1', 0.01) #NFkB nuclear import

    #A20 and IkBa synthesis Parameters
    Parameter('c4', 0.5) #A20, IkBa transformation rate
    Parameter('c5', 0.0005) #A20 degredation rate
    Parameter('c1', 0.1) #inducible A20 mRNA synthesis
    Parameter('G', 0) #initial status of A20 promoter
    Parameter('c3', 0.00075) #A20 and IkBa mRNA deg rate
    Parameter('q1', 4e-7) #NFkB attaching @ A20 and IkBa site
    Parameter('q2', 1e-6) #IkBa inducible detaching from A20, IkBa site

    #added parameters
    Parameter('a1_mult_kv', a1.value*kv.value) #for volume IkBa association NFkB

    #parameters for initial conditions
    Parameter('IKKKa_0', 0)
    Parameter('IKKn_0', 2e5)
    Parameter('IKKa_0', 0)
    Parameter('IKKi_0', 0)
    Parameter('IkBa_p_0', 0)
    Parameter('IkBapc_NFkBpc_0', 0)
    Parameter('NFkB_c_0', 0)
    Parameter('NFkBn_0', 1)
    Parameter('A20_0', 10000)
    Parameter('A20t_0', 10)
    Parameter('IkBa_0', 0.14*100000)
    Parameter('IkBan_0', 0.06*100000)
    Parameter('IkBat_0', 10)
    Parameter('IkBa_NFkB_0', 100000)
    Parameter('IkBan_NFkBn_0', 0)
    Parameter('TNFR1a_0', 0)
    Parameter('A20_on_0', 0)
    Parameter('IkBa_on_0', 0)
    Parameter('TNF_ext_0', 0.1) #1e-1
    Parameter('IKKKn_0', KN.value)
    Parameter('IKKii_0', KNN.value-IKKn_0.value)
    Parameter('TNFR1i_0', M.value)
    Parameter('A20_off_0', AN.value)
    Parameter('ANa_0', ANa.value)

    alias_model_components()

def declare_initial_conditions():
    Initial(IKKKa(), IKKKa_0)
    Initial(IKKn(), IKKn_0)
    Initial(IKKa(), IKKa_0)
    Initial(IKKi(), IKKi_0)
    Initial(IkBa_p(), IkBa_p_0)
    Initial(IkBapc_NFkBpc(), IkBapc_NFkBpc_0)
    Initial(NFkB_c(), NFkB_c_0)
    Initial(NFkBn(), NFkBn_0)
    Initial(A20(), A20_0)
    Initial(A20t(), A20t_0)
    Initial(IkBa(), IkBa_0)
    Initial(IkBan(), IkBan_0)
    Initial(IkBat(), IkBat_0)
    Initial(IkBa_NFkB(), IkBa_NFkB_0)
    Initial(IkBan_NFkBn(), IkBan_NFkBn_0)
    Initial(TNFR1a(), TNFR1a_0)
    Initial(A20_on(), A20_on_0)
    Initial(IkBa_on(), IkBa_on_0)
    Initial(TNF_ext(), TNF_ext_0) #1e-1
    Initial(IKKKn(), IKKKn_0)
    Initial(IKKii(), IKKii_0)
    Initial(TNFR1i(), TNFR1i_0)
    Initial(A20_off(), A20_off_0)
    Initial(IkBa_off(), ANa_0)

def declare_observables():
    Observable('IKKKa_obs', IKKKa()) #
    Observable('IKKKn_obs', IKKKn())
    Observable('IKKa_obs', IKKa()) #
    Observable('IKKi_obs', IKKi()) #
    Observable('IKKii_obs', IKKii()) #
    Observable('TNF_ext_obs', TNF_ext()) #
    Observable('IkBa_p_obs', IkBa_p()) #
    Observable('IkBapc_NFkBpc_obs', IkBapc_NFkBpc()) #
    Observable('NFkB_c_obs', NFkB_c())
    Observable('NFkBn_obs', NFkBn())
    Observable('TNFR1a_obs', TNFR1a())
    Observable('TNFR1i_obs', TNFR1i())
    Observable('A20_off_obs', A20_off())
    Observable('IkBa_off_obs', IkBa_off())
    Observable('A20_on_obs', A20_on())
    Observable('IkBa_on_obs', IkBa_on())
    Observable('A20t_obs', A20t())
    Observable('IkBa_obs', IkBa())
    Observable('IkBan_obs', IkBan())
    Observable('IkBat_obs', IkBat())
    Observable('IkBa_NFkB_obs', IkBa_NFkB())
    Observable('IkBan_NFkBn_obs', IkBan_NFkBn())
    Observable('IKKn_obs', IKKn()) #
    Observable('A20_obs', A20())

def declare_functions():
    Expression('keff', sympify("ka*ka20/(ka20+A20_obs)")) #10000 #michaelis menten


#Declaring rules
def declare_rules():
    Rule('IKKKa_to_IKKKn', IKKKa() >> IKKKn(), ki) #IKKKa creates KN (total amount of IKKK kinase molecules)
    Rule('IKKKa_and_IKKn', IKKKa() + IKKKa() + IKKn() >> IKKKa() + IKKKa() + IKKa(), k1) #proportion of IKKKa changes IKKn to IKKa
    Rule('IKKa_to_IKKi', IKKa() >> IKKi(), k3) #IKKa active to IKKi inactive
    Rule('IKKa_and_A20', IKKa() + A20() >> A20() + IKKi(), k3_div_k2) #Exp1) #A20 mediated IKKa to IKKi
    Rule('IKKi_to_IKKii', IKKi() >> IKKii(), k4) #
    Rule('IkBa_p_deg', IkBa_p() >> None, tp)
    Rule('IkBapc_NFkBpc_to_NFkB_c', IkBapc_NFkBpc() >> NFkB_c(), tp)
    Rule('IkBa_NFkB_to_NFkB_c', IkBa_NFkB() >> NFkB_c(), c6a)
    Rule('NFkB_c_to_NFkBn', NFkB_c() >> NFkBn(), i1)
    Rule('NFkB_c_and_IkBa', NFkB_c() + IkBa() >> IkBa_NFkB(), a1)
    Rule('NFkBn_and_IkBan', NFkBn() + IkBan() >> IkBan_NFkBn(), a1_mult_kv)
    Rule('A20t_to_A20', A20t() >> A20t() + A20(), c4)
    Rule('A20_deg', A20() >> None, c5)
    Rule('TNFR1a_and_IKKKn', TNFR1a() + IKKKn() >> TNFR1a() + IKKKa(), keff)
    Rule('A20_on_to_A20t', A20_on() >> A20_on() + A20t(), c1)
    Rule('A20t_deg', A20t() >> None, c3)
    Rule('IKKa_and_IkBa', IKKa() + IkBa() >> IKKa() + IkBa_p(), a2)
    Rule('IkBat_to_IkBa', IkBat() >> IkBat() + IkBa(), c4)
    Rule('IkBa_deg', IkBa() >> None, c5a)
    Rule('IkBa_to_IkBan', IkBa() >> IkBan(), i1a)
    Rule('IkBan_to_IkBa', IkBan() >> IkBa(), e1a)
    Rule('IkBa_on_to_IkBat', IkBa_on() >> IkBa_on() + IkBat(), c1a)
    Rule('IkBat_deg', IkBat() >> None, c3)
    Rule('IKKa_and_IkBa_NFkB', IKKa() + IkBa_NFkB() >> IKKa() + IkBapc_NFkBpc(), a3)
    Rule('IkBan_NFkBn_to_IkBa_NFkB', IkBan_NFkBn() >> IkBa_NFkB(), e2a)
    Rule('TNFR1a_to_TNFR1i', TNFR1a() >> TNFR1i(), kf)
    Rule('IkBan_and_A20_on', IkBan() + A20_on() >> IkBan() + A20_off(), q2)
    Rule('IkBan_and_IkBa_on', IkBan() + IkBa_on() >> IkBan() + IkBa_off(), q2)
    Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg)
    Rule('IKKii_to_IKKn', IKKii() >> IKKn(), k4)
    Rule('TNF_ext_and_TNFR1i', TNF_ext() + TNFR1i() >> TNF_ext() + TNFR1a(), kb)
    Rule('NFkBn_and_A20_off', NFkBn() + A20_off() >> NFkBn() + A20_on(), q1)
    Rule('NFkBn_and_IkBa_off', NFkBn() + IkBa_off() >> NFkBn() + IkBa_on(), q1)
