from pysb import *
# from pysb.util import alias_model_components
# from earm import shared
#
# from earm import lopez_modules as lopez
# from earm import albeck_modules as albeck
import Necroptosis_module as irvin
import merge
import pickle
import numpy as np
import pylab as p
from numtools import calibratortools as ct
from numtools import simulator_1_0 as sim

from irvin_anrm_model import model
from pysb.integrate  import odesolve


# -----Monomers-----
def compile_monomers():
    Model('m')
    # From irvin_modules
    irvin.TNFa_to_ComplexI_Monomers()
    irvin.ComplexII_to_Bid_Monomers()
    irvin.NFkB_Activation_and_Signaling_monomers()
    # irvin.Bid_Hypothesis_monomers()
    irvin.Momomers_zVad_to_C8()

    # From lopez_modules
    # lopez.momp_monomers()
    #
    # # From albeck_modules
    # albeck.apaf1_to_parp_monomers()
    return m.monomers


revised_monomers = {'Bid': (['bf', 'state'], {'state': ['U', 'T', 'M', 'po4']}),
                    'C6': (['bf1', 'bf2', 'state'], {'state': ['pro', 'A']}),
                    'PARP': (['bf', 'state'], {'state': ['U', 'C', 'A']}),
                    'C3': (['bf', 'state'], {'state': ['pro', 'A', 'ub', 'I']})}

monomer_edits = merge.Edit_Monomers(compile_monomers(), revised_monomers)
merged_monomers = monomer_edits.merged_monomers

Model('model')
model.monomers = merged_monomers

irvin.TNFa_to_ComplexI_Initials()
irvin.TNFa_to_ComplexI()
irvin.CompI_TRADD_RIP1_Dissociation()
"""Hypothesis 1: FADD transiently localizes to TNFR1 to retrieve RIP1
    Hypothesis 2: FADD replaces TRADD in TRADD:RIP1
    Hypothesis 3: FADD binds TRADD in TRADD:RIP1"""
irvin.CompII_Hypothesis_1_FADD_CompI_interaction()
irvin.CompII_Hypothesis_2_FADD_displaces_TRADD()
irvin.CompII_Hypothesis_3_FADD_binds_TRADD()
irvin.ComplexII_to_Bid_Initials()
irvin.ComplexIIa_Assembly()
irvin.ComplexIIb_to_MLKL()
irvin.RIP1_truncation_ComplexII()
irvin.C8_catalyzed_truncations()
irvin.NFkB_Activation_and_Signaling_Initials()
irvin.NFkB_Activation_and_Signaling()

irvin.Bid_Hypothesis_initials()
"""Hypothesis 0: Bid-po4 does not occur
    Hypothesis 1: Bid-po4 recruits proC8 and cFlip_L and mediates RIP1 and CYLD truncation."""
irvin.Bid_Hypothesis_2()  # Bid mediated inhibition of necrosis (revised hypothesis)
irvin.Bid_proC8_cleaves_substrates()  # Bid mediated inhibition of necrosis (revised hypothesis)

irvin.C3_inhibits_MLKL()
irvin.Initials_zVad_to_C8()
irvin.zVad_to_C8()
irvin.observables()

# From lopez_modules
lopez.declare_initial_conditions()
lopez.translocate_tBid_Bax_BclxL()
lopez.tBid_activates_Bax_and_Bak()
lopez.tBid_binds_all_anti_apoptotics()
lopez.sensitizers_bind_anti_apoptotics()
lopez.effectors_bind_anti_apoptotics()
lopez.lopez_pore_formation(do_pore_transport=True)

# From irvin_modules
irvin.pore_to_parp()


#-----------Calibrated Parameters-----------------------
position = pickle.load(open('irvin_anrm_model_fitted_params.pkl'))

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,20000,1000)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-3
sims.atol = 1e-6

solve = sim.Solver(sims)
solve.run()

ic_params  = model.parameters_initial_conditions()

#-----------Apoptotic Conditions--------------------------
# apoptotic_conditions = ct.initial_conditions(['Bak_0', 'Bax_0', 'Bid_0', 'zVad_0', 'TNFa_0'], [0.2e5, 40165, 12044, 0, 1500], ic_params)
# y_apop = solve.simulate(position, observables=True, initial_conc = apoptotic_conditions)
# yout1 = ct.extract_records(y_apop, ['Obs_cPARP', 'Obs_MLKL', 'ComplexI', 'ComplexI_ub', 'TRADD_RIP1', 'RIP1_FADD', 'RIP1_Trunc', 'RIP1_po4'])

#-----------Necroptotic Conditions--------------------------
necroptotic_conditions = ct.initial_conditions(['Bak_0', 'Bax_0', 'Bid_0', 'zVad_0'], [0.2e5, 40165, 12044, 9.6e6], ic_params)
#20uM zVad == 9.6e6 zVad per cell for a cell volume of 8e-13L
y_necr = solve.simulate(position, observables=True, initial_conc = necroptotic_conditions)
yout2 = ct.extract_records(y_necr, ['Obs_cPARP', 'Obs_MLKL'])

#-----------Experimental Conditions--------------------------
apoptotic_conditions = ct.initial_conditions(['Bak_0', 'Bax_0', 'Bid_0', 'zVad_0', 'TNFa_0'], [0, 0, 0, 0, 1500], ic_params)
#20uM zVad == 9.6e6 zVad per cell for a cell volume of 8e-13L
y_necr = solve.simulate(position, observables=True, initial_conc = necroptotic_conditions)
yout3 = ct.extract_records(y_necr, ['Obs_cPARP', 'Obs_MLKL', 'RIP1_FADD'])

p.figure('Condition')
p.ion()
#p.plot(sims.tspan, yout3[:,2], label = 'RIP1:FADD')
p.plot(sims.tspan, yout3[:,1], label = 'MLKL')
p.legend()
p.show()