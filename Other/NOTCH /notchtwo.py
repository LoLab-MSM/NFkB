from pysb.core import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import ScipyOdeSimulator as SOS
import pandas as pd

Model()

Monomer('L_ex', ['r_in'])
Monomer('R_in', ['l_ex', 'ysc'])
Monomer('YSC', ['lr'])
Monomer('ICD', ['csl', 'l'], {'l': ['cyto', 'nuc']})
Monomer('CSL', ['icd', 'gh', 'gsl'])
Monomer('gHes1', ['ic', 'ph'])
Monomer('mHes1', ['l'], {'l': ['cyto', 'nuc']})
Monomer('pHES1', ['gh', 'ga1', 'gn2'])
Monomer('Endo_on')
Monomer('Endo_off')
Monomer('Endo_off_L_ex')
Monomer('Endo_off_R_ex')
Monomer('gSlug', ['ic'])
Monomer('mSlug', ['l'], {'l': ['nuc', 'cyto']})
Monomer('pSLUG')
Monomer('gAscl1', ['ph'])
Monomer('mAscl1', ['l'], {'l': ['nuc', 'cyto']})
Monomer('pASCL1')
Monomer('gNgn2', ['ph'])
Monomer('mNgn2', ['l'], {'l': ['nuc', 'cyto']})
Monomer('pNGN2', ['gl_in'])
Monomer('gL_in', ['n2'])
Monomer('mL_in', ['l'], {'l': ['nuc', 'cyto']})
Monomer('pL_in', ['r_ex'])
Monomer('R_ex', ['l_in', 'ysc'])

# Defining Observables
Observable('ICD_obs', ICD(csl=None, l='nuc'))
Observable('CSL_on', ICD(csl=1, l='nuc') % CSL(icd=1, gh=None, gsl=None))
Observable('Ex_L', L_ex(r_in=None))
Observable('In_R', R_in(l_ex=None, ysc = None))
Observable('ExL_InR', L_ex(r_in=1) % R_in(l_ex=1))
Observable('LR_on', R_in(l_ex=1, ysc = None) % L_ex(r_in=1))
Observable('CSL_off', CSL(icd=None, gh=None, gsl=None))
Observable('gHes1_Hes1', gHes1(ph=1) % pHES1(gh=1))
Observable('mHes1_Hes1', mHes1(l = 'nuc'))
Observable('protein_HES1', pHES1(gh=None, ga1=None, gn2=None))
Observable('Slug_levels', pSLUG())
Observable('Ascl1_levels', pASCL1())
Observable('In_L', pL_in(r_ex=None))
Observable('protNGN2', pNGN2(gl_in=None))

# Defining initial conditions and parameters
Parameter('L_ex_0', 100)
Parameter('R_in_0', 200)
Parameter('YSC_0', 10)
Parameter('CSL_0', 50)
Parameter('gHes1_0', 1)
Parameter('gAscl1_0', 1)
# Parameter('pASCL1_0', 200)
Parameter('gNgn2_0', 1)
# Parameter('pNGN2_0', 20)
Parameter('gSlug_0', 1)
# Parameter('gL_in_0', 1)
Parameter('pL_in_0', 5)

Initial(L_ex(r_in=None), L_ex_0)
Initial(R_in(l_ex=None, ysc=None), R_in_0)
Initial(YSC(lr=None), YSC_0)
Initial(CSL(icd=None, gh=None, gsl=None), CSL_0)
Initial(gHes1(ic=None, ph=None), gHes1_0)
Initial(gSlug(ic=None), gSlug_0)
Initial(gAscl1(ph=None), gAscl1_0)
# Initial(pASCL1(), pASCL1_0)
Initial(gNgn2(ph=None), gNgn2_0)
# Initial(pNGN2(gl_in=None), pNGN2_0)
# Initial(gL_in(n2=None), gL_in_0)
Initial(pL_in(r_ex=None), pL_in_0)


#synthesis of Hes1
# Parameter('a', 4.5)
# Parameter('b', 0.23)
# Parameter('c', 0.23)
# Parameter('k', 33)
# Parameter('pcrit', 40)
# Parameter('hill', 2)


#synthesis of Hes1
# Parameter('a', 0.1)
# Parameter('b', 0.0)
# Parameter('c', 0.15)
# Parameter('k', 0.015)
# Parameter('pcrit', 40)
# Parameter('hill', 2)
# Parameter('gp_kf', 0.01)
# Parameter('gp_kr', 1.0)
# Parameter('gm_kf', 0.1)
# Parameter('mp_kf', 0.01)
#
# Expression('m_exp', b/(1.0 + (protein_HES1/25.0)**2))
# # Expression('p_exp', a/(0.031+ 0.031*protein_HES1**2))
# Expression('p_exp', a/(0.031+ 0.031*protein_HES1**2))
#
# # Rule('g_to_p', gHes1(ic = None, ph = None) + pHES1(gh = None, ga1 = None, gn2 = None) <> gHes1(ic = None, ph = 1) % pHES1(gh = 1, ga1 = None, gn2 = None), gp_kf, gp_kr)
# Rule('g_to_m', None >>  mHes1(), m_exp)
# Rule('m_to_p', mHes1() >> pHES1() + mHes1(), a)
# Rule('hes1_mrna_deg', mHes1() >> None, c)
# Rule('hes1_deg', pHES1() >> None, k)
# Expression('hes1_exp', k/((protein_HES1)**hill + 1))
# Expression('hes1mrna_exp', a*mHes1_Hes1)
# Rule('hes1_mrna', None >> mHes1(l = 'nuc'), hes1_exp)
# Rule('hes1_prot', None >> pHES1(gh = None, ga1 = None, gn2 = None), hes1mrna_exp)


# Define perturbations by their time points



# perturbation_params = {
#     15: {'b': 1.0},
#     # 30: {'a': 0.1}, #
#     # 40: {'c': 0.15},  # 0.15
#     50: {'k': 0.051}, # 0.051 50 #180
#
# }

# Expression('hes1_exp', k/((protein_HES1/pcrit)**hill + 1))
#
# Rule('hes1_mrna', None >> mHes1(l = 'nuc'), hes1_exp)
# Rule('hes1_prot', None >> pHES1(gh = None, ga1 = None, gn2 = None), a)
# Rule('hes1_mrna_deg', mHes1(l = 'nuc') >> None, c)
# Rule('hes1_deg', pHES1(gh = None, ga1 = None, gn2 = None) >> None, c)

# Notch Ligand binding and unbinding to Notch Receptor; Gamma SC
Parameter('kLRf', 1e-2)
Parameter('kLRr', 1e-4)
Parameter('kLRYSCf', 1e-2)
Parameter('kLRYSCr', 1e-5)
Parameter('kYSCICD', 1e-2)

Rule('LR_binding', L_ex(r_in=None) + R_in(l_ex=None, ysc = None) <> L_ex(r_in=1) % R_in(l_ex=1, ysc = None), kLRf, kLRr)
Rule('YSC_LR_binding', YSC(lr=None) + L_ex(r_in=1) % R_in(l_ex=1, ysc = None) <>
     YSC(lr=2) % L_ex(r_in=1) % R_in(l_ex=1, ysc = 2), kLRYSCf, kLRYSCr)
Rule('ICD_formation', YSC(lr=2) % L_ex(r_in=1) % R_in(l_ex=1, ysc = 2) >>
     ICD(csl=None, l='cyto') , kYSCICD)

# Activation of CSL
Parameter('kICD_c_n', 1e-1)
Parameter('kICDCSLf', 1e-4)
Parameter('kICDCSLr', 1e-2)
Rule('ICD_cyto_to_nucleus', ICD(csl = None, l='cyto') >> ICD(csl = None, l='nuc'), kICD_c_n)
Rule('ICD_CSL_complex', ICD(csl=None, l='nuc') + CSL(icd=None) <> ICD(csl=1, l='nuc') % CSL(icd=1), kICDCSLf, kICDCSLr)



Parameter('a', 0.1)
Parameter('b', 1.0)
Parameter('c', 0.15)
Parameter('k', 0.051)
Parameter('pcrit', 40)
Parameter('hill', 2)
Parameter('gp_kf', 0.01)
Parameter('gp_kr', 1.0)
Parameter('gm_kf', 0.1)
Parameter('mp_kf', 0.01)

Expression('m_exp', b/(1.0 + (protein_HES1/25.0)**2))
Expression('p_exp', a/(0.031+ 0.031*protein_HES1**2))
# Expression('hes1_exp', k/((protein_HES1/pcrit)**hill + 1))

Rule('hes1_mrna', ICD(csl=1, l='nuc') % CSL(icd=1) >> mHes1(l = 'nuc') ICD(csl=1, l='nuc') % CSL(icd=1), m_exp)
Rule('hes1_prot', mHes1(l = 'nuc') >> pHES1(gh = None, ga1 = None, gn2 = None) + mHes1(l = 'nuc'), a)
Rule('hes1_mrna_deg', mHes1(l = 'nuc') >> None, c)
Rule('hes1_deg', pHES1(gh = None, ga1 = None, gn2 = None) >> None, k)

#
# Expression('m_exp', b/(1.0 + (protein_HES1/25.0)**2))
# Expression('p_exp', a/(0.031+ 0.031*protein_HES1**2))

# Rule('g_to_p', gHes1(ic = None, ph = None) + pHES1(gh = None, ga1 = None, gn2 = None) <> gHes1(ic = None, ph = 1) % pHES1(gh = 1, ga1 = None, gn2 = None), gp_kf, gp_kr)
# Rule('g_to_m', None >>  mHes1(), m_exp)
# Rule('m_to_p', mHes1() >> pHES1() + mHes1(), a)
# Rule('hes1_mrna_deg', mHes1() >> None, c)
# Rule('hes1_deg', pHES1() >> None, k)
# Expression('hes1_exp', k/((protein_HES1)**hill + 1))
# Expression('hes1mrna_exp', a*mHes1_Hes1)
# Rule('hes1_mrna', None >> mHes1(l = 'nuc'), hes1_exp)
# Rule('hes1_prot', None >> pHES1(gh = None, ga1 = None, gn2 = None), hes1mrna_exp)


# Define perturbations by their time points
# perturbation_params = {
#     0: {'b': 1.0},
    # 30: {'a': 0.1}, #
    # 40: {'c': 0.15}, # 0.15
    # 50: {'k': 0.051}, # 0.051 50 #180
# }





# Parameter('kfhes1', 0.02)
# Parameter('hes1mrna', 3e-3)
# Parameter('n', 3)
# Parameter('hes1_deg', 1e-2)
# Expression('khes1', CSL_on**n/(kfhes1**n + CSL_on**n))
# # Rule('mrna_hes1', None >> gHes1(ph = None, ic = None), hes1mrna)
# Rule('synth_hes1', None >> pHES1(gh = None, ga1 = None, gn2 = None), khes1)
# Rule('protein_HES1_decay', pHES1() >> None, hes1_deg)
# HEs1 sysnthesis

# Parameter('kgICf', 1e-3)
# Parameter('kgICr', 1e-3)
# Parameter('kmH_n_c', 1e-1)

Parameter('ktrcxn', 1e-1)
Parameter('ktrlxn', 1e-2)

# Hes1 self inhibition
Parameter('kgpHf', 1)
Parameter('kgpHr', 1e-3)

# Rule('pHES1_inhibits_gHes1_1', gHes1(ph=None) + pHES1(gh=None) <> gHes1(ph=1) % pHES1(gh=1), kgpHf, kgpHr)

# SLUG synthesis
Parameter('kfslug', 0.02)
Parameter('n', 3)
Expression('kslug', CSL_on**n/(kfslug**n + CSL_on**n))
Rule('synth_slug', None >> pSLUG(), kslug)


# ASCL1 synthesis

Parameter('kfascl1', 0.02)
Parameter('n2', 3)
Expression('kascl1', 1.0/(1.0 + protein_HES1**n2))
Rule('synth_ascl1', None >> pASCL1(), kascl1)


# NGN2 Synthesis
Parameter('kfngn2', 0.02)
# Parameter('n', 3)
Expression('kngn2', 1.0/(1.0 + protein_HES1**n2))
Rule('synth_ngn2', None >> pNGN2(gl_in=None), kngn2)


# Intrinsic Ligand Synthesis
Parameter('kflin', 0.02)
Parameter('nl', 2)
Expression('kLin', protNGN2**n/(kflin**n + protNGN2**n))
Rule('synth_Lin', None >> pL_in(r_ex=None), kLin)


# Decay of species
Parameter('d_L', 1e-2)
Parameter('d_R', 1e-2)
Parameter('d_ICD', 1e-2)
Parameter('d_mH', 1e-3)
Parameter('d_pH', 1e-5)
Parameter('d_protein', 0.01)
Parameter('d_mRNA', 1e-3)

Rule('Ligand_decay', L_ex(r_in=None) >> None, d_L)
Rule('Receptor_decay', R_in(l_ex=None, ysc = None) >> None, d_R)
Rule('ICD_decay', ICD(l='cyto') >> None, d_ICD)
# Rule('mRNA_Hes1_decay', mHes1(l='cyto') >> None, d_mH)
# Rule('protein_HES1_decay', pHES1() >> None, d_pH)
Rule('pSlug_decay', pSLUG() >> None, d_protein)
# Rule('mSlug_decay', mSlug() >> None, d_mRNA)
Rule('pAscl1_decay', pASCL1() >> None, d_protein)
# Rule('mAscl1_decay', mAscl1() >> None, d_mRNA)
Rule('pNgn2_decay', pNGN2() >> None, d_protein)
# Rule('mNGN2_decay', mNgn2() >> None, d_mRNA)
Rule('pL_in_decay', pL_in() >> None, d_protein)
# Rule('mL_in_decay', mL_in() >> None, d_mRNA)

# Recyling of Ligand and Receptor in Endosomes

# Parameter('rR_in', 1e-2)
# Parameter('rL_ex', 1e-1)
# Parameter('rR_ex', 1e-12)
# Parameter('kEndosplit', 1e-1)
#
# Rule('Recycle_R_in', Endo_on() >> R_in(l_ex=None, ysc=None), rR_in)  # rate is medium
#
# Rule('Endo_off_split', Endo_off() >> Endo_off_L_ex() + Endo_off_R_ex(), kEndosplit)
# Rule('Recycle_Ligand_Endo_L_to_L', Endo_off_L_ex() >> L_ex(r_in=None), rL_ex)  # rate is v high
# Rule('Recycle_Receptor_Endo_L_to_R', Endo_off_R_ex() + ICD(l='cyto') >> R_ex(l_in=None, ysc=None),
#      rR_ex)  # rate is v low

# species_dict = {0:'L',
#                1:'R',
#                2:'YSC',
#                3:'CSL',
#                4:'gHES1_off',
#                5:'__source()',
#                6:'LR',
#                7:'__sink()',
#                8:'LRYSC',
#                9:'ICDc',
#               10:'ICDn',
#               11:'ICDnCSL',
#               12:'gHes1ICDCSL_off',
#               13:'CSLgHES1_off',
#               14:'gHes1ICDCSL_on',
#               15:'CSLgHES1_on',
#               16:'mHES1n_off',
#               17:'mHES1c_off',
#               18:'mHES1c_on',
#               19:'pHES1',
#               20:'gHES1pHES1_off',
#               21:'ICDnCSLgHES1pHES1_off',
#               22:'CSLgHES1pHES1_off'
# }

# tspan = np.linspace(0, 720, 7201)
# sim1 = SOS(model, tspan)
# sim_result = sim1.run()
#


tspan = np.linspace(0, 360, 3601)
sim = ScipyOdeSimulator(model, tspan = tspan)
sim_result = sim.run()
# df = run_sim_with_perturbations(sim, tspan, perturbation_params)



#
# for i,sp in enumerate(model.species):
#     print i,":",sp
#
# print(len(model.species))
# print(len((model.parameters)))
# print(len(model.odes))
# print(len(model.rules))
# # for  j,ode in enumerate(model.odes):
# #    for i in range(len(model.species)):
# #       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
# #    print j,":",ode
#
# # print('Active CSL complex levels')
# # print(sim_result.observables(['CSL_on']))
# print('EL')
# print(sim_result.observables['Ex_L'])
# # print('mRNA_Hes1')
# # print(sim_result.observables['mRNA_Hes1'])
# # print('protein_HES1')
# print(sim_result.observables['protein_HES1'])
# print('Slug_levels')
# print(sim_result.observables['Slug_levels'])
# print('Ascl1_levels')
# print(sim_result.observables['Ascl1_levels'])


# plt.figure(figsize=(15, 10))
# plt.subplot(231)
# plt.plot(tspan, df['Ex_L'], label='Extrinsic Ligand')
# plt.plot(tspan, df['In_R'], label='Intrinsic Receptor')
# plt.plot(tspan, df['ExL_InR'], label='exLinR')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)
#
# plt.subplot(232)
# plt.plot(tspan, df['ICD_obs'], label='Nuclear ICD')
# plt.plot(tspan, df['CSL_on'], label='ICDCSL')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)
#
# plt.subplot(233)
# plt.plot(tspan, df['mHes1_Hes1'], label='mHes1_Hes1')
# plt.plot(tspan, df['protein_HES1'], label='HES1')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)
#
# plt.subplot(234)
# plt.plot(tspan, df['Slug_levels'], label='SLUG')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)
#
# plt.subplot(235)
# plt.plot(tspan, df['In_L'], label='Intrinsic Ligand')
# plt.plot(tspan, df['protNGN2'], label='NGN2 protien')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)
#
# plt.subplot(236)
# plt.plot(tspan, df['Ascl1_levels'], label='ASCL1')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)

plt.figure(figsize=(15, 10))
plt.subplot(231)
plt.plot(tspan/60, sim_result.observables['Ex_L'], label='Extrinsic Ligand')
plt.plot(tspan/60, sim_result.observables['In_R'], label='Intrinsic Receptor')
plt.plot(tspan/60, sim_result.observables['ExL_InR'], label='exLinR')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(232)
plt.plot(tspan/60, sim_result.observables['ICD_obs'], label='Nuclear ICD')
plt.plot(tspan/60, sim_result.observables['CSL_on'], label='ICDCSL')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(233)
plt.plot(tspan/60, sim_result.observables['mHes1_Hes1'], label='mHes1_Hes1')
plt.plot(tspan/60, sim_result.observables['protein_HES1'], label='HES1')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(234)
plt.plot(tspan/60, sim_result.observables['Slug_levels'], label='SLUG')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(235)
plt.plot(tspan/60, sim_result.observables['In_L'], label='Intrinsic Ligand')
plt.plot(tspan/60, sim_result.observables['protNGN2'], label='NGN2 protien')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(236)
plt.plot(tspan/60, sim_result.observables['Ascl1_levels'], label='ASCL1')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)



plt.show()