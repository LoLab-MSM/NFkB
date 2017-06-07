from pysb.core import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import ScipyOdeSimulator as SOS

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
Observable('In_R', R_in(l_ex=None))
Observable('LR_on', R_in(l_ex=1) % L_ex(r_in=1))
Observable('CSL_off', CSL(icd=None, gh=None, gsl=None))
Observable('mRNA_Hes1', mHes1(l='cyto'))
Observable('protein_HES1', pHES1(gh=None, ga1=None, gn2=None))
Observable('Slug_levels', pSLUG())
Observable('Ascl1_levels', pASCL1())
Observable('In_L', pL_in(r_ex=None))
Observable('protNGN2', pNGN2(gl_in=None))

# Defining initial conditions and parameters
Parameter('L_ex_0', 100)
Parameter('R_in_0', 500)
Parameter('YSC_0', 10)
Parameter('CSL_0', 10)
Parameter('gHes1_0', 1)
Parameter('gAscl1_0', 1)
Parameter('pASCL1_0', 200)
Parameter('gNgn2_0', 1)
Parameter('pNGN2_0', 20)
Parameter('gSlug_0', 1)
Parameter('gL_in_0', 1)
Parameter('pL_in_0', 5)

Initial(L_ex(r_in=None), L_ex_0)
Initial(R_in(l_ex=None, ysc=None), R_in_0)
Initial(YSC(lr=None), YSC_0)
Initial(CSL(icd=None, gh=None, gsl=None), CSL_0)
Initial(gHes1(ic=None, ph=None), gHes1_0)
Initial(gSlug(ic=None), gSlug_0)
Initial(gAscl1(ph=None), gAscl1_0)
Initial(pASCL1(), pASCL1_0)
Initial(gNgn2(ph=None), gNgn2_0)
Initial(pNGN2(gl_in=None), pNGN2_0)
Initial(gL_in(n2=None), gL_in_0)
Initial(pL_in(r_ex=None), pL_in_0)

# Notch Ligand binding and unbinding to Notch Receptor; Gamma SC

Parameter('kLRf', 1e-2)
Parameter('kLRr', 1e-2)
Parameter('kLRYSCf', 1e-3)
Parameter('kLRYSCr', 1e-3)
Parameter('kYSCICD', 1e-2)

Rule('LR_binding', L_ex(r_in=None) + R_in(l_ex=None) <> L_ex(r_in=1) % R_in(l_ex=1), kLRf, kLRr)
Rule('YSC_LR_binding', YSC(lr=None) + R_in(l_ex=ANY, ysc=None) <> YSC(lr=1) % R_in(l_ex=ANY, ysc=1), kLRYSCf, kLRYSCr)
Rule('ICD_formation', R_in(ysc=ANY) >> ICD(csl=None, l='cyto') + YSC(lr=None) + Endo_on() + Endo_off(), kYSCICD)

# Activation of CSL

Parameter('kICD_c_n', 1e-1)
Parameter('kICDCSLf', 1e-4)
Parameter('kICDCSLr', 1e-4)
Rule('ICD_cyto_to_nucleus', ICD(l='cyto') >> ICD(l='nuc'), kICD_c_n)
Rule('ICD_CSL_complex', ICD(csl=None, l='nuc') + CSL(icd=None) <> ICD(csl=1, l='nuc') % CSL(icd=1), kICDCSLf, kICDCSLr)

# HEs1 sysnthesis

Parameter('kgICf', 1e-3)
Parameter('kgICr', 1e-3)
Parameter('kmH_n_c', 1e-1)

Parameter('ktrcxn', 1e-1)
Parameter('ktrlxn', 1e-2)

Rule('gHes1_binds_ICDCSL', gHes1(ic=None) + CSL(icd=ANY, gh=None) <> gHes1(ic=1) % CSL(icd=ANY, gh=1), kgICf, kgICr)
Rule('gHes1_IC_complex_to_mRNA_Hes1', gHes1(ic=ANY, ph=None) >> mHes1(l='nuc') + gHes1(ic=ANY, ph=None), ktrcxn)
Rule('mRNA_Hes1_nuc_to_cyto', mHes1(l='nuc') >> mHes1(l='cyto'), kmH_n_c)
Rule('mRNA_Hes1_to_Protein_HES1', mHes1(l='cyto') >> pHES1(gh=None, ga1=None, gn2=None) + mHes1(l='cyto'), ktrlxn)

# Hes1 self inhibition

Parameter('kgpHf', 1e-3)
Parameter('kgpHr', 1e-3)

Rule('pHES1_inhibits_gHes1_1', gHes1(ph=None) + pHES1(gh=None) <> gHes1(ph=1) % pHES1(gh=1), kgpHf, kgpHr)

# SLUG synthesis

Parameter('kICDCSLSlugf', 1e-3)
Parameter('kICDCSLSlugr', 1e-5)
Parameter('kSlug_n_to_c', 1e-1)

Rule('geneSlug_binds_ICDCSL', gSlug(ic=None) + CSL(icd=ANY, gsl=None) <> gSlug(ic=1) % CSL(icd=ANY, gsl=1),
     kICDCSLSlugf, kICDCSLSlugr)
Rule('gSlug_ICDCSL_complex_forms_mRNA', gSlug(ic=ANY) >> gSlug(ic=ANY) + mSlug(l='nuc'), ktrcxn)
Rule('mSLug_to_cyto', mSlug(l='nuc') >> mSlug(l='cyto'), kSlug_n_to_c)
Rule('mRNASlug_to_protein', mSlug(l='cyto') >> mSlug(l='cyto') + pSLUG(), ktrlxn)

# ASCL1 synthesis

Parameter('kHes1Ascl1f', 1e-3)
Parameter('kHes1Ascl1r', 1e-5)
Parameter('kmAscl1_n_to_c', 1e-1)

Rule('geneAscl1_bind_proteinHES1', gAscl1(ph=None) + pHES1(ga1=None) <> gAscl1(ph=1) % pHES1(ga1=1), kHes1Ascl1f,
     kHes1Ascl1r)
Rule('geneAscl1_to_mRNA_basal', gAscl1(ph=None) >> gAscl1(ph=None) + mAscl1(l='nuc'), ktrcxn)
Rule('mAscl1_n_to_c', mAscl1(l='nuc') >> mAscl1(l='cyto'), kmAscl1_n_to_c)
Rule('mAscl1_to_protein', mAscl1(l='cyto') >> mAscl1(l='cyto') + pASCL1(), ktrlxn)

# NGN2 Synthesis

Parameter('kHes1Ngn2f', 1e-3)
Parameter('kHes1Ngn2r', 1e-5)
Parameter('kmNgn2_n_to_c', 1e-1)

Rule('geneNgn2_bind_Hes1', gNgn2(ph=None) + pHES1(gn2=None) <> gNgn2(ph=1) % pHES1(gn2=1), kHes1Ngn2f, kHes1Ngn2r)
Rule('geneNgn2_to_mRNA', gNgn2(ph=None) >> mNgn2(l='nuc') + gNgn2(ph=None), ktrcxn)
Rule('mNgn2_n_to_c', mNgn2(l='nuc') >> mNgn2(l='cyto'), kmNgn2_n_to_c)
Rule('mNgn2_to_protein', mNgn2(l='cyto') >> pNGN2(gl_in=None), ktrlxn)

# Intrinsic Ligand Synthesis

Parameter('kNgn2L_inf', 1e-4)
Parameter('kNgn2L_inr', 1e-4)
Parameter('kmL_in_n_to_c', 1e-1)

Rule('gLig_in_binds_to_Ngn2', gL_in(n2=None) + pNGN2(gl_in=None) <> gL_in(n2=1) % pNGN2(gl_in=1), kNgn2L_inf,
     kNgn2L_inr)
Rule('gLinNgn2complex_to_mRNA', gL_in(n2=ANY) >> gL_in(n2=ANY) + mL_in(l='nuc'), ktrcxn)
Rule('mLin_n_to_c', mL_in(l='nuc') >> mL_in(l='cyto'), kmL_in_n_to_c)
Rule('mL_in_to_protein', mL_in(l='cyto') >> pL_in(r_ex=None) + mL_in(l='cyto'), ktrlxn)

# Decay of species
Parameter('d_L', 1e-5)
Parameter('d_R', 1e-5)
Parameter('d_ICD', 1e-5)
Parameter('d_mH', 1e-3)
Parameter('d_pH', 1e-5)
Parameter('d_protein', 1e-5)
Parameter('d_mRNA', 1e-3)

Rule('Ligand_decay', L_ex(r_in=None) >> None, d_L)
Rule('Receptor_decay', R_in(l_ex=None) >> None, d_R)
Rule('ICD_decay', ICD(l='cyto') >> None, d_ICD)
Rule('mRNA_Hes1_decay', mHes1(l='cyto') >> None, d_mH)
Rule('protein_HES1_decay', pHES1() >> None, d_pH)
Rule('pSlug_decay', pSLUG() >> None, d_protein)
Rule('mSlug_decay', mSlug() >> None, d_mRNA)
Rule('pAscl1_decay', pASCL1() >> None, d_protein)
Rule('mAscl1_decay', mAscl1() >> None, d_mRNA)
Rule('pNgn2_decay', pNGN2() >> None, d_protein)
Rule('mNGN2_decay', mNgn2() >> None, d_mRNA)
Rule('pL_in_decay', pL_in() >> None, d_protein)
Rule('mL_in_decay', mL_in() >> None, d_mRNA)

# Recyling of Ligand and Receptor in Endosomes

Parameter('rR_in', 1e-2)
Parameter('rL_ex', 1e-1)
Parameter('rR_ex', 1e-12)
Parameter('kEndosplit', 1e-1)

Rule('Recycle_R_in', Endo_on() >> R_in(l_ex=None, ysc=None), rR_in)  # rate is medium

Rule('Endo_off_split', Endo_off() >> Endo_off_L_ex() + Endo_off_R_ex(), kEndosplit)
Rule('Recycle_Ligand_Endo_L_to_L', Endo_off_L_ex() >> L_ex(r_in=None), rL_ex)  # rate is v high
Rule('Recycle_Receptor_Endo_L_to_R', Endo_off_R_ex() + ICD(l='cyto') >> R_ex(l_in=None, ysc=None),
     rR_ex)  # rate is v low

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

tspan = np.linspace(0, 10000000, 1001)
sim1 = SOS(model, tspan)
sim_result = sim1.run()

# for  j,ode in enumerate(model.odes):
#    for i in range(len(model.species)):
#       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
#    print j,":",ode

# print('Active CSL complex levels')
# print(sim_result.observables(['CSL_on']))
print('mRNA_Hes1')
print(sim_result.observables['mRNA_Hes1'])
print('protein_HES1')
print(sim_result.observables['protein_HES1'])
print('Slug_levels')
print(sim_result.observables['Slug_levels'])
print('Ascl1_levels')
print(sim_result.observables['Ascl1_levels'])

plt.figure(figsize=(5, 15))
plt.subplot(511)
plt.plot(tspan, sim_result.observables['Ex_L'], label='Extrinsic Ligand')
plt.plot(tspan, sim_result.observables['In_R'], label='Intrinsic Receptor')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(512)
plt.plot(tspan, sim_result.observables['ICD_obs'], label='Nuclear ICD')
plt.plot(tspan, sim_result.observables['CSL_on'], label='Active CSL')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(513)
plt.plot(tspan, sim_result.observables['mRNA_Hes1'], label='mRNA Hes1')
plt.plot(tspan, sim_result.observables['protein_HES1'], label='HES1')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(514)
plt.plot(tspan, sim_result.observables['Slug_levels'], label='SLUG')
plt.plot(tspan, sim_result.observables['Ascl1_levels'], label='ASCL1')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(515)
plt.plot(tspan, sim_result.observables['In_L'], label='Intrinsic Ligand')
plt.plot(tspan, sim_result.observables['protNGN2'], label='NGN2 protien')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)
plt.show()