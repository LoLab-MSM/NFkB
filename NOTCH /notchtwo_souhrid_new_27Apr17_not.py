from pysb.core import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import ScipyOdeSimulator as SOS
from pysb.bng import generate_equations


Model ()


Monomer('L',['r','loc'], {'loc':['cell_1','cell_2']})
Monomer('R_ex', ['l','r_tm'])
Monomer('R_tm', ['r_ex', 'r_icd','ysc'])
Monomer('R_ICD', ['r_tm','csl'])
Monomer('YSC', ['lr'])
Monomer('L_endo')
Monomer('CSL', ['icd'])
Monomer('mHes1')
Monomer('pHES1')
Monomer('pSLUG')
Monomer('pASCL1')
Monomer('pNGN2')

#print model.monomers

#Compartments

Parameter('VolECM', 1000)
Parameter('VappCM', 100)
Parameter('VCyto1', 400)
Parameter('VCyto2', 400)
Parameter('VappNM', 40)
Parameter('VNuc1', 50)
Parameter('VNuc2', 50)

Compartment('ECM', None, 3, VolECM)
Compartment('CellMem1',ECM, 2, VappCM)
Compartment('CellMem2',ECM, 2, VappCM)
Compartment('Cyto1',CellMem1,3,VCyto1)
Compartment('Cyto2',CellMem2,3,VCyto2)
Compartment('NM1',Cyto1, 2, VappNM)
Compartment('NM2',Cyto2,2,VappNM)
Compartment('Nuc1',NM1,3,VNuc1)
Compartment('Nuc2',NM2,3,VNuc2)

 
#print model.compartments

#Defining Observables

Observable('L_1', L(r=ANY,loc='cell_1')**ECM)
Observable('L_2', L(r=ANY,loc='cell_2')**ECM)
        
Observable('R_1', R_tm(r_ex=ANY,r_icd=ANY)**CellMem1)
Observable('R_2', R_tm(r_ex=ANY,r_icd=ANY)**CellMem2)

Observable('R_ICD_1', R_ICD()**Nuc1)
Observable('R_ICD_2', R_ICD()**Nuc2)

Observable('CSL_on_1', CSL(icd = ANY)**Nuc1)          
Observable('CSL_on_2', CSL(icd = ANY)**Nuc2)    
       
Observable('mRNA_Hes1_1', mHes1()**Cyto1)
Observable('mRNA_Hes1_2', mHes1()**Cyto2)

Observable('prot_Hes1_1', pHES1()**Nuc1)
Observable('prot_Hes1_2', pHES1()**Nuc2)

Observable('Slug_levels_1', pSLUG()**Cyto1)
Observable('Slug_levels_2', pSLUG()**Cyto2)

Observable('Ascl1_levels_1', pASCL1()**Cyto1)
Observable('Ascl1_levels_2', pASCL1()**Cyto2)

Observable('NGN2_levels_1', pNGN2()**Nuc1)
Observable('NGN2_levels_2', pNGN2()**Nuc2)

print model.observables

#Defining initial conditions and parameters

Parameter('L_1_0',100)
Parameter('L_2_0',10)
Parameter('R_1_0', 200)
Parameter('R_2_0', 200)
Parameter('YSC_0', 10)
Parameter('CSL_0', 50)
Parameter('pASCL1_0', 200)
Parameter('pNGN2_0', 100)

Parameter('gSlug_0', 1)
Parameter('gL_0', 1)
Parameter('gHes1_0', 1)
Parameter('gAscl1_0', 1)
Parameter('gNgn2_0', 1)

Initial(L(r=None,loc='cell_1')**ECM, L_1_0)
Initial(L(r=None,loc='cell_2')**ECM, L_2_0)
Initial(R_ex(l=None,r_tm=1)**ECM % R_tm(r_ex=1,ysc=None,r_icd=2)**CellMem1 % R_ICD(r_tm=2,csl=None)**Cyto1, R_1_0)
Initial(R_ex(l=None,r_tm=1)**ECM % R_tm(r_ex=1,ysc=None,r_icd=2)**CellMem2 % R_ICD(r_tm=2,csl=None)**Cyto2, R_2_0)
Initial(YSC(lr = None)**CellMem1, YSC_0)
Initial(YSC(lr = None)**CellMem2, YSC_0)
Initial(CSL(icd = None)**Nuc1, CSL_0)
Initial(CSL(icd = None)**Nuc2, CSL_0)
Initial(pASCL1()**Cyto1, pASCL1_0)
Initial(pASCL1()**Cyto2, pASCL1_0)
Initial(pNGN2()**Cyto1, pNGN2_0)
Initial(pNGN2()**Cyto2, pNGN2_0)


#for ic in model.initial_conditions:
#    print ic


#Notch Ligand binding and unbinding to Notch Receptor; Gamma SC

Parameter('kRform', 1e-1)
Parameter('kLform', 1e-1)
Parameter('kLRf', 1e-2)
Parameter('kLRr', 1e-2)
Parameter('kLRYSCf', 1e-3)
Parameter('kLRYSCr', 1e-3)  
Parameter('kYSCICD', 1e-2)

#Ligand of one cell binds receptor of another

Rule('L_1_binds_R_2', L(r=None,loc='cell_1') **ECM + R_ex(l=None, r_tm=1)**ECM % R_tm(r_ex=1,r_icd=ANY,ysc=None)**CellMem2 <> 
         L(r=2,loc='cell_1') ** ECM %  R_ex(l=2, r_tm=1) ** ECM % R_tm(r_ex=1, r_icd=ANY, ysc=None) ** CellMem2 , kLRf,kLRr)


Rule('L_2_binds_R_1', L(r=None,loc='cell_2') **ECM + R_ex(l=None, r_tm=1)**ECM % R_tm(r_ex=1,r_icd=ANY,ysc=None)**CellMem1 <> 
         L(r=2,loc='cell_2') ** ECM %  R_ex(l=2, r_tm=1) ** ECM % R_tm(r_ex=1, r_icd=ANY, ysc=None) ** CellMem1 , kLRf,kLRr)


Rule('YSC_binding', YSC(lr=None) + R_tm(r_ex=1,ysc=None,r_icd=ANY) % R_ex(l=ANY, r_tm=1) <>
    YSC(lr=2) % R_tm(r_ex=1,ysc=2,r_icd=ANY) % R_ex(l=ANY,r_tm=1), kLRYSCf,kLRYSCr)


Rule('ICD_formation_1',L(r=1, loc='cell_2') ** ECM % R_ICD(r_tm=2, csl=None) ** Cyto1 % R_ex(l=1, r_tm=3) ** ECM % R_tm(r_ex=3, r_icd=2, ysc=4) ** CellMem1 % YSC(lr=4) ** CellMem1
     >> R_ICD(csl=None,r_tm=None)**Cyto1 + YSC(lr=None)**CellMem1 + L_endo()**Cyto2, kYSCICD)
Rule('ICD_formation_2',L(r=1, loc='cell_1') ** ECM % R_ICD(r_tm=2, csl=None) ** Cyto2 % R_ex(l=1, r_tm=3) ** ECM % R_tm(r_ex=3, r_icd=2, ysc=4) ** CellMem2 % YSC(lr=4) ** CellMem2
     >> R_ICD(csl=None,r_tm=None)**Cyto2 + YSC(lr=None)**CellMem2 + L_endo()**Cyto1, kYSCICD)

# Activation of CSL

Parameter('kICD_c_n', 1e-1)
Parameter('kICDCSLf', 1e-4)
Parameter('kICDCSLr', 1e-4)
Rule('ICD_cyto_to_nucleus_1', R_ICD(csl=None,r_tm=None)**Cyto1 >> R_ICD(csl=None,r_tm=None)**Nuc1, kICD_c_n)
Rule('ICD_cyto_to_nucleus_2', R_ICD(csl=None,r_tm=None)**Cyto2 >> R_ICD(csl=None,r_tm=None)**Nuc2, kICD_c_n)

Rule('ICD_CSL_complex', R_ICD(csl = None,r_tm=None) + CSL(icd = None) <> R_ICD(csl = 1,r_tm=None)%CSL(icd = 1), kICDCSLf, kICDCSLr)

##synthesis of Hes1

Parameter('ktrcxn', 1e-1)
Parameter('ktrlxn', 1e-2)
Parameter('kmrnadeg', 1e-3)
Parameter('km_n_to_c', 1e-1)
Parameter('ktf_c_to_n', 1e-1)
Parameter('Kprotsynth', ktrcxn.value*ktrlxn.value*km_n_to_c.value/kmrnadeg.value)
Parameter('kgICf', 1e-2)
Parameter('kgICr', 2e-4)
Parameter('kgpHf', 1e-2)
Parameter('kgpHr', 2e-4)
Parameter('kgICpHf', 1e-2)
Parameter('kgICpHr', 2e-4)
Parameter('kgpHICf', 1e-2)
Parameter('kgpHICr',2e-4)
Parameter('hes1_deg', 1e-3)


Expression('kHes1_1', 
           ktrcxn * (kgICpHr/kgICpHf) * CSL_on_1/
           (
            (kgICr/kgICf) * (kgICpHr/kgICpHf) +
            (kgICpHr/kgICpHf) * CSL_on_1 +
            (kgpHICr/kgpHICf) * prot_Hes1_1   +
            prot_Hes1_1*CSL_on_1
            ))

Expression('kHes1_2', 
           ktrcxn * (kgICpHr/kgICpHf) * CSL_on_2/
           (
            (kgICr/kgICf) * (kgICpHr/kgICpHf) +
            (kgICpHr/kgICpHf) * CSL_on_2 +
            (kgpHICr/kgpHICf) * prot_Hes1_2   +
            prot_Hes1_2*CSL_on_2
            ))

Rule('mHes1_synth_1', None >> mHes1()**Cyto1, kHes1_1)
Rule('mHes1_synth_2', None >> mHes1()**Cyto2, kHes1_2)

Parameter('km_to_p', ktrlxn.value* km_n_to_c.value * ktf_c_to_n.value)

Rule('pHes1_synth_1', mHes1()**Cyto1 >> pHES1()**Nuc1,km_to_p)
Rule('pHes1_synth_2', mHes1()**Cyto2 >> pHES1()**Nuc2,km_to_p)

Rule('mHes1_decay', mHes1() >> None, kmrnadeg)
Rule('pHES1_decay', pHES1() >> None, hes1_deg)


## SLUG synthesis

Parameter('kgSlICf', 1e-2)
Parameter('kgSlICr', 2e-4)

Parameter('n', 3)
Expression('kSlug_1', Kprotsynth*CSL_on_1**n/((kgSlICr/kgSlICf)**n + CSL_on_1**n))
Expression('kSlug_2', Kprotsynth*CSL_on_2**n/((kgSlICr/kgSlICf)**n + CSL_on_2**n))
#
Rule('synth_Slug_1', None >> pSLUG()**Cyto1, kSlug_1)
Rule('synth_Slug_2', None >> pSLUG()**Cyto2, kSlug_2)

# ASCL1 synthesis

Parameter('kgAs1pHf', 1e-2)
Parameter('kgAs1pHr', 2e-4)


Expression('kAscl1_1', Kprotsynth*(kgAs1pHr/kgAs1pHf)**n/((kgAs1pHr/kgAs1pHf)**n + prot_Hes1_1**n))
Expression('kAscl1_2', Kprotsynth/((kgAs1pHr/kgAs1pHf)**n + prot_Hes1_2**n))
#
Rule('synth_Ascl1_1', None >> pASCL1()**Cyto1, kAscl1_1)
Rule('synth_Ascl1_2', None >> pASCL1()**Cyto2, kAscl1_2)



# NGN2 Synthesis
Parameter('kgNg2pHf', 1e-2)
Parameter('kgNg2pHr', 2e-4)
#
Expression('kNgn2_1', Kprotsynth*ktf_c_to_n/((kgNg2pHr/kgNg2pHf)**n + prot_Hes1_1**n))
Expression('kNgn2_2', Kprotsynth*ktf_c_to_n/((kgNg2pHr/kgNg2pHf)**n + prot_Hes1_2**n))
#
Rule('synth_Ngn2_1', None >> pNGN2()**Nuc1, kNgn2_1)
Rule('synth_Ngn2_2', None >> pNGN2()**Nuc2, kNgn2_2)
#

#Ligand Synthesis

Parameter('kgLNg2f', 1e-2)
Parameter('kgLNg2r', 2e-4)
#
Expression('kLig_1', Kprotsynth*NGN2_levels_1**n/((kgLNg2r/kgLNg2f)**n + NGN2_levels_1**n))
Expression('kLig_2', Kprotsynth*NGN2_levels_2**n/((kgLNg2r/kgLNg2f)**n + NGN2_levels_2**n))
#
Rule('synth_Lig_1', None >> L(r=None,loc='cell_1')**ECM, kLig_1)
Rule('synth_Lig_2', None >> L(r=None,loc='cell_2')**ECM, kLig_2)
#


## Decay of species

Parameter('d_L', 1e-2)
Parameter('d_R', 1e-4)
Parameter('d_ICD', 1e-1)
Parameter('d_mH', 1e-3)
Parameter('d_pH', 1e-5)
Parameter('d_protein', 0.01)
Parameter('d_mRNA', 1e-3)
#
Rule('Ligand_decay', L() >> None, d_L)
Rule('Receptor_decay', (R_ex(r_tm=1) % R_tm(r_ex=1,r_icd=2) % R_ICD(r_tm=2)) >> None, d_R)
Rule('ICD_decay', R_ICD(r_tm=None) >> None, d_ICD)
Rule('pSlug_decay', pSLUG() >> None, d_protein)
Rule('pAscl1_decay', pASCL1() >> None, d_protein)
Rule('pNgn2_decay', pNGN2() >> None, d_protein)

## Recyling of Ligand and Receptor in Endosomes

Parameter('kRsynth', 1)
Parameter('krL', 1)

Rule('R_synthesis_1', None >> R_ex(l=None,r_tm=1)**ECM % R_tm(r_ex=1,ysc=None,r_icd=2)**CellMem1 % R_ICD(r_tm=2,csl=None)**Cyto1, kRsynth)
Rule('R_synthesis_2', None >> R_ex(l=None,r_tm=1)**ECM % R_tm(r_ex=1,ysc=None,r_icd=2)**CellMem2 % R_ICD(r_tm=2,csl=None)**Cyto2, kRsynth)
#
Rule('L_recycle_1', L_endo()**Cyto1 >> L(r=None,loc='cell_1')**ECM,krL)
Rule('L_recycle_2', L_endo()**Cyto2 >> L(r=None,loc='cell_2')**ECM,krL)
#
## species_dict = {0:'L',
##                1:'R',
##                2:'YSC',
##                3:'CSL',
##                4:'gHES1_off',
##                5:'__source()',
##                6:'LR',
##                7:'__sink()',
##                8:'LRYSC',
##                9:'ICDc',
##               10:'ICDn',
##               11:'ICDnCSL',
##               12:'gHes1ICDCSL_off',
##               13:'CSLgHES1_off',
##               14:'gHes1ICDCSL_on',
##               15:'CSLgHES1_on',
##               16:'mHES1n_off',
##               17:'mHES1c_off',
##               18:'mHES1c_on',
##               19:'pHES1',
##               20:'gHES1pHES1_off',
##               21:'ICDnCSLgHES1pHES1_off',
##               22:'CSLgHES1pHES1_off'
## }
#

#print model.rules
#generate_equations(model, verbose=True)
#print
for sp in model.species:
    print sp
print 
#for rx in model.reactions:
#    print rx
#
#from pysb.generator.bng import BngGenerator
#print BngGenerator(model).get_content()
#
#
tspan = np.linspace(0, 720, 7201)
sim1 = SOS(model, tspan)
sim_result = sim1.run()

##
##for i,sp in enumerate(model.species):
##    print i,":",sp
#
## for  j,ode in enumerate(model.odes):
##    for i in range(len(model.species)):
##       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
##    print j,":",ode
#
print('Active CSL complex levels Cell 1')
print(sim_result.observables(['CSL_on_1']))
print('Active CSL complex levels Cell 2')
print(sim_result.observables(['CSL_on_2']))
print('Ligand Cell 1')
print(sim_result.observables['L_1'])
print('Ligand Cell 2')
print(sim_result.observables['L_2'])
print('Receptor Cell 1')
print(sim_result.observables['R_1'])
print('Receptor Cell 2')
print(sim_result.observables['R_2'])
print('mRNA_Hes1 Cell 1')
print(sim_result.observables['mRNA_Hes1_1'])
print('mRNA_Hes1 Cell 2')
print(sim_result.observables['mRNA_Hes1_2'])
print('protein_HES1, Cell 1')
print(sim_result.observables['prot_Hes1_1'])
print('protein_HES1, Cell 2')
print(sim_result.observables['prot_Hes1_2'])
print('Slug_levels, Cell 1')
print(sim_result.observables['Slug_levels_1'])
print('Slug_levels, Cell 2')
print(sim_result.observables['Slug_levels_2'])
print('Ascl1_levels, Cell 1')
print(sim_result.observables['Ascl1_levels_1'])
print('Ascl1_levels, Cell 2')
print(sim_result.observables['Ascl1_levels_2'])

plt.figure(figsize=(15,5))
plt.subplot(611)
plt.plot(tspan/60, sim_result.observables['L_1'], label='Ligand Cell 1')
plt.plot(tspan/60, sim_result.observables['L_2'], label='Ligand Cell 2')
plt.plot(tspan/60, sim_result.observables['R_1'], label='Receptor Cell 1')
plt.plot(tspan/60, sim_result.observables['R_2'], label='Receptor Cell 2')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(612)
plt.plot(tspan/60, sim_result.observables['R_ICD_1'], label='Nuclear ICD Cell 1')
plt.plot(tspan/60, sim_result.observables['R_ICD_2'], label='Nuclear ICD Cell 2')
plt.plot(tspan/60, sim_result.observables['CSL_on_1'], label='Activated CSL Cell 1')
plt.plot(tspan/60, sim_result.observables['CSL_on_2'], label='Activated CSL Cell 2')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(613)

plt.plot(tspan/60, sim_result.observables['prot_Hes1_1'], label='pHes1 Cell 1')
plt.plot(tspan/60, sim_result.observables['prot_Hes1_2'], label='pHes1 Cell 2')
plt.plot(tspan/60, sim_result.observables['mRNA_Hes1_1'], label='mHes1 Cell 1')
plt.plot(tspan/60, sim_result.observables['mRNA_Hes1_2'], label='mHes1 Cell 2')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(614)
plt.plot(tspan/60, sim_result.observables['Slug_levels_1'], label='SLUG Cell 1')
plt.plot(tspan/60, sim_result.observables['Slug_levels_2'], label='SLUG Cell 2')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(615)
plt.plot(tspan/60, sim_result.observables['Ascl1_levels_1'], label='Ascl1 Cell 1')
plt.plot(tspan/60, sim_result.observables['Ascl1_levels_2'], label='Ascl1 Cell 2')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)

plt.subplot(616)
plt.plot(tspan/60, sim_result.observables['NGN2_levels_1'], label='NGN2 Cell 1')
plt.plot(tspan/60, sim_result.observables['NGN2_levels_2'], label='NGN2 Cell 2')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
plt.legend(loc=0)


plt.show()
