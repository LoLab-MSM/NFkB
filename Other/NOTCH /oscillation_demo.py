# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 18:17:49 2017

@author: Souhrid
"""

from pysb.core import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import ScipyOdeSimulator as SOS
from pysb.bng import generate_equations

Model()

Monomer('gene', ['p1', 'p2'])
Monomer('mrna',['e'])
Monomer('prot',['g','e'])
Monomer('enzp', ['p'])
Monomer('enzm', ['m'])


Parameter('enzp_0',10)
Parameter('enzm_0',10)
Parameter('gene_0', 1)
Initial(gene(p1=None,p2=None), gene_0)
Initial(enzp(p=None),enzp_0)
Initial(enzm(m=None),enzm_0)

Observable('mRNA_obs', mrna())
Observable('Protein_obs', prot())

Parameter('kgpf', 1)
Parameter('kgpr', 4) # 1 #8

Rule('gene_binds_protein_1', gene(p1=None) + prot(g=None,e=None) <>
     gene(p1=1) % prot(g=1,e=None), kgpf,kgpr)

# Rule('gene_binds_protein_2', gene(p2=None) + prot(g=None,e=None) <>
#      gene(p2=1) % prot(g=1,e=None), kgpf,kgpr)

Parameter('ktrcxn', 100)
Rule('mRNA_synth', gene(p1=None,p2=None) >> mrna(e=None) + gene(p1=None,p2=None), ktrcxn)

Parameter('ktrlxn', 0.03) # 1
Rule('Protein_synth', mrna(e=None) >> mrna(e=None) + prot(g=None,e=None), ktrlxn)

Parameter('kmenzf', 1)
Parameter('kmenzr', 1)
Parameter('kpenzf', 1)
Parameter('kpenzr', 1)
Parameter('kmdeg', 0.15) # 1
Parameter('kpdeg', 1)

Rule('m_deg_complex', mrna(e=None) + enzm(m=None) <> mrna(e=1) % enzm(m=1), kmenzf,kmenzr)
Rule('p_deg_complex', prot(g=None,e=None) + enzp(p=None) <> prot(g=None,e=1) % enzp(p=1), kpenzf,kpenzr)

Rule('mRNA_deg',mrna(e=1) % enzm(m=1) >> enzm(m=None), kmdeg)

Rule('protein_deg', prot(e=1) % enzp(p=1) >> enzp(p=None), kpdeg)

generate_equations(model, verbose=True)
print
for sp in model.species:
    print sp

tspan = np.linspace(0, 5000, 1001)
sim1 = SOS(model, tspan, verbose=True)
sim_result = sim1.run()

plt.plot(tspan, sim_result.observables['mRNA_obs'], lw=2, label='mRNA')
plt.plot(tspan, sim_result.observables['Protein_obs'], lw=2, label='Protein')
plt.xlabel("Time")
plt.ylabel("Molecules")
plt.legend(loc=0)
plt.show()