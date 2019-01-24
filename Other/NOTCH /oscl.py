from pysb import *
from pylab import *
from pysb import bng
from pysb.integrate import odesolve
from pysb.bng import generate_network
from pysb.bng import generate_equations
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import random
# import numpy as np
# from pysb.kappa import contact_map, set_kappa_path, influence_map
# # from pysb.tools.render_reactions import run
# import pygraphviz as pyg
import model2 as m2
import pandas as pd
from pysb.core import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import ScipyOdeSimulator as SOS

Model()

#Oct4 and Cdx2 monomers for promoter, and gene
#
# Monomer('oct4p', ['oct4gene'])
# Monomer('cdx2p', ['cdx2gene'])
# Monomer('cdx2_gene', ['p', 'cdx2', 'oct4', 'switch'], {'switch': ('on', 'off')})
# Monomer('cdx2_gene', ['p', 'cdx2', 'oct4', 'switch'], {'switch': ('on', 'off')})

# Monomer('x_gene', ['p','switch'], {'switch': ['on', 'off']})
# Monomer('x_promoter', ['gene'])
# Monomer('x_mrna')
# Monomer('x')

# no_promoter = False
# if no_promoter == True:
#     Monomer('x_gene', ['p', 'switch'], {'switch': ['on', 'off']})
#     # Monomer('x_promoter', ['gene'])
#     Monomer('x_mrna')
#     Monomer('x')
#     #initial conditions
#     Parameter('dna_off_switch_x', 0)
#     Initial(x_gene(p = None, switch = 'off'), dna_off_switch_x)
#     Parameter('dna_on_switch_x', 1)
#     Initial(x_gene(p = None, switch = 'on'), dna_on_switch_x)
#     Parameter('mrna_initial', 20)
#     Initial(x_mrna(), mrna_initial)
#     Parameter('protein_initial', 200)
#     Initial(x(), protein_initial)
#
#     #Observables for cdx2 and oct4
#     Observable('x_obs', x())
#     Observable('x_off_obs', x_gene(p = None, switch = 'off'))
#     Observable('x_on_obs', x_gene(p = None, switch = 'on'))
#     Observable('x_mrna_obs', x_mrna())
#     # Observable('oct4_mrna_obs', oct4_mrna())
#     # Observable('cdx2_mrna_obs', cdx2_mrna())
#
#
#     #reversible activation and deactivation of dna for oct4 and cdx2
#     # dna_oct4(off) + dna_oct4(on) <> dna_oct4(on)
#     # dna_cdx2(off) + dna_cdx2(on) <> dna_cdx2(on)
#     e = random.randint(0,10)
#     f = random.randint(0,10)
#     # Parameter('dna_switch_on', e) # .30)
#     # Parameter('dna_switch_off', f) # .10)
#     Parameter('dna_switch_on', .1)
#     Expression('prop_switch', dna_switch_on) #*x_off_obs)
#     Parameter('dna_switch_off', .5)
#     Rule('dna_switch_x', x_gene(p = None, switch = 'off') <> x_gene(p = None, switch = 'on'), prop_switch, dna_switch_off)
#
#     #dna activation of mrna for oct4 and cdx2
#     # dna(on) + oct4_mrna >> oct4_mrna
#     # dna(on) + cdx2_mrna >> cdx2_mrna
#     # oct4_mrna >> None
#     # cdx2_mrna >> None
#
#     a = random.randint(0,10)
#     b = random.randint(0,1)
#     # print(a)
#     # print(b)
#     Parameter('mrna_x', 1)
#     # print(mrna_x)
#     Parameter('mrna_x_deg', .001)
#     # Rule('synthesize_oct4_mrna', None >> oct4_mrna(), mrna_oct4)
#     # Rule('synthesize_cdx2_mrna', None >> cdx2_mrna(), mrna_cdx2)
#     Rule('dna_to_mrna_x', x_gene(p = None, switch = 'on') >> x_mrna() + x_gene(p = None, switch = 'on') , mrna_x)
#     Rule('deg_x_mrna', x_mrna() >> None, mrna_x_deg)
#
#     # mrna to protein reactions and degredation
#     # oct4_mrna + oct4 >> oct4
#     # cdx2_mrna + cdx2 >> cdx2
#     # oct4 >> None
#     # cdx2 >> None
#     # c = random.randint(0,10)
#     # d = random.randint(0,1)
#     # print(c)
#     # print(d)
#     Parameter('mrna_x_protein', 1)
#     Parameter('x_deg', .001)
#     Rule('x_mrna_to_protein', x_mrna() >> x() + x_mrna(), mrna_x_protein)
#     Rule('deg_x', x() >> None, x_deg)

dnabind = True
if dnabind == True:
    Monomer('x_gene', ['p', 'switch'], {'switch': ['on', 'off']})
    Monomer('x_promoter', ['gene'])
    Monomer('x_mrna')
    Monomer('x')
    Parameter('dna_off_switch_x', 0)
    Initial(x_gene(p = None, switch = 'off'), dna_off_switch_x)
    Parameter('dna_on_switch_x', 1)
    Initial(x_gene(p = None, switch = 'on'), dna_on_switch_x)
    Parameter('mrna_initial', 100)
    Initial(x_mrna(), mrna_initial)
    Parameter('protein_initial', 250)
    Initial(x(), protein_initial)
    Parameter('promoter_initial', .4)
    Initial( x_promoter(gene= 1) % x_gene(p=1, switch = 'on'), promoter_initial)

    #Observables for cdx2 and oct4
    Observable('x_obs', x())
    Observable('x_off_obs', x_gene(p = None, switch = 'off'))
    Observable('x_on_obs', x_gene(p = None, switch = 'on'))
    Observable('x_mrna_obs', x_mrna())
    # Observable('x_gene_obs', x_gene(p = None, switch = 'on'))
    Observable('x_promoter_bind_obs', x_promoter(gene= 1) % x_gene(p=1, switch = 'on'))

    Parameter('dna_switch_on', 80)
    Parameter('dna_switch_off', 10)
    Rule('dna_switch_x', x_gene(p = None, switch = 'off') <> x_gene(p = None, switch = 'on'), dna_switch_on, dna_switch_off)

    Parameter('x_promoter_bind_dna', .01)
    Parameter('x_promoter_unbind_dna', .5)
    Parameter('x_promoter_dna_mrna', .002)
    Expression('protein_induced', x_promoter_bind_dna*x_obs)
    Expression('protein_deduced', x_promoter_unbind_dna)
    Rule('x_dna_switch_protein', x_gene(p = None, switch = 'off')  <> x_gene(p = 1, switch = 'on') % x_promoter(gene = 1), protein_induced, protein_deduced)
    Rule('xdna_to_mrna', x_gene(p = 1, switch = 'on') % x_promoter(gene = 1) >> x_mrna() + x_gene(p = 1, switch = 'on') % x_promoter(gene = 1), x_promoter_dna_mrna)


    Parameter('mrna_x', .002)
    Rule('dna_to_mrna_x', x_gene(p = None, switch = 'on') >> x_mrna() + x_gene(p = None, switch = 'on') , mrna_x)
    Parameter('mrna_x_deg', .0001)
    Rule('deg_x_mrna', x_mrna() >> None, mrna_x_deg)
    Parameter('mrna_x_protein', .002)
    Parameter('x_deg', .0001)
    Rule('x_mrna_to_protein', x_mrna() >> x() + x_mrna(), mrna_x_protein)
    Rule('deg_x', x() >> None, x_deg)


# generate_network(model)
# generate_equations(model)

tspan = np.linspace(0, 240, 2401)
# sim1 = SOS(model, tspan)
# sim_result = sim1.run()


# ssa_sim = bng.run_ssa(model, t_end=36000, step=36001)
# ssa_sim = bng.run_ssa(model, t_end=72000, step=72001)
ssa_sim = bng.run_ssa(model, t_end=180000, step=180001)
# tspan = np.linspace(0, 720, 721)
# phase = odesolve(model,tspan,verbose=True)


last_conc = [ssa_sim[x][-1] for x in ['__s%d' % i for i in np.arange(len(model.species))]]
concs = [ssa_sim[x][:] for x in ['__s%d' % i for i in np.arange(len(model.species))]]
df = pd.DataFrame(concs)
df['mean'] = df.mean(axis=1)
print(df['mean'])
print("mean 1")
# print(df['mean'][3])
# df.to_csv('grn_data.csv')
#
# print(last_conc)
#
# for i,sp in enumerate(model.species):
#     print i,":",sp

# plt.figure()
# # plt.figure(figsize=(15, 10))
# # plt.subplot(121)
# plt.plot(tspan/60, sim_result.observables['x_mrna_obs'], label='x_mrna_obs')
# plt.plot(tspan/60, sim_result.observables['x_obs'], label='x_obs')
# # plt.plot(tspan/60, sim_result.observables['prot_gene'], label='bound gene')
# plt.xlabel("Time (in hr)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)
# plt.show()
hr = True
if hr == True:
    figure()
    plt.plot((ssa_sim['time']/60)/60, ssa_sim['x_obs'], color='g', label='x_ssa')
    plt.plot(df['mean'][3], color = 'r', linestyle = '-')
    # plt.plot(ssa_sim['time'], ssa_sim['x_off_obs'], color='r', label = 'x_off_ssa')
    # plt.plot(t, ode_sim['x_obs'], color='g', label = 'x_ode')
    # plt.xlabel("Time (in min)", fontsize=16)
    plt.xlabel("Time (in hr)", fontsize=16)
    plt.ylabel("Molecules per Cell", fontsize=16)
    # plt.ylim(ymin = -.5, ymax = 250)
    plt.legend(loc=0)
    # plt.yscale('log', nonposy='clip')

    figure()
    plt.plot((ssa_sim['time']/60)/60, ssa_sim['x_off_obs'], color='r', label='x_off_ssa')
    plt.plot((ssa_sim['time']/60)/60, ssa_sim['x_on_obs'], color='b', label='x_on_ssa')
    plt.plot((ssa_sim['time']/60)/60, ssa_sim['x_promoter_bind_obs'], color='c', label='x_promoter_bind_ssa')
    # plt.plot(t, ode_sim['x_off_obs'], color='g', label = 'x_off_ode')
    # plt.xlabel("Time (in min)", fontsize=16)
    plt.xlabel("Time (in hr)", fontsize=16)
    plt.ylabel("Molecules per Cell", fontsize=16)
    plt.ylim(ymin = -.5, ymax = 2.1)
    # plt.xlim()
    plt.legend(loc=0)
    # plt.yscale('log', nonposy='clip')

    figure()
    plt.plot((ssa_sim['time']/60)/60, ssa_sim['x_mrna_obs'], color='y', label='x_mrna_ssa')
    # plt.plot(t, ode_sim['x_on_obs'], color='g', label = 'x_on_ode')
    # plt.xlabel("Time (in min)", fontsize=16)
    plt.xlabel("Time (in hr)", fontsize=16)
    plt.ylabel("Molecules per Cell", fontsize=16)
    # plt.ylim(ymin = 15, ymax = 22)
    plt.legend(loc=0)
    # plt.yscale('log', nonposy='clip')

    # figure()
    # # plt.plot(ssa_sim['time'], ssa_sim['x_gene_obs'], color='p', label='x_gene_ssa')
    # plt.plot(ssa_sim['time']/60, ssa_sim['x_promoter_bind_obs'], color='c', label='x_promoter_bind_ssa')
    # # plt.plot(t, ode_sim['x_on_obs'], color='g', label = 'x_on_ode')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = -5, ymax = 5)
    # plt.legend(loc=0)
    # plt.yscale('log', nonposy='clip')

    plt.show()
#
# min = False
# if min == True:
#     figure()
#     plt.plot(ssa_sim['time']/60, ssa_sim['x_obs'], color='g', label='x_ssa')
#     # plt.plot(ssa_sim['time'], ssa_sim['x_off_obs'], color='r', label = 'x_off_ssa')
#     # plt.plot(t, ode_sim['x_obs'], color='g', label = 'x_ode')
#     # plt.xlabel("Time (in min)", fontsize=16)
#     plt.xlabel("Time (in min)", fontsize=16)
#     plt.ylabel("Molecules per Cell", fontsize=16)
#     # plt.ylim(ymin = -.5, ymax = 250)
#     plt.legend(loc=0)
#     # plt.yscale('log', nonposy='clip')
#
#     figure()
#     plt.plot(ssa_sim['time']/60, ssa_sim['x_off_obs'], color='r', label='x_off_ssa')
#     plt.plot(ssa_sim['time']/60, ssa_sim['x_on_obs'], color='b', label='x_on_ssa')
#     plt.plot(ssa_sim['time']/60, ssa_sim['x_promoter_bind_obs'], color='c', label='x_promoter_bind_ssa')
#     # plt.plot(t, ode_sim['x_off_obs'], color='g', label = 'x_off_ode')
#     # plt.xlabel("Time (in min)", fontsize=16)
#     plt.xlabel("Time (in min)", fontsize=16)
#     plt.ylabel("Molecules per Cell", fontsize=16)
#     plt.ylim(ymin = -.5, ymax = 1.5)
#     # plt.xlim()
#     plt.legend(loc=0)
#     # plt.yscale('log', nonposy='clip')
#
#     figure()
#     plt.plot(ssa_sim['time']/60, ssa_sim['x_mrna_obs'], color='y', label='x_mrna_ssa')
#     # plt.plot(t, ode_sim['x_on_obs'], color='g', label = 'x_on_ode')
#     # plt.xlabel("Time (in min)", fontsize=16)
#     plt.xlabel("Time (in min)", fontsize=16)
#     plt.ylabel("Molecules per Cell", fontsize=16)
#     # plt.ylim(ymin = 15, ymax = 22)
#     plt.legend(loc=0)
#     # plt.yscale('log', nonposy='clip')
#
#     # figure()
#     # # plt.plot(ssa_sim['time'], ssa_sim['x_gene_obs'], color='p', label='x_gene_ssa')
#     # plt.plot(ssa_sim['time']/60, ssa_sim['x_promoter_bind_obs'], color='c', label='x_promoter_bind_ssa')
#     # # plt.plot(t, ode_sim['x_on_obs'], color='g', label = 'x_on_ode')
#     # plt.xlabel("Time (in min)", fontsize=16)
#     # plt.ylabel("Concentration", fontsize=16)
#     # plt.ylim(ymin = -5, ymax = 5)
#     # plt.legend(loc=0)
#     # plt.yscale('log', nonposy='clip')
#
#     plt.show()
