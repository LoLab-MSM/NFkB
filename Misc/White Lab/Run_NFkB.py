__author__ = 'geena'
from Model_NFkB import model as m1
# from NFkB_Unstruct import model as m2

from pysb.integrate import odesolve
import matplotlib.pyplot as plt
import numpy as np
from pysb.core import *
from pysb.bng import generate_network
from pysb.bng import generate_equations

from pysb.kappa import contact_map, set_kappa_path, influence_map
from pysb.tools.render_reactions import run
import pygraphviz as pyg


Model ()

#Running simulation
generate_network(model)
generate_equations(model)
time = np.linspace(1, 18000, 1801)
x = odesolve(m1, time, verbose = True)
# y = run_ssa(m1, time, verbose = True)
# y = odesolve(m2, time, verbose = True)
# colors = np.linspace(0,0.9,8)
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '0.5']

set_kappa_path('/Users/geenaildefonso/Projects/KaSim')
x = contact_map(model)
x.draw('contact_map.pdf', format='pdf', prog='dot')

x = run(model)
g = pyg.AGraph(x)
g.draw('render_reactions.pdf', format='pdf', prog='dot')


#plotting figures of TNF, TNFR1, IkBa (in N and C), and IKK (active)
plt.figure()
# for i,sp in enumerate([2, 12, 4, 5]):
for i,sp in enumerate([2]):
        plt.plot(time/60, x["__s%d"%sp], color = str(colors[i]), lw = 3, label=str(m1.species[sp])) #@TODO not correct?
        # plt.plot(time/60, y["__s%d"%sp], "--", color = str(colors[i]), lw = 3, label=str(m2.species[sp]))
        plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig("Figure 2", format= "png")
# plt.yscale('log')
# plt.show()

plt.figure()
# for i,sp in enumerate([15, 19]):
for i,sp in enumerate([10,11]):
        plt.plot(time/60, x["__s%d"%sp], color = str(colors[i]), lw = 3, label=str(m1.species[sp])) #@TODO not correct?
        # plt.plot(time/60, y["__s%d"%sp], "--", color = str(colors[i]), lw = 3, label=str(m2.species[sp]))
        plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig("Figure 3", format= "png")
# plt.yscale('log')
# plt.show()

plt.figure()
for i,sp in enumerate([18]):
        plt.plot(time/60, x["__s%d"%sp], color = str(colors[i]), lw = 3, label=str(m1.species[sp])) #@TODO not correct?
        # plt.plot(time/60, y["__s%d"%sp], "--", color = str(colors[i]), lw = 3, label=str(m2.species[sp]))
        plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig("Figure 4", format= "png")
# plt.yscale('log')
# plt.show()


plt.figure()
for i,sp in enumerate([15]):
        plt.plot(time/60, x["__s%d"%sp], color = str(colors[i]), lw = 3, label=str(m1.species[sp])) #@TODO not correct?
        # plt.plot(time/60, y["__s%d"%sp], "--", color = str(colors[i]), lw = 3, label=str(m2.species[sp]))
        plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig("Figure 5", format= "png")
# plt.yscale('log')

plt.figure()
for i,sp in enumerate([6,7]):
        plt.plot(time/60, x["__s%d"%sp], color = str(colors[i]), lw = 3, label=str(m1.species[sp])) #@TODO not correct?
        # plt.plot(time/60, y["__s%d"%sp], "--", color = str(colors[i]), lw = 3, label=str(m2.species[sp]))
        plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig("Figure 5", format= "png")
# plt.yscale('log')
# plt.show()


