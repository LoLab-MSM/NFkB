__author__ = 'geena'
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
from Model_NFkB import model as m1
# from NFkB_Unstruct import model as m2

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
from numpy import array

Model ()

#Running simulation
time = np.linspace(1, 18000, 1801)
x = odesolve(m1, time, verbose = True)
# y = run_ssa(m1, time, verbose = True)
# y = odesolve(m2, time, verbose = True)
# colors = np.linspace(0,0.9,8)
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '0.5']


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
plt.show()


