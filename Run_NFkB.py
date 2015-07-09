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
from pysb.util import alias_model_components

import Model_NFkB


Model ()

Model_NFkB.declare_monomers()
Model_NFkB.declare_parameters()
alias_model_components()
Model_NFkB.declare_initial_conditions()
Model_NFkB.declare_observables()
Model_NFkB.declare_functions()
Model_NFkB.declare_rules()


#generate_equations(model, verbose=True)

# for i,sp in enumerate(model.species):
#     print i,":",sp
# print
# for i,rxn in enumerate(model.reactions):
#     print i,":",rxn
# print
# for i,ode in enumerate(model.odes):
#     print i,":",ode

# Simulate the model
time = np.linspace(0, 18000, 1801)

# ODE simulation
#@TODO fix array problem
# param_values = np.array([p.value for p in model.parameters])
# param_values[18]=0.1
plt.figure(1)
x = odesolve(model, time, verbose=True) #integrator='lsoda',
for obs in ["IKKKa_obs", "IKKii_obs", "IKKa_obs", "IKKi_obs", "IkBa_p_obs", "IkBan_NFkBn_obs"]:
    plt.plot(time, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    # re.match(r"(\w+)_obs", obs).group()
    #plt.plot(time/60, x[obs.name], label=obs.name)
    plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

# plt.figure()
# for i in range(len(model.species)):
#     plt.plot(time/60, x["__s%d" %i])
plt.figure(2)
x = odesolve(model, time, verbose=True)
for obs in ["TNFR1a_obs", "A20_on_obs"]:
    plt.plot(time, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.figure(3)
x = odesolve(model, time, verbose=True)
for obs in ["A20_obs", "A20_off_obs", "IkBa_obs", "IkBan_obs"]:
    plt.plot(time, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentrations", fontsize=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

x = odesolve(model, time, verbose=True)
plt.figure(4)
plt.plot(time/60, x["NFkBn_obs"], label=NFkBn_obs)
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
plt.legend()

plt.figure(5)
plt.plot(time/60, x["IkBa_off_obs"], label=IkBa_off_obs)
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
plt.legend()

plt.figure(6)
plt.plot(time/60, x["IkBa_on_obs"], label=IkBa_on_obs)
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
plt.legend()

plt.show()
# # SSA simulation
# x = run_ssa(model, time, verbose=True)
# plt.figure('SSA')
# plt.plot(time, x['B_obs'], label=B_obs, lw=2)
# plt.xlabel('Time (seconds)')
# plt.ylabel('Amount of B_obs')
# plt.legend(loc=0)

# plt.figure('ODE')
# plt.plot(time, x['IKKKa_obs'], label=IKKKa_obs, lw=2)
# plt.xlabel('Time (seconds)')
# plt.ylabel('Amount of IKKKa_obs')
# plt.legend(loc=0)
#plt.show()
