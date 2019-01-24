from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import numpy as np
import matplotlib.pyplot as plt

Model ()

Monomer('R', ['state'], {'state': ['off', 'on']})
Monomer('P', ['state'], {'state': ['off', 'on']})


Parameter('Ri_0', 1)
Parameter('R_0', 1)
Parameter('Pi_0', 1)
Parameter('P_0', 1)

Initial(R(state = 'off'), Ri_0)
Initial(R(state = 'on'), R_0)
Initial(P(state = 'off'), Pi_0)
Initial(P(state = 'on'), Pi_0)

Observable('R_obs', R(state = 'on'))
Observable('Ri_obs', R(state = 'off'))
Observable('P_obs', P(state = 'on'))
Observable('Pi_obs', P(state = 'off'))


Parameter('degr', 1)

# Expression('R_exp', Ri_obs/(2.0 + 2.0*(P_obs)**4))
# Expression('P_exp', Pi_obs/(2.0 + 2.0*(R_obs)**4))


Expression('R_exp', 1.0/(2.0 + 2.0*(P_obs)**4))
Expression('P_exp', 1.0/(2.0 + 2.0*(R_obs)**4))

Rule('Ri_R', R(state = 'off') >> R(state ='on'), R_exp)
Rule('Pi_P', P(state = 'off') >> P(state ='on'), P_exp)

Rule('R_Ri', R(state = 'on') >> R(state ='off'), degr)
Rule('P_Pi', P(state = 'on') >> P(state ='off'), degr)


tspan = np.linspace(0,10,101)
sim = ScipyOdeSimulator(model, tspan = tspan)
sim_result = sim.run()

print(model.species)
species_dict = {
    0: 'Ri',
    1: 'R',
    2: 'Pi',
    3: 'P'
}

for  j,ode in enumerate(model.odes):
    for i in range(len(model.species)):
       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
    print j,":",ode

print(model.expressions)
print(model.parameters)
quit()

# print('R_obs')
# print(sim_result.observables['R_obs'])
# print('Ri_obs')
# print(sim_result.observables['Ri_obs'])
# print('P_obs')
# print(sim_result.observables['P_obs'])
# print('Pi_obs')
# print(sim_result.observables['Pi_obs'])

plt.figure(figsize = (10,5))
plt.subplot(121)
plt.plot(sim_result.observables['R_obs'], sim_result.observables['Ri_obs'])
# plt.plot(tspan, sim_result.observables['Ri_obs'], label = 'Ri_obs')
# plt.plot(tspan, sim_result.observables['P_obs'], label = 'P_obs')
# plt.plot(tspan, sim_result.observables['Pi_obs'], label = 'Pi_obs')
plt.xlabel("R_active", fontsize=16)
plt.ylabel("R_inactive", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)


plt.subplot(122)
plt.plot(sim_result.observables['P_obs'], sim_result.observables['Pi_obs'])
# plt.plot(tspan, sim_result.observables['Ri_obs'], label = 'Ri_obs')
# plt.plot(tspan, sim_result.observables['P_obs'], label = 'P_obs')
# plt.plot(tspan, sim_result.observables['Pi_obs'], label = 'Pi_obs')
plt.xlabel("P_active", fontsize=16)
plt.ylabel("P_inactive", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan, sim_result.observables['R_obs'], label = 'R_obs')
# plt.plot(tspan, sim_result.observables['Ri_obs'], label = 'Ri_obs')
# plt.plot(tspan, sim_result.observables['P_obs'], label = 'P_obs')
# plt.plot(tspan, sim_result.observables['Pi_obs'], label = 'Pi_obs')
# plt.xlabel("Time (in seconds)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0, ymax =1)
# plt.legend(loc=0)
plt.show()

# plt.figure(figsize = (15,10))
# # plt.figure()
# plt.subplot(231)
# plt.plot(tspan/60, simu_result.observables['TNF_obs'], marker = '*',label = 'TNF')
# plt.plot(tspan/60, sim_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Concentration", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(232)
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'],marker = '*',label = 'TNFR')
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Concentration", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.tight_layout()
# plt.show()

