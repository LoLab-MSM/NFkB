from pysb import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from pysb.logging import setup_logger

setup_logger(level = 10)



Model()

ncells = 101
r = np.random.random(ncells)
rxnpat = []
fexps = []
obs = model.observables

Parameter('c1', 50)
Parameter('c2', 100)
for n in range(ncells):
    mon = Monomer('N%d' %n)
    par = Parameter('N%d_0' %n, r[n])
    monpat = MonomerPattern(mon, {}, None)
    cpxpat = ComplexPattern([monpat], None)
    Initial(cpxpat, par)
    rxnpat.append(ReactionPattern([cpxpat]))
    Observable('N%d_obs' %n, rxnpat[-1])
    gexp = Expression('g%dexp' %n, 1.0 / (1 + c2 * obs[n] ** 2))
    fexps.append(Expression('f%dexp' %n, gexp**2/(c1 + gexp**2)))

Parameter('kdeg', 1)
for n in range(ncells):
    # degradation rule
    rxpr = RuleExpression(rxnpat[n], ReactionPattern([]), False)
    Rule('N%d_deg' % n, rxpr, kdeg)
    #synthesis rules
    rxpr = RuleExpression(ReactionPattern([]), rxnpat[n], False)
    if n == 0:
        Rule('N%d_synth_right' % n, rxpr, fexps[n + 1])
    elif n == ncells - 1:
        Rule('N%d_synth_left' % n, rxpr, fexps[n - 1])
    else:
        Rule('N%d_synth_right' % n, rxpr, fexps[n + 1])
        Rule('N%d_synth_left' % n, rxpr, fexps[n - 1])

tspan = np.linspace(0, 30, 31)
# x = odesolve(model,tspan, verbose=False)
sim = ScipyOdeSimulator(model, tspan=tspan)
sim_result = sim.run()
df = sim_result.dataframe
# print(df['N1_obs'].iloc[-1])
# print(df['N2_obs'].iloc[-1])
# quit()
k = 0

for n in range(2, ncells):
    k = k + abs(2.0 * \
                (df['N%d_obs' % n].iloc[-1] - df['N%d_obs' % (n - 1)].iloc[-1]) / \
                (df['N%d_obs' % n].iloc[-1] + df['N%d_obs' % (n - 1)].iloc[-1]))
    # kbar[-1] = kbar[-1] + abs(k)
k = k / (ncells - 1)
print(k)
#
# l = k
# for n in range(0,ncells):
#     plt.plot(tspan, df['N%d_obs' % n].iloc[:], lw =2)
# plt.xlabel('Time', fontsize=16)
# plt.ylabel('Ni Cells', fontsize=16)
# plt.title('Kbar = %d' %l)
# # plt.legend(loc = 0)
# plt.show()

# c1pars = np.linspace(0.01, 100, 5)
# c2pars = np.linspace(0.01, 100, 5)
#
# c1v = []
# c2v = []
# for i in c1pars:
#     for j in c2pars:
#         c1v.append(i)
#         c2v.append(j)

#
# sim = ScipyOdeSimulator(model, tspan = tspan)
# # sim_result = sim.run(param_values= {'c1': c1v, 'c2': c2v})
# sim_result = sim.run()
# print('it finished')
# df = sim_result.dataframe
# # w = odesolve(model,tspan, param_values= {'c1': c1v, 'c2': c2v},verbose=False)
# # print 'hello'
# # print [w[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]]
# # quit()
# kbar = []
# k = []
# #
# # for i in c1pars:
# #     for j in c2pars:
#         # print (i,j)
# kbar.append(0.0)
# for n in range(2,ncells):
#     print n
#     # x = odesolve(model,tspan, param_values= {'c1': c1v, 'c2': c2v},verbose=False)
#     k = 2.0 * \
#         (df['N%d_obs' %n].iloc[-1] - df['N%d_obs' %(n-1)].iloc[-1]) / \
#         (df['N%d_obs' %n].iloc[-1] + df['N%d_obs' %(n-1)].iloc[-1])
#     # print(df['N%d_obs' %n].iloc[-1])
#     # print()
#     # print(df['N%d_obs' %(n-1)].iloc[-1])
#     # print(k)
#     # k = 2.0 * \
#     #     ([w[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]] - sim_result.observables[n-1][-1]) / \
#     #     ([w[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]] + sim_result.observables[n-1][-1])
#     kbar[-1] = kbar[-1] + abs(k)
#     # print kbar[-1]
#     # print()
# kbar[-1] = (1.0/(ncells -1)) * kbar[-1]
# print('final kbar')
# print(kbar[-1])

# for obs in model.observables:
#     plt.plot(tspan, x[obs.name], lw =2, label = obs.name)
# # plt.legend(loc = 0)
# plt.show()