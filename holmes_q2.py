from pysb import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from pysb.logging import setup_logger
import pickle as p
import pandas as pd
import seaborn as sns

setup_logger(level = 10)

Model()

ncells = 101
r = np.random.random(ncells)
rxnpat = []
fexps = []
obs = model.observables

Parameter('c1', 0.1)
Parameter('c2', 10)
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
c1pars = np.linspace(0.02, 99, 5)
c2pars = np.linspace(0.02, 99, 5)

c1v = []
c2v = []
for i in c1pars:
    for j in c2pars:
        c1v.append(i)
        c2v.append(j)


sim = ScipyOdeSimulator(model, tspan = tspan)
sim_result = sim.run(param_values= {'c1': c1v, 'c2': c2v})
df = sim_result.dataframe
print(df.loc[0]['N1_obs'].iloc[:])
# print(df.loc[0]['N2_obs'].iloc[-1])
# print(df.loc[0]['N3_obs'].iloc[-1])
print(df.loc[1]['N1_obs'].iloc[:])
print(df.loc[2]['N1_obs'].iloc[:])
quit()

kbar = []
# for i in c1pars:
#     for j in c2pars:
#         print('ij')
#         print(i,j)

for g in range(0,2499):
    # print(g)
    kbar.append(0.0)
    for n in range(2, ncells):
        # print(n)
        kbar[-1] = kbar[-1] + abs(2.0 * \
            (df.loc[g]['N%d_obs' %n].iloc[-1] - df.loc[g]['N%d_obs' %(n-1)].iloc[-1]) / \
            (df.loc[g]['N%d_obs' %n].iloc[-1] + df.loc[g]['N%d_obs' %(n-1)].iloc[-1]))
    kbar[-1] = kbar[-1]/ (ncells - 1)
    # print('k')
    # print(kbar[-1])
print(len(kbar))
print(kbar)

X = c1v[:-1]
Y = c2v[:-1]
Z = kbar
data = pd.DataFrame({'c1': X, 'c2': Y, 'kbar': Z})
data_pivoted = data.pivot("c1", "c2", "kbar")
ax = sns.heatmap(data_pivoted)
plt.show()


















# # return k/ (ncells - 1)
# #         print('kbarafkn')
# #         print(kbar[-1])
# #         print(kbar)
# # print('final kbar')
# # print(len(kbar))
# # print(kbar)
# print('DONE')




# for i in c1pars:
#     for j in c2pars:
#         kbar.append(0.0)
#         for n in range(2, ncells):
#             k = 2.0 * \
#                 (df['N%d_obs' %n].iloc[-1] - df['N%d_obs' %(n-1)].iloc[-1]) / \
#                 (df['N%d_obs' %n].iloc[-1] + df['N%d_obs' %(n-1)].iloc[-1])
#             kbar[-1] = kbar[-1] + abs(k)
#         kbar[-1] = kbar[-1]/ (ncells - 1)
#         print('kbarafkn')
#         print(kbar)
# #             # l.append(kbar[-1])
# # print('kbar2')
# # print(kbar[-1])
# # l.append(kbar)
# print('final kbar')
# print(len(kbar))
# print(kbar)
# # print('all kbars')
# # print(final)
# #
# print('DONE')
f = open("c1", "w")
f.write("\n".join(map(lambda x: str(x), c1v[:-1])))
f.close()

f = open("c2", "w")
f.write("\n".join(map(lambda x: str(x), c2v[:-1])))
f.close()

f = open("kbar", "w")
f.write("\n".join(map(lambda x: str(x), kbar)))
f.close()



# with open('kbars.txt', 'wb') as fp:
#     p.dump((final, fp))



# for item in final:
#   thefile.write("%s\n" % item)

# for obs in model.observables:
#     plt.plot(tspan, x[obs.name], lw =2, label = obs.name)
# # plt.legend(loc = 0)
# plt.show()