from pysb import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator

Model ()

Monomer('V')
Monomer('W')


Parameter('V_0', 5)
Parameter('W_0', 3)

Initial(V(), V_0)
Initial(W(), W_0)


Observable('V_obs', V())
Observable('W_obs', W())

Parameter('deg_w', 5e-7)
Parameter('v_w', 1)
Parameter('v_w2', .005)
Expression('v_exp', V_obs*(V_obs - 0.1)*(1.0-V_obs))
Expression('v_exp_deg', v_w*W_obs)
Expression('w_exp', v_w2*V_obs)

Rule('v_synth', None >> V(), v_exp)
Rule('v_deg', V() >> None, v_w)
# Rule('w_synth', None >> W(), w_exp)
# Rule('v_deg', V() >> None, v_exp_deg)
Rule('w_synth', None >> W(), w_exp)

Rule('w_deg', W() >> None, deg_w)


# tspan = np.linspace(0,10,101)
#
# V_init =  np.linspace(0.01,10,10)
# W_init =  np.linspace(0.01,10,10)
#
# x = list(V_init)
# y = list(W_init)
#
#

# import sympy as sy
# L = sy.symbols('L')
# M = sy.Matrix([[-1 - L, .439, -.421, 0],
#                [1, -.439 - L, .421, 0],
#                [-.421, 0, -1- L, .439],
#                [.421, 0, 1, -.439 - L]])
# print 'Characteristic Polynomial: P(L) = ', M.det()
# print 'Eigenvalues: ', sy.solve(M.det())



#
# c1v = []
# c2v = []
# for i in x:
#     for j in y:
#         c1v.append(i)
#         c2v.append(j)
# print(c1v)
# print()
# print(c2v)
# quit()
#
# sim = ScipyOdeSimulator(model, tspan = tspan)
# sim_result = sim.run(param_values= {'V_0': x, 'W_0': y})

tspan = np.linspace(0, 10, 101)
sim1 = ScipyOdeSimulator(model, tspan)
sim2 =ScipyOdeSimulator(model, tspan)
sim3 = ScipyOdeSimulator(model, tspan)
sim4 =ScipyOdeSimulator(model, tspan)

L1 = sim1.run(param_values={'V_0':1, 'W_0':4})
L2 = sim2.run(param_values={'V_0':3, 'W_0':2.5})
L3 = sim3.run(param_values={'V_0':5, 'W_0':3})
L4 = sim4.run(param_values={'V_0':4, 'W_0':1})


print(model.species)
print(model.parameters)
print(model.expressions)

for i,ode in enumerate(model.odes):
    print i,":",ode

quit()

print(L1.observables['V_obs'][-1], L1.observables['W_obs'][-1])
print(L2.observables['V_obs'][-1])
print(L3.observables['V_obs'][-1])
print(L4.observables['V_obs'][-1])
print()
print(L1.observables['W_obs'][-1])
print(L2.observables['W_obs'][-1])
print(L3.observables['W_obs'][-1])
print(L4.observables['W_obs'][-1])
# plt.figure(figsize = (15,10))
# plt.subplot(231)
plt.figure()
plt.plot(L1.observables['V_obs'], L1.observables['W_obs'], label = 'x1y4', )
plt.plot(L2.observables['V_obs'], L2.observables['W_obs'], label = 'x3y2.5')
plt.plot(L3.observables['V_obs'], L3.observables['W_obs'], label = 'x5y3')
plt.plot(L4.observables['V_obs'], L4.observables['W_obs'], label = 'x4y1')
# plt.plot(L1.observables['x_obs'][0], L1.observables['y_obs'][0])
plt.xlabel("X", fontsize=10)
plt.ylabel("Y", fontsize=10)
plt.ylim(ymin = -3, ymax = 7)
plt.xlim(xmin = -3, xmax = 7)
plt.legend(loc=0)
plt.show()



df = sim_result.dataframe
print(df.shape)

# print('V')
# print(sim_result.observables['V_obs'])
# print('W')
# print(sim_result.observables['W_obs'])
#
plt.figure()
plt.plot(df['V_obs'].iloc[:], df['W_obs'].iloc[:])
# # plt.plot(tspan, sim_result.observables['Ri_obs'], label = 'Ri_obs')
# # plt.plot(tspan, sim_result.observables['P_obs'], label = 'P_obs')
# # plt.plot(tspan, sim_result.observables['Pi_obs'], label = 'Pi_obs')
plt.xlabel("V", fontsize=16)
plt.ylabel("W", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
plt.show()
# plt.plot(tspan, sim_result.observables['W_obs'], label = 'W')
# # plt.plot(tspan, sim_result.observables['Ri_obs'], label = 'Ri_obs')
# # plt.plot(tspan, sim_result.observables['P_obs'], label = 'P_obs')
# # plt.plot(tspan, sim_result.observables['Pi_obs'], label = 'Pi_obs')
# plt.xlabel("W", fontsize=16)
# plt.ylabel("time", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.plot(tspan, sim_result.observables['V_obs'], label = 'V')
# # plt.plot(tspan, sim_result.observables['Ri_obs'], label = 'Ri_obs')
# # plt.plot(tspan, sim_result.observables['P_obs'], label = 'P_obs')
# # plt.plot(tspan, sim_result.observables['Pi_obs'], label = 'Pi_obs')
# plt.xlabel("V", fontsize=16)
# plt.ylabel("time", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.show()