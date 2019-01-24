from pysb import *
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import scipy.interpolate
#from pysb.simulator import ScipyOdeSimulator
import numpy as np
# import holmes_model_two as m2
# import holmes_three as m3
# import holmes_four as m4
from pysb.util import alias_model_components
from matplotlib import *

Model()



#Making X and Y monomers

Monomer('X')
Monomer('Y')


Parameter('x_0', 1)
Parameter('y_0', 1)

Initial(X(), x_0)
Initial(Y(), y_0)

Observable('x_obs', X())
Observable('y_obs', Y())


Parameter('b', 0)
Parameter('a', 2)
Parameter('I', 0.5)
Parameter('n', 5)
Parameter('d', 1)

Expression('x_n', a*x_obs**n/(1+x_obs**n))
Expression('y_i', I/(1+y_obs**n))

Expression('y_n', a*y_obs**n/(1+y_obs**n))
Expression('x_i', I/(1+x_obs**n))

Rule('synth_x_b', None >> X(), b)
Rule('synth_x_obs', None >> X(), x_n)
Rule('synth_x_I', None >> X(), y_i)
Rule('x_deg', X() >> None, d)


Rule('synth_y_b', None >> Y(), b)
Rule('synth_y_obs', None >> Y(), y_n)
Rule('synth_y_I', None >> Y(), x_i)
Rule('y_deg', Y() >> None, d)


tspan = np.linspace(0, 50, 60)
sim1 = ScipyOdeSimulator(model, tspan)
sim2 =ScipyOdeSimulator(model, tspan)
sim3 = ScipyOdeSimulator(model, tspan)
sim4 =ScipyOdeSimulator(model, tspan)

L1 = sim1.run(param_values={'x_0':1, 'y_0':20})
L2 = sim2.run(param_values={'x_0':3, 'y_0':2.5})
L3 = sim3.run(param_values={'x_0':25, 'y_0':3})
L4 = sim4.run(param_values={'x_0':1.7, 'y_0':1.7})


# plt.figure(figsize = (15,10))
# plt.subplot(231)
plt.figure()
plt.plot(L1.observables['x_obs'], L1.observables['y_obs'], label = 'x1y20', )
plt.plot(L2.observables['x_obs'], L2.observables['y_obs'], label = 'x3y2.5')
plt.plot(L3.observables['x_obs'], L3.observables['y_obs'], label = 'x25y3')
plt.plot(L4.observables['x_obs'], L4.observables['y_obs'], label = 'x1.7y1.7')
# plt.plot(L1.observables['x_obs'][0], L1.observables['y_obs'][0])
plt.xlabel("X", fontsize=10)
plt.ylabel("Y", fontsize=10)
plt.ylim(ymin = -3, ymax = 5)
plt.xlim(xmin = -3, xmax = 5)
plt.legend(loc=0)
plt.show()


# # plt.figure()
# plt.subplot(232)
# plt.plot(tspan/60, simulation_result.observables['TNFRM_obs'], label = 'TNFRM_obs')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Concentration", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(233)
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], label = 'TNFR_obs')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Concentration", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(234)
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], label = 'IKKa_obs')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Concentration", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(235)
# plt.plot(tspan/60, simulation_result.observables['obs_A20t'], label = 'A20_obs')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Concentration", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(236)
# plt.plot(tspan/60, simulation_result.observables['IkBe_obs'], label = 'IkBe_obs')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Concentration", fontsize=10)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.tight_layout()
# plt.show()



# print(np.transpose(simulation_result.observables['x_obs']))
# print(L['x_obs'])
# print()
# print(L['y_obs'])

# a = [L.observables['x_obs'], L.observables['y_obs']]
# numpy.savetxt("booooo.csv", a, delimiter=",")


# plt.figure()
# plt.plot(tspan, simulation_result.observables['x_obs'], label = 'x')
# plt.plot(tspan, simulation_result.observables['y_obs'], label = 'y')
# # plt.plot(simulation_result2.observables['x_obs'], simulation_result2.observables['y_obs'])
# # plt.plot(simulation_result3.observables['x_obs'], simulation_result3.observables['y_obs'])
# # plt.plot(simulation_result4.observables['x_obs'], simulation_result4.observables['y_obs'])
# plt.xlabel("time", fontsize=16)
# plt.ylabel("conc", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(simulation_result.observables['x_obs'], simulation_result.observables['y_obs'])
# # plt.plot(simulation_result2.observables['x_obs'], simulation_result2.observables['y_obs'])
# # plt.plot(simulation_result3.observables['x_obs'], simulation_result3.observables['y_obs'])
# # plt.plot(simulation_result4.observables['x_obs'], simulation_result4.observables['y_obs'])
# plt.xlabel("x", fontsize=16)
# plt.ylabel("y", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# plt.show()



