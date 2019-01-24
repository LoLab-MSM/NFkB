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


Model()


#Making X and Y monomers
Monomer('X')
Monomer('Y')

Parameter('x_0', 10)
Parameter('y_0', 1)

Initial(X(), x_0)
Initial(Y(), y_0)

Observable('x_obs', X())
Observable('y_obs', Y())

Parameter('b', .05)
Parameter('a', 1.58)
Parameter('c', 50)
Parameter('n', 2)
Parameter('d', 12)
Parameter('h', 1)
Parameter('si')

Expression('x_n', d*a*(1 + c*x_obs**n)/(1+x_obs**n + si*y_obs**n))
# Expression('y_n', d*a*b*(1 + c*x_obs**n)/(1+x_obs**n))


# Expression('x_n', d*a*(-1 - c*(x_obs)**n)/(-1-(x_obs)**n - y_obs**n))
Expression('y_n', d*a*b*(1 + c*x_obs**n)/(1+x_obs**n))
# Expression('y_i', I/(1+y_obs**n))
# Expression('y_n', a*y_obs**n/(1+y_obs**n))
# Expression('x_i', I/(1+x_obs**n))


#normal system
Rule('synth_x_b', X() >> None, d)
Rule('synth_x_obs', None >> X(), x_n)

Rule('synth_y_obs', None >> Y(), y_n)
Rule('y_deg', Y() >> None, h)


#for time reversed system
# Rule('synth_x_b', None >> X(), d)
# Rule('synth_x_obs', X() >> None, x_n)
#
# Rule('synth_y_obs', Y() >> None, y_n)
# Rule('y_deg', None >> Y(), h)

tspan = np.linspace(0, 120, 400)
sim = ScipyOdeSimulator(model, tspan = tspan)
simulation_result = sim.run()

print(simulation_result.observables['x_obs'])
print()
print(simulation_result.observables['y_obs'])
# x1 = simulation_result.observables['x_obs']
# y1 = simulation_result.observables['y_obs']
# # print(simulation_result.species)
# df = simulation_result.dataframe
# data = df.loc[0:, ['x_obs', 'y_obs']]
# data.to_csv('model_data.csv', sep = '\t')


# raw_data = {'xlow_yhigh': [df.loc[0:, ['x_obs']]], }
# for i,ode in enumerate(model.odes):
#     print i,":",ode

# # for  j,ode in enumerate(model.odes):
# #     for i in range(len(model.species)):
# #        ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
# #     print j,":",ode

# print(model.parameters)
# print(model.species)


plt.figure(figsize= (15,5))
plt.subplot(121)
plt.plot(simulation_result.observables['x_obs'],simulation_result.observables['y_obs'], label = 'x')
# plt.plot(tspan, simulation_result.observables['y_obs'], label = 'y')
plt.xlabel("x", fontsize=16)
plt.ylabel("y", fontsize=16)
plt.ylim(ymin = -2)
plt.xlim(xmin = -2)
plt.legend(loc=0)
# plt.show()

# plt.figure()
plt.subplot(122)
plt.plot(tspan, simulation_result.observables['x_obs'], label = 'x')
plt.plot(tspan, simulation_result.observables['y_obs'], label = 'y')
plt.xlabel("time", fontsize=16)
plt.ylabel("conc", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
plt.show()