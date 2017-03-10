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


Model ()


Monomer('X')
Monomer('Y')

Parameter('kf', 0.5)
Parameter('x_0', 1)
Parameter('y_0', 2)

Initial(X(), x_0)
Initial(Y(), y_0)

Rule('x_to_xx', X() >> X() + X(), kf)
Rule('y_deg', Y() >> None, kf)
Rule('xy_to_yy', X() + Y() >> Y() + Y(), kf)
Observable('x_obs', X())
Observable('y_obs', Y())

tspan = np.linspace(0,100,1001)
w = odesolve(model, tspan,verbose=True)
#
# print(w['x_obs'])
# print()
# print(w['y_obs'])



plt.figure()
plt.plot(tspan, w['x_obs'], lw = 2, color = 'b',label = 'x')
plt.plot(tspan, w['y_obs'],lw = 2, color = 'r',label = 'y')
plt.legend(loc = 0)


plt.figure()
plt.plot(w['y_obs'], w['x_obs'], lw = 2)
# plt.plot(w['x_obs'], w['y_obs'], lw = 2, color = 'r')
plt.xlim(xmax = 5, xmin = 0)
plt.ylim(ymax = 3.5, ymin =0)
plt.xlabel('y obs')
plt.ylabel('x obs')
plt.show()
