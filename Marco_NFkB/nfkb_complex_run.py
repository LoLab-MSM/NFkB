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
import numpy as np
from nfkb_complex import model

tspan = np.linspace(0, 720, 721)
sim = ScipyOdeSimulator(model, tspan = tspan)
simulation_result = sim.run()

plt.figure(figsize = (15,10))
# plt.figure()
plt.subplot(231)
plt.plot(tspan/60, simulation_result.observables['TNF_obs'], marker = '*',label = 'TNF')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(232)
plt.plot(tspan/60, simulation_result.observables['TNFR_obs'],marker = '*',label = 'TNFR')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(233)
plt.plot(tspan/60, simulation_result.observables['IKKa_obs'],marker = '*', label = 'IKKa')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(234)
plt.plot(tspan/60, simulation_result.observables['NFkBn_obs'], marker = '*', label = 'NFkBn')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(235)
plt.plot(tspan/60, simulation_result.observables['obs_A20t'],marker = '*', label = 'A20')

plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(236)
plt.plot(tspan/60, simulation_result.observables['IkBa_obs'], marker = '*',label = 'IkBa')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Concentration", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.tight_layout()
plt.show()