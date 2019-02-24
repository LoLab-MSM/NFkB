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
from pysb.simulator.bng import BngSimulator
from nfkb_complex import model

tspan = np.linspace(0, 720, 721)
sim = BngSimulator(model, tspan = tspan)
simulation_result = sim.run(method='ode')


plt.figure(figsize = (15,10))
# plt.figure()
plt.subplot(231)
plt.plot(tspan/60, simulation_result.observables['TNF_obs'], label = 'TNF')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(232)
plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], label = 'TNFR')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(233)
plt.plot(tspan/60, simulation_result.observables['IKKa_obs'],label = 'IKKa')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(234)
plt.plot(tspan/60, simulation_result.observables['NFkBn_obs'], label = 'NFkBn')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(235)
plt.plot(tspan/60, simulation_result.observables['obs_A20t'], label = 'A20')

plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(236)
plt.plot(tspan/60, simulation_result.observables['IkBa_obs'],label = 'IkBa')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
plt.savefig('bng_nfkb')
plt.tight_layout()
plt.show()