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
from pysb.simulator import ScipyOdeSimulator
import numpy as np


def run_sim_with_perturbations(sim, tspan, perturbations):
    perturbation_times = sorted(perturbations.keys())

    for t in perturbation_times:
        if not any(np.isclose(t, tspan)):
            raise ValueError('Perturbation time %f is not defined in tspan'
                             % t)

    # Create the tspans for each simulation stage
    tspans = np.split(tspan, [np.abs(tspan-value).argmin() for value in
                              perturbation_times])
    for i in range(len(tspans) - 1):
        tspans[i] = np.append(tspans[i], tspans[i + 1][0])

    # Do an initial simulation with the above setup until the point at which
    # we want to change a parameter (perturbation_time)
    res = sim.run(tspan=tspans[0])

    df_out = res.dataframe

    for i, t in enumerate(perturbation_times, 1):
        res = sim.run(initials=res.species[-1],
                      param_values=perturbation_params[t],
                      tspan=tspans[i])

        df_out = pd.concat([df_out, res.dataframe.iloc[1:]])

    return df_out



time = tuple(range(0, 721))
# print(time)
data = pd.read_csv('/Users/geenaildefonso/data_ikk.csv')
IKK_flux = tuple(data.loc[:, "IKK_flux"])
# print(IKK.shape)
# print(IKK)




Model ()


Monomer('IKK')
Parameter('IKK_0', 0.0002) #TRADD-TRAF-RIP
Initial(IKK(), IKK_0)



Parameter('on', 0)

# for i in tuple(range(721)):
perturbation_params = {
    1: {'on': [w for w in IKK_flux]}
}

Rule('none_IKK', None >> IKK(), on)

tspan = np.linspace(0, 720, 721)
# print(len(tspan))
# print(tspan)
sim = ScipyOdeSimulator(model)
# simulation_result = sim.run()

df = run_sim_with_perturbations(sim, tspan, perturbation_params)
print(df.shape)
print(df)
