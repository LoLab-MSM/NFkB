from pysb import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import scipy.interpolate
import pandas as pd
#from pysb.simulator import ScipyOdeSimulator


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


Model ()


# Declaring Species

Monomer('S')
Monomer('M')
Monomer('CD')
Monomer('E', ['r'])
Monomer('CE')
Monomer('R', ['e','state'], {'state': ['unmod','phos']})

#Initial Conditions

Parameter('Rb_0', 0.55)
Initial(R(e = None, state = 'phos'), Rb_0)
Parameter('S_0', 10)
Initial(S(), S_0)

#Observables

Observable('S_obs', S())
Observable('M_obs', M())
Observable('E_obs', E(r = None))
Observable('CD_obs', CD())
Observable('CE_obs', CD())
Observable('RE_obs', E(r = 1)%R(e = 1, state = 'unmod'))
Observable('R_obs', R(e = None, state = 'unmod'))
Observable('Rp_obs', R(e = None, state = 'phos'))



perturbation_params = {1: {'S_0': 0.3}}


Parameter('km', 1.0)
Parameter('Ks', 0.5)
Parameter('kCDS', 0.45)

Expression('myc_synth', km*S_obs/(Ks + S_obs))
Rule('myc_synth_s', None >> M(), myc_synth)
Expression('cycd_synth', kCDS*S_obs/(Ks + S_obs))
Rule('cycd_synth_s', None >> CD(), cycd_synth)


Parameter('kE', 0.4)
Parameter('KM', 0.15)
Parameter('KE', 0.15)
Parameter('kb', 0.003)
Expression('e2f_synth', kE*(M_obs/(KM + M_obs))*(E_obs/(KE + E_obs)) + kb*M_obs/(KM + M_obs))
Rule('e2f_synth_myc', None >> E(r = None), e2f_synth)


Parameter('kCE', 0.35)
Parameter('kCD', 0.03)
Parameter('kR', 0.18)
Expression('cyce_synth', kCE*E_obs/(KE + E_obs))
Expression('cycd_synth_myc', kCD*M_obs/(KM + M_obs))
Rule('cyce_synth_e2f', None >> CE(), cyce_synth)
Rule('cycd_synth_by_myc', None >> CD(), cycd_synth_myc)
Rule('rb_synth', None >> R(e = None, state = 'unmod'), kR)


Parameter('kp1', 18)
Parameter('kp2', 18)
Parameter('KCD', 0.92)
Parameter('KCE', 0.92)
Parameter('kRE', 180)
Parameter('kDP', 3.6)
Parameter('KRP', 0.01)
Expression('e2f_rp', kp1*CD_obs*RE_obs/(KCD + RE_obs) + kp2*CE_obs*RE_obs/(KCE + RE_obs))
Expression('Rp_R_exp', kp1*CD_obs*R_obs/(KCD + R_obs) + kp2*CE_obs*R_obs/(KCE + R_obs))
Expression('RE_exp', kRE*R_obs*E_obs)
Expression('Rp_to_R_exp', kDP*Rp_obs/(KRP + Rp_obs))
Rule('E_bind_R', E(r = None) + R(e = None, state = 'unmod') >> E(r = 1)%R(e = 1, state = 'unmod'), RE_exp)
Rule('RE_to_E_Rp', E(r = 1)%R(e = 1, state = 'unmod') >> E(r = None) + R(e = None, state = 'phos'), e2f_rp)
Rule('R_to_Rp', R(e = None, state = 'unmod') >> R(e = None, state = 'phos'), Rp_R_exp)
Rule('Rp_to_R', R(e = None, state = 'phos') >> R(e = None, state = 'unmod'), Rp_to_R_exp)


#Degradation Reactions
Parameter('dM', 0.7)
Parameter('dE', 0.25)
Parameter('dCE', 1.5)
Parameter('dCD', 1.5)
Parameter('dR', 0.06)
Parameter('dRP', 0.06)
Parameter('dRE', 0.03)


Rule('myc_deg', M() >> None, dM)
Rule('e2f_deg', E(r = None) >> None, dE)
Rule('cyce_deg', CE() >> None, dCE)
Rule('cycd_deg', CD() >> None, dCD)
Rule('rb_deg', R(e = None, state = 'unmod') >> None, dR)
Rule('rbp_deg', R(e = None, state = 'phos') >> None, dRP)
Rule('rb32f_deg', E(r = 1)%R(e = 1, state = 'unmod') >> None, dRE)

tspan = np.linspace(0, 3000, 30001)
sim = ScipyOdeSimulator(model)
# simulation_result = sim.run()

df = run_sim_with_perturbations(sim, tspan, perturbation_params)

print(model.species)

plt.figure(figsize= (15,5))
plt.subplot(121)
plt.plot(tspan/60, df['__s4'], label = 'CycD')
# plt.plot(tspan, simulation_result.observables['y_obs'], label = 'y')
plt.xlabel("time(hours)", fontsize=16)
plt.ylabel("Concentration uM", fontsize=16)
# plt.ylim(ymin = -2)
# plt.xlim(xmin = -2)
plt.legend(loc=0)
# plt.show()

# plt.figure()
plt.subplot(122)
plt.plot(tspan/60, df['__s5'], label = 'E2F')
# plt.plot(tspan, simulation_result.observables['y_obs'], label = 'y')
plt.xlabel("time(hours)", fontsize=16)
plt.ylabel("Concentration uM", fontsize=16)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
plt.show()


# plt.figure(figsize= (15,5))
# plt.subplot(121)
# plt.plot(tspan/60, simulation_result.observables['E_obs'], label = 'E2F')
# # plt.plot(tspan, simulation_result.observables['y_obs'], label = 'y')
# plt.xlabel("time(hours)", fontsize=16)
# plt.ylabel("Concentration uM", fontsize=16)
# # plt.ylim(ymin = -2)
# # plt.xlim(xmin = -2)
# plt.legend(loc=0)
# # plt.show()
#
# # plt.figure()
# plt.subplot(122)
# plt.plot(tspan/60, simulation_result.observables['CD_obs'], label = 'CycD')
# # plt.plot(tspan, simulation_result.observables['y_obs'], label = 'y')
# plt.xlabel("time(hours)", fontsize=16)
# plt.ylabel("Concentration uM", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# plt.show()