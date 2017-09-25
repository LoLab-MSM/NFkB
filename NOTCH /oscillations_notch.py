from pysb.core import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import ScipyOdeSimulator as SOS


from pysb import Model, Monomer, Initial, Parameter, Observable, Rule
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import pandas as pd


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


Model()
# Monomer('ICD', ['csl', 'l'], {'l': ['cyto', 'nuc']})
# Monomer('CSL', ['icd', 'gh', 'gsl'])
# Monomer('gHes1', ['ic', 'ph'])
Monomer('mHes1')
Monomer('pHES1')

# Parameter('gHes1_0', 2)
# Initial(gHes1(ic=None, ph=None), gHes1_0)
# Parameter('mHes1_0', 20)
# Initial(mHes1(l = 'nuc'), mHes1_0)
Parameter('pHes1_0', 0.0)
Initial(pHES1(), pHes1_0)

# Defining Observables
# Observable('ICD_obs', ICD(csl=None, l='nuc'))
# Observable('CSL_on', ICD(csl=1, l='nuc') % CSL(icd=1, gh=None, gsl=None))
# Observable('Ex_L', L_ex(r_in=None))
# Observable('In_R', R_in(l_ex=None, ysc = None))
# Observable('ExL_InR', L_ex(r_in=1) % R_in(l_ex=1))
# Observable('LR_on', R_in(l_ex=1, ysc = None) % L_ex(r_in=1))
# Observable('CSL_off', CSL(icd=None, gh=None, gsl=None))
Observable('mHes1_Hes1', mHes1())
Observable('protein_HES1', pHES1())
# Observable('prot_gene',gHes1(ic = None, ph = 1) % pHES1(gh = 1, ga1 = None, gn2 = None))


#synthesis of Hes1
Parameter('a', 0.1)
Parameter('b', 0.0)
Parameter('c', 0.15)
Parameter('k', 0.015)
Parameter('pcrit', 40)
Parameter('hill', 2)
Parameter('gp_kf', 0.01)
Parameter('gp_kr', 1.0)
Parameter('gm_kf', 0.1)
Parameter('mp_kf', 0.01)

Expression('m_exp', b/(1.0 + (protein_HES1/25.0)**2))
# Expression('p_exp', a/(0.031+ 0.031*protein_HES1**2))
Expression('p_exp', a/(0.031+ 0.031*protein_HES1**2))

# Rule('g_to_p', gHes1(ic = None, ph = None) + pHES1(gh = None, ga1 = None, gn2 = None) <> gHes1(ic = None, ph = 1) % pHES1(gh = 1, ga1 = None, gn2 = None), gp_kf, gp_kr)
Rule('g_to_m', None >>  mHes1(), m_exp)
Rule('m_to_p', mHes1() >> pHES1() + mHes1(), a)
Rule('hes1_mrna_deg', mHes1() >> None, c)
Rule('hes1_deg', pHES1() >> None, k)
# Expression('hes1_exp', k/((protein_HES1)**hill + 1))
# Expression('hes1mrna_exp', a*mHes1_Hes1)
# Rule('hes1_mrna', None >> mHes1(l = 'nuc'), hes1_exp)
# Rule('hes1_prot', None >> pHES1(gh = None, ga1 = None, gn2 = None), hes1mrna_exp)


# Define perturbations by their time points
perturbation_params = {
    15: {'b': 1.0},
    # 30: {'a': 0.1}, #
    # 40: {'c': 0.15},  # 0.15
    50: {'k': 0.051}, # 0.051 50 #180

}

# End define perturbation

# Define simulations
tspan = np.linspace(0, 120, 1201)
sim = ScipyOdeSimulator(model)
df = run_sim_with_perturbations(sim, tspan, perturbation_params)
print(df)


# quit()
#
# sim1 = SOS(model, tspan)
# sim_result = sim1.run()
#
# print('phes1')
# print(sim_result.observables['protein_HES1'])
# print('mhes1')
# print(sim_result.observables['mHes1_Hes1'])
#
# for i,sp in enumerate(model.species):
#     print i,":",sp
#
# print(len(model.species))
# print(len((model.parameters)))
# print(len(model.odes))
# print(len(model.rules))
# for  j,ode in enumerate(model.odes):
#    for i in range(len(model.species)):
#       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
#    print j,":",ode

plt.figure()
# plt.figure(figsize=(15, 10))
# plt.subplot(121)
plt.plot(tspan, df['mHes1_Hes1'], label='mHes1_Hes1')
plt.plot(tspan, df['protein_HES1'], label='protein_HES1')
# plt.plot(tspan, sim_result.observables['prot_gene'], label='bound gene')
plt.xlabel("Time (in min)", fontsize=10)
plt.ylabel("Molecules per Cell", fontsize=10)
# plt.ylim(ymax = 21)
plt.legend(loc=0)
#
# plt.subplot(122)
# plt.plot(tspan/60, sim_result.observables['protein_HES1'], label='protein_HES1')
# plt.plot(tspan/60, sim_result.observables['CSL_on'], label='ICDCSL')
# plt.xlabel("Time (in min)", fontsize=10)
# plt.ylabel("Molecules per Cell", fontsize=10)
# plt.legend(loc=0)

plt.show()