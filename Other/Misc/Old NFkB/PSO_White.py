from simplepso.pso import PSO
from Model_NFkB import model
from pysb.integrate import Solver
import numpy as np
import os

exp_data = np.loadtxt('/Users/geenaildefonso/Downloads/nfkbn_traces.csv', delimiter=',', skiprows=1)

exp_avg = exp_data[:,1]
exp_var = exp_data[:,2] **2
exp_time = exp_data[:,0]
equil_time = exp_time[0:25]
exp_time_pts = exp_time[25::]

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters_rules()])

param_log = np.log10(param_values)
bounds_radius = 2


solver1 = Solver(model, equil_time)
solver2 = Solver(model, exp_time_pts)

def likelihood(params):
    params_to_sub = 10**params
    for p, v in zip(model.parameters_rules(), params_to_sub):
        p.value = v

    model.parameters['TNF_ext_0'].value = 0
    solver1.run()
    pre_equil_nnfkb = solver1.yobs['Nuclear_NFkBn']
    last_conc = [solver1.y[:,x][-1] for x in np.arange(len(model.species))]

    solver2.run(y0 = last_conc)
    post_equil_nnfkb = solver2.yobs['Nuclear_NFkBn']
    full_nnfkb = np.concatenate((pre_equil_nnfkb, post_equil_nnfkb))

    error =  np.sum((exp_avg - full_nnfkb) ** 2 / (2 * exp_var))

    return error,

pso = PSO(save_sampled=True)
pso.set_cost_function(likelihood)
pso.set_start_position(param_log)
pso.set_bounds(2)
pso.set_speed(-.25, .25)
values = pso.run(100, 5000)

ranked_params = pso.return_ranked_populations()
np.savetxt('nfkbn_pso_fitness_vals_100_5000_2.txt', values)
np.savetxt('nfkbn_best_vals_100_5000_2.txt', pso.best)
np.savetxt('nfkbn_ranked_pop_vals_100_5000_2.txt', ranked_params)

