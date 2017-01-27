# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014
@author: james
"""
import numpy as np
import scipy.interpolate

try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False
    pass
import os
from pysb.integrate import Solver
import sys
sys.path.append('/Users/geenaildefonso/Projects/ParticleSwarmOptimization')
from exported_nnrm import model
from simplepso.pso import PSO


#defining observable names from model to be used for calibration
#data_names and var_names from csv data file for all time points
obs_names = ['IkBa_obs', 'IkBb_obs', 'IkBe_obs', 'IkBd_obs']
data_names = ['ikba_avg', 'ikbb_avg','ikbe_avg','ikbd_avg']
# print(data_names)
var_names = ['ikba_var', 'ikbb_var','ikbe_var','ikbd_var']
# print(var_names)
# Total starting amounts of proteins in obs_names, for normalizing simulations
obs_totals = [model.parameters['IkBa_0'].value, model.parameters['IkBb_0'].value, model.parameters['IkBe_0'].value, model.parameters['IkBd_0'].value]
exp_data = np.genfromtxt('/Users/geenaildefonso/nfkbn_data.csv', delimiter=',', names=True)

ntimes = len(exp_data['Time'])
tmul = 10
tspan = np.linspace(exp_data['Time'][0], exp_data['Time'][-1],
                    (ntimes - 1) * tmul + 1)

solver = Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6)
rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
# print(len(param_values))
rate_mask = np.array([i for i,p in enumerate(model.parameters) if p in rate_params])
# print(rate_mask)
# print(len(rate_mask))
# print(len(model.parameters_rules()))
# print(len(model.parameters))
k_ids = [p.value for p in model.parameters_rules()]
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
bounds_radius = 2

print(len(model.parameters))

def display(position):
    exp_obs_norm = exp_data[data_names].view(float).reshape(len(exp_data), -1).T
    print(exp_obs_norm)
    var_norm = exp_data[var_names].view(float).reshape(len(exp_data), -1).T
    print(var_norm)
    std_norm = var_norm ** 0.5
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    obs_names_disp = obs_names
    sim_obs = solver.yobs[obs_names_disp].view(float).reshape(len(solver.yobs), -1)
    totals = obs_totals
    sim_obs_norm = (sim_obs / totals).T
    colors = ('r', 'b', 'g', 'y')
    for exp, exp_err, sim, c, name, in zip(exp_obs_norm, std_norm, sim_obs_norm, colors, data_names):
        plt.plot(exp_data['Time'], exp, color=c, marker='.', linestyle=':', label=name)
        plt.errorbar(exp_data['Time'], exp, yerr=exp_err, ecolor=c, elinewidth=0.5, capsize=0)
        plt.plot(tspan, sim, color=c, label=name)
        plt.legend()
    plt.savefig('nfkb_trained.png')
    plt.show()

def likelihood(position):
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values=param_values)
    for obs_name, data_name, var_name, obs_total in \
            zip(obs_names, data_names, var_names, obs_totals):
        ysim = solver.yobs[obs_name][::tmul]
        ysim_norm = ysim / obs_total
        ydata = exp_data[data_name]
        yvar = exp_data[var_name]
        if obs_name == 'IkBa_obs':
            e1 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
        elif obs_name == 'IkBb_obs':
            e2 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
        elif obs_name == 'IkBe_obs':
            e3 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
        else:
            e4 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
    print("error")
    error = e1 + e2 + e3 + e4

    return error,

def run_example():
    # display(position=xnominal)
    # quit()
    pso_fn = PSO(save_sampled=False, verbose=True, num_proc=4)
    pso_fn.set_cost_function(likelihood)
    pso_fn.set_start_position(xnominal)
    pso_fn.set_bounds(2)
    pso_fn.set_speed(-.25, .25)
    pso_fn.run(20, 100) #particles in swarm and iterations
    display(pso_fn.best)

if __name__ == '__main__':
    run_example()
