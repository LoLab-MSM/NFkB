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
# obs_names = ['mBid', 'cPARP']
# data_names = ['norm_ICRP', 'norm_ECRP']
# var_names = ['nrm_var_ICRP', 'nrm_var_ECRP']

obs_names = ['IkBa_obs', 'IkBb_obs', 'IkBe_obs', 'IkBd_obs']
data_names = ['ikba_rep1', 'ikbb_rep1','ikbe_rep1','ikbd_rep1']
var_names = ['ikba_rep2', 'ikbb_rep2','ikbe_rep2','ikbd_rep2']
# Total starting amounts of proteins in obs_names, for normalizing simulations
obs_totals = [model.parameters['IkBa_0'].value, model.parameters['IkBb_0'].value, model.parameters['IkBe_0'].value, model.parameters['IkBd_0'].value]
exp_data = np.genfromtxt('/Users/geenaildefonso/nfkbn_data.csv', delimiter=',', names=True)
# directory = os.path.dirname(__file__)
# data_path = os.path.join(directory, 'data',
#                          'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
# exp_data = np.genfromtxt(data_path, delimiter=',', names=True)

# Model observable corresponding to the IMS-RP reporter (MOMP timing)
# momp_obs = 'aSmac'
# # Mean and variance of Td (delay time) and Ts (switching time) of MOMP, and
# # yfinal (the last value of the IMS-RP trajectory)
# momp_obs_total = model.parameters['Smac_0'].value
# momp_data = np.array([9810.0, 180.0, momp_obs_total])
# momp_var = np.array([7245000.0, 3600.0, 1e4])

ntimes = len(exp_data['Time'])
tmul = 10
tspan = np.linspace(exp_data['Time'][0], exp_data['Time'][-1],
                    (ntimes - 1) * tmul + 1)

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
k_ids = [p.value for p in model.parameters_rules()]
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
bounds_radius = 2
solver = Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6, )

def display(position):
    exp_obs_norm = exp_data[data_names].view(float).reshape(len(exp_data), -1).T
    var_norm = exp_data[var_names].view(float).reshape(len(exp_data), -1).T
    std_norm = var_norm ** 0.5
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    obs_names_disp = obs_names #+ ['aSmac']
    sim_obs = solver.yobs[obs_names_disp].view(float).reshape(len(solver.yobs), -1)
    totals = obs_totals #+ [momp_obs_total]
    sim_obs_norm = (sim_obs / totals).T
    colors = ('r', 'b', 'g', 'o')
    for exp, exp_err, sim, c in zip(exp_obs_norm, std_norm, sim_obs_norm, colors):
        plt.plot(exp_data['Time'], exp, color=c, marker='.', linestyle=':')
        plt.errorbar(exp_data['Time'], exp, yerr=exp_err, ecolor=c, elinewidth=0.5, capsize=0)
        plt.plot(solver.tspan, sim, color=c)
    plt.plot(solver.tspan, sim_obs_norm[2], color='g')
    # plt.vlines(momp_data[0], -0.05, 1.05, color='g', linestyle=':')
    plt.savefig('nfkb_trained.png')
    plt.show()

def likelihood(position):
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
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
    # ysim_momp = solver.yobs[momp_obs]
    # if np.nanmax(ysim_momp) == 0:
    #     print 'No aSmac!'
    #     ysim_momp_norm = ysim_momp
    #     t10 = 0
    #     t90 = 0
    # else:
    #     ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
    #     st, sc, sk = scipy.interpolate.splrep(tspan, ysim_momp_norm)
    #     try:
    #         t10 = scipy.interpolate.sproot((st, sc - 0.10, sk))[0]
    #         t90 = scipy.interpolate.sproot((st, sc - 0.90, sk))[0]
    #     except IndexError:
    #         t10 = 0
    #         t90 = 0
    # td = (t10 + t90) / 2
    # ts = t90 - t10
    # yfinal = ysim_momp[-1]
    # momp_sim = [td, ts, yfinal]
    # e3 = np.sum((momp_data - momp_sim) ** 2 / (2 * momp_var)) / 3
    error = e1 + e2 + e3 + e4

    return error,

def run_example():
    pso_fn = PSO(save_sampled=False, verbose=True, num_proc=4)
    pso_fn.set_cost_function(likelihood)
    pso_fn.set_start_position(xnominal)
    pso_fn.set_bounds(2)
    pso_fn.set_speed(-.25, .25)
    pso_fn.run(20, 100)
    display(pso_fn.best)

if __name__ == '__main__':
    run_example()
#
#
# if '__main__' == __name__:
#     run_example()



# from simplepso.pso import PSO
# from Model_NFkB import model
# from pysb.integrate import Solver
# import numpy as np
# import os
#
# exp_data = np.loadtxt('/Users/geenaildefonso/Downloads/nfkbn_traces.csv', delimiter=',', skiprows=1)
#
# exp_avg = exp_data[:,1]
# exp_var = exp_data[:,2] **2
# exp_time = exp_data[:,0]
# equil_time = exp_time[0:25]
# exp_time_pts = exp_time[25::]
#
# rate_params = model.parameters_rules()
# param_values = np.array([p.value for p in model.parameters_rules()])
#
# param_log = np.log10(param_values)
# bounds_radius = 2
#
#
# solver1 = Solver(model, equil_time)
# solver2 = Solver(model, exp_time_pts)
#
# def likelihood(params):
#     params_to_sub = 10**params
#     for p, v in zip(model.parameters_rules(), params_to_sub):
#         p.value = v
#
#     model.parameters['TNF_ext_0'].value = 0
#     solver1.run()
#     pre_equil_nnfkb = solver1.yobs['Nuclear_NFkBn']
#     last_conc = [solver1.y[:,x][-1] for x in np.arange(len(model.species))]
#
#     solver2.run(y0 = last_conc)
#     post_equil_nnfkb = solver2.yobs['Nuclear_NFkBn']
#     full_nnfkb = np.concatenate((pre_equil_nnfkb, post_equil_nnfkb))
#
#     error =  np.sum((exp_avg - full_nnfkb) ** 2 / (2 * exp_var))
#
#     return error,
#
# pso = PSO(save_sampled=True)
# pso.set_cost_function(likelihood)
# pso.set_start_position(param_log)
# pso.set_bounds(2)
# pso.set_speed(-.25, .25)
# values = pso.run(100, 5000)
#
# ranked_params = pso.return_ranked_populations()
# np.savetxt('nfkbn_pso_fitness_vals_100_5000_2.txt', values)
# np.savetxt('nfkbn_best_vals_100_5000_2.txt', pso.best)
# np.savetxt('nfkbn_ranked_pop_vals_100_5000_2.txt', ranked_params)