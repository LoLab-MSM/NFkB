# -*- coding: utf-8 -*-


import matplotlib, os, sys
import numpy as np
import pandas as pd
# sys.path.append('/Users/geenaildefonso/Projects/ParticleSwarmOptimization')
from simplepso.pso import PSO

r = os.system('python -c "import matplotlib.pyplot as plt;plt.figure()"')
if r != 0:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    show = False
else:
    import matplotlib.pyplot as plt

    show = True

# from exported_nnrm import model
from calibration_nfkb import model
from pysb.integrate import Solver
# from pysb.simulator.scipyode import ScipyOdeSimulator
show = False
# you can do this in pycharm by right clicking on folder and
# do "Mark directory as" and do "Source"


# defining observable names from model to be used for calibration
# data_names and var_names from csv data file for all time points
obs_names = ['IkBa_obs', 'IkBb_obs', 'IkBe_obs', 'IkBd_obs']
data_names = ['ikba_avg', 'ikbb_avg', 'ikbe_avg', 'ikbd_avg']
var_names = ['ikba_var', 'ikbb_var', 'ikbe_var', 'ikbd_var']

# Total starting amounts of proteins in obs_names, for normalizing simulations
obs_totals = [model.parameters['IkBa_0'].value,
              model.parameters['IkBb_0'].value,
              model.parameters['IkBe_0'].value,
              model.parameters['IkBd_0'].value]

# handle experimental data
exp_data = pd.read_csv('nfkbn_data_jp.csv', delimiter='\t', )
# all_points = [0.0, 10., 30., 60., 90., 120., 180., 300., 480., 960., 1440.]

# this is if you want to train to subsets of the points
# exclude_last_points = all_points[:-4]
# exp_data = exp_data[exp_data['Time'].isin(exclude_last_points)]

exp_obs_norm = exp_data[data_names].as_matrix().T
var_norm = exp_data[var_names].as_matrix().T
std_norm = var_norm ** 0.5

tspan = np.array(exp_data['Time'])

# following 3 lines can be deleted once you see what it was doing below
tmul = 10
ntimes = len(tspan)
old_tspan = np.linspace(tspan[0], tspan[-1], (ntimes - 1) * tmul + 1)[::tmul]

t_to_plot = tspan
t_exp = tspan
# before they would not match up! So when you compare two vectors, they were
# not aligned, thus give crazy errors! Only first and last time point aligned
print("original sim time {}".format(old_tspan))
print("algined sim time {}".format(t_to_plot))
print("exp time {}".format(t_exp))

# deprecated solver
solver = Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6)

# this is the solver you want to move towards
# the integrate.solver is being deprecated

# solver = ScipyOdeSimulator(model=model, tspan=tspan, integrator='vode',
#                            integrator_options={'rtol':1e-6,
#                                                'atol':1e-6})

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])

rate_mask = np.array(
    [i for i, p in enumerate(model.parameters) if p in rate_params])

k_ids = [p.value for p in model.parameters_rules()]
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
bounds_radius = 2


def display_jp(position, save_name='nfkb_untrained'):
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)

    obs1 = solver.yobs['IkBa_obs']
    obs1 = obs1 / model.parameters['IkBa_0'].value

    obs2 = solver.yobs['IkBb_obs']
    obs2 = obs2 / model.parameters['IkBb_0'].value

    obs3 = solver.yobs['IkBe_obs']
    obs3 /= model.parameters['IkBe_0'].value

    obs4 = solver.yobs['IkBd_obs']
    obs4 /= model.parameters['IkBd_0'].value

    # create subplots of each to see what is going on

    plt.figure(figsize=(10, 8))
    plt.subplot(411)
    plt.plot(t_to_plot, obs1, 'b-o', label='IkBa_obs')
    plt.errorbar(t_exp, exp_obs_norm[0], yerr=std_norm[0], ecolor='r',
                 color='r', elinewidth=0.5, capsize=0, label=data_names[0])
    plt.legend(loc=0)

    plt.subplot(412)
    plt.plot(t_to_plot, obs2, 'b-o', label='IkBb_obs')
    plt.errorbar(t_exp, exp_obs_norm[1], yerr=std_norm[1], ecolor='r',
                 color='r', elinewidth=0.5, capsize=0, label=data_names[1])
    plt.legend(loc=0)

    plt.subplot(413)
    plt.plot(t_to_plot, obs3, 'b-o', label='IkBe_obs')
    plt.errorbar(t_exp, exp_obs_norm[2], yerr=std_norm[2], ecolor='r',
                 color='r', elinewidth=0.5, capsize=0, label=data_names[2])
    plt.legend(loc=0)

    plt.subplot(414)
    plt.plot(t_to_plot, obs4, 'b-o', label='IkBd_obs')
    plt.errorbar(t_exp, exp_obs_norm[3], yerr=std_norm[3], ecolor='r',
                 color='r', elinewidth=0.5, capsize=0, label=data_names[3])
    plt.legend(loc=0)

    plt.tight_layout()
    plt.savefig('{}.png'.format(save_name), bbox_inches='tight')
    if show:
        plt.show()
    plt.close()
    observable = [obs1, obs2, obs3, obs4]

    error = calculate_error(observable)

    print('')
    for n, name in enumerate(obs_names):
        print("{} error is {}".format(name, error[n]))


def calculate_error(obs_list):
    error = []
    for ysim_norm, data_name, var_name in zip(obs_list, data_names, var_names):
        ydata = exp_data[data_name]
        yvar = exp_data[var_name]
        error.append(
            np.sum((ydata[1:] - ysim_norm[1:]) ** 2 / (2 * yvar[1:])) / len(
                ydata[1:]))

    return error


def likelihood(position):
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values=param_values)
    obs1 = solver.yobs['IkBa_obs']
    obs1 /= model.parameters['IkBa_0'].value

    obs2 = solver.yobs['IkBb_obs']
    obs2 /= model.parameters['IkBb_0'].value

    obs3 = solver.yobs['IkBe_obs']
    obs3 /= model.parameters['IkBe_0'].value

    obs4 = solver.yobs['IkBd_obs']
    obs4 /= model.parameters['IkBd_0'].value

    observable = [obs1, obs2, obs3, obs4]
    error = calculate_error(observable)
    error = np.array(error).sum()
    return error,


if __name__ == '__main__':

    display_jp(position=xnominal)
    # quit()

    out_dir = 'PSO_results'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # create 100 best fits parameter sets
    for i in range(100):
        pso_fn = PSO(save_sampled=False, verbose=False, num_proc=16)
        print("running {}".format(i))
        pso_fn.set_cost_function(likelihood)
        pso_fn.set_start_position(xnominal)
        pso_fn.set_bounds(3)
        pso_fn.set_speed(-.25, .25)
        pso_fn.run(100, 300)  # particles in swarm and iterations
        print("Iteration = {} : best value = {}".format(i, pso_fn.best.fitness.values[0]))
        csv_out = os.path.join(out_dir, 'best_fit_{}.csv'.format(i))
        fig_out = os.path.join(out_dir, 'best_fit_{}'.format(i))
        np.savetxt(csv_out, pso_fn.best)
        display_jp(pso_fn.best, fig_out)
