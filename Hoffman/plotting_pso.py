import os
# from model_nfkb import model
import csv
import numpy as np

def listdir_fullpath(d):
    """
    Gets a list of path files in directory
    :param d: path to directory
    :return: a list of path of files in directory
    """
    return [os.path.join(d, f) for f in os.listdir(d)]


def list_pars_infile(f, new_path=None):
    """

    :param f: File that contain paths to parameter set values
    :param new_path: parameter paths of f may be different to where they are in the local computer
    :return: list of parameter paths
    """
    par_sets = pd.read_csv(f, names=['parameters'])['parameters'].tolist()
    if new_path:
        par_sets = [w.replace(w.rsplit('/', 1)[0], new_path) for w in par_sets]
    return par_sets

def read_pars(par_path):
    """
    Reads parameter file
    :param par_path: path to parameter file
    :return: Return a list of parameter values from csv file
    """
    f = open(par_path)
    data = csv.reader(f)
    param = [float(d[0]) for d in data]
    return param


def read_all_pars(pars_path, new_path=None):
    """
    Reads all pars in file or directory
    :param new_path:
    :param pars_path: Parameter file or directory path
    :return: DataFrame with all the parameters
    """
    if type(pars_path) is list:
        par_sets = pars_path
    elif os.path.isfile(pars_path):
        par_sets = list_pars_infile(pars_path, new_path)
    elif os.path.isdir(pars_path):
        par_sets = listdir_fullpath(pars_path)
    else:
        raise Exception("Not valid parameter file or path")
    print(par_sets[0])
    pars_0 = read_pars(par_sets[0])
    all_pars = np.zeros((len(par_sets), len(pars_0)))
    all_pars[0] = pars_0
    for idx in range(1, len(par_sets)):
        all_pars[idx] = read_pars(par_sets[idx])
    return all_pars

x = read_all_pars('/Users/geenaildefonso/Projects/NFkB/Hoffman/PSO_results')
x = 10 ** x

exp_data = pd.read_csv('nfkbn_data_jp.csv', delimiter='\t', )
# all_points = [0.0, 10., 30., 60., 90., 120., 180., 300., 480., 960., 1440.]

# this is if you want to train to subsets of the points
# exclude_last_points = all_points[:-4]
# exp_data = exp_data[exp_data['Time'].isin(exclude_last_points)]

# exp_obs_norm = exp_data[data_names].as_matrix().T
# var_norm = exp_data[var_names].as_matrix().T
# std_norm = var_norm ** 0.5

tspan = np.array(exp_data['Time'])

# following 3 lines can be deleted once you see what it was doing below
tmul = 10
ntimes = len(tspan)
old_tspan = np.linspace(tspan[0], tspan[-1], (ntimes - 1) * tmul + 1)[::tmul]

t_to_plot = tspan
t_exp = tspan

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


# def display_jp(position, save_name='nfkb_untrained'):
#     Y = np.copy(position)
for col in x:
    # param_values[rate_mask] = 10 ** Y
    solver.run(col)

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

