# import pydelay and numpy and pylab
import numpy as np
import pylab as pl
from pydelay import dde23

# define the equations
eqns = {
    'x' : '0.5 * x*(1 - y(t - tau1))',
    'y' : '-0.5 * y*(1 - x(t - tau2))'
    }

#define the parameters
params = {
    'tau1': 15,
    'tau2': 15
    }

# Initialise the solver
dde = dde23(eqns=eqns, params=params)

# set the simulation parameters
# (solve from t=0 to t=1000 and limit the maximum step size to 1.0)
dde.set_sim_params(tfinal=30, dtmax=0.1)

# set the history of to the constant function 0.5 (using a python lambda function)
histfunc = {
    'x': lambda t: 0.5,
    'y': lambda t: 0.5
    }
dde.hist_from_funcs(histfunc, 51)

# run the simulator
dde.run()

# Make a plot of x(t) vs x(t-tau):
# Sample the solution twice with a stepsize of dt=0.1:

# once in the interval [515, 1000]
sol1x = dde.sample(515, 1000, 0.1)
x1 = sol1x['x']
sol1y = dde.sample(515, 1000, 0.1)
y1 = sol1y['y']

# and once between [500, 1000-15]
sol2x = dde.sample(500, 1000-15, 0.1)
x2 = sol2x['x']
sol2y = dde.sample(500, 1000-15, 0.1)
y2 = sol2y['x']

pl.figure()
pl.plot(x1, y1)
pl.xlabel('$x(t)$')
pl.ylabel('$y(t)$')

pl.figure()
pl.plot(x2, y2)
pl.xlabel('$x(t - 15)$')
pl.ylabel('$y(t - 15)$')
pl.show()