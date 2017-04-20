from pylab import *
from ddeint import ddeint
import matplotlib as plt

# def model(Y,t,d):
#     # print(tuple(Y(t)))
#     x,y = tuple(Y(t))
#     print(x,y)
#     xd,yd = tuple(Y(t-d))
#     return array([0.5*x-0.5*x*y, -0.5*y + 0.5*xd*yd])
#
# g = lambda t : array([1,2])
# tt = linspace(2,30,20000)
#
# for d in [0]:
#     yy = ddeint(model,g,tt,fargs=(d,))
#     # WE PLOT X AGAINST Y
#     plot(yy[:,0], yy[:,1], lw=2, label='delay = %.01f'%d)
#
# legend()


def model(Y,t,d):
    x,y = Y(t)
    xd,yd = Y(t-d)
    return array([0.5*x*(1-yd), -0.5*y*(1-xd)])

g = lambda t : array([1,2])
tt = linspace(2,30,20000)

for d in [0, 0.2]:
    yy = ddeint(model,g,tt,fargs=(d,))
    # WE PLOT X AGAINST Y
    plot(yy[:,0], yy[:,1], lw=2, label='delay = %.01f'%d)
legend() # display the legend
show()