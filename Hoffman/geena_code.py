# import numpy as np
# x = np.searchsorted([1,2,3,4,5], 3.5)
# print(x)
# quit()
# dummy data
# import pandas as pd
# import numpy as np
# # inner_dict_values = [5*i for i in range(10)]
# # print(inner_dict_values)
# # construct dictionary with format:
# # {int: {'some_key': float} }
#
#
# # time = tuple(range(0, 721))
# # print(time)
# data = pd.read_csv('/Users/geenaildefonso/data_ikk.csv')
# IKK_flux = tuple(data.loc[:, "IKK_flux"])
# print(len(IKK_flux))
# print(np.transpose(IKK_flux))
# # print(IKK.shape)
# # print(IKK)
#
#
# my_dict = {}
#
# # you have ten values, so loop using range
# for i in range(721):
# 	my_dict[i] = {'some_key': IKK_flux[i]}
#
# print("dictionary sorted\n")
# print(sorted(my_dict.items(), key=lambda x: x[0]))

# -----------------------------
# convert outer key values into inner key values
# {int: {'some_key': float} } --> {float: {'some_key': int} }


# make new dict
# my_new_dict = {}
#
# for keys in my_dict:
# 	my_new_key = my_dict[keys]['some_key']
# 	my_new_dict[my_new_key] = {'some_key': keys}
#
# print("\n\nnew dictionary sorted\n")
# print(sorted(my_new_dict.items(), key=lambda x: x[0]))



from pylab import *
from ddeint import ddeint
import matplotlib as plt
# model = lambda Y, t:  -Y( t-3*cos( Y(t) )**2 )
# tt = linspace(0, 30, 2000)
# yy = ddeint(model, lambda t:1, tt)
# print('plottin...')
# plot(tt, yy)

def model(Y,t,d):
    # print(tuple(Y(t)))
    x,y = tuple(Y(t))
    print(x,y)
    xd,yd = tuple(Y(t-d))
    return array([0.5*x-0.5*x*y, -0.5*y + 0.5*xd*yd])

g = lambda t : array([1,2])
tt = linspace(2,30,20000)

for d in [0]:
    print('d', d)
    yy = ddeint(model,g,tt,fargs=(d,))
    # WE PLOT X AGAINST Y
    plot(yy[:,0], yy[:,1], lw=2, label='delay = %.01f'%d)
#
legend()


# from pylab import *
# from ddeint import ddeint
#
# model = lambda Y,t : Y(t - 3*pi/2) # Model
# tt = linspace(0,50,10000) # Time start, time end, nb of points/steps
# g=sin # Expression of Y(t) before the integration interval
# yy = ddeint(model,g,tt) # Solving
#
# # PLOTTING
# plot(tt,yy,c='r',label="$y(t)$")
# plot(tt,sin(tt),c='b',label="$sin(t)$")
# set_ylim(ymax=2) # make room for the legend
# legend()