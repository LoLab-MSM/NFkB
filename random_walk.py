# import numpy
#
# # Populate the lattice
# lattice = numpy.concatenate([numpy.ones(90), numpy.zeros(10)])
# numpy.random.shuffle(lattice)
#
# # Intialize problem
# in_trap = False
# steps = 0
# pos = int(numpy.random.randint(0,len(lattice),1))
# history = []
#
# while in_trap == False:
#     # Step of -1 is backward, 1 is forward
#     step = numpy.random.permutation([-1,1])[0]
#
#     # Check position for edges and fix if required
#     if pos + step > len(lattice) - 1:
#         pos = 0
#     elif pos + step < 0:
#         pos = len(lattice) - 1
#     else:
#         pos += step
#
#     # Keep track of random walk
#     history.append(pos)
#
#     # Check if it's a trap
#     if lattice[pos] == 0:
#         in_trap = True
#
#     # If not, continue
#     steps += 1
#
#
# print steps
# print history
# print lattice
# import plotly.plotly as py
# from plotly.tools import FigureFactory as FF
# import plotly.graph_objs as go
#
# import numpy as np
# import pandas as pd
# import scipy
# import random
#
# x = [0]
#
# for j in range(100):
#     step_x = random.randint(0, 1)
#     if step_x == 1:
#         x.append(x[j] + 1 + 0.05 * np.random.normal())
#     else:
#         x.append(x[j] - 1 + 0.05 * np.random.normal())
#
# y = [0.05 * np.random.normal() for j in range(len(x))]
#
# trace1 = go.Scatter(
#     x=x,
#     y=y,
#     mode='markers',
#     name='Random Walk in 1D',
#     marker=dict(
#         color=[i for i in range(len(x))],
#         size=7,
#         colorscale=[[0, 'rgb(178,10,28)'], [0.50, 'rgb(245,160,105)'],
#                     [0.66, 'rgb(245,195,157)'], [1, 'rgb(220,220,220)']],
#         showscale=True,
#     )
# )
#
# layout = go.Layout(
#     yaxis=dict(
#         range=[-1, 1]
#     )
# )
#
# data = [trace1]
# fig = go.Figure(data=data, layout=layout)
# py.iplot(fig)

# import random
# import math
# def rw3(n,tries):
#     s = 0
#     for m in range(1,tries+1):
#         x = 0
#         y = 0
#         z = 0
#         pi = math.pi
#         for step in range(1,n+1):
#             t = random.uniform(0,2*pi)
#             f = random.uniform(0,2*pi)
#             p = 1
#             x += p*math.sin(f)*math.cos(t)
#             y += p*math.sin(f)*math.sin(t)
#             z += p*math.cos(f)
#         s += (x**2+y**2+z**2)**.5
#     return s/tries
# life = 42
# while life:
#     n = int(input("Please enter the number of steps: "))
#     tries = int(input("How many times should I perform the experiment? "))
#     print()
#     print(rw3(n,tries))
#     print()


# def two_dimensional_random_walk():
#     steps = 0 # Steps counter for understand how many steps that our drunken man take
#     grid_size = 11 # Grid size variable,
#     # Creating Two dimensional array by using lists
#     times = [0] * grid_size
#     for i in range(0,grid_size):
#         times[i] = [0] * grid_size
#     # Initial variables to start in the middle of grid
#     x = 5
#     y = 5
#     # Tuples to get directions and decide where to go
#     moves = [(1,0, "right"), (0,1, "up"), (-1,0, "left"), (0,-1, "down")]
#     # My loop for evaluate the steps
#     while True:
#         dx, dy, position = moves[randint(0,3)] # By using randint I could make decision randomly
#         x += dx
#         y += dy
#         print("He moved", position)
#         try:
#             times[x][y] += 1 # And here is, how many times have he stood on each square
#             steps += 1
#         except IndexError: # The exit of loop
#             break
#     # My print function which answers these questions (How long will it be until he reaeches the end of the sidewalk, and how many times will he have stood on each square)
#     for i in range(0,11):
#         for j in range(0,11):
#             print("He took {steps} steps until he reaches end of the sidewalk.".format(steps = steps),  "He stood on {1}x{2} square at {0} times".format(times[i][j], i+1,j+1) )





from random import randint

def drunken_man():
    steps = 0
    times = [0] * 11
    left_move = 0
    right_move = 0
    i = 0
    while left_move < 5 or right_move < 5:
        value = randint(0,1)
        times[5] = 1
        if value == 1:
            steps += 1
            # print("He moved left")
            left_move += 1
            if right_move > 0:
                right_move -= 1
            if left_move == 1:
                times[4] += 1
            elif left_move == 2:
                times[3] += 1
            elif left_move == 3:
                times[2] += 1
            elif left_move == 4:
                times[1] += 1
            #elif left_move == 5:
                #times[0] += 1
        elif value == 0:
            steps += 1
            # print("He moved right")
            right_move += 1
            if left_move > 0:
                left_move -= 1
            if right_move == 1:
                times[6] += 1
            elif right_move == 2:
                times[7] += 1
            elif right_move == 3:
                times[8] += 1
            elif right_move == 4:
                times[9] += 1
            #elif right_move == 5:
                #times[10] += 1
        times[i] += 1
    for i in range(1,100):
        print("He took {steps} steps until he reaches end of the sidewalk.".format(steps = steps),  "He stood on {1} square at {0} times".format(times[i], i) )

def main():
    drunken_man()

    return 0
if __name__ == '__main__':
    main()