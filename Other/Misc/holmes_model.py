import math
from matplotlib import pyplot as plt
import random

B = random.uniform(0, 0.5)
A = random.uniform(0,0.5)

print(B)
print(A)

p = (1/(B - A))*(math.exp(-0.5/B) - math.exp(-0.5/A))

print(p)

# plt.figure()
# plt.plot(p)
# plt.show()