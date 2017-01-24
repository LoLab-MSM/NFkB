from pysb import *
import nnrm_modules as g
from pysb.bng import generate_equations, generate_network
import numpy as np


Model('model')

g.TNF_ComplexI_Monomers()
g.TNF_ComplexI_Initials()
g.TNF_ComplexI_Assembly()
g.TNF_ComplexI_K63_Polyubiquitylation()
g.RIP1ub_complex_assembly()
g. ComplexI_PCD_Reactions()
g.ComplexI_Amplication_Reactions()
g.ComplexI_Amplication_Reactions_2_Monomers()
g.ComplexI_Amplication_Reactions_2_Initials()
# g.ComplexI_Amplication_Reactions_2()
g.survival()
g.TNF_ComplexI_Observables()


generate_equations(model)
generate_network(model)

print("printing species length")
print(len(model.species))
