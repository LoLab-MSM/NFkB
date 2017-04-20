from pysb import *
from pysb import *
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import scipy.interpolate
#from pysb.simulator import ScipyOdeSimulator
import numpy as np
# import holmes_model_two as m2
# import holmes_three as m3
# import holmes_four as m4
from pysb.util import alias_model_components
from matplotlib import *


Model ()


Monomer('L', ['r'])
Monomer('R', ['l'])
Monomer('ySEC', ['lr'])
Monomer('ICD', ['csl'])
Monomer('CSL', ['icd'])


Observable('ICD_obs', ICD(csl = None))
Observable('ICDCSL_obs', ICD(csl = 1)%CSL(icd = 1))


Parameter('kflr', 1e-5)
Parameter('krlr', 1e-5)
Parameter('k1f_y', 1e-5)
Parameter('k1r_y', 1e-5)
Parameter('k2_y', 1e-5)
Parameter('kf_icd', 1e-5)
Parameter('kr_icd', 1e-5)
Parameter('L_0', 10000)
Parameter('R_0', 5000)
Parameter('ySEC_0', 500)
Parameter('CSL_0', 1000)


Initial(L(r = None), L_0)
Initial(R(l = None), R_0)
Initial(ySEC(lr = None), ySEC_0)
Initial(CSL(icd = None),CSL_0)

