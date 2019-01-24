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

Model ()


Monomer('A20')
Monomer('C1', ['tnf', 'state'], {'state': ['a', 'i']})
Monomer('A20t')
Monomer('NFkB', ['ikb', 'S'], {'S': ['C', 'N']})

# def initial nfkb_conditions():
Parameter('NFkB_0', 0.125) #Nuclear Factor-kappaB
Initial(NFkB(ikb=None, S='N'), NFkB_0)
Observable('NFkBn_free', NFkB(ikb=None, S='N'))

Parameter('A20_mRNA', 2e-6)
Parameter('A20n', 0.4)
Parameter('A20_mRNA_c_deg', 0.035)
Parameter('a1d_c_deg', 0.36)
Parameter('A20_synth', 0.25)
Parameter('A20_deg', 0.0029)
Parameter('hill', 3)

Observable('obs_A20t', A20t())

Expression('A20t_NFkBn', A20n*(NFkBn_free)**(hill))
Expression('A20_synthesis', A20_synth*obs_A20t)

Rule('A20t_synth', None >> A20t(), A20_mRNA)
Rule('A20t_mediated_nfkbn', None >> A20t(), A20t_NFkBn)
Rule('A20t_deg', A20t() >> None, A20_mRNA_c_deg)

Rule('synth_A20', None >> A20(), A20_synthesis)
Rule('deg_A20', A20() >> None, A20_deg)

Parameter('C1_f_A20', 1000.0)
Rule('C1_i_A20', C1(tnf = None, state = 'a') + A20() >> C1(tnf = None, state = 'i') + A20(), C1_f_A20)