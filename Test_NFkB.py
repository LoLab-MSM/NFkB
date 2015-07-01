
from pysb import *
Model ()

from pysb.integrate import odesolve
from pylab import linspace, plot, xlabel, ylabel, show
import pylab as pl
#from scipy.integrate import ode


#Declaring the monomers
Monomer('s1')
Monomer('s2')
Monomer('s3')
Monomer('s4')
Monomer('s5')
Monomer('s6')
Monomer('s7')
Monomer('s8')
Monomer('s9')
Monomer('s10')
Monomer('s11')
Monomer('s12')
Monomer('s13')
Monomer('s14')
Monomer('s15')
Monomer('s16')
Monomer('s17')
Monomer('s18')
Monomer('s19')
Monomer('s20')
Monomer('s21')
Monomer('s22')
Monomer('s23')
Monomer('s24')

#Declaring initial conditions
Initial(s1(), Parameter('s1_0'))
Initial(s2(), Parameter('s2_0', 2e5))
Initial(s3(), Parameter('s3_0'))
Initial(s4(), Parameter('s4_0'))
Initial(s5(), Parameter('s5_0'))
Initial(s6(), Parameter('s6_0'))
Initial(s7(), Parameter('s7_0'))
Initial(s8(), Parameter('s8_0', 1))
Initial(s9(), Parameter('s9_0', 10000))
Initial(s10(), Parameter('s10_0', 10))
Initial(s11(), Parameter('s11_0', .14*2e-5))
Initial(s12(), Parameter('s12_0', .06*2e-5))
Initial(s13(), Parameter('s13_0', 10))
Initial(s14(), Parameter('s14_0', 2e5))
Initial(s15(), Parameter('s15_0'))
Initial(s16(), Parameter('s16_0', 2000))
Initial(s17(), Parameter('s17_0'))
Initial(s18(), Parameter('s18_0'))
Initial(s19(), Parameter('s19_0', 1e-1))
Initial(s20(), Parameter('s20_0', 10^5))
Initial(s21(), Parameter('s21_0', 1e-9))
Initial(s22(), Parameter('s22_0', 2000))
Initial(s23(), Parameter('s23_0', 2))
Initial(s24(), Parameter('s24_0', 2))

#Declaring parameters
Parameter('kb', .000012)
Parameter('kf', 1.2e-3)
Parameter('Tdeg', 2e-4)
Parameter('ka',2e-05)
Parameter('Kn', 10^5)
Parameter('ka20', 10^5)
Parameter('ki', .01)
Parameter('B', 0)
Parameter('KNN', 2e-05)
Parameter('k1', 2*6e-10)
Parameter('k3', .002)
Parameter('k2', 10000)
Parameter('a2', 10^-7)
Parameter('tp', .01)
Parameter('i1a', 2e-3)
Parameter('e1a', 5e-3)
Parameter('a3', 5e-7)
Parameter('e2a', 5e-2)
Parameter('a1', 5e-7)
Parameter('i1', .01)
Parameter('c6a', .00002)
Parameter('kv', 5)
Parameter('c4', .5)
Parameter('c5', .0005)
Parameter('c1', .01)
Parameter('G', 0)
Parameter('c3', .00075)
Parameter('c5a', .0001)
Parameter('k4', 10^-3)
Parameter('q1', 4e-7)
Parameter('q2', 10^-6)
Parameter('c1a', 1*.1)

#Declaring expression
Expression('keff', ka*ka20/(ka20+10000))
Expression('Exp1', k3/k2)
Expression('Exp2', a1*kv)
#Expression('k1', 2*k_1)

#Declaring rules
Rule('s1_to_s20', s1() >> s20(), ki)
Rule('s1_and_s2', s1() + s1() + s2() >> s1() + s1() + s3(), k1)
Rule('s3_to_s4', s3() >> s4(), k3)
Rule('s3_and_s9', s3() + s9() >> s9() + s4(), Exp1)
Rule('s4_to_s21', s4() >> s21(), k4)
Rule('s5_deg', s5() >> None, tp)
Rule('s6_to_s7', s6() >> s7(), tp)
Rule('s14_to_s7', s14() >> s7(), c6a)
Rule('s7_to_s8', s7() >> s8(), i1)
Rule('s7_and_s11', s7() + s11() >> s14(), a1)
Rule('s8_and_s12', s8() + s12() >> s15(), Exp2)
Rule('s10_to_s9', s10() >> s10() + s9(), c4)
Rule('s9_deg', s9() >> None, c5)
Rule('s16_and_s20', s16() + s20() >> s16() + s1(), keff)
Rule('s17_to_s10', s17() >> s17() + s10(), c1)
Rule('s10_deg', s10() >> None, c3)
Rule('s3_and_s11', s3() + s11() >> s3() + s5(), a2)
Rule('s13_to_s11', s13() >> s13() + s11(), c4)
Rule('s11_deg', s11() >> None, c5a)
Rule('s11_to_s12', s11() >> s12(), i1a)
Rule('s12_to_s11', s12() >> s11(), e1a)
Rule('s18_to_s13', s18() >> s18() + s13(), c1a)
Rule('s13_deg', s13() >> None, c3)
Rule('s3_and_s14', s3() + s14() >> s3() + s6(), a3)
Rule('s15_to_s14', s15() >> s14(), e2a)
Rule('s16_to_s22', s16() >> s22(), kf)
Rule('s12_and_s17', s12() + s17() >> s12() + s23(), q2)
Rule('s12_and_s18', s12() + s18() >> s12() + s24(), q2)
Rule('s19_deg', s19() >> None, Tdeg)
Rule('s21_to_2', s21() >> s2(), k4)
Rule('s19_and_s22', s19() + s22() >> s19() + s16(), kb)
Rule('s8_and_s23', s8() + s23() >> s8() + s17(), q1)
Rule('s8_and_s24', s8() + s24() >> s8() + s18(), q1)

#Declaring observables
myobs= 's1_obs'
Observable(myobs, s1())
# Observable('s2_obs', s2())
# Observable('s3_obs', s3())
# Observable('s4_obs', s4())
# Observable('s5_obs', s5())
# Observable('s6_obs', s6())
# Observable('s7_obs', s7())
# Observable('s8_obs', s8())
# Observable('s9_obs', s9())
# Observable('s10_obs', s10())
# Observable('s11_obs', s11())
# Observable('s12_obs', s12())
# Observable('s13_obs', s13())
# Observable('s14_obs', s14())
# Observable('s15_obs', s15())
# Observable('s16_obs', s16())
# Observable('s17_obs', s17())
# Observable('s18_obs', s18())
# Observable('s19_obs', s19())
# Observable('s20_obs', s20())
# Observable('s21_obs', s21())
# Observable('s22_obs', s22())
# Observable('s23_obs', s23())
# Observable('s24_obs', s24())
#Observable('s8_and_s15', s8()%s15())

#if __name__ == '__main__':
    # Simulate the model through 40 seconds
# time = linspace(0, 200, 201)
# print "Simulating..."

from pysb.bng import generate_equations

generate_equations(model, verbose=True)

for i,sp in enumerate(model.species):
    print i,":",sp
print
for i,rxn in enumerate(model.reactions):
    print i,":",rxn
print
for i,ode in enumerate(model.odes):
    print i,":",ode

# @TODO fix sumulation errors
#if __name__ == '__main__':
    # Simulate the model through 40 seconds
    time = pl.linspace(0, 1000, 10000)
    print "Simulating..."
    x = odesolve(model, time, integrator='lsoda')
    # Plot the trajectory of LR
    plot(time, x[myobs])
    xlabel('Time (seconds)')
    ylabel('Amount of ' +myobs)
    show()