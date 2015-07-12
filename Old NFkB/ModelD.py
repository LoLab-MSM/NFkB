__author__ = 'geena'
from pysb import *
Model ()
from pysb.macros import catalyze
from pysb.macros import synthesize
from pysb.macros import degrade
dy = (t,y,Ga,G,GR,B)

#TNFa binds to TNFR1 to activate receptor
Monomer('TNFa', ['b'])
Monomer('TNFR1', ['b'])

#full activation of IKKK requires its phosphorylation at two serine residues
# s177, s181. Rate of IKK activation proportional to (IKKKa)^2
Monomer('IKK', ['s177', 's181', 'k'], {'s177': ['u', 'p'], 's181': ['u', 'p']})

#by nuclear lysate (NL) active IKKKa molecules activate IKKK molecules,
#transforming IKKn to IKKa
Monomer('IKKKa', ['s'], {'s': ['u', 'p']})
Monomer('IKKa', ['s'])

#active IKKa phosphorylates IkBa molecules leading to ubiquitination and degradation
Monomer('IkBa', ['s'], {'s': ['u', 'p']})

Monomer('PIkBac', ['s'])
Monomer('PIkBa|NFkBc',['s'])
Monomer('NFkB', ['n', 'c'])
Monomer('NFkBc', ['s'])
Monomer('NFkB_nn', ['s'])
Monomer('A20', ['s'])
Monomer('A20t', ['s'])
Monomer('IkBa', ['s', 'k'], {'s': ['u', 'p']})
Monomer('IkBan', ['s'])
Monomer('IkBat', ['s'])
Monomer('IkBa|NFkBc', ['n', 'c'])
Monomer('IkBan|NFkB_nn', ['s'])
Monomer('IKKi',['s'], {'s': ['u', 'p']})

#TNFa binds to TNFR1 to activate receptor then inactivates and degrades
Parameter('kb', 1.2*10^-5)
Parameter('kf', 1.2*10^-3)
Parameter('TNFa_deg', 2*10^-4)
Rule('TNFa_bind_TNFR1', TNFa(b=None) % TNFR1(b=None), (kb,kf))
degrade('TNFa', TNFa_deg)

#IKKK in active state (IKKKa) (ODE 1)
Parameter('ka',2e-05)
Parameter('Kn', 10^5)
Parameter('Ka20', 10^5)
Parameter('ki', .01)
Parameter('B', 0)
Parameter('A20', 10000)
Parameter('RIF')
Parameter('IKKKa', 0)
Rule('RIP|TRAF2', RIF(b=None), TRAF2(b=None), kf)
Rule('IKKK_bind_B', IKKK(b=None) + B(b=None) >> IKKK(b=1)%B(b=1), Ka)
degrade("IKKK", ki)
#Expression('aIKKK', ka*B*(kn-IKKKa)*ka20/(ka20+A20)-ki*IKKKa)

#IKK in natural state IKKn (ODE 2)
Parameter('IKKn', 2e-05)
Parameter('KNN', 2e-05)
Parameter('kf1', 6*10^-10) #kf= k1
Expression('nIKKn', -(IKKKa)^2 * kcop1*IKKn + k4*(KNN-IKKn-IKKa-IKKi))
catalyze(Ikka, 's', IKK(s = 'u'), 's', IKK(s = 'p'), kf1)


#IKK in the active state IKKa (ODE3)
Parameter('k3', .002)
Parameter('k2', 10000)
Expression('IKKa_rate', (ka)^2)
Rule('')
degrade('IKKa', k2)
Expression('aIKK', (IKKKa)^2*k1*IKKn-k3*IKKa*(k2+A20))

#IKK in inactive state IKKi (ODE4)
Expression('iIKKi', k3*IKKa*(k2+A20)/k2-k4*IKKi)

#Phosphorylation- Ikba(Ikbap) (ODE5)
Parameter('a2', 10^-7)
Parameter('tp', .01)
Parameter('i1a', 2*10^-3)
Parameter('e1a', 5*10^-3)
Parameter('IKKKa', .14*IkBa|NFkBc)
catalyze(IKKa,'s', IkBa(s = 'u'), 's', IkBap(s = 'p'), a2 )
degrade(IkBap(),tp)
Rule('trans_IkBa_IkBap', IkBa(b=None) <> IkBap(b=None), (i1a,e1a))


#Phosphorylated IkBa complex to NFkB (ODE6)
Parameter('a3', 5*10^-7)
Parameter('IkBa|NFkBc', 10^5)
Parameter('e2a', 5*10^-2)
Rule('IkBa_NFkB', IkBa(b=None) + NFkB(b=None), k1)
catalyze(IKKa(), 'b', IkBa_NFkB(y='u'), IkBa_NFkB(y='p'), 'b', a3)
degrade(IkBa_NFkB, tp)
Rule('trans_IkBa_NFkB', IkBa_NFkB(b=None) <> IkBa_NFkB(b=None), e2a)


#Free Cytoplasmic NFkB (ODE7)
Parameter('a1', 5*10^-7)
Parameter('i1', .01)
Parameter('c6a', .00002)
degrade(IkBa(), c6a)
#Expression('fNFkB',c6a*IkBa|NFkB_c-a1*NFkBc*IkBa +tp*IkBa|NFkBp)
Rule('trans_NFkB_Nuc', NFkB(b=None) <> NFkB(b=None), i1)

#Free Nuclear nFkB (ODE8)
Parameter('kv', 5)
Parameter('IkBan',.06*10^-5)
Expression('fNFkBn', -a1*kv*IkBan*NFkBn)

#A20 Protein (ODE 9)
Parameter('c4', .5)
Parameter('A20t', 10)
Parameter('c5', .0005)
synthesize(A20t(), c4)
degrade(A20(), c5)

#A20 Transcript (ODE10)
Parameter('c1', .01)
Parameter('G', 0)
Parameter('c3', .00075)
Rule('Synthesize_NFkB')
Rule('degrad_A20t')

#Free Cytoplasmic IkBa protein (ODE11)
Parameter('a2', 10^-7)
Parameter('IkBa', .14*10^-5)
Parameter('IkBat', 10)
Parameter('c5a', .0001)
Parameter('i1a', .002)
Parameter('e1a', .005)
catalyze(IKKa. IkBa, 'u', 'p', a2)




#Declaring the monomers
Monomer('y1')
Monomer('y2')
Monomer('y3')
Monomer('y4')
Monomer('y5')
Monomer('y6')
Monomer('y7')
Monomer('y8')
Monomer('y9')
Monomer('y10')
Monomer('y11')
Monomer('y12')
Monomer('y13')
Monomer('y14')
Monomer('y15')
Monomer('y16')
Monomer('y17')
Monomer('y18')
Monomer('y19')
Monomer('y20')
Monomer('y21')
Monomer('y22')
Monomer('y23')
Monomer('y24')

#Declaring parameters
Parameter('kb', .000012)
Parameter('kf', 1.2*10^-3)
Parameter('Tdeg', 2*10^-4)
Parameter('ka',2e-05)
Parameter('Kn', 10^5)
Parameter('Ka20', 10^5)
Parameter('ki', .01)
Parameter('B', 0)
Parameter('KNN', 2e-05)
Parameter('k1', 6*10^-10)
Parameter('k3', .002)
Parameter('k2', 10000)
Parameter('a2', 10^-7)
Parameter('tp', .01)
Parameter('i1a', 2*10^-3)
Parameter('e1a', 5*10^-3)
Parameter('a3', 5*10^-7)
Parameter('e2a', 5*10^-2)
Parameter('a1', 5*10^-7)
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
Parameter('q1', 4*10^-7)
Parameter('q2', 10^-6)
Parameter('y14', 2*10^5)
Parameter('y2', 2*10^5)
Parameter('y11', .14*(2*10^-5))
Parameter('y12', .06*(2*10^-5))
Parameter('y13', 10)
Parameter('y10', 10)
Parameter('y9', 10000)
Parameter('y8', 1)

#Declaring expression
Expression('keff', ka*ka20/(ka20+y9))

#Declaring initial conditions
Initial(y14(), y14_0)
Initial(y2(), y2_0)
Initial(y13(), y13_0)
Initial(y11(), y11_0)
Initial(y12(), y12_0)
Initial(y10(), y10_0)
Initial(y8(), y8_0)
Initial(y9(), y9_0)

#Declaring rules
Rule('y1_to_y20', y1() >> y20(), ki)
Rule('y1_and_y2', 2*y1() + y2() >> 2*y1() + y3(), k1)
Rule('y3_to_y4', y3() >> y4(), k3)
Rule('y3_and_y9', y3() + y9() >> y9() + y4(), k3/k2)
Rule('y4_to_y21', y4() >> y21(), k4)
Rule('y5_deg', y5() >> None, tp)
Rule('y6_to_y7', y6() >> y7(), tp)
Rule('create_y14', None >> y14(), c6a)
Rule('y7_to_y8', y7() >> y8(), i1)
Rule('y7_and_y11', y7() + y11() >> y14(), a1)
Rule('y8_and_y12', y8() + y12() >> y15(), a1*kv)
Rule('y10_to_y9', y10() >> y10() + y9(), c4)
Rule('y9_deg', y9() >> None, c5)
Rule('y16_and_y20', y16() + y20() >> y16() + y1(), keff)
Rule('y17_to_y10', y17() >> y17() + y10(), c1)
Rule('y10_deg', y10() >> None, c3)
Rule('y3_and_y11', y3() + y11() >> y3() + y5(), a2)
Rule('y13_to_y11', y13() >> y13() + y11(), c4)
Rule('y11_deg', y11() >> None, c5a)
Rule('y11_to_y12', y11() >> y12(), i1a)
Rule('y12_to_y11', y12() >> y11(), e1a)
Rule('y18_to_y13', y18() >> y18() + y13(), c1a)
Rule('y13_deg', y13() >> None, c3)
Rule('y14_to_y7', y14() >> y7(), c6a)
Rule('y3_and_y14', y3() + y14() >> y3() + y6(), a3)
Rule('y15_to_y14', y15() >> y14(), e2a)
Rule('y16_to_y22', y16() >> y22(), kf)
Rule('y12_and_y17', y12() + y17() >> y12() + y23(), q2)
Rule('y12_and_y18', y12() + y18() >> y12() + y24(), q2)
Rule('y19_deg', y19() >> None, Tdeg)
Rule('y21_to_2', y21() >> y2(), k4)
Rule('y19_and_y22', y19() + y22() >> y19() + y16(), kb)
Rule('y8_and_y23', y8() + y23() >> y8() + y17(), q1)
Rule('y8_and_y24', y8() + y24() >> y8() + y18(), q1)

#Declaring observables
Observable('y1', y1())

#Simulation
if __name__ == '__main__':
    # Simulate the model through 100 seconds
    time = linspace(0, 100, 101)
    print "Simulating..."
    x = odesolve(model, time, verbose=true)
    # Plot the trajectory of LR
    plot(time, x['y1'])
    xlabel('Time (seconds)')
    ylabel('Amount of y1')
    show()