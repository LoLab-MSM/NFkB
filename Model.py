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


