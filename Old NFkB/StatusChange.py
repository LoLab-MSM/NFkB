kb,1.2*10e-5    #default 1.2*10e-5 - receptor activation rate
kf,1.2*10e-3    #default 1.2*10e-3 - receptor inactivation rate
q1,4*10e-7      #default 4*10e-7 - NF-kB ataching at A20 and IkBa site                                     
q2,10e-6        #default 10e-6 - IkBa inducible detaching from A20 and IkBa site
Tdeg,7.7*10e-4  #TNF loss
KN,10e5         #default 10e5 - total number of IKKK kinase molecules, Assumption
KNN,2*10e5      #default 2*10e5 - total number of IKK kinase molecules, Assumption
ka,10e-5        #default 10e-5 - IKKK kinase activation rate (at most 1/s), Assumption                        
ki,0.01         #default 0.01 - IKKK kinase inactivation rate, Assumption
c1,0.1          #inducible A20 mRNA synthesis
c3,0.00075      #default 0.00075 - A20 and IkBa mRNA degradation rate
c4,0.5          #default 0.5 - A20 and IkBa translation rate, FIT
c5,0.0005       #default 0.0005 - A20 degradation rate, FIT
ka20,10e5       #default 10000 - A20 TNFR1 block, FIT                                                      
k2,10000        #default 10000 - IKKa inactivation caused by A20, FIT                                         
k1,6*10e-10     #default 6*10e-10 IKKn activation caused by active IKKK, FIT
k3,0.002        #default 0.002 - IKKa inactivation, FIT                                                          
k4,0.001        #default 0.001 - IKKii transformation, FIT
c1a,0.1         #inducible IkBa mRNA synthesis 
a1,5*10e-7      #default 5*10e-7 - IkBA*NFkB association, Assumption
a2,10e-7        #default 10e-7 - IkBa phosphoryation due to action of IKKa, FIT
a3,5*10e-7      #default 5*10e-7 - (IkBa|NFkb) phosphorylation due to action of IKKa, (at most 0.01/s) FIT
tp,0.01         #default 0.01 - degradation of phospho-IkBa and phospho-IkBa complexed to NF-kB, FIT  
c5a,0.0001      #default 0.0001 - IkBa degradation rate
c6a,0.00002     #default 0.00002 - spontaneous (IkBa|NFkB) degradation of IkBa  complexed to NF-kB
i1,0.01         #default 0.01 - NFkB nuclear import, FIT 
e2a,0.05        #default 0.05 - (IkBa|NFkB) nuclear export, FIT 
i1a,0.002       #default 0.002 - IkBa nuclear import, FIT
e1a,0.005       #default 0.005 - IkBa nuclear export, FIT 