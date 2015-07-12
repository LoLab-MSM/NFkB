__author__ = 'geena'
import math
import numpy as np
import random
from scipy.integrate import ode

from Model import dy



#initial conditions
TNF_dose = 0.1 # TNF dose
Ana = 2 #IkBa alleles
An = 2 #A20 alleles
ANR = 2 #reporter gene alleles

#Reporter genes
IkBa_RG = 0
A20_RG = 0
TNF_RG = 0

#simulation time points
IC_rand = 10*3600 #10h randomization of initial conditions
Equil_wait = 10*3600 #10h equilibrium waiting time
TNF_intro = 50*60 #TNF introduced into the system
TNFw_1 = 50*60 # length of first simulation
TNFe_1 = 100*60 # length of first break
TNFw_2 = 5*60 # length of second simulation
TNFe_2 = 100*60 # length of second break
TNFw_3 = 5*60 # length of third simulation
TNFe_3 = 100*60 # length of third break

tt = 1000 #time for ODE solving

YYY = 0 # matrix of average, all variables y0(i)(t)
NFKB = 0 #total nuclear NF-kB

#Status of IkBa, A20, and TNF reporter genes
GGa = 0
GG = 0
GGT = 0
GGR = 0
Bb = 0 # number of active receptors
MM = 0
NFF = 0
for i=1: N  #beginning main loop

    i #cell number

    #Randomizations of total TNF receptors and NFkB levels
    NFO = 10^5
    NF1 = 1/math.sqrt(2)
    NF2 = -1/4
    MO = 5000
    M1 = math.sqrt(2)
    M2 = -1

    NF = round(NFO*math.exp(NF2 + random.gauss(0,1)*NF1)) # for the NFkB level

    while NF > 10*NFO:
        NF = round(NFO*math.exp(NF2*NF1))
        break

    #NF = NFO #uncomment to remove the extrinsic noise

    M = round(MO*math.exp(M2+random.gauss(0,1)*M1))

    while M > 10*MO:
        M = round(MO*math.exp(M2+random.gauss(0,1)*M1))
        break

    #M = M0 #uncomment to remove the extrinsic noise
    AB = 1 # A20 (on or off)
    y0=np.zeros(1,19)     #initial conditions set to zero and next:

    y0[14]=NF          #NF-kB is given in cytoplasmic complex(IkBa|NFkB) (total NF-kB kept constant), standard = 10^5
    y0[2]=2*10^5       #initial IKKn, total IKK kept constant
    y0[11]=0.14*y0[14] #free cytoplasmic IkBa protein
    y0[12]=0.06*y0[14] #free nuclear IkBa protein
    y0[13]=10          #IkBa mRNA
    y0[10]=10          #10 A20 mRNA
    y0[9]=10000        #10000 A20 protein


    y0[10]=AB*y0[10]
    y0[9]=AB*y0[9]

    Ga=0              #initial status of IkBa promoter
    G=0              #initial status of A20 promoter
    GT=0              #initial status of TNF promoter
    GR=0              #initial status of reporter gene promoter
    B=0              #initial number of active receptors
    yy0=y0            #initial conditions y0(i)
    dt = 10
    ##################################################
    # first step- Randomization of initial conditions

    real_time = 0 #simulated time
    phase = round(random.uniform(0,1)*IC_rand/dt)*dt
    span=[0,dt,tt]                #time for which the solution is derived to find the switching time, tt=1h

    while (real_time<phase):
        [TNF_intro,Y0]=ode(dy,span,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]               #TNF level
        Gax=Ga
        Gx=G
        GRx=GR
        Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Y_act,Yin,M)   #function determining the change of gene status, calls statuschange

        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        real_time=real_time+tc
        break

    if (real_time>phase):
        nn=(real_time-phase)/dt
        yy0=Y0[mk-nn,:]
        a=Gax
        G=Gx
        GR=GRx
        B=Bx                         #status before the last change it occured outside of the time interval
    else:
        print("0")

    ##################################################
    # 0-step waiting for "equilibrium"

    real_time = 0

    while (real_time<Equil_wait):
        [TNF_intro,Y0]=ode(dy,tspan,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]               #TNF level
        Gax=Ga
        Gx=G
        GRx=GR
        Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Y_act,Yin,M)   #function determining the change of gene status, calls statuschange

        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        real_time=real_time+tc
        break

    if (real_time>phase):
        nn=(real_time-phase)/dt
        yy0=Y0[mk-nn,:]
        a=Gax
        G=Gx
        GR=GRx                     #status before the last change it occured outside of the time interval
        B=Bx
    else:
        print("0")


    #####################################################
    # 1 step- still no TNF

    real_time=0
    ga=[Ga]
    g=[G]
    gR=[GR]                  #saves activity of IkBa A20 reporter genes
    bb=[B]                                 #saves number of active receptors
    Y=yy0                                  #variables where single cell run is stored
    T=np.zeros(1,1)

    while (real_time<TNF_intro):
        [TNF_intro,Y0]=math.[dy,tspan,yy0,[],Ga,G,GR,B]
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]                #TNF level
        Gax=Ga; Gx=G; GRx=GR; Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Y_act,Yin,M); #function determining the change of gene status, call statuschange
        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        Y=[Y,Y0[2:mk,:]]           #rows from 2 do mk, all columns
        T=[T,TNF_intro[2:mk]+real_time]
        ga=[ga,Gax*np.ones(mk-1,1)]
        g=[g,Gx*np.ones(mk-1,1)]
        gR=[gR,GRx*np.ones(mk-1,1)]
        bb=[bb,Bx*np.ones(mk-1,1)]
        real_time=real_time+tc
        break


    nn=(real_time-TNF_intro)/dt
    x=a.shape(Y)
    Y=Y[1:(x(1)-nn),:]
    T=T[1:(x(1)-nn)]

    Y0[mk-nn,16]=TNF_dose        #setting TNF ON for the next step

    yy0=Y0[mk-nn,:]
    x1=len(ga)
    ga=ga[1:x1-nn]; g=g[1:x1-nn]; gR=gR[1:x1-nn]; bb=bb[1:x1-nn]
    Ga=Gax; G=Gx; GR=GRx; B=Bx

    ##################################################################
    # 2 step TNF on 1 time


    real_time=0

    while (real_time<TNF_intro):
        [TNF_intro,Y0]=(dy,tspan,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]                #TNF level
        Gax=Ga
        Gx=G
        GTx=GT
        GRx=GR
        Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Y_act,Yin,M) #function determining the change of gene status, call statuschange
        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        Y=[Y,Y0[2:mk,:]]           #rows from 2 do mk, all columns
        T=[T,TNF_intro[2:mk]+real_time]
        ga=[ga,Gax*np.ones(mk-1,1)]
        g=[g,Gx*np.ones(mk-1,1)]
        gR=[gR,GRx*np.ones(mk-1,1)]
        bb=[bb,Bx*np.ones(mk-1,1)]
        real_time=real_time+tc
        break


    nn=(real_time-TNF_intro)/dt
    x=a.shape(Y)
    Y=Y[1:(x(1)-nn),:]
    T=T[1:(x(1)-nn)]

    Y0[mk-nn,16]=TNF_dose        #setting TNF ON for the next step

    yy0=Y0[mk-nn,:]
    x1=len(ga)
    ga=ga[1:x1-nn]
    g=g[1:x1-nn]
    gR=gR[1:x1-nn]
    bb=bb[1:x1-nn]
    Ga=Ga; xG=Gx; GR=GR; xB=Bx

    ####### 3 step TNF washed out 1 time ###########
    #########################################


    real_time=TNF_intro+TNFw_1

    while (real_time<TNF_intro+TNFw_1+TNFe_1):
        [TNF_intro,Y0]=ode23tb(dy,tspan,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]                #TNF level
        Gax=Ga
        Gx=G
        GRx=GR
        Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) #function determining the change of gene status, call statuschange
        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        Y=[Y,Y0[2:mk,:]]           #rows from 2 do mk, all columns
        T=[T,TNF_intro[2:mk]+real_time]
        ga=[ga,Gax*np.ones(mk-1,1)]
        g=[g,Gx*np.ones(mk-1,1)]
        gR=[gR,GRx*np.ones(mk-1,1)]
        bb=[bb,Bx*np.ones(mk-1,1)]
        real_time=real_time+tc
        break


    nn=(real_time-TNF_intro-TNFw_1-TNFe_1)/dt
    x=a.shape(Y)
    Y=Y[1:(x(1)-nn),:]
    T=T[1:(x(1)-nn)]

    Y0[mk-nn,16]=TNF_dose        #setting TNF ON for the next step

    yy0=Y0[mk-nn,:]
    x1=len(ga)
    ga=ga[1:x1-nn]
    g=g[1:x1-nn]
    gR=gR[1:x1-nn]
    bb=bb[1:x1-nn]
    Ga=Gax
    G=Gx
    GR=GRx
    B=Bx



        #########################################
        ####### 4 step TNF on for 2 time ###############
        #########################################


    real_time=TNF_intro+TNFw_1+TNFe_1

    while (real_time<TNF_intro+TNFw_1+TNFe_1+TNFw_2):
        [TNF_intro,Y0]=ode23tb(dy,tspan,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]                #TNF level
        Gax=Ga; Gx=G; GRx=GR; Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) #function determining the change of gene status, call statuschange
        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        Y=[Y,Y0[2:mk,:]]           #rows from 2 do mk, all columns
        T=[T,TNF_intro[2:mk]+real_time]
        ga=[ga,Gax*np.ones(mk-1,1)]
        g=[g,Gx*np.ones(mk-1,1)]
        gR=[gR,GRx*np.ones(mk-1,1)]
        bb=[bb,Bx*np.ones(mk-1,1)]
        real_time=real_time+tc
        break


    nn=(real_time-TNF_intro-TNFw_1-TNFe_1-TNFw_2)/dt
    x=a.shape(Y)
    Y=Y[1:(x(1)-nn),:]
    T=T[1:(x(1)-nn)]

    Y0[mk-nn,16]=0        #setting TNF OFF for the next step

    yy0=Y0[mk-nn,:]
    x1=len(ga)
    ga=ga[1:x1-nn]; g=g[1:x1-nn]; gR=gR[1:x1-nn]; bb=bb[1:x1-nn]
    Ga=Gax; G=Gx; GR=GRx; B=Bx

    #########################################
    ####### 5 step TNF washed out 2 time ###########
    #########################################


    real_time=TNF_intro+TNFw_1+TNFe_1+TNFw_2

    while (real_time<TNF_intro+TNFw_1+TNFe_1+TNFw_2+TNFe_2):
        [TNF_intro,Y0]=ode23tb(dy,tspan,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]                #TNF level
        Gax=Ga; Gx=G; GRx=GR; Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) #function determining the change of gene status, call statuschange
        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        Y=[Y,Y0[2:mk,:]]           #rows from 2 do mk, all columns
        T=[T,TNF_intro[2:mk]+real_time]
        ga=[ga,Gax*np.ones(mk-1,1)]
        g=[g,Gx*np.ones(mk-1,1)]
        gR=[gR,GRx*np.ones(mk-1,1)]
        bb=[bb,Bx*np.ones(mk-1,1)]
        real_time=real_time+tc
        break


    nn=(real_time-TNF_intro-TNFw_1-TNFe_1-TNFw_2-TNFe_2)/dt
    x=a.shape(Y)
    Y=Y[1:(x(1)-nn),:]
    T=T[1:(x(1)-nn)]

    Y0[mk-nn,16]=TNF_dose        #setting TNF ON for the next step

    yy0=Y0[mk-nn,:]
    x1=len(ga)
    ga=ga[1:x1-nn]; g=g[1:x1-nn]; gR=gR[1:x1-nn]; bb=bb[1:x1-nn]
    Ga=Gax; G=Gx; GR=GRx; B=Bx


    #########################################
    ####### 6 step TNF on for the 3 time ###########
    #########################################


    real_time=TNF_intro+TNFw_1+TNFe_1+TNFw_2+TNFe_2

    while (real_time<TNF_intro+TNFw_1+TNFe_1+TNFw_2+TNFe_2+TNFw_3):
        [TNF_intro,Y0]=ode23tb(dy,tspan,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]                #TNF level
        Gax=Ga; Gx=G; GRx=GR; Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) #function determining the change of gene status, call statuschange
        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        Y=[Y,Y0[2:mk,:]]           #rows from 2 do mk, all columns
        T=[T,TNF_intro[2:mk]+real_time]
        ga=[ga,Gax*np.ones(mk-1,1)]
        g=[g,Gx*np.ones(mk-1,1)]
        gR=[gR,GRx*np.ones(mk-1,1)]
        bb=[bb,Bx*np.ones(mk-1,1)]
        real_time=real_time+tc
        break


    nn=(real_time-TNF_intro-TNFw_1-TNFe_1-TNFw_2-TNFe_2-TNFw_3)/dt
    x=a.shape(Y)
    Y=Y[1:(x(1)-nn),:]
    T=T[1:(x(1)-nn)]

    Y0[mk-nn,16]=0        #setting TNF ON for the next step

    yy0=Y0[mk-nn,:]
    x1=len(ga)
    ga=ga[1:x1-nn]; g=g[1:x1-nn]; gR=gR[1:x1-nn]; bb=bb[1:x1-nn]
    Ga=Gax; G=Gx; GR=GRx; B=Bx

        #########################################
        ####### 7 step TNF washed out 3 time ###########
        #########################################


    real_time=TNF_intro+TNFw_1+TNFe_1+TNFw_2+TNFe_2+TNFw_3;


    while (real_time<TNF_intro+TNFw_1+TNFe_1+TNFw_2+TNFe_2+TNFw_3):
        [TNF_intro,Y0]=ode23tb(dy,tspan,yy0,[],Ga,G,GR,B)
        Y_act=Y0[:,8]               #amount of NF-kBn
        Yin=Y0[:,12]               #amount of IkBan
        TR=Y0[:,16]                #TNF level
        Gax=Ga; Gx=G; GRx=GR; Bx=B
        [mk,Ga,G,GR,B]=StatusChange(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) #function determining the change of gene status, call statuschange
        tc=TNF_intro(mk)                  #time when the status changes
        yy0=Y0[mk,:]               #transfer of initial conditions to the next iteration
        Y=[Y,Y0[2:mk,:]]           #rows from 2 do mk, all columns
        T=[T,TNF_intro[2:mk]+real_time]
        ga=[ga,Gax*np.ones(mk-1,1)]
        g=[g,Gx*np.ones(mk-1,1)]
        gR=[gR,GRx*np.ones(mk-1,1)]
        bb=[bb,Bx*np.ones(mk-1,1)]
        real_time=real_time+tc
        break


    nn=(real_time-TNF_intro-TNFw_1-TNFe_1-TNFw_2-TNFe_2-TNFw_3-TNFe_3)/dt;
    x=a.shape(Y)
    Y=Y[1:(x(1)-nn),:]
    T=T[1:(x(1)-nn)]

    Y0[mk-nn,16]=TNF_dose        #setting TNF ON for the next step

    yy0=Y0[mk-nn,:]
    x1=len(ga)
    ga=ga[1:x1-nn]; g=g[1:x1-nn]; gR=gR[1:x1-nn]; bb=bb[1:x1-nn]
    Ga=Gax; G=Gx; GR=GRx; B=Bx

    #########################################################################


    YYY=YYY+Y
    GGa=GGa+ga
    GG=GG+g
    GGR=GGR+gR
    Bb=Bb+bb

    MM(i)=M
    NFF(i)=NF

    ##############################
    ####### DATA FOR PLOTS #######
    ##############################


    XXX(i,:,:)=Y(:,:)
    XB(i,:)=bb(:
    XG(i,:)=g(:)
    XGa(i,:)=ga(:)
    XGR(i,:)=gR(:)  # data for all cells

print(i)                # end of the Main LOOP


    YYY=YYY/N                      #average over population
    GGa=GGa/N
    GG=GG/N
    GGR=GGR/N
    Bb=Bb/N
    T=(T-TNF_intro)/60



