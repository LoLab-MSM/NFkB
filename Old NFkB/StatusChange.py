__author__ = 'geena'

import random
import math
#########################################################################
####### Changes the values of the discrete variables              #######
### First time to the next reaction is determined, then the reaction ####
#########################################################################

def StatusChange(AN,ANa,ANR,TRx,Gax,Gx,GRx,Bx,Yact,Yin,M):
        
    Ga=Gax; G=Gx; GR=GRx ;B=Bx

    #   Yact- amount of NFkBn -Y(:,8)
    #   Yin - amount of IKBa  -Y(:,12)
    #   mk index (time) of gene status change

    ro=(ANa-Gax)*q1*Yact+Gax*(q2*Yin)+(AN-Gx)*q1*Yact+Gx*(q2*Yin)+(ANR-GRx)*q1r*Yact+GRx*(q2r*Yin+q2rr)+(M-Bx)*(kb*TRx)+Bx*kf #total propensity function

    roint=dt*cumtrapz(ro)        # propensity function integrated
    fd=1-math.exp(-roint)             # Distribution of the switching time

    r=rand
    if (fd(len(fd))<r):

        mk=len(fd)
        break

    if (fd(len(fd))>=r):

    a=abs(fd-r)
    mk=find((a-min(a))==0)         # mk = index (time) of next reaction


    ######################################################
    ####### Determining which reaction takes place #######
    ######################################################

    p1a=(ANa-Gax)*q1*Yact(mk)  # risk of NF-kB association to IKBa site at time mk
    p2a=Gax*(q2*Yin(mk))       # risk of  NF-kB dissociation from IkBa at time mk

    p1=(AN-Gx)*q1*Yact(mk)     # risk of NF-kB association to A20 site
    p2=Gx*(q2*Yin(mk))         # risk of  NF-kB dissociation from A20 site

    p1r=(ANR-GRx)*q1r*Yact(mk) # risk of NF-kB association to reporter gene
    p2r=GRx*(q2r*Yin(mk)+q2rr) # risk of  NF-kB dissociation from reporter gene

    p3=(M-Bx)*kb*TRx(mk)  # risk of TNFR1-TNF binding
    p4=Bx*kf             # rist of TNFR1 inactivation

    ss=(p1a+p2a+p1+p2+p1r+p2r+p3+p4)
    p1a=p1a/ss; p2a=p2a/ss
    p1=p1/ss; p2=p2/ss
    p1r=p1r/ss; p2r=p2r/ss
    p3=p3/ss; p4=p4/ss

    rnumber= random.gauss(0,1)

    if (rnumber<p1a):
        Ga=Ga+1  #IKBa activates

    if (rnumber>=p1a)&(rnumber<p1a+p2a):
        Ga=Ga-1  #IKBa inactivates

    if (rnumber>=p1a+p2a)&(rnumber<p1a+p2a+p1):
        G=G+1     #A20 activates

    if (rnumber>=p1a+p2a+p1)&(rnumber<p1a+p2a+p1+p2):
        G=G-1     #A20 inactivates

    if (rnumber>=p1a+p2a+p1+p2)&(rnumber<p1a+p2a+p1+p2+p1r):
        GR=GR+1   #reporter gene activates

    if (rnumber>=p1a+p2a+p1+p2+p1r)&(rnumber<p1a+p2a+p1+p2+p1r+p2r):
        GR=GR-1   #reporter gene inactivates

    if (rnumber>=p1a+p2a+p1+p2+p1r+p2r)&(rnumber<p1a+p2a+p1+p2+p1r+p2r+p3):
        B=B+1    #receptor activation

    if (rnumber>=p1a+p2a+p1+p2+p1r+p2r+p3)&(rnumber<p1a+p2a+p1+p2+p1r+p2r+p3+p4):
        B=B-1    #receptor deactivation
