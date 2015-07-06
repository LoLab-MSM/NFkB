# NFkB
                                                              
Function includes system of ODEs describing       
NF-kB regulatory pathway.                         
Present molecules are coded as follows:           
                                                 
y(1)   IKKKa active                               
y(2)   IKKn   neutral                                 
y(3)   IKKa   active                              
y(4)   IKKi   inactive                            
y(5)   phospho-IkBa cytoplasmic                   
y(6)   (phospho-IkBa|NFkB) cytoplasmic              
y(7)   NFkB  cytoplasmic                          
y(8)   NFkBn  nuclear                             
y(9)   A20                                        
y(10)  A20t                                       
y(11)  IkBa                                       
y(12)  IkBan                                      
y(13)  IkBat                                            
y(14)  (IkBa|NFkB) cytoplasmic                     
y(15)  (IkBan|NFkBn) nuclear                        
y(16)  Active receptors                           
y(17)  A20 gene state                                   
y(18)  IkBa gene state                            
y(19)  extracellular TNF                           

#added ODEs for model
y(20) total # IKKK kinase molecules (KN)
y(21) total # IKK kinase molecules (KNN)
y(22) average number active receptors (M)
y(23) A20 alleles (AN)
y(24) IkBa alleles (Ana)

     
#ODEs not used in model
y(20) Reporter gene state                         
y(21) Reporter transcript                         

function dy=ModelD(t,y,Ga,G,M,AN,ANa,ANR)

[NF0,NF1,NF2,M0,M1,M2,k4,ka20,AB,kv,q1,q2,c1,c3,c4,c5,k1,k2,k3,a1,a2,a3,c1a,c5a,c6a,i1,i1a,e1a,e2a,dt,tp,KN,KNN,ka,ki,kb,kf,Tdeg,q1r,q2r,q2rr,c1r,c1rr,c3r]=ParametersD
 ################################################################
 
dy=zeros(21,1)

dy(1)=ka*y(16)*(KN-y(1))* ka20/(ka20+y(9))-ki*y(1)                          #active IKKK kinase 
dy(2)=-y(1)^2*k1*y(2)+k4*(KNN-y(2)-y(3)-y(4))                                 #neutral IKK   
dy(3)=y(1)^2*k1*y(2)-k3*y(3)*(k2+y(9))/k2                                     #free active IKK                                                            
dy(4)=k3*y(3)*(k2+y(9))/k2-k4*y(4)                                            #inactive IKK   
dy(5)=a2*y(3)*y(11)-tp*y(5)                                                   #Phospo-IkBa cytoplasmic 
dy(6)=a3*y(3)*y(14)-tp*y(6)                                                   #cytoplasmic (phospho-IkBa|NF-kB) 
dy(7)=c6a*y(14)-a1*y(7)*y(11)+tp*y(6)-i1*y(7)                                 #free cytoplasmic NFkB
dy(8)=i1*y(7)-a1*kv*y(12)*y(8)                                                #free nuclear NFkB
dy(9)=c4*y(10)-c5*y(9)                                                        #cytoplasmic A20
dy(10)=c1*y(17)-c3*y(10)                                                          #A20 transcript
dy(11)=-a2*y(3)*y(11)-a1*y(11)*y(7)+c4*y(13)-c5a*y(11)-i1a*y(11)+e1a*y(12)    #free cytoplasmic IkBa
dy(12)=-a1*kv*y(12)*y(8)+i1a*y(11)-e1a*y(12)                                  #free nuclear IkBan
dy(13)=c1a*y(18)-c3*y(13)                                                        #IkBa transcript
dy(14)=a1*y(11)*y(7)-c6a*y(14)-a3*y(3)*y(14)+e2a*y(15)                        #cytoplasmic (IkBa|NFkB) complex
dy(15)=a1*kv*y(12)*y(8)-e2a*y(15)                                             #nuclear (IkBa|NFkB) complex
dy(16)=kb*y(19)*(M-y(16))-kf*y(16)                                            #Active receptors
dy(17)=q1*y(8)*(AN-y(17))-q2*y(12)*y(17)                                      #A20 gene state
dy(18)=q1*y(8)*(ANa-y(18))-q2*y(12)*y(18)                                     #IkBa gene state
dy(19)=-Tdeg*y(19)                                                            #extracellular TNF

#added ODES for model 
dy(20)=-ka*y(16)*y(20)*

#ODEs not used in model
dy(20)=q1r*y(8)*(ANR-y(20))-(q2rr+q2r*y(12))*y(20)                      # Reporter gene state 
dy(21)=c1r*y(20)-c3r*y(21)                                              # Reporter transcript  

dy(1)=ka*y(16)*(KN-y(1))* ka20/(ka20*y(9))-ki*y(1)  
Let (KN-y(1)) = y(20)
Let ka*ka20/(ka20+y(9) = Exp1
dy(1)=Exp1*y(16)*y(20)-ki*y(1)                         
#active IKKK kinase 

dy(2)=-y(1)^2*k1*y(2)+k4*(KNN-y(2)-y(3)-y(4))
Let (KNN-y(2)-y(3)-y(4)) = y(21)
dy(2)=-y(1)^2*k1*y(2)+k4*y(21)                                 
#neutral IKK   

dy(3)=y(1)^2*k1*y(2)-k3*y(3)*(k2+y(9))/k2                                    
#free active IKK                                                
                                    
dy(4)=k3*y(3)*(k2+y(9))/k2-k4*y(4)                                            
#inactive IKK   

dy(5)=a2*y(3)*y(11)-tp*y(5)                                                   
#Phospo-IkBa cytoplasmic 

dy(6)=a3*y(3)*y(14)-tp*y(6)                                                   
#cytoplasmic (phospho-IkBa|NF-kB) 

dy(7)=c6a*y(14)-a1*y(7)*y(11)+tp*y(6)-i1*y(7)                                 
#free cytoplasmic NfkB

dy(8)=i1*y(7)-a1*kv*y(12)*y(8)                                                
#free nuclear NfkB

dy(9)=c4*y(10)-c5*y(9)                                                        
#cytoplasmic A20

dy(10)=c1*y(17)-c3*y(10)                                                          
#A20 transcript

dy(11)=-a2*y(3)*y(11)-a1*y(11)*y(7)+c4*y(13)-c5a*y(11)-i1a*y(11)+e1a*y(12)    
#free cytoplasmic IkBa

dy(12)=-a1*kv*y(12)*y(8)+i1a*y(11)-e1a*y(12)                                  
#free nuclear IkBan

dy(13)=c1a*y(18)-c3*y(13)                                                        
#IkBa transcript

dy(14)=a1*y(11)*y(7)-c6a*y(14)-a3*y(3)*y(14)+e2a*y(15)                        
#cytoplasmic (IkBa|NFkB) complex

dy(15)=a1*kv*y(12)*y(8)-e2a*y(15)                                             
#nuclear (IkBa|NFkB) complex
 
dy(16)=kb*y(19)*(M-y(16))-kf*y(16)
Let (M-y(16)) = y(22)
dy(16)=kn*y(19)*y(22)-kf*y(16)                                            
#Active receptors

dy(17)=q1*y(8)*(AN-y(17))-q2*y(12)*y(17)
Let (AN-y(17))=y(23)
dy(17)=q1*y(8)*y(23)-q2*y(12)*y(17)                                     
#A20 gene state

dy(18)=q1*y(8)*(ANa-y(18))-q2*y(12)*y(18)
let (Ana-y(18))=y(24)
dy(18)=q1*y(8)*y(24)-q2*y(12)*y(18)                                     
#IkBa gene state

dy(19)=-Tdeg*y(19)                                                            
#extracellular TNF
 
 
#Omitted ODEs
dy(20)=q1r*y(8)*(ANR-y(20))-(q2rr+q2r*y(12))*y(20)                      # Reporter gene state 
dy(21)=c1r*y(20)-c3r*y(21)                                              # Reporter transcript  

#Added ODEs
dy(20)=-Exp1*y(16)*y(20) + ki*y(1)
#total # IKKK kinase molecules (KN)

dy(21) = k4*y(4)-k4*y(21)
#total # IKK kinase molecules (KNN)

dy(22)=-kb*y(19)*y(22) + kf*y(16)
#average number active receptors (M)

dy(23)=-q1*y(8)*y(23) + q2*y(12)*y(17)
#A20 alleles (AN)

dy(24)=q1*y(8)*y(24) + q2*y(12)*y(18)
#IkBa alleles (Ana)
