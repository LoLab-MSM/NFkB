from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
import nfkb_modules

Model ()

#Declaration of monomers
nfkb_modules.ikba_and_mRNA_monomers()
nfkb_modules.ikbb_and_mRNA_monomers()
nfkb_modules.ikbe_and_mRNA_monomers()
nfkb_modules.ikbd_and_mRNA_monomers()
print('hello1')
nfkb_modules.nfkb_and_ikk_monomers()
print('hello2')
# nfkb_modules.a20_mrna_to_a20()
# nfkb_modules.ligand_to_receptor_monomers()
# nfkb_modules.complex_monomers()
nfkb_modules.kinase_monomers()
# nfkb_modules.initial_conditions()
nfkb_modules.observables()

nfkb_modules.ikb_mrna_to_ikb()
#IkB(a,b,e) association and dissociation from IKK2 and IkBd association and dissociation from IKK2
nfkb_modules.ikb_assoc_diss_nfkb()

# #IkB and NFkB cellular localization reactions
nfkb_modules.ikb_nfkb_localization()

# IkB Protein Degradation Reactions
nfkb_modules.ikb_deg_reactions()

#IKK-mediated IkB degradation reactions
nfkb_modules.ikb_ikk_mediated_deg()

#A20 mRNA and Protein Synthesis and Degradation Reactions
nfkb_modules.a20_reactions()

#IKK Activation Module
#TNF-Independent Complex 1 Activity Reactions
nfkb_modules.tnf_independent_to_c1()

#TNF-Dependent Complex 1 Activity Reactions
nfkb_modules.tnf_dependent_to_c1()

#IKKK (TAB1/2-TAK1 complex) Activity Reactions
nfkb_modules.ikkk_to_ikk_complex()

#Dictionary to substitute in species names to match matlab files
nfkb_modules.species_dict()

#phase 2 of model using Solver
tspan = np.linspace(0, 300, 301)
x = odesolve(model,tspan,verbose=True)

#updating the species names in the odes
for  j,ode in enumerate(model.odes):
    for i in range(len(model.species)):
       ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
    print j,":",ode





# #equilibrium phase 1 of model using odesolve
# # time_equil = np.linspace(-8000, 0, 8001)
# # equil = odesolve(model, time_equil, verbose=True)
#
# # plt.figure(1)
# # plt.plot(time_equil/60., equil['NFkBn_free']*1000, label= NFkBn_free.name)
# # plt.xlabel("Time (in hours)", fontsize=16)
# # plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = 0, ymax =100)
# # plt.legend(loc=0)
#
# # last_conc = [equil[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]]
# # print(last_conc)#setting previous conc species equal to none
#


#
# set_kappa_path('/Users/geenaildefonso/Projects/KaSim')
# x = contact_map(model)
# x.draw('contact_map.pdf', format='pdf', prog='dot')
#
# x = run(model)
# g = pyg.AGraph(x)
# g.draw('render_reactions.pdf', format='pdf', prog='dot')

# #phase 2 of model using Solver
tspan = np.linspace(0, 720, 721)
phase = odesolve(model,tspan,verbose=True)
last_conc = [phase[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]]
print(last_conc)

#updating the species names in the odes
# for  j,ode in enumerate(model.odes):
#     for i in range(len(model.species)):
#        ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
#     print j,":",ode

with PdfPages('multipage_pdf.pdf') as pdf:
    plt.figure()
    plt.plot(tspan, phase['IKKK_off_obs'], label = IKKK_off_obs.name)
    # plt.plot(tspan, pandas_df['IKKK_off'][0:721], label = 'IKKK_off')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0, ymax = 0.5)
    plt.legend(loc=0)
    plt.title('IKKK_off')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # plt.show()

    plt.figure()
    plt.plot(tspan, phase['IKKK_obs'], label = IKKK_obs.name)
    # plt.plot(tspan, pandas_df['IKKK'][0:721], label = 'IKKK')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IKKK')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IKK_i_obs'], label = IKK_i_obs.name)
    # plt.plot(tspan, pandas_df['IKK_i'][0:721], label = 'IKK_i_off')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IKK_i')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IKK_obs'], label = IKK_obs.name)
    # plt.plot(tspan, pandas_df['IKK'][0:721], label = 'IKK')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IKK')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    #
    plt.figure()
    plt.plot(tspan, phase['IKK_off_obs'], label = IKK_off_obs.name)
    # plt.plot(tspan, pandas_df['IKK_off'][0:721], label = 'IKK_off')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IKK_off')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['TNFR_obs'], label = TNFR_obs.name)
    # plt.plot(tspan, pandas_df['tnfr'][0:721], label = 'TNFR')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0, ymax = 0.00000008)
    plt.legend(loc=0)
    plt.title('tnfr')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['NFkB_obs'], label = NFkB_obs.name)
    # plt.plot(tspan, pandas_df['NFkB'][0:721], label = 'NFkB')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0, ymax = 0.00000008)
    plt.legend(loc=0)
    plt.title('tnfr')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['NFkBn_obs'], label = NFkBn_obs.name)
    # plt.plot(tspan, pandas_df['NFkBn'][0:721], label = 'NFkBn')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0, ymax = 0.00000008)
    plt.legend(loc=0)
    plt.title('tnfr')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    #
    # # plt.figure()
    # # plt.plot(tspan, phase['TNFR_obs'], label = TNFR_obs.name)
    # # plt.plot(tspan, pandas_df['TNFR values'][0:721], label = 'TNFR')
    # # plt.xlabel("Time (in min)", fontsize=16)
    # # plt.ylabel("Concentration", fontsize=16)
    # # plt.ylim(ymin = 0, ymax = 0.00000008)
    # # plt.legend(loc=0)
    # # plt.show()
    # #
    # plt.figure()
    # plt.plot(tspan, phase['tnfrtnf_obs'], label = TNFRtnf_obs.name)
    # plt.plot(tspan, pandas_df['tnfrtnf'][2:723], label = 'TNFRtnf')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0)
    # plt.legend(loc=0)
    #
    # plt.figure()
    # plt.plot(tspan, phase['TNFRtnf_obs'], label = TNFRtnf_obs.name)
    # plt.plot(tspan, pandas_df['TNFRtnf values'][2:723], label = 'TNFRtnf')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0)
    # plt.legend(loc=0)


    plt.figure()
    plt.plot(tspan, phase['tnfrm_obs'], label = tnfrm_obs.name)
    # plt.plot(tspan, pandas_df['tnfrm'][0:721], label = 'tnfrm')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('tnfrm')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
# plt.show()
#       print(model.reactions)
#
# #
# # # for i,ode in enumerate(model.odes):
# # #     print i,":", ode
# # #
# for i,sp in enumerate(model.species):
#     print i,":", sp
#
#
# for i,pr in enumerate(model.parameters):
#     print i,":", pr
#
# # #
#
#

    plt.figure()
    plt.plot(tspan, phase['IkBat_obs'], label = IkBat_obs.name)
    # plt.plot(tspan, pandas_df['IkBat'][0:721], label = 'IkBat')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBat')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBbt_obs'], label = IkBbt_obs.name)
    # plt.plot(tspan, pandas_df['IkBbt'][0:721], label = 'IkBbt')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBbt')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBet_obs'], label = IkBet_obs.name)
    # plt.plot(tspan, pandas_df['IkBet'][0:721], label = 'IkBet')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBet')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBdt_obs'], label = IkBdt_obs.name)
    # plt.plot(tspan, pandas_df['IkBdt'][0:721], label = 'IkBdt')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBdt')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBa_obs'], label = IkBa_obs.name)
    # plt.plot(tspan, pandas_df['IkBa'][0:721], label = 'IkBa')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBa')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBb_obs'], label = IkBb_obs.name)
    # plt.plot(tspan, pandas_df['IkBb'][0:721], label = 'IkBa')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBb')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBe_obs'], label = IkBe_obs.name)
    # plt.plot(tspan, pandas_df['IkBe'][0:721], label = 'IkBe')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBe')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBd_obs'], label = IkBd_obs.name)
    # plt.plot(tspan, pandas_df['IkBd'][0:721], label = 'IkBd')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBd')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    #
    # #

    #
    # #
    # plt.figure()
    # plt.plot(tspan, phase['IKKK_off_obs'], label = IKKK_off_obs.name)
    # # plt.plot(tspan, pandas_df['IKKK_off values'][0:721], label = 'IKKK_off')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0, ymax = 0.5)
    # plt.legend(loc=0)
    # # plt.show()
    #
    #
    # plt.figure()
    # plt.plot(tspan, phase['IKK_i_obs'], label = IKK_i_obs.name)
    # # plt.plot(tspan, pandas_df['IKK_i values'][0:721], label = 'IKK_i_off')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0)
    # plt.legend(loc=0)
    #
    # plt.figure()
    # plt.plot(tspan, phase['IKK_obs'], label = IKK_obs.name)
    # # plt.plot(tspan, pandas_df['IKK values'][0:721], label = 'IKK')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0)
    # plt.legend(loc=0)
    #


    plt.figure()
    plt.plot(tspan, phase['IkBaNFkB_obs'], label = IkBaNFkB_obs.name)
    # plt.plot(tspan, pandas_df['IkBaNFkB'][0:721], label = 'IkBaNFkB')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0, ymax = 0.5)
    plt.legend(loc=0)
    plt.title('IkBaNFkB')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # plt.show()
    #
    plt.figure()
    plt.plot(tspan, phase['IkBbNFkB_obs'], label = IkBbNFkB_obs.name)
    # plt.plot(tspan, pandas_df['IkBbNFkB'][0:721], label = 'IkBbNFkB')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBbNFkB')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBeNFkB_obs'], label = IkBeNFkB_obs.name)
    # plt.plot(tspan, pandas_df['IkBeNFkB'][0:721], label = 'IkBebNFkB')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBeNFkB')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBdNFkB_obs'], label = IkBdNFkB_obs.name)
    # plt.plot(tspan, pandas_df['IkBdNFkB'][0:721], label = 'IkBdNFkB')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBdNFkB')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


    plt.figure()
    plt.plot(tspan, phase['IkBaNFkBn_obs'], label = IkBaNFkBn_obs.name)
    # plt.plot(tspan, pandas_df['IkBaNFkBn'][0:721], label = 'IkBaNFkBn')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0, ymax = 0.5)
    plt.legend(loc=0)
    plt.title('IkBaNFkBn')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    # plt.show()
    #
    plt.figure()
    plt.plot(tspan, phase['IkBbNFkBn_obs'], label = IkBbNFkBn_obs.name)
    # plt.plot(tspan, pandas_df['IkBbNFkBn'][0:721], label = 'IkBbNFkBn')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBbNFkBn')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBeNFkBn_obs'], label = IkBeNFkB_obs.name)
    # plt.plot(tspan, pandas_df['IkBeNFkBn'][0:721], label = 'IkBebNFkBn')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBeNFkBn')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['IkBdNFkBn_obs'], label = IkBdNFkBn_obs.name)
    # plt.plot(tspan, pandas_df['IkBdNFkBn'][0:721], label = 'IkBdNFkBn')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('IkBdNFkBn')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    # plt.figure()
    # plt.plot(tspan, phase['IkBa_obs'], label = IkBa_obs.name)
    # plt.plot(tspan, pandas_df['IkBa values'][0:721], label = 'IkBa')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0, ymax = 0.5)
    # plt.legend(loc=0)
    # # plt.show()
    # #
    # plt.figure()
    # plt.plot(tspan, phase['IkBb_obs'], label = IkBb_obs.name)
    # plt.plot(tspan, pandas_df['IkBb values'][0:721], label = 'IkBb')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0)
    # plt.legend(loc=0)
    #
    # plt.figure()
    # plt.plot(tspan, phase['IkBe_obs'], label = IkBe_obs.name)
    # plt.plot(tspan, pandas_df['IkBe values'][0:721], label = 'IkBe')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0)
    # plt.legend(loc=0)
    #
    # plt.figure()
    # plt.plot(tspan, phase['IkBd_obs'], label = IkBd_obs.name)
    # plt.plot(tspan, pandas_df['IkBd values'][0:721], label = 'IkBd')
    # plt.xlabel("Time (in min)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.ylim(ymin = 0)
    # plt.legend(loc=0)



    plt.figure()
    plt.plot(tspan, phase['TNF_obs'], label = TNF_obs.name)
    # plt.plot(tspan, pandas_df['TNF'][0:721], label = 'TNF')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('TNF')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['A20_obs'], label = A20_obs.name)
    # plt.plot(tspan, pandas_df['A20'][0:721], label = 'A20')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('A20')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['A20t_obs'], label = A20t_obs.name)
    # plt.plot(tspan, pandas_df['A20t'][0:721], label = 'A20t')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('A20t')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['C1_obs'], label = C1_obs.name)
    # plt.plot(tspan, pandas_df['C1'][0:721], label = 'C1')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('C1')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['C1_off_obs'], label = C1_off_obs.name)
    # plt.plot(tspan, pandas_df['C1_off'][0:721], label = 'C1_off')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('C1_off')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['C1tnf_obs'], label = C1tnf_obs.name)
    # plt.plot(tspan, pandas_df['C1tnf'][0:721], label = 'C1tnf')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('C1tnf')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    plt.plot(tspan, phase['C1tnf_off_obs'], label = C1tnf_off_obs.name)
    # plt.plot(tspan, pandas_df['C1tnf_off'][0:721], label = 'C1tnf_off')
    plt.xlabel("Time (in min)", fontsize=16)
    plt.ylabel("Concentration", fontsize=16)
    plt.ylim(ymin = 0)
    plt.legend(loc=0)
    plt.title('C1tnf_off')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

# plt.show()
# plt.show()
# plt.figure()
# plt.plot(tspan, x['IkBbt_obs'], label = IkBbt_obs.name)
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan, x['IkBet_obs'], label = IkBet_obs.name)
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0)
# plt.legend(loc=0)
# # plt.show()
#
# plt.figure()
# plt.plot(tspan, x['IkBdt_obs'], label = IkBdt_obs.name)
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0)
# plt.legend(loc=0)
# plt.show()
#
# plt.figure()
# plt.plot(tspan, x['NFkBn_free'], label = NFkBn_free.name)
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0, ymax =0.125)
# plt.legend(loc=0)
#
#
# plt.figure()
# plt.plot(tspan, x['NFkBn_obs'], label = NFkBn_obs.name)
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan, x['NFkBn_bound'], label = NFkBn_bound.name)
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0, ymax = .0012)
# plt.legend(loc=0)
# #
# plt.figure()
# plt.plot(tspan, x['IkBat_obs'], label = IkBat_obs.name)
# plt.xlabel("Time (in min)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0, ymax =.016)
# plt.legend(loc=0)
# plt.show()
# plt.figure(1)
# plt.plot(tspan/60, x['NFkBc_free'], label = NFkBc_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# plt.show()


# plt.figure(1)
# plt.plot(tspan/60, x['NFkBn_free'], label = NFkBn_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure(1)
# plt.plot(tspan/60, x['NFkBn_free'], label = NFkBn_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure(1)
# plt.plot(tspan/60, x['NFkBn_bound'], label = NFkBn_bound.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan/60, x['IKK_obs'], label = IKK.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.figure()
# plt.plot(tspan/60, x['A20_obs'], label = A20.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# #
# # ikk2_vals = []
# nsims = len(tspan) - 1 #simulation time is 720
# # yobs = np.empty((nsims, len(model.observables)))
# rows = len(tspan) # nSims
# cols = len(model.observables) #length of observables
# yobs = np.empty((rows, cols)) #creating empty array of R and C
#
# assert model.observables[1].name == 'IKK2_obs'
# assert str(model.species[1]) == "IKK2(ikb=None, S='C')"
#
# # print(model.parameters['IKK2_0'].value)
# #using imported IKK2 values during phase 2 simulation
# for i in range(nsims):
#     print(i) #printing each simulation iteration
#     last_conc[1] = ikk2[i]
#     solver.tspan = [tspan[i], tspan[i+1]] #simulating in 2 step intervals
#     solver.run(y0 = last_conc) #setting initial concentrations to previous simulation value (equil phase 1)
#     # print(last_conc)
#     if i == 0:
#         yobs[0, :] = solver.yobs_view[0]
#     yobs[i + 1 , :] = solver.yobs_view[1] #taking each simulation index and all obs and keeping the last obs in last time point of simulation
#     last_conc = solver.y[1,:] # running solver taking first row and all simulation species


#plotting NFkBn free
# plt.figure(1)
# plt.plot(tspan/60., yobs[:,0]*1000, label= NFkBn_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0, ymax =100)
# plt.legend(loc=0)

# #plotting NFkBc free
# plt.figure(2)
# plt.plot(tspan/60., yobs[:,2], label= NFkBc_free.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)
#
# #plotting IkBa unbound cytosol
# plt.figure(2)
# plt.plot(tspan/60., yobs[:,10], label= IkBa_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)

plt.figure()
plt.plot(tspan/60, x['A20_obs'], label = A20_obs.name)
plt.xlabel("Time (in hours)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
# plt.ylim(ymin = 0, ymax =0.04)
plt.legend(loc=0)
plt.show()





