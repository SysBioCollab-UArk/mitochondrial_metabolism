from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

# Hplus initial? source?
# oxygen for complexes I-III vs complex IV
# How to reflect mitochondrial membrane potential (electrical & chemical)?
# NAD+/NADH ratio
# FAD/(FAD + NAD(P)H) ratio
# initial FAD amount
# How to make RET conditional?
# ADP & Pi as substrates
# Can the model be used to create an equation relating the redox ratio and ROS?

Model()

# Monomers
Monomer('Complex_I', ['fmn', 'q'])
Monomer('FMN', ['cI', 'e_', 'Hp'], {'e_': ['_0', '_1', '_2'], 'Hp': ['_0', '_1', '_2']})
#Monomer('NADH')
Monomer('NADplus', ['e_','Hp'], {'e_':['_0', '_2'], 'Hp': ['_0', '_1']})
Monomer('Hplus', ['loc'], {'loc': ['mat','ims']})
#Monomer('e_', ['loc'], {'loc': ['mat', 'mem']})
Monomer('O2', ['e_'], {'e_': ['_0', '_1']})
#Monomer('O2_') # superoxide
Monomer('Q',['b','e_','Hp'], {'e_': ['_0', '_1', '_2'], 'Hp': ['_0', '_1', '_2']})
#Monomer('Q_Complex_I', ['cI', 'e_'], {'e_': ['_0', '_1', '_2']})
#Monomer('QH2_Complex_I')
#Monomer('Q_Complex_II', ['cII', 'e_'], {'e_': ['_0', '_1', '_2']})
#Monomer('QH2_Complex_II')
Monomer('Complex_II', ['fad', 'q'])
Monomer('FAD', ['cII', 'e_', 'Hp'], {'e_': ['_0', '_1', '_2'], 'Hp': ['_0', '_1', '_2']})
Monomer('Succinate')
Monomer('Fumarate')
Monomer('Complex_III', ['qo', 'qi'])


# Initial conditions
Parameter('CI_0', 100)
Parameter('NADH_0', 1e4)
Parameter('O2_0', 100)
Parameter('CII_0', 100)
Parameter('SUCC_0', 100)

#Initial(Complex_I(fmn=1, q=2) % FMN(e_='_0', cI=1) % Q_Complex_I(e_='_0', cI=2), CI_0)
Initial(Complex_I(fmn=1, q=2) % FMN(e_='_0', cI=1) % Q(e_='_0', b=2), CI_0)
#Initial(NADH(), NADH_0)
Initial(NADplus(e_='_2', Hp='_1'), NADH_0)
Initial(O2(e_='_0'), O2_0)
#Initial(Complex_II(fad=1, q=2) % FAD(e_='_0', cII=1) % Q_Complex_II(e_='_0', cII=2), CII_0)
Initial(Complex_II(fad=1, q=2) % FAD(e_='_0', cII=1) % Q(e_='_0', b=2), CII_0)
Initial(Succinate(), SUCC_0)


# Rules
Parameter('k_NADH_ox', 1)
Parameter('k_superox_FMN', 1)
Parameter('k_FMN_Q_reduc', 1)
Parameter('k_CI_Q_QH2', 1)
Parameter('k_SUCC_ox', 1)
Parameter('k_FADH2_ox', 1)
Parameter('k_FADH_ox', 1)
Parameter('k_FAD_Q_reduc', 1)
Parameter('k_CII_Q_QH2', 1)

#Rule('NADH_oxidation',
#    NADH() + FMN(e_='_0') >> NADplus() + Hplus(loc='mat') + FMN(e_='_2'), k_NADH_ox)

Rule('NADH_oxidation',
     NADplus(e_='_2', Hp='_1') + FMN(e_='_0', Hp='_0') + Hplus(loc='mat') >> NADplus(e_='_0', Hp='_0') + FMN(e_='_2', Hp='_2'), k_NADH_ox)

Rule('superoxide_FMN_2_1',
     O2(e_='_0') + FMN(e_='_2', Hp='_2') >> O2(e_='_1') + FMN(e_='_1', Hp='_1') + Hplus(loc='ims'), k_superox_FMN)

Rule('superoxide_FMN_1_0',
     O2(e_='_0') + FMN(e_='_1', Hp='_1') >> O2(e_='_1') + FMN(e_='_0', Hp='_0') + Hplus(loc='ims'), k_superox_FMN)

# Rule('FMN_Q_reduction_2_0',
#      FMN(e_='_2') + Q_Complex_I(e_='_0') >> FMN(e_='_1') + Q_Complex_I(e_='_1'), k_FMN_Q_reduc)
#
# Rule('FMN_Q_reduction_1_0',
#      FMN(e_='_1') + Q_Complex_I(e_='_0') >> FMN(e_='_0') + Q_Complex_I(e_='_1'), k_FMN_Q_reduc)
#
# Rule('FMN_Q_reduction_2_1',
#      FMN(e_='_2') + Q_Complex_I(e_='_1') >> FMN(e_='_1') + Q_Complex_I(e_='_2'), k_FMN_Q_reduc)
#
# Rule('FMN_Q_reduction_1_1',
#      FMN(e_='_1') + Q_Complex_I(e_='_1') >> FMN(e_='_0') + Q_Complex_I(e_='_2'), k_FMN_Q_reduc)

Rule('FMN_Q_reduction_2_0',
     FMN(e_='_2', Hp='_2') % Q(e_='_0', Hp='_0') >> FMN(e_='_1', Hp='_1') % Q(e_='_1', Hp='_1'), k_FMN_Q_reduc)

Rule('FMN_Q_reduction_1_0',
     FMN(e_='_1', Hp='_1') % Q(e_='_0', Hp='_0') >> FMN(e_='_0', Hp='_0') % Q(e_='_1', Hp='_1'), k_FMN_Q_reduc)

Rule('FMN_Q_reduction_2_1',
     FMN(e_='_2', Hp='_2') % Q(e_='_1', Hp='_1') >> FMN(e_='_1', Hp='_1') % Q(e_='_2', Hp='_2'), k_FMN_Q_reduc)

Rule('FMN_Q_reduction_1_1',
     FMN(e_='_1', Hp='_1') % Q(e_='_1', Hp='_1') >> FMN(e_='_0', Hp='_0') % Q(e_='_2', Hp='_2'), k_FMN_Q_reduc)

Rule('superoxide_Q_2_1',
    O2(e_='_0') + Complex_I() % Q(e_='_2', Hp='_2') >> O2(e_='_1') + Complex_I() % Q(e_='_1', Hp='_1') + Hplus(loc='ims'), k_superox_Q_CI)

Rule('superoxide_Q_1_0',
    O2(e_='_0') + Complex_I() % Q(e_='_1', Hp='_1') >> O2(e_='_1') + Complex_I() % Q(e_='_0', Hp='_0') + Hplus(loc='ims'), k_superox_Q_CI)



# Rule('CI_Q_reduction_QH2',
#      Q_Complex_I(e_='_2') + Hplus() + Hplus() >> QH2_Complex_I(), k_CI_Q_QH2)

#Rule('CI_Q_reduction_QH2',
#    FMN() % Q(e_='_2', Hp='_0') + Hplus(loc='mat') + Hplus(loc='mat') >> FMN() % Q(e_='_2', Hp='_2'), k_CI_Q_QH2)

Rule('Succinate_oxidation',
     Succinate() + FAD(e_='_0', Hp='_0') >> Fumarate() + FAD(e_='_2', Hp='_2'), k_SUCC_ox)

#Rule('superoxide_FAD_2_1',
#    O2() + FAD(e_='_2') >> O2_() + FAD(e_='_1'), k_FADH2_ox)

#Rule('superoxide_FAD_1_0',
#     O2() + FAD(e_='_1') >> O2_() + FAD(e_='_0'), k_FADH_ox)

Rule('FAD_Q_reduction_2_0',
     FAD(e_='_2', Hp='_2') % Q(e_='_0', Hp='_0') >> FAD(e_='_1', Hp='_1') % Q(e_='_1', Hp='_1'), k_FAD_Q_reduc)

Rule('FAD_Q_reduction_1_0',
     FAD(e_='_1', Hp='_1') % Q(e_='_0', Hp='_0') >> FAD(e_='_0', Hp='_0') % Q(e_='_1', Hp='_1'), k_FAD_Q_reduc)

Rule('FAD_Q_reduction_2_1',
     FAD(e_='_2', Hp='_2') % Q(e_='_1', Hp='_1') >> FAD(e_='_1', Hp='_1') % Q(e_='_2', Hp='_2'), k_FAD_Q_reduc)

Rule('FAD_Q_reduction_1_1',
     FAD(e_='_1', Hp='_1') % Q(e_='_1', Hp='_1') >> FAD(e_='_0', Hp='_0') % Q(e_='_2', Hp='_2'), k_FAD_Q_reduc)

#Rule('CII_Q_reduction_QH2',
#     Q_Complex_II(e_='_2') + Hplus() + Hplus() >> QH2_Complex_II(), k_CII_Q_QH2)

Rule('CI_QH2_CIII',
     Complex_I() % Q(e_='_2', Hp='_2') + Complex_III() >> Complex_I() + Q(e_='_2', Hp='_2') % Complex_III(), k_Q_CItoCIII)

Rule('CII_QH2_CIII',
     Complex_II() % Q(e_='_2', Hp='_2') + Complex_III() >> Complex_II() + Q(e_='_2', Hp='_2') % Complex_III(), k_Q_CIItoCIII)

# Observables
# Complex I
Observable('FMN_e0', FMN(e_='_0'))
Observable('FMN_e1', FMN(e_='_1'))
Observable('FMN_e2', FMN(e_='_2'))
Observable('QI_e0', Q_Complex_I(e_='_0'))
Observable('QI_e1', Q_Complex_I(e_='_1'))
Observable('QI_e2', Q_Complex_I(e_='_2'))
Observable('QIH2_tot', QH2_Complex_I())

# Complex II
Observable('Succinate_to_fumarate', Succinate())
Observable('FAD_e0', FAD(e_='_0'))
Observable('FAD_e1', FAD(e_='_1'))
Observable('FAD_e2', FAD(e_='_2'))
Observable('QII_e0', Q_Complex_II(e_='_0'))
Observable('QII_e1', Q_Complex_II(e_='_1'))
Observable('QII_e2', Q_Complex_II(e_='_2'))
Observable('QIIH2_tot', QH2_Complex_II())

# ROS generation
Observable('O2_tot', O2())
Observable('SuperOx', O2_())

Observable('NADplus_tot', NADplus())
Observable('NADH_tot', NADH())
Expression('redox_ratio', NADplus_tot/NADH_tot)
Observable('FAD_ox', FAD(e_='_0', Hp='_0'))
Expression('redox_ratio_lab', FAD_ox/(FAD_ox + NADH_tot))

# from sympy import Piecewise
# Parameter('thres', 10)
# Observable('X_tot', X())
# Expression('kf', Piecewise( (1,X_tot<thres), (0,True) ) )
# Expression('kr', Piecewise( (1,X_tot>thres), (0,True) ) )
# Rule('A_to_B', A() | B(), kf, kr)

# Simulation commands + plotting

tspan = np.linspace(0, 0.1, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

plt.figure()
for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc=(0.6, 0.05))
plt.tight_layout()

plt.figure()
for obs in [o for o in model.observables if o.name not in ['O2_tot', 'SuperOx', 'Succinate_to_fumarate', 'FAD_e0', 'FAD_e1',
                                                           'FAD_e2', 'QII_e0', 'QII_e1', 'QII_e2', 'QIIH2_tot']]:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc=(0.5,0.55), ncol=2)
plt.tight_layout()

plt.figure()
for obs in [o for o in model.observables if o.name not in ['O2_tot', 'SuperOx', 'FMN_e0', 'FMN_e1', 'FMN_e2', 'QI_e0', 'QI_e1', 'QI_e2', 'QIH2_tot']]:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best')
plt.tight_layout()

plt.figure()
for obs in [o for o in model.observables if o.name in ['O2_tot', 'SuperOx']]:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best')
plt.tight_layout()

plt.show()
