#!/usr/bin/env python
# -*- coding: utf-8 -*-

from basf2 import *
from modularAnalysis import *
from simulation import add_simulation
from reconstruction import add_reconstruction

myMain = create_path()

# generation of 100 events according to the specified DECAY table
# Y(4S) -> Btag- Bsig+
# Btag- -> D0 pi-; D0 -> K- pi+
# Bsig+ -> pi0 e+ nu_e
#generateEvents(394, 'Upsilon4S.dec'
#		  , myMain)

evtinfosetter = register_module('EventInfoSetter')

evtinfosetter.param('expList', [0])
evtinfosetter.param('runList', [1])
evtinfosetter.param('evtNumList', [100])  # number of events to process

# create geometry
gearbox = register_module('Gearbox')
geometry = register_module('Geometry')

# simulate only inner tracking detectors:
geometry.param('Components', ['MagneticField', 'BeamPipe', 'PXD', 'SVD','CDC'])

evtgeninput = register_module('EvtGenInput')
evtgeninput.param('ParentParticle', "vpho")
evtgeninput.param('boost2LAB', True)
evtgeninput.param('userDECFile',"Upsilon4S.dec")

myMain.add_module(evtinfosetter)

myMain.add_module(gearbox)
myMain.add_module(geometry)
myMain.add_module(evtgeninput)

# simulation
add_simulation(myMain)

# reconstruction
add_reconstruction(myMain)

# do the analysis
loadReconstructedParticles(myMain)

selectParticle('K-', -321, True, myMain)
selectParticle('K+',  321, True, myMain)



# create final state particle lists
selectParticle('K-', -321, [], True, myMain)
selectParticle('K+',  321, [], True, myMain)
#selectParticle('pi+', 211, [], True, myMain)
selectParticle('pi-', -211, [], True, myMain)
#selectParticle('pi0', 111, [], True, myMain)
#selectParticle('e+', -11, [], True, myMain)
#selectParticle('g', 22, [], True, myMain)

# reconstruct D0 -> K+ K- decay (and c.c.)
makeParticle(
    'D0',
    421,
    ['K+', 'K-'],
    1.800,
    1.900,
    True,
    myMain,
    )
    
# reconstruct D(*)- -> D0 pi- decay (no c.c. yet)
makeParticle(
    'D*-',
    413,
    ['D0', 'pi-'],
    1.800,
    2.200,
    True,
    myMain,
    )

# define what should be dumped to ntuple for Btag
toolsD = ['MCTruth', '^D0 -> K+ K-']
toolsD += ['DeltaEMbc', '^D0 -> K+ K-']
toolsD += ['ROEMultiplicities', '^D0 -> K+ K-']
toolsD += ['M','^D0 -> K+ K-']
toolsD += ['p','^D0 -> K+ K-']

# define what should be dumped to ntuple for Btag
toolsDstar = ['MCTruth', '^D*- -> D0 pi-']
toolsDstar += ['DeltaEMbc', '^D*- -> D0 pi-']
toolsDstar += ['ROEMultiplicities', '^D*- -> D0 pi-']
toolsDstar += ['M','^D*- -> D0 pi-']
toolsDstar += ['p','^D*- -> D0 pi-']

ntupleTree('Dzero', 'D0', toolsD, myMain)
ntupleTree('Dstar', 'D*-', toolsDstar, myMain)


process(myMain)
print statistics

