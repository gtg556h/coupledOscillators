import numpy as np
import pdb
import matplotlib.pyplot as plt
import oscillatorLib
import numpy.random

filename = "leaky_nonPhaseDependent_manyOscillators_v2"
title = r"$\dot{c}(c)=0$, $\partial{k_{cell}}/\partial{c}\neq 0$, oscillator pair"

maxTime = 1000 #100
dt = 0.01
t = np.arange(0, maxTime, dt)
x1 = -1
y1 = 0
x2 = 1
y2 =0

nCells = 8

### OSCILLATOR PARAMETERS:
staticRate = 1/1.78  # Computed as 1/((1/.429) - 0.55), where .429 is the frequency of cells from 20160516
leakRate = 0.4 #0.2
randomStd = 0.01
peakEpsilon = 1
epsilon_floor = 0.8
APLength = 0.5#.5  # measured from contractions on 20160516.  This is the visible part only...
contractionDelay = 0.00  # pretty random i suppose
peakForce = 1
peakCouplingRate = 1/(nCells-1)
c_initial0 = 0.7
c_initial1 = 0
sensitivityWinType = 0
sensitivityMean = 0.3
sensitivityStd = 0.1
cellTitle = "cardiacOscillator"

sensitivityWinParam = {'sensitivityWinType' : sensitivityWinType, 'sensitivityMean' : sensitivityMean, 'sensitivityStd' : sensitivityStd}

title = "integrate and fire oscillator"

#strainFunctionParameters = {'funcType' : funcType, 'minStrain' : minStrain, 'maxStrain' : maxStrain, 'phi0' : phi0}

c_initial = numpy.random.rand(nCells)
staticRate = (1 + 0.01*numpy.random.rand(nCells))*staticRate

######################################
# Simulate coupled

cells = []
for i in range(0, nCells):
    cells.append(oscillatorLib.oscillator(x1, y1, dt, t, staticRate[i], leakRate, randomStd, peakEpsilon, epsilon_floor, APLength, contractionDelay, peakForce, c_initial[i], peakCouplingRate, sensitivityWinParam, title))

for j in range(1, t.size):
    forces = 0
    for i in range(0, nCells):
        forces = forces + cells[i].f[j-1]
    for i in range(0, nCells):
        cells[i].epsilon[j] = forces - cells[i].f[j-1]
        cells[i].step(j, cells[i].epsilon[j])

for i in range(0, nCells):
    cells[i].phaseGen()

# dtheta_coupled = oscillatorLib.relativePhase(c0.theta, c1.theta)

######################################
# Simulate uncoupled

ucells = []
peakCouplingRate = 0
for i in range(0, nCells):
    ucells.append(oscillatorLib.oscillator(x1, y1, dt, t, staticRate[i], leakRate, randomStd, peakEpsilon, epsilon_floor, APLength, contractionDelay, peakForce, c_initial[i], peakCouplingRate, sensitivityWinParam, title))

for j in range(1, t.size):
    forces = 0
    for i in range(0, nCells):
        forces = forces + ucells[i].f[j-1]
    for i in range(0, nCells):
        ucells[i].epsilon[j] = forces - ucells[i].f[j-1]
        ucells[i].step(j, ucells[i].epsilon[j])

for i in range(0, nCells):
    ucells[i].phaseGen()

# dtheta_coupled = oscillatorLib.relativePhase(c0.theta, c1.theta)

if 1:
    f0 = 0.5
    f1 = 0.75
    title = r"$\dot{c}(c)\neq0$, $\partial{k_{cell}}/\partial{c}= 0$, many oscillators"
    oscillatorLib.plotManyOscillators(cells, ucells, f0, f1, title, filename)



