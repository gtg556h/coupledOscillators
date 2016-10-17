import numpy as np
import pdb
import matplotlib.pyplot as plt
import oscillatorLib


filename = "leaky_nonPhaseDependent_oscillatorPair"
title = r"$\dot{c}(c)<0$, $\partial{k_{cell}}/\partial{c}= 0$, oscillator pair"

maxTime = 100
dt = 0.001
t = np.arange(0, maxTime, dt)
x1 = -1
y1 = 0
x2 = 1
y2 =0

### OSCILLATOR PARAMETERS:
staticRate = 1/1.78  # Computed as 1/((1/.429) - 0.55), where .429 is the frequency of cells from 20160516
leakRate = 0.2
randomStd = 0.01
peakEpsilon = 1
epsilon_floor = 0.8
APLength = 0.5  # measured from contractions on 20160516.  This is the visible part only...
contractionDelay = 0.005  # pretty random i suppose
peakForce = 1
peakCouplingRate = 4
c_initial0 = 0.7
c_initial1 = 0
sensitivityWinType = 0
sensitivityMean = 0.6
sensitivityStd = 0.1
cellTitle = "cardiacOscillator"

sensitivityWinParam = {'sensitivityWinType' : sensitivityWinType, 'sensitivityMean' : sensitivityMean, 'sensitivityStd' : sensitivityStd}

#strainFunctionParameters = {'funcType' : funcType, 'minStrain' : minStrain, 'maxStrain' : maxStrain, 'phi0' : phi0}


######################################
# Simulate coupled

staticRate0 = staticRate
staticRate1 = 1.1*staticRate
peakCouplingRate0 = peakCouplingRate
peakCouplingRate1 = peakCouplingRate 

c0 = oscillatorLib.oscillator(x1, y1, dt, t, staticRate0, leakRate, randomStd, peakEpsilon, epsilon_floor, APLength, contractionDelay, peakForce, c_initial0, peakCouplingRate0, sensitivityWinParam, title)

c1 = oscillatorLib.oscillator(x2, y2, dt, t, staticRate1, leakRate, randomStd, peakEpsilon, epsilon_floor, APLength, contractionDelay, peakForce, c_initial1, peakCouplingRate1, sensitivityWinParam, title)

c0.epsilon = np.zeros_like(t)
c1.epsilon = np.zeros_like(t)
for i in range(1,c0.t.size):
    c0.epsilon[i] = c1.f[i-1]
    c1.epsilon[i] = c0.f[i-1]
    c0.step(i, c0.epsilon[i])
    c1.step(i, c1.epsilon[i])

c0.phaseGen()
c1.phaseGen()
dtheta_coupled = oscillatorLib.relativePhase(c0.theta,c1.theta)


######################################
# Simulate uncoupled

staticRate0 = staticRate
staticRate1 = 1.1*staticRate
peakCouplingRate0 = 0
peakCouplingRate1 = 0

u0 = oscillatorLib.oscillator(x1, y1, dt, t, staticRate0, leakRate, randomStd, peakEpsilon, epsilon_floor, APLength, contractionDelay, peakForce, c_initial0, peakCouplingRate0, sensitivityWinParam, title)

u1 = oscillatorLib.oscillator(x2, y2, dt, t, staticRate1, leakRate, randomStd, peakEpsilon, epsilon_floor, APLength, contractionDelay, peakForce, c_initial1, peakCouplingRate1, sensitivityWinParam, title)

u0.epsilon = np.zeros_like(t)
u1.epsilon = np.zeros_like(t)
for i in range(1,u0.t.size):
    u0.epsilon[i] = u1.f[i-1]
    u1.epsilon[i] = u0.f[i-1]
    u0.step(i, u0.epsilon[i])
    u1.step(i, u1.epsilon[i])

u0.phaseGen()
u1.phaseGen()
dtheta_uncoupled = oscillatorLib.relativePhase(u0.theta, u1.theta)


if 1:
    f1 = 0.5
    f2 = 0.75

    oscillatorLib.plotOscillatorPair(c0, c1, u0, u1, f1, f2, title, filename)






