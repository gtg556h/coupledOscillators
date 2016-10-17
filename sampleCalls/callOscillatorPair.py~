import numpy as np
import pdb
import matplotlib.pyplot as plt
import oscillatorLib




maxTime = 50
dt = 0.001
t = np.arange(0, maxTime, dt)
x = 0
y = 0

### OSCILLATOR PARAMETERS:
staticRate = 1/1.78  # Computed as 1/((1/.429) - 0.55), where .429 is the frequency of cells from 20160516
leakRate = 0 #0.2
randomStd = 0.1
peakEpsilon = 1
epsilon_floor = 0.8
APLength = 0.5  # measured from contractions on 20160516.  This is the visible part only...
contractionDelay = 0.005  # pretty random i suppose
peakForce = 1
peakCouplingRate = 10
c0 = 0
sensitivityWinType = 1
sensitivityMean = 0.6
sensitivityStd = 0.1
cellTitle = "cardiacOscillator"

sensitivityWinParam = {'sensitivityWinType' : sensitivityWinType, 'sensitivityMean' : sensitivityMean, 'sensitivityStd' : sensitivityStd}

title = "integrate and fire oscillator"

#strainFunctionParameters = {'funcType' : funcType, 'minStrain' : minStrain, 'maxStrain' : maxStrain, 'phi0' : phi0}


o1 = oscillatorLib.oscillator(x, y, dt, t, staticRate, leakRate, randomStd, peakEpsilon, epsilon_floor, APLength, contractionDelay, peakForce, c0, peakCouplingRate, sensitivityWinParam, title)


epsilon = np.zeros_like(t)
epsilon = np.sin(3*t)+1
for i in range(0,o1.t.size):
    o1.step(i, epsilon[i])


if 1:
    plt.plot(o1.t, o1.c, o1.t, epsilon)
    plt.show()


