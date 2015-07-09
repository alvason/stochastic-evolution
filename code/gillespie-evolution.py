# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Stochastic evolution
# https://github.com/alvason/stochastic-evolution
# 
# ### Evolutionary insights provided by stochastic tools

# <codecell>

'''
author: Alvason Zhenhua Li
date:   07/07/2015
'''
%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
import time
import os
dir_path = '/Users/al/Desktop/GitHub/stochastic-evolution/figure'
file_name = 'gillespie-evolution'

import alva_machinery_probability as alva

AlvaFontSize = 23
AlvaFigSize = (16, 6)
numberingFig = 0

# <codecell>

'''uniform randomness distribution'''
total_event = int(1000)
gInput = np.arange(total_event)
randomSeed = np.random.uniform(0, 1, total_event)

meanP = 0.5
sumP = 0
for i in range(total_event):
    sumP = sumP + (meanP - randomSeed[i])**2
deviationP = (sumP/total_event)**(1.0/2)

totalLevel = int(total_event/10)
category = alva.AlvaLevel(randomSeed, totalLevel, False)
gLevel = category[0]
gLevel_int = gLevel.astype(int)
numberLevel = category[1]
#print ('level =', gLevel)
#print ('level_int =', gLevel_int)

# plotting
figure_name = ''
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
plot1 = figure.add_subplot(1, 2, 1)
plot1.plot(gInput, randomSeed, color = 'gray', marker = 'o', label = 'data')
plot1.plot(gInput, alva.AlvaMinMax(randomSeed), color = 'red', marker = 'o', label = 'minMaxListing')
if total_event < 100:
    plot1.set_xticks(gInput, minor = True) 
    plot1.set_yticks(randomSeed, minor = True)
    plot1.grid(True, which = 'minor')
else:
    plot1.grid(True, which = 'major')
plt.title(r'$ Exponential \ (mean = {:1.3f},\ deviation = {:1.3f}) $'.format(meanP, deviationP), fontsize = AlvaFontSize)
plt.xlabel(r'$ event-input $', fontsize = AlvaFontSize)
plt.ylabel(r'$ output $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))
plt.xticks(fontsize = AlvaFontSize*0.6)
plt.yticks(fontsize = AlvaFontSize*0.6) 

plot2 = figure.add_subplot(1, 2, 2)
#plot2.plot(exp_D, gLevel_int, color = 'blue', marker = 'o', label = 'Exponential Distribution') 
plot2.plot(numberLevel, gLevel, color = 'red', marker = 'o', label = 'category')
if totalLevel < 100:
    plot2.set_xticks(numberLevel, minor = True) 
    plot2.set_yticks(gLevel, minor = True)
    plot2.grid(True, which = 'minor')
else:
    plot2.grid(True, which = 'major')
plt.title(r'$ Exponential \ (events = {:},\ levels = {:}) $'.format(total_event, totalLevel)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ event/level $', fontsize = AlvaFontSize)
plt.ylabel(r'$ level-range $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))
plt.xticks(fontsize = AlvaFontSize*0.6)
plt.yticks(fontsize = AlvaFontSize*0.6) 

figure.tight_layout()
plt.savefig(save_figure, dpi = 30)
plt.show()

# <codecell>

''' starting from one infected '''
# setting parameter
timeUnit = 'day'
if timeUnit == 'day':
    day = 1
    year = 365 
elif timeUnit == 'year':
    year = 1
    day = float(1)/365 
    
# one random number per event
total_event = 300
rate0 = 1.0/10

def Susceptible(total_event, rate0):
    dt = 1.0/10
    gT = np.arange(total_event)*dt
    gS = np.zeros(total_event)
    gS[0] = 100
    randomSeed = np.random.uniform(0, 1, total_event)
    for i in range(total_event - 1):
        if randomSeed[i] < gS[i] * rate0 * dt:
            gS[i + 1] = gS[i] - 1
        else: gS[i + 1] = gS[i]
    return (gT, gS)
 
def Suscept(total_event, rate0):
    gT = np.zeros(total_event)
    gS = np.zeros(total_event)
    gS[0] = 100
    randomSeed = np.random.uniform(0, 1, total_event)
    for i in range(total_event - 1):
        if gS[i] > 1:
            dt = (1.0/(gS[i]*rate0)) * np.log(1.0/randomSeed[i])
            gT[i + 1] = gT[i] + dt 
            gS[i + 1] = gS[i] - 1
        else: 
            dt = 0.0
            gT[i + 1] = gT[i] + dt 
            gS[i + 1] = gS[i]
    return (gT, gS)

# plotting
figure_name = ''
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
plot1 = figure.add_subplot(1, 2, 1)
plot1.plot(Susceptible(total_event, rate0)[0], Susceptible(total_event, rate0)[1], label = 'way 0', drawstyle = 'steps')
plot1.plot(Susceptible(total_event, rate0)[0], Susceptible(total_event, rate0)[1], label = 'way 1', drawstyle = 'steps')
plot1.plot(Susceptible(total_event, rate0)[0], Susceptible(total_event, rate0)[1], label = 'way 2', drawstyle = 'steps')

plot1.grid(True)
plt.title(r'$ (events = {:}, \ rate = {:1.3f}) $'.format(total_event, rate0), fontsize = AlvaFontSize)
plt.xlabel(r'$ time \ ({:})$'.format(timeUnit), fontsize = AlvaFontSize)
plt.ylabel(r'$ Population \ of \ Susceptible $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.3))
plt.xticks(fontsize = AlvaFontSize*0.6)
plt.yticks(fontsize = AlvaFontSize*0.6) 

plot2 = figure.add_subplot(1, 2, 2)
plot2.plot(Suscept(total_event, rate0)[0], Suscept(total_event, rate0)[1], label = 'way 0', drawstyle = 'steps')
plot2.plot(Suscept(total_event, rate0)[0], Suscept(total_event, rate0)[1], label = 'way 1', drawstyle = 'steps')
plot2.plot(Suscept(total_event, rate0)[0], Suscept(total_event, rate0)[1], label = 'way 2', drawstyle = 'steps')

plot2.grid(True)
plt.title(r'$ (events = {:}, \ rate = {:1.3f}) $'.format(total_event, rate0), fontsize = AlvaFontSize)
plt.xlabel(r'$ time \ ({:})$'.format(timeUnit), fontsize = AlvaFontSize)
plt.ylabel(r'$ Population \ of \ Susceptible $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.3))
plt.xticks(fontsize = AlvaFontSize*0.6)
plt.yticks(fontsize = AlvaFontSize*0.6) 

figure.tight_layout()
plt.savefig(save_figure, dpi = 30)
plt.show()

# <codecell>

''' starting from one infected '''
# setting parameter
timeUnit = 'day'
if timeUnit == 'day':
    day = 1
    year = 365 
elif timeUnit == 'year':
    year = 1
    day = float(1)/365 
    
# two random numbers per event
total_event = 300
rate0 = 1.0/10
rate1 = 10.0

def Infected(total_event, rate0, rate1):
    gT = np.zeros(total_event)
    gI = np.zeros(total_event)
    randomSeed_0 = np.random.uniform(0, 1, total_event)
    randomSeed_1 = np.random.uniform(0, 1, total_event)
    for i in range(total_event - 1):
        rr = gI[i]*rate0 + rate1
        dt = (1.0/rr) * np.log(1.0/randomSeed_0[i])
        gT[i + 1] = gT[i] + dt 
        if randomSeed_1[i] < rate1/rr:
            gI[i + 1] = gI[i] + 1
        else: gI[i + 1] = gI[i] - 1
    return (gT, gI)
        
# plotting
figure_name = '-eventSIR'
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
plot1 = figure.add_subplot(1, 2, 1)
plot1.plot(Infected(total_event, rate0, rate1)[0], Infected(total_event, rate0, rate1)[1], label = 'way 0', drawstyle = 'steps')
plot1.plot(Infected(total_event, rate0, rate1)[0], Infected(total_event, rate0, rate1)[1], label = 'way 1', drawstyle = 'steps')
plot1.plot(Infected(total_event, rate0, rate1)[0], Infected(total_event, rate0, rate1)[1], label = 'way 2', drawstyle = 'steps')
plot1.grid(True)
plt.title(r'$ (events = {:}, \ rate = {:1.3f}) $'.format(total_event, rate0), fontsize = AlvaFontSize)
plt.xlabel(r'$ time \ ({:}) $'.format(timeUnit), fontsize = AlvaFontSize)
plt.ylabel(r'$ Population \ of \ Infectious $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.3))
plt.xticks(fontsize = AlvaFontSize*0.6)
plt.yticks(fontsize = AlvaFontSize*0.6) 

figure.tight_layout()
plt.savefig(save_figure, dpi = 300)
plt.show()

# <codecell>


