
# coding: utf-8

# # Stochastic infectious pulse
# https://github.com/alvason/stochastic-infectious-pulse
# 
# ### Stochastic version for evolutionary insights

# In[7]:

'''
author: Alvason Zhenhua Li
date:   07/07/2015
'''
get_ipython().magic(u'matplotlib inline')

import numpy as np
import matplotlib.pyplot as plt
import time
import os
dir_path = '/Users/al/Desktop/GitHub/stochastic-infectious-pulse/figure'
file_name = 'gillespie-evolution'

import alva_machinery_probability as alva

AlvaFontSize = 23
AlvaFigSize = (16, 6)
numberingFig = 0


# In[8]:

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


# In[9]:

# algorithm for stochastic evolution
# plotting
figure_name = '-equation'
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)

totalStep = 4
numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize=(9, 6))
plt.axis('off')
plt.title(r'$ Stochastic-evolution $',fontsize = AlvaFontSize)
plt.text(0, 3.0/totalStep, r'$ 1. \ event = random.uniform(0, 1, totalEvent) $', fontsize = 1.2*AlvaFontSize)
plt.text(0, 2.0/totalStep, r'$ 2. \ event = exp[-S(t) \beta \Delta t] $'
         , fontsize = 1.2*AlvaFontSize)
plt.text(0, 1.0/totalStep, r'$ 3. \ \Longrightarrow \Delta t = \frac{1}{S(t)\beta} (-ln[event]) = \frac{1}{S(t)\beta} ln[\frac{1}{event}] $'
         , fontsize = 1.2*AlvaFontSize)
plt.text(0, 0.0/totalStep, r'$ 4. \ \Longrightarrow S(t + \Delta t) = S(t) - 1 $'
         , fontsize = 1.2*AlvaFontSize)
plt.savefig(save_figure, dpi = 100)
plt.show()


# In[10]:

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
rate_0 = 3.0/365
initial_value = 99

def Susceptible(total_event, rate_0, initial_value):
    dt = 1.0/10
    gT = np.arange(total_event)*dt
    gS = np.zeros(total_event)
    gS[0] = initial_value
    randomSeed = np.random.uniform(0, 1, total_event)
    for i in range(total_event - 1):
        if randomSeed[i] < gS[i] * rate_0 * dt:
            gS[i + 1] = gS[i] - 1
        else: gS[i + 1] = gS[i]
    return (gT, gS)
 
def Suscept(total_event, rate_0, initial_value):
    gT = np.zeros(total_event)
    gS = np.zeros(total_event)
    gS[0] = initial_value
    randomSeed = np.random.uniform(0, 1, total_event)
    for i in range(total_event - 1):
        # ending condition
        if gS[i] > 1:
            dt = (1.0/(gS[i]*rate_0)) * np.log(1.0/randomSeed[i])
            gT[i + 1] = gT[i] + dt 
            gS[i + 1] = gS[i] - 1
        else: 
            dt = 0.0
            gT[i + 1] = gT[i] + dt 
            gS[i + 1] = gS[i]
    return (gT, gS)

total_approach = 10
approach_T = np.zeros([total_approach, total_event])
approach_S = np.zeros([total_approach, total_event])
for i in range(total_approach):
    approach_T[i] = Suscept(total_event, rate_0, initial_value)[0]
    approach_S[i] = Suscept(total_event, rate_0, initial_value)[1]
    
# plotting
figure_name = ''
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
plot1 = figure.add_subplot(1, 2, 1)
plot1.plot(Susceptible(total_event, rate_0, initial_value)[0], Susceptible(total_event, rate_0, initial_value)[1]
           , label = 'way 0', drawstyle = 'steps')
plot1.grid(True)
plt.title(r'$ (events = {:}, \ rate = {:1.3f}) $'.format(total_event, rate_0), fontsize = AlvaFontSize)
plt.xlabel(r'$ time \ ({:})$'.format(timeUnit), fontsize = AlvaFontSize)
plt.ylabel(r'$ Population \ of \ S(t) $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.3))
plt.xticks(fontsize = AlvaFontSize*0.6)
plt.yticks(fontsize = AlvaFontSize*0.6) 

plot2 = figure.add_subplot(1, 2, 2)
for i in range(total_approach):
    plot2.plot(approach_T[i], approach_S[i], drawstyle = 'steps')
plot2.grid(True)
plt.title(r'$ (events = {:}, \ rate = {:1.3f}) $'.format(total_event, rate_0), fontsize = AlvaFontSize)
plt.xlabel(r'$ time \ ({:})$'.format(timeUnit), fontsize = AlvaFontSize)
plt.ylabel(r'$ Population \ of \ S(t) $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.3))
plt.xticks(fontsize = AlvaFontSize*0.6)
plt.yticks(fontsize = AlvaFontSize*0.6) 

figure.tight_layout()
plt.savefig(save_figure, dpi = 10)
plt.show()


# In[11]:

approach_S[0]


# In[14]:

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
total_event = 3000
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
plt.savefig(save_figure, dpi = 100)
plt.show()


# In[ ]:



