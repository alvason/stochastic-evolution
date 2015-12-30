
# coding: utf-8

# # Stochastic evolution
# https://github.com/alvason/stochastic-evolution
# 
# ### Evolutionary insights provided by stochastic tools

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
dir_path = '/Users/al/Desktop/GitHub/stochastic-evolution/figure'
file_name = 'gillespie-evolution'

import alva_machinery_probability as alva

AlvaFontSize = 23
AlvaFigSize = (16, 8)
numberingFig = 0

numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize=(12, 3))
plt.axis('off')
plt.title(r'$ Susceptible-Infectious-Recovered \ equation $', fontsize = AlvaFontSize)
plt.text(0, 2.0/3, r'$ \frac{\partial S(t)}{\partial t} =          -\beta S(t)I(t) +\mu N -\mu S(t)$', fontsize = 1.2*AlvaFontSize)
plt.text(0, 1.0/3, r'$ \frac{\partial I(t)}{\partial t} =          +\beta S(t)I(t) - \gamma I(t) -\mu I(t) $', fontsize = 1.2*AlvaFontSize)
plt.text(0, 0.0/3, r'$ \frac{\partial R(t)}{\partial t} =          +\gamma I(t) - \mu R(t) $', fontsize = 1.2*AlvaFontSize)
plt.show()


# In[23]:

''' starting from one infected '''
# setting parameter
timeUnit = 'day'
if timeUnit == 'day':
    day = 1
    year = 365 
elif timeUnit == 'year':
    year = 1
    day = float(1)/365 
    
initial_N = 100
initial_S = 98
initial_I = 2
initial_R = initial_N - initial_S - initial_I

reprodNum = float(1.5) # basic reproductive number R0: one infected person will transmit to 1.8 person 
recovRate = float(1)/(4*day) # 4 days per period ==> rate/year = 365/4
inOutRate = float(1)/(30*year) # birth rate per year
infecRate = reprodNum*(recovRate + inOutRate)/1 # per year, per person, per total-population

# initial boundary condition
minT = float(0*day)
maxT = float(40*day)

# stochastic evolution way
total_way = int(3)
total_step = int(3000)
gTT = np.zeros([total_way, total_step]) 
gNN = np.zeros([total_way, total_step])
gSS = np.zeros([total_way, total_step]) 
gII = np.zeros([total_way, total_step]) 
gRR = np.zeros([total_way, total_step]) 
gT = np.zeros([total_step]) 
gN = np.zeros([total_step])
gS = np.zeros([total_step]) 
gI = np.zeros([total_step]) 
gR = np.zeros([total_step]) 
for i in range(total_way):   
    j = int(0)
    # intialized
    gT[j] = minT
    gN[j] = initial_N
    gS[j] = initial_S
    gI[j] = initial_I
    gR[j] = initial_R 
    # storing
    gTT[i, j] = gT[j]
    gNN[i, j] = gN[j]
    gSS[i, j] = gS[j]
    gII[i, j] = gI[j]
    gRR[i, j] = gR[j] 
    # all possible events
    event_N = inOutRate*gN[j]
    event_S = inOutRate*gS[j]
    event_SI = infecRate*gS[j]*gI[j]/gN[j]
    event_I = inOutRate*gI[j]
    event_IR = recovRate*gI[j]
    event_R = inOutRate*gR[j]
    while (gT[j] < maxT):
        event_all = event_N + event_S + event_SI + event_I + event_IR + event_R          
        dt = -np.log(np.random.random())/event_all 
        # S in event
        if np.random.random() < (event_N/event_all):                      
            gN[j] = gN[j] + 1 
            gS[j] = gS[j] + 1 
            event_N = inOutRate*gN[j]
            event_S = inOutRate*gS[j]
        # S out event    
        elif np.random.random() < ((event_N + event_S)/event_all):   
            gN[j] = gN[j] - 1 
            gS[j] = gS[j] - 1 
            event_N = inOutRate*gN[j]
            event_S = inOutRate*gS[j]
        # infected event
        elif np.random.random() < ((event_N + event_S + event_SI)/event_all):   
            gS[j] = gS[j] - 1 
            gI[j] = gI[j] + 1 
            event_S = inOutRate*gS[j]
            event_I = inOutRate*gI[j]
            event_IR = recovRate*gI[j]
        # I out event
        elif np.random.random() < ((event_N + event_S + event_SI + event_I)/event_all):   
            gN[j] = gN[j] - 1 
            gI[j] = gI[j] - 1 
            event_I = inOutRate*gI[j]
            event_IR = recovRate*gI[j]
        # recovered event
        elif np.random.random() < ((event_N + event_S + event_SI + event_I + event_IR)/event_all):    
            gI[j] = gI[j] - 1
            gR[j] = gR[j] + 1
            event_I = inOutRate*gI[j]
            event_IR = recovRate*gI[j]
            event_R = inOutRate*gR[j]
        # R out event
        else:
            gN[j] = gN[j] - 1 
            gR[j] = gR[j] - 1
            event_R = inOutRate*gR[j] 
        # update SI-rate after all
        event_SI = infecRate*gS[j]*gI[j]/gN[j]
        # next step is based on current step
        j = j + 1
        gT[j] = gT[j - 1] + dt 
        gN[j] = gN[j - 1]
        gS[j] = gS[j - 1]
        gI[j] = gI[j - 1]
        gR[j] = gR[j - 1]
    # set the value of remaining steps = value of the last step (for ending)
    gT[j:] = gT[j]
    gN[j:] = gN[j]
    gS[j:] = gS[j]
    gI[j:] = gI[j]
    gR[j:] = gR[j]
    # storing
    gTT[i] = gT
    gNN[i] = gN
    gSS[i] = gS
    gII[i] = gI
    gRR[i] = gR 

# plotting
figure_name = '-SIR'
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
for i in range(total_way):
    plt.plot(gTT[i], gSS[i], drawstyle = 'steps', label = r'$ S_{:}(t) $'.format(i), linewidth = (1 + i)
             , color = 'blue', alpha = float(0.5 + i/total_way))
    plt.plot(gTT[i], gII[i], drawstyle = 'steps', label = r'$ I_{:}(t) $'.format(i), linewidth = (1 + i)
             , color = 'red', alpha = float(0.5 + i/total_way))
    plt.plot(gTT[i], gRR[i], drawstyle = 'steps', label = r'$ R_{:}(t) $'.format(i), linewidth = (1 + i)
             , color = 'green', alpha = float(0.5 + i/total_way))
    plt.plot(gTT[i], (gSS[i] + gII[i] + gRR[i]), drawstyle = 'steps', label = r'$ N_{:}(t) $'.format(i)
             , linewidth = (1 + i), color = 'black', alpha = float(0.5 + i/total_way))    
plt.grid(True)
plt.title(r'$ Stochastic \ Susceptible-Infected-Recovered $', fontsize = AlvaFontSize)
plt.xlabel(r'$ time \ ({:})$'.format(timeUnit), fontsize = AlvaFontSize)
plt.ylabel(r'$ Population $', fontsize = AlvaFontSize)
plt.ylim(0, initial_N*1.1)
plt.legend(loc = (1,0))
plt.text(maxT, initial_N*6.0/6, r'$ R_0 = %f $'%(reprodNum), fontsize = AlvaFontSize)
plt.text(maxT, initial_N*5.0/6, r'$ \gamma = %f $'%(recovRate), fontsize = AlvaFontSize)
plt.text(maxT, initial_N*4.0/6, r'$ \beta = %f $'%(infecRate), fontsize = AlvaFontSize)
plt.text(maxT, initial_N*3.0/6, r'$ \mu = %f $'%(inOutRate), fontsize = AlvaFontSize)
plt.xticks(fontsize = AlvaFontSize*0.7)
plt.yticks(fontsize = AlvaFontSize*0.7) 
figure.tight_layout()
plt.savefig(save_figure, dpi = 100)
plt.show()


# In[ ]:



