
# coding: utf-8

# # Stochastic infectious pulse
# https://github.com/alvason/stochastic-infectious-pulse
# 
# ### Stochastic version for evolutionary insights

# In[1]:

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
file_name = 'simple-stochastic-growth'

import alva_machinery_probability as alva

AlvaFontSize = 23
AlvaFigSize = (16, 8)
numberingFig = 0


# In[2]:

# algorithm for stochastic evolution
figure_name = '-Gillespie-algorithm'
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)
text_list = [r'$ Gillespie-algorithm: $'
             , r'$ 1. \ initialize \ the \ number \ of \ each \ group: \ S(t=0), I(t=0), R(t=0) $'
             , r'$ 2. \ compute \ the \ probability \ of \ each \ possible \ event_i \ at \ the \ moment \ \bf{t} $'
             , r'$ 3. \ randomly \ select \ event_{next} \
                   \ according \ to \ random{[0,1)} < \frac{\sum_{k=1}^{next}event_{k}}{\sum_{i=1}^{all} event_i} $'
             , r'$ 4. \ update \ the \ number \ of \ corresponding \ group $'
             , r'$ 5. \ compute \ \Delta t = \frac{-log_{e}(event_{next})}{\sum_{i}^{} event_i} $'
             , r'$ \ (according \ to \ probability-density-function: \ Pr(t < event_{next} < t+\Delta t) = \
                   exp(-\Delta t \sum_{i}^{} event_i )) $'
             , r'$ 7. \ update \ t = t + \Delta t $'
             , r'$ 6. \ go \ to \ step-2 $'
            ]
total_list = np.size(text_list)
numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize=(total_list, total_list*1.5))
plt.axis('off')
for i in range(total_list):
    plt.text(0, (total_list - float(i))/total_list
             , text_list[i].replace('\\\n', '')
             , fontsize = 1.2*AlvaFontSize)
plt.savefig(save_figure, dpi = 100)
plt.show()


# In[19]:

''' define simple stochastic function '''
def simpleStochasticGrowth(total_step, minT, maxT, initial_G, inRate, noRate):
    # intialized
    gT = np.zeros([total_step]) 
    gG = np.zeros([total_step]) 
    j = int(0)
    gT[j] = minT
    gG[j] = initial_G
    # all possible events
    event_in = inRate*gG[j]
    event_no = noRate*gG[j]
    # configuration table
    eventRate_updateNum = np.array([[event_in, +1]
                                  , [event_no, 0]])
    ###
    while (gT[j] < maxT):       
        # randomly choose event
        if np.random.random() < (eventRate_updateNum[0:1, 0].sum() / eventRate_updateNum[:, 0].sum()):
            k = 0
        else:
            k = 1
        # update number of section
        gG[j] = gG[j] + eventRate_updateNum[k, 1]
        # update event_rate
        event_in = inRate*gG[j]
        event_no = noRate*gG[j]
        # configuration table
        eventRate_updateNum = np.array([[event_in, +1]
                                      , [event_no, 0]])
        dt = -np.log(np.random.random()) / eventRate_updateNum[:, 0].sum()
        dt = -np.log(np.random.random()) / eventRate_updateNum[:, 0].sum() * 1
        # next step is based on current step
        gT[j + 1] = gT[j] + dt 
        gG[j + 1] = gG[j]
        j = j + 1
    # set the value of remaining steps = value of the last step (for ending)
    gT[j:] = gT[j]
    gG[j:] = gG[j]
    ###
    return(gT, gG)

''' starting from one infected '''
# setting parameter
timeUnit = 'day'
if timeUnit == 'day':
    day = 1
    year = 365 
elif timeUnit == 'year':
    year = 1
    day = float(1)/365 
    
total_N = 1000
initial_G = 1
# set parameter
inRate = float(1)/(day) # birth rate per year
noRate = float(1)/(day) # birth rate per year
# initial boundary condition
minT = float(0*day)
maxT = float(3*day)

total_step = int(maxT*total_N)
# stochastic evolution way
total_way = int(4)
gTT = np.zeros([total_way, total_step]) 
gGG = np.zeros([total_way, total_step]) 


for i in range(total_way):
    evolution = simpleStochasticGrowth(total_step, minT, maxT, initial_G, inRate, noRate)
    gTT[i] = evolution[0]
    gGG[i] = evolution[1]
# plotting
figure_name = '-growth'
file_suffix = '.png'
save_figure = os.path.join(dir_path, file_name + figure_name + file_suffix)

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = (7, 7))
for i in range(total_way):
    plt.plot(gTT[i], gGG[i], label = r'$ G_{:}(t) $'.format(i), drawstyle = 'steps') 
plt.plot(gTT[0], initial_G*np.exp(inRate*gTT[0]), linewidth = 6.0
         , color = 'black', alpha = 0.3, label = r'$ Natural \ Growth $') 
plt.grid(True)
plt.title(r'$ Simple \ Stochastic \ Growth $', fontsize = AlvaFontSize)
plt.xlabel(r'$ time \ ({:})$'.format(timeUnit), fontsize = AlvaFontSize)
plt.ylabel(r'$ Population $', fontsize = AlvaFontSize)
plt.legend(loc = (1,0))
plt.text(maxT, 1.0, r'$ \mu = {:} $'.format(inRate), fontsize = AlvaFontSize)
plt.xticks(fontsize = AlvaFontSize*0.7)
plt.yticks(fontsize = AlvaFontSize*0.7) 
figure.tight_layout()
plt.savefig(save_figure, dpi = 100, bbox_inches='tight')
plt.show()


# In[ ]:



