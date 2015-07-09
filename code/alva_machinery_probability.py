# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Probability-insighter
# https://github.com/alvason/probability-insighter
# 
# ### home-made machinery for insight into randomness

# <codecell>

'''
author: Alvason Zhenhua Li
date:   04/16/2015

Home-made machinery for sorting a list from min-max
'''
import numpy as np

# min-max sorting
def AlvaMinMax(data):
    totalDataPoint = np.size(data)
    minMaxListing = np.zeros(totalDataPoint)   
    for i in range(totalDataPoint):
        # searching the minimum in current array
        jj = 0 
        minMaxListing[i] = data[jj] # suppose the 1st element [0] of current data-list is the minimum
        for j in range(totalDataPoint - i):
            if data[j] < minMaxListing[i]: 
                minMaxListing[i] = data[j]
                jj = j # recording the position of selected element
        # reducing the size of searching zone (removing the minmum from current array)
        data = np.delete(data, jj)
    return (minMaxListing)

'''
author: Alvason Zhenhua Li
date:   04/16/2015

Home-made machinery for leveling a list by using min-max way
'''
# leveling by using min-max way
def AlvaLevel(data, totalLevel, normalization = True):
    totalDataPoint = np.size(data)
    minMaxListing = AlvaMinMax(data)
    # searching minimum and maximum values
    minValue = minMaxListing[0]
    maxValue = minMaxListing[-1]
    spacingValue = np.linspace(minValue, maxValue, num = totalLevel + 1, retstep = True)        
    gLevel = np.delete(spacingValue[0], 0)
    # catogerizing the level set
    # initialize the levelspace by a 'null' space
    levelSpace = np.zeros([2])
    numberLevel = np.zeros([totalLevel])
    jj = 0 # counting the checked number
    for i in range(totalLevel): 
        n = 0 # counting the number in each level
        for j in range(jj, totalDataPoint):
            if minMaxListing[j] <= gLevel[i]: 
                levelSpace = np.vstack((levelSpace, [i, minMaxListing[j]]))
                n = n + 1
        numberLevel[i] = n
        jj = jj + n
    # delete the inital 'null' space
    levelSpace = np.delete(levelSpace, 0, 0) 
    if normalization == True:
        numberLevel = numberLevel/AlvaMinMax(numberLevel)[-1]
    return (gLevel, numberLevel, levelSpace)

