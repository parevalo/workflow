#!/usr/env/python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

plt.style.use('ggplot')

# Read area files and create necessary variables for plotting
area = pd.read_csv("area_ha.csv", sep=',', header=0, index_col=0)
area_lower = pd.read_csv("area_lower.csv", sep=',', header=0, index_col=0)
area_upper = pd.read_csv("area_upper.csv", sep=',', header=0, index_col=0)
years = range(2002, 2017)
strata_names = ("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                "Stable pasture-cropland", "Stable regrowth", "Stable water", "Forest to pasture", 
                "Forest to regrowth", "All others to regrowth", "Loss of regrowth")

# Margin of error in percentage= area_ci / area
margin_error = (area_upper - area)/area *100


# In[24]:

# Fnc to change the ticks to the desired color
def axiscolors(axis, color):
    for tl in axis.get_yticklabels():
        tl.set_color(color)

# Plots
for i in range(len(area.columns)):
    # Create the plots
    fig, ax1 = plt.subplots(1)
    
    ax1.plot(years, area.iloc[:,i], '-o', color='b')
    ax1.fill_between(years,area_lower.iloc[:,i], area_upper.iloc[:,i], alpha=0.5)
    ax1.set_ylim(ymin=min(area_lower.min()), ymax=area_upper.iloc[:,i].max()*1.1)
    ax1.set_xlabel('Years')
    ax1.set_ylabel('Area [ha]', color='b')
    ax1.get_yaxis().set_major_formatter(tkr.FuncFormatter(lambda x, p: format(int(x), ',')))
    axiscolors(ax1, 'b')    

    ax2 = ax1.twinx()
    ax2.plot(years, margin_error.iloc[:,i], 'r-')
    ax2.set_ylabel('Margin of error [%]', color='r')
    ax2.set_ylim(ymin=0, ymax=margin_error.iloc[:,i].max()*1.1)
    axiscolors(ax2, 'r')    

    # Save fig
    fig.savefig(strata_names[i]+"_python", dpi=300, bbox_inches='tight')


