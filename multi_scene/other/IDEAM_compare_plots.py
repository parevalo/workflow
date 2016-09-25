__author__ = 'paulo'
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

os.chdir("/home/paulo/Downloads")
plt.style.use('ggplot')

# Store tables in a dictionary
areas={}

for i, file in enumerate(glob.glob('comparison_*.csv')):
    areas['comp_' + file[11:16]] = pd.read_csv(file, sep=',', header=0, index_col=0)

# Concatenate and transpose to get all areas in a single dataframe
areas_merged = pd.concat(areas, axis=1).transpose()
areas_merged.drop(areas_merged.columns[[3,10,11,12,13]], axis=1, inplace=True) # Drop unnecessary columns
areas_merged.index = areas_merged.index.droplevel(1) # Drop unnecessary index, could be prevented from input file
periods = ['03-05', '05-07', '07-09', '09-11', '11-13', '13-14', '14-15']
classes = ["Other changes", 'Forest agreement', "Deforestation agreement", "Regrowth agreement", "Stable non-forest agreement",
           "Deforestation (IDEAM) - Forest (BU)", "Deforestation (IDEAM) - Forest to secondary/secondary to others",
           "Stable non-forest (IDEAM) - Forest (BU)", "Stable non-forest (IDEAM) - Stable secondary (BU"]

areas_merged = **2 / 100**2
# Plot the areas! Keep in mind that there are still some disagreements in the total data area (0 to 14) that are
# probably caused by the inconsistent NoData areas in the IDEAM files.
for i in range(len(areas_merged-1)): # We don't want to plot column 14
    fig, ax1 = plt.subplots(1)
    ax1.plot(areas_merged.iloc[:,i], '-o') # We need to index it by location and not by name
    ax1.grid(True)
    ax1.set_xticklabels(periods)
    ax1.set_xlabel('Periods')
    ax1.set_ylabel('Area [ha]')
    # Format axis as hectares
    ax1.get_yaxis().set_major_formatter(tkr.FuncFormatter(lambda x, p: format(int(x*30**2 / 100**2), ',')))
    plt.title(classes[i])
    fig.savefig(classes[i], dpi=300, bbox_inches='tight')

#for k in areas:
#    print k, areas[k][0:14].sum()

