# plots.py

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pydataset import data

plt.ioff()

# Read in the data set.
crime = pd.read_csv("crime_data.txt", header=1, index_col=0)

# Basic plot w/ matplotlib
plt.plot(crime["Population"], label="Population")
plt.xlabel("Year")
plt.xlim(min(crime.index), max(crime.index))
plt.legend(loc="best")
plt.savefig("population.pdf")
plt.clf(); plt.close()

# Equivalent plot w/ Pandas
crime.plot(y="Population")
plt.clf(); plt.close()

# Other basic plots w/ matplotlib
plt.plot(crime["Population"], crime["Burglary"])
plt.savefig("pltBurglary.pdf")
plt.clf(); plt.close()

crime.plot(x="Population", y="Burglary")
plt.savefig("dfBurglary.pdf")
plt.clf(); plt.close()

crime.plot(subplots=True, layout=(4,3), linewidth=3, style='--', legend=False)
plt.savefig("subplots.pdf")
plt.clf(); plt.close()

# Histograms
crime.hist()
plt.savefig("hist1.pdf")
plt.clf(); plt.close()

crime.hist(column=['Population', 'Total', 'Property'], bins=15)
plt.savefig("hist2.pdf")
plt.clf(); plt.close()

crime[["Total", "Property"]].plot(kind="hist", alpha=.75)
plt.savefig("hist3.pdf")
plt.clf(); plt.close()

crime[["Total", "Property"]].plot(kind="hist", stacked=True,
                                        orientation="horizontal")
plt.savefig("hist4.pdf")
plt.clf(); plt.close()

# Bar Plots
crime.iloc[-5:][["Aggravated-assault", "Violent", "Burglary"]].plot(kind='bar')
# plt.savefig("bar1.pdf")
# plt.show()
plt.clf(); plt.close()

crime.iloc[-5:][["Aggravated-assault", "Violent", "Burglary"]].plot(kind='bar',
                                                    legend=False, stacked=True)
# plt.savefig("bar2.pdf")
# plt.show()
plt.clf(); plt.close()

crime[["Robbery", "Aggravated-assault", "Vehicle-Theft"]].plot(kind='box')
plt.savefig("box1.pdf")
plt.clf(); plt.close()
crime[["Robbery", "Aggravated-assault", "Vehicle-Theft"]].plot(kind='box',
                                                        vert=False, rot=30)
plt.savefig("box2.pdf")
plt.clf(); plt.close()
plt.show()

############# Middle section ############
hec = data("HairEyeColor")
X = np.unique(hec["Hair"], return_inverse=True)
Y = np.unique(hec["Eye"], return_inverse=True)
hec["Hair"] = X[1]
hec["Eye"] = Y[1]
hec.plot(kind="scatter", x="Hair", y="Eye", s=hec["Freq"]*20)
plt.xticks([0,1,2,3], X[0])
plt.yticks([0,1,2,3], Y[0])
plt.savefig("HairEyeColorscatter.png")

##########Plots for Data Visualization Section###########
data('Icecream').plot(kind='scatter', x='temp', y='cons')
plt.savefig('Nolabels.png')

msleep = data('msleep')
msleep.plot(y='sleep_total', title='Mammalian Sleep Data', legend=False)
plt.xlabel('Animal Index')
plt.ylabel('Sleep in Hours')
plt.savefig('Msleep1.png')

vore = msleep.groupby('vore')
means = vore.mean()
errors = vore.std()
means.loc[:,['sleep_total', 'sleep_rem','sleep_cycle']].plot(kind='bar', yerr=errors, title='Mean Mammallian Sleep Data')
plt.xlabel('Mammal diet classification (vore)')
plt.ylabel('Hours')
plt.tight_layout()
plt.savefig('MeanMammal.png')

Diamonds = data('diamonds')
DiaColor = Diamonds.groupby('color')
Ddiamond = DiaColor.get_group('D')
Jdiamond = DiaColor.get_group('J')
Ddiamond.plot(kind='hexbin', x='carat', y='price', gridsize=10, title='D Color Quality Diamonds: Best', cmap='viridis')
plt.ylabel('Price in US dollars')
plt.xlabel('Diamond Weight in Carats')
plt.tight_layout()
plt.savefig('DiamondD.png')
Jdiamond.plot(kind='hexbin', x='carat', y='price', gridsize=10, title='J Color Quality Diamonds: Worst', cmap='viridis')
plt.ylabel('Price in US dollars')
plt.xlabel('Diamond Weight in Carats')
plt.tight_layout()
plt.savefig('DiamondJ.png')