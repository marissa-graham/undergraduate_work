import pandas as pd
from pydataset import data
from matplotlib import pyplot as plt
import numpy as np
%matplotlib inline


# Problem 1
print "Problem 1:"

nottem_data = data('nottem')
VADeaths_data = data('VADeaths')
Arbuthnot_data = data('Arbuthnot')

#data('nottem', show_doc=True)
print 'Nottem:'
print 'This data set contains the average air temperatures at Nottingham Castle for 20 years.'
print 'Since this data set consists of a single series, a line plot is best.'
nottem_data.plot(y='nottem',legend=False)
plt.show()

#data('VADeaths', show_doc=True)
print 'VA Deaths:'
print 'This data set contains death rates in Virginia in 1940, cross-classified \
by age group, gender, and urban/rural.'
print 'Since this data set contains multiple categories, a bar plot is more effective.'
VADeaths_data.iloc[:][['Rural Male','Rural Female','Urban Male','Urban Female']].plot(kind='bar')
plt.show()

#data('Arbuthnot', show_doc=True)
#Arbuthnot_data.plot(subplots=True, layout=(7,1))
print 'Arbuthnot:'
print 'This data set contains information on male and female birth ratios in London from 1629-1710,\
as well as deaths from plague and total mortality over the same time period.'
print 'We plot each category against the year. Since there are a large number of years, a bar plot \
becomes unwieldy, so we simply use line plots. Related categories are plotted on the same chart.'
Arbuthnot_data.iloc[:][['Males','Females']].plot()
plt.show()
Arbuthnot_data.plot(x='Year',y='Ratio')
Arbuthnot_data.plot(x='Year',y='Total')
Arbuthnot_data.iloc[:][['Plague','Mortality']].plot()
plt.show()

# Problem 2
print "Problem 2:"

trees_data = data('trees')
# Girth, Height, Volume

road_data = data('road')
# state, deaths, drivers, popden, rural, temp, fuel

birthdeathrates_data = data('birthdeathrates')
# birth, death (over 69 observations)


print 'Tree data:'
print 'In this case, a box plot is more appropriate because it gives us a rough idea of\
the structure of each category of data, as well as how they compare to each other.'
trees_data[['Girth','Height','Volume']].plot(kind='hist',alpha=.75)
trees_data[['Girth','Height','Volume']].plot(kind='box')
plt.show()

print 'Road data:'
road_data[['deaths','popden','drivers','rural', 'temp','fuel']].plot(kind='hist',stacked=True)
road_data[['drivers','rural', 'temp','fuel']].plot(kind='box')
plt.show()
road_data[['deaths','popden']].plot(kind='box')
plt.show()

print 'Birth and death rates data:'
birthdeathrates_data[['birth','death']].plot(kind='hist',stacked=True,alpha=.75)
birthdeathrates_data[['birth','death']].plot(kind='box')

# Problem 3 is very open ended and easy enough to produce.

# Problem 5 is similar to the example given in the lab.