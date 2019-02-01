# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import collections as col
import random as r
import csv
import sys
import itertools

def hungergames(filename):
    
    events = list() #put events in a list
    with open('events.txt', 'r') as event:
        for word in event:
            events.append(word.strip())
    
    male = range(1,13)
    r.shuffle(male)
    female = range(1,13)
    r.shuffle(female)
    maletributes = list()
    femaletributes = list()
    Tribute = col.namedtuple('Tribute', 'name, gender, district')
    with open('tributes.csv', 'r') as trib:
        csv_reader = csv.reader(trib, delimiter = ',')
        rowtracker = 0
        for row in csv_reader:
            boy = Tribute(row[0], 'boy', male[rowtracker])
            girl = Tribute(row[1], 'girl', female[rowtracker])
            maletributes.append(boy)
            femaletributes.append(girl)
            rowtracker += 1
        
    tributes = maletributes + femaletributes
    
    days=itertools.count(start=1,step=1)
    threshold = 0.5
    while len(tributes) > 1:
        daynumber = days.next()
        with open(filename, 'a') as f:
            f.write("Day " + str(daynumber) + ":\n\n")
        currentsurvivors = len(tributes)
        dead = list()
        for i in range(currentsurvivors):
            death = r.random()
            if death < threshold:
                dead.append(tributes[i])
                with open(filename, 'a') as f:
                    f.write(tributes[i].name + " experienced " + r.choice(events) + " and died\n")
            else:
                with open(filename, 'a') as f:
                    f.write(tributes[i].name + " experienced " + r.choice(events) + " and survived\n")
        for i in range(len(dead)):
            tributes.remove(dead[i])
        toll = currentsurvivors - len(tributes)
        with open(filename, 'a') as f:
            f.write("\nEnd of day " + str(daynumber) + ": cannon fired " + str(toll) + " times\n\n")
      
    if len(tributes) == 0:
        with open(filename, 'a') as f:
            f.write("This hunger games left no survivors\n")
    if len(tributes) == 1:
        with open(filename, 'a') as f:
            f.write("The victor of this hunger games is the " + tributes[0].gender + " from District " + str(tributes[0].district) + ": " + tributes[0].name + "\n")

if __name__ == '__main__':
    from sys import argv
    filename = sys.argv[1]
    hungergames(filename)

# <codecell>


