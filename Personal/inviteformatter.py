import csv

names = []
addresses = []
citystatezips = []
people = []
with open('responses.csv') as mylist:
    for line in mylist:
        people.append(line.strip())
        
for i in xrange(0, len(people)-3, 3):
    if len(people[i]) < 4:
        people.pop(i)
    names.append(people[i])
    addresses.append(people[i+1])
    citystatezips.append(people[i+2])
    
for i in xrange(len(names)):
    citystatezips[i] = citystatezips[i].replace('Utah', 'UT')
    names[i] = names[i].replace('"','')
    names[i] = names[i].replace('family', 'Family')
    addresses[i] = addresses[i].replace('Dr', 'Drive')
    addresses[i] = addresses[i].replace('Driveive', 'Drive')
    addresses[i] = addresses[i].replace('Cir', 'Circle')
    addresses[i] = addresses[i].replace('Circlecle', 'Circle')
    addresses[i] = addresses[i].replace('Ct', 'Court')
    addresses[i] = addresses[i].replace('.','')
    addresses[i] = addresses[i].replace('Rd','Road')
    addresses[i] = addresses[i].replace('NE', 'Northeast')
    addresses[i] = addresses[i].replace('Ln','Lane')
    addresses[i] = addresses[i].replace(' W ',' West ')
    addresses[i] = addresses[i].replace(' E ',' East ')
    addresses[i] = addresses[i].replace(' N ',' North ')
    addresses[i] = addresses[i].replace(' S ',' South ')
    if addresses[i][-1] == 'W':
        addresses[i] = addresses[i][:-1] + 'West'
    if addresses[i][-1] == 'S':
        addresses[i] = addresses[i][:-1] + 'South'
    if addresses[i][-1] == 'E':
        addresses[i] = addresses[i][:-1] + 'East'
    if addresses[i][-1] == 'N':
        addresses[i] = addresses[i][:-1] + 'North'
    addresses[i] = addresses[i].replace('Sgate', 'Southgate')
    addresses[i] = addresses[i].replace('St','Street')
    addresses[i] = addresses[i].replace('Streetreet', 'Street')
    addresses[i] = addresses[i].replace('Streetate', 'State')
    addresses[i] = addresses[i].replace('Streetonewood', 'Stonewood')
    addresses[i] = addresses[i].replace('Ct','Court')
    addresses[i] = addresses[i].replace(' w ', ' West ')
    addresses[i] = addresses[i].replace(' s', ' South ')
    addresses[i] = addresses[i].replace(' n ', ' North ')
    addresses[i] = addresses[i].replace(' e ', ' East ')
    addresses[i] = addresses[i].replace('P O', 'P.O.')
    citystatezips[i] = citystatezips[i].replace('"','')
    citystatezips[i] = citystatezips[i].replace(',','')
    citystatezips[i] = citystatezips[i].replace('.','')
    citystatezips[i] = citystatezips[i].replace('  ',' ')
    citystatezips[i] = citystatezips[i].replace('St ', 'St. ')

cities = []
states = []
zips = []
for i in xrange(len(citystatezips)):
    states.append(citystatezips[i][-9:-6])
    states[i] = states[i].replace(' ','')
    zips.append(citystatezips[i][-6:])
    zips[i] = zips[i].replace(' ','')
    cities.append(citystatezips[i][0:-9])
    if cities[i][-1] == ' ':
        cities[i] = cities[i][0:-1]

with open('Invitation Addresses.csv','w') as output:
    outwriter = csv.writer(output, lineterminator = '\n')
    for i in xrange(len(zips)):
        outwriter.writerow([names[i], addresses[i], cities[i], states[i], zips[i]])
    print 'Done!'




