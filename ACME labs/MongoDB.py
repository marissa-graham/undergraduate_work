import pymongo,json
import re

conn = pymongo.MongoClient()
db = conn.mydb
collection = db.rest

# Dropped : 2nd and 3rd parts of 3, third part of 4

# Problem 1
with open('restaurants.json') as data_file:
    content = data_file.readlines()
    
for c in content:
    collection.insert(json.loads(c))

# Problem 2
with open('mylansbistro.json') as new_data:
	new_content = new_data.readlines()

for c in new_content:
	collection.insert(json.loads(c))

print "Restaurants that close at 1800:"
A = collection.find({'closing_time':'1800'})
for thing in A:
	print thing['name']

# Problem 3-4
print "Number of restaurants in Manhattan:"
print collection.find({'borough':'Manhattan'}).count()

print "Number of restaurants with 'grill' in their name:"
grills = collection.find({'name': re.compile('grill')})
for grill in grills:
	print grill['name']
	grill['name'] = grill['name'.replace(re.compile('Grill'),'Magical Fire Table')
	print grill['name']

has_id = collection.find({'restaurant_id':'$exists'})
for restaurant in has_id:
	id = int(restaurant['restaurant_id'])
	id += 1000
	restaurant['restaurant_id'] = str(id)

conn.drop_database('mydb')