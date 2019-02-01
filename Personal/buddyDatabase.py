import csv
import numpy as np
import geopy

def get_gender(gender):
	"""
	-1 is prefer not to answer
	0  is female
	1  is male to female transgender
	2  is genderqueer
	3  is female to male transgender
	4  is male
	"""
	if gender.startswith("Prefer"): 
		return -1
	elif gender == "Female": 
		return 0		
	elif gender.startswith("Male to"): 
		return 1
	elif gender.startswith("Genderqueer"): 
		return 2
	elif gender.startswith("Female to"): 
		return 3
	elif gender == "Male": 
		return 4

def get_pair_gender(ans):
	"""True: pair with the same gender identity."""
	if len(ans) > 1:
		return True
	else:
		return False

def get_height(height):
	""" 0: < 5'
		1: 5'0-5'3
		2: 5'4-5'7
		3: 5'8-6'0
		4: 6'0+
	"""
	if height.startswith('< 5'):
		return 0
	elif height.startswith("5'0"):
		return 1
	elif height.startswith("5'4"):
		return 2
	elif height.startswith("5'8"):
		return 3
	elif height.startswith("6'0"):
		return 4

def get_fitness_goals(goals):
	"""
	0: Better health
	1: Mental health
	2: Get strong
	3: Look better naked
	4: Lose weight
	5: Gain weight
	6: Better endurance
	7: Training for competitions
	8: Be better at sport/activity
	"""
	goals = goals.replace('(relaxation, anxiety/depression/ED management, etc.)','')
	goals = goals.split(',')
	out = np.zeros(9,dtype=bool)
	while len(goals) > 0:
		if goals[-1].startswith("Better health"):
			out[0] = True
		elif goals[-1].startswith("Mental health"):
			out[1] = True
		elif goals[-1].startswith("Get strong"):
			out[2] = True
		elif goals[-1].startswith("Look better"):
			out[3] = True
		elif goals[-1].startswith("Lose weight"):
			out[4] = True
		elif goals[-1].startswith("Gain weight"):
			out[5] = True
		elif goals[-1].startswith("Better endurance"):
			out[6] = True
		elif goals[-1].startswith("Training for"):
			out[7] = True
		elif goals[-1].startswith("Be better at"):
			out[8] = True
		goals = goals[:-1]
	return out

def get_competitions(comps):
	"""
	0: Powerlifting
	1: Olympic weightlifting
	2: 5/10k
	3: Half/Full Marathons
	4: Triathlons
	5: Bodybuilding
	6: Strongman
	7: Dance/performance
	"""
	comps = comps.split(',')
	out = np.zeros(9,dtype=bool)
	while len(comps) > 0:
		if "Powerlifting" in comps[-1]:
			out[0] = True
		elif "Olympic" in comps[-1]:
			out[1] = True
		elif "5/10K" in comps[-1]:
			out[2] = True
		elif "Half/Full" in comps[-1]:
			out[3] = True
		elif "Triathlons" in comps[-1]:
			out[4] = True
		elif "Bodybuilding" in comps[-1]:
			out[5] = True
		elif "Strongman" in comps[-1]:
			out[6] = True
		elif "Dance/performance" in comps[-1]:
			out[7] = True
		elif len(comps[-1]) > 1:
			out[8] = True
		comps = comps[:-1]
	return out

def get_workout_freq(freq):
	"""
	0: Workout 1-2 times per week
	1: Workout 3-4 times per week
	2: Workout 5-6 times per week
	3: Workout 7+ times per week
	"""
	if freq.startswith('1-2'):
		return 0
	elif freq.startswith('3-4'):
		return 1
	elif freq.startswith('5-6'):
		return 2
	elif freq.startswith('7+'):
		return 3

def get_training_age(age):
	"""
	0: Just starting
	1: 1-6 months
	2: 6-12 months
	3: 1-2 yrs
	4: 2-4 yrs
	5: 4+ yrs
	"""
	if age.startswith("Just"):
		return 0
	elif age.startswith("1-6"):
		return 1
	elif age.startswith("6-12"):
		return 2
	elif age.startswith("1-2"):
		return 3
	elif age.startswith("2-4"):
		return 4
	elif age.startswith("4+"):
		return 5

def get_athlete_type(types):
	"""
	0: Running
	1: General lifting
	2: Powerlifting
	3: Olympic weightlifting
	4: Bodybuilding
	5: Yoga
	6: Sports
	7: Biking
	8: Crossfit
	9: Rock climbing
	10: Dance/performance arts
	11: Strength group classes
	12: Cardio group classes
	13: Outdoor adventuring
	14: Combat/martial arts
	15: Other
	"""
	types = types.replace('dance, pole, aerial','dance/pole/aerial')
	types = types.replace('TRX, BodyPump','TRX/BodyPump')
	types = types.replace('Zumba, spinning','Zumba/spinning')
	types = types.replace('biking, kayaking, hiking','biking/kayaking/hiking')
	types = types.split(',')
	out = np.zeros(16,dtype=bool)
	while len(types) > 0:
		if "Running" in types[-1]:
			out[0] = True
		elif "General lifting" in types[-1]:
			out[1] = True
		elif "Powerlifting" in types[-1]:
			out[2] = True
		elif "Olympic" in types[-1]:
			out[3] = True
		elif "Bodybuilding" in types[-1]:
			out[4] = True
		elif "Yoga" in types[-1]:
			out[5] = True
		elif "Sports" in types[-1]:
			out[6] = True
		elif "Biking" in types[-1]:
			out[7] = True
		elif "Crossfit" in types[-1]:
			out[8] = True
		elif "Rock Climbing" in types[-1]:
			out[9] = True
		elif "Dance/Performance" in types[-1]:
			out[10] = True
		elif "Strength group" in types[-1]:
			out[11] = True
		elif "Cardio group" in types[-1]:
			out[12] = True
		elif "Outdoor Adventuring" in types[-1]:
			out[13] = True
		elif "Combat/martial arts" in types[-1]:
			out[14] = True
		elif len(types[-1]) > 1:
			out[15] = True
		types = types[:-1]
	return out

def get_diet_type(diet):
	"""
	0: Nothing specific, not focused on diet.
	1: Just trying to eat better
	2: Calorie counting
	3: IIFYM(If it fits your macros)/Flexible dieting
	4: Paleo/primal
	5: Low-carb/Atkins
	6: Keto
	7: Intermittent Fasting
	8: Clean eating
	9: Celiac/Gluten-free
	10: IBD/IBS
	11: Vegetarian
	12: Vegan
	13: Intuitive eating
	"""
	diet = diet.replace('Nothing specific, not','Not')
	diet = diet.split(',')
	out = np.zeros(14,dtype=bool)
	while len(diet) > 0:
		if "Not focused" in diet[-1]:
			out[0] = True
		elif "Just trying" in diet[-1]:
			out[1] = True
		elif "Calorie counting" in diet[-1]:
			out[2] = True
		elif "IIFYM" in diet[-1]:
			out[3] = True
		elif "Paleo" in diet[-1]:
			out[4] = True
		elif "Low-carb" in diet[-1]:
			out[5] = True
		elif "Keto" in diet[-1]:
			out[6] = True
		elif "Intermittent" in diet[-1]:
			out[7] = True
		elif "Clean" in diet[-1]:
			out[8] = True
		elif "Celiac" in diet[-1]:
			out[9] = True
		elif "IBD/IBS" in diet[-1]:
			out[10] = True
		elif "Vegetarian" in diet[-1]:
			out[11] = True
		elif "Vegan" in diet[-1]:
			out[12] = True
		elif "Intuitive" in diet[-1]:
			out[13] = True
		diet = diet[:-1]
	return out

def get_weight_goal(goal):
	"""
	-5: Lose more than 100lbs
	-4: Lose 50-100lbs
	-3: Lost 25-49lbs
	-2: Lose 10-25lbs
	-1: Lose 10lbs or less
	0: Happy to maintain
	1: Gain 5-10lbs
	2: Gain 11-25lbs
	3: Gain 25-49lbs
	4: Gain more than 50lbs
	"""
	if goal.startswith("Gain more than"):
		return 4
	elif goal.startswith("Gain 25"):
		return 3
	elif goal.startswith("Gain 11"):
		return 2
	elif goal.startswith("Gain 5"):
		return 1
	elif goal.startswith("Happy"):
		return 0
	elif goal.startswith("Lose 10lbs"):
		return -1
	elif goal.startswith("Lose 10-25"):
		return -2
	elif goal.startswith("Lost 25-49"):
		return -3
	elif goal.startswith("Lose 50"):
		return -4
	elif goal.startswith("Lose more"):
		return -5

def get_special(special):
	"""
	0: I'm currently pregnant
	1: I had a baby in the last 6 months
	2: I'm recovering from an eating disorder
	3: I have PCOS
	4: I'm getting in shape for an upcoming event
	5: I'm injured/recovering from an injury
	6: I'm diabetic
	"""
	special = special.split(',')
	out = np.zeros(7,dtype=bool)
	while len(special) > 0:
		if "currently pregnant" in special[-1]:
			out[0] = True
		elif "had a baby" in special[-1]:
			out[1] = True
		elif "eating disorder" in special[-1]:
			out[2] = True
		elif "PCOS" in special[-1]:
			out[3] = True
		elif "getting in shape" in special[-1]:
			out[4] = True
		elif "injured" in special[-1]:
			out[5] = True
		elif "diabetic" in special[-1]:
			out[6] = True
		special = special[:-1]
	return out

def get_buddy_prefs(prefs):
	"""
	0: Someone to text/chat with about all kinds of random fitness stuff.
	1: Someone to keep me accountable and motivated.
	2: Someone to give me advice.
	3: Someone to give advice to.
	4: Someone to work out with.
	5: Someone to talk about non-fitness stuff with.
	"""
	prefs = prefs.split(',')
	out = np.zeros(6,dtype=bool)
	while len(prefs) > 0:
		if "text" in prefs[-1]:
			out[0] = True
		elif "accountable" in prefs[-1]:
			out[1] = True
		elif "give me advice" in prefs[-1]:
			out[2] = True
		elif "give advice" in prefs[-1]:
			out[3] = True
		elif "work out with" in prefs[-1]:
			out[4] = True
		elif "talk about non-fitness stuff" in prefs[-1]:
			out[5] = True
		prefs = prefs[:-1]
	return out

def get_personality(ptype):
	"""
	-1: None of the above/Other
	0: Introvert
	1: Ambivert (situation-dependent)
	2: Extrovert
	"""
	if ptype.startswith("Extrovert"):
		return 2
	elif ptype.startswith("Ambivert"):
		return 1
	elif ptype.startswith("Introvert"):
		return 0
	else:
		return -1

def get_comm_freq(comm_freq):
	"""
	0: Once a week
	1: Several times a week
	2: Every day
	3: Several times a day
	"""
	if comm_freq.startswith("Once a week"):
		return 0
	elif comm_freq.startswith("Several times a week"):
		return 1
	elif comm_freq.startswith("Every day"):
		return 2
	elif comm_freq.startswith("Several times a day"):
		return 3
	else:
		return None

def get_contact(contact):
	"""
	0: Reddit Mail
	1: Email
	2: Text
	3: Facebook
	4: Twitter
	5: Instagram
	6: MyFitnessPal
	7: Fitocracy
	8: Snapchat
	"""
	contact = contact.split(',')
	out = np.zeros(9,dtype=bool)
	while len(contact) > 0:
		if "Reddit Mail" in contact[-1]:
			out[0] = True
		elif "Email" in contact[-1]:
			out[1] = True
		elif "Text" in contact[-1]:
			out[2] = True
		elif "Facebook" in contact[-1]:
			out[3] = True
		elif "Twitter" in contact[-1]:
			out[4] = True
		elif "Instagram" in contact[-1]:
			out[5] = True
		elif "MyFitnessPal" in contact[-1]:
			out[6] = True
		elif "Fitocracy" in contact[-1]:
			out[7] = True
		elif "Snapchat" in contact[-1]:
			out[8] = True
		contact = contact[:-1]
	return out

def get_local_pref(pref):
	"""
	0: I could not care less.
	1: Same timezone is enough for me.
	2: Important. I want an in-person buddy.
	"""
	if pref.startswith("Important"):
		return 2
	elif pref.startswith("Same timezone"):
		return 1
	else:
		return 0

def get_priorities(priorities):
	"""
	0: Location (ignore?)
	1: Age
	2: How long we've been working out
	3: Frequency of working out
	4: Similar fitness goals
	5: Similar fitness interests
	6: Similar weight loss/gain goals
	7: Training for the same type of competition
	8: Height
	9: Dietary preferences
	10: Personality type (introvert/extrovert)
	11: Special circumstances as listed above
	"""
	pris = priorities.split(',')
	out = np.zeros(12,dtype=bool)
	while len(pris) > 0:
		if "Location" in pris[-1]:
			out[0] = True
		elif "Age" in pris[-1]:
			out[1] = True
		elif "How long" in pris[-1]:
			out[2] = True
		elif "Frequency of working out" in pris[-1]:
			out[3] = True
		elif "Similar fitness goals" in pris[-1]:
			out[4] = True
		elif "Similar fitness interests" in pris[-1]:
			out[5] = True
		elif "Similar weight loss" in pris[-1]:
			out[6] = True
		elif "Training for the same" in pris[-1]:
			out[7] = True
		elif "Height" in pris[-1]:
			out[8] = True
		elif "Dietary" in pris[-1]:
			out[9] = True
		elif "Personality" in pris[-1]:
			out[10] = True
		elif "Special" in pris[-1]:
			out[11] = True
		pris = pris[:-1]
	return out

def get_time_zone(time_zone):
	if len(time_zone) > 1:
		zone = time_zone[3:6]
		zone = zone.replace(':','')
		return int(zone)
	else:
		return None

def get_group_share(share_opt):
	if share_opt.startswith("Yes"):
		return True
	else:
		return False

def get_top_five(top_five):
	if top_five == "Yes":
		return True
	else:
		return False

class Person(object):
	def __init__(self, attrs):
		"""
		attrs is the list of responses to the survey.

		userID: Reddit username (string) 
		age: Age (int)
		gender: Gender identity (int)
		pair_gender: Pair with same gender (pair_gender, bool)
		height: Height (int)
		fitness_goals: Fitness goals (array<int>) 
		competitions: Fitness competitions (array<int>)
		workout_freq: Workout frequency (int)
		training_age: training age (int)
		athlete_type: exercise types (array<int>)
		diet_type: diet type (array<int>)
		attrs[13]: food allergies (none found)
		weight_goal: weight goals (int)
		special: special circumstances (array<int>)
		buddy_prefs: buddy preferences (array<int>)
		ptype: introvert/extrovert (int)
		comm_freq: communication frequency (int)
		contact: contact type (int)
		local: local_buddy (int)
		priorities: buddy priorities (array<int> (or string?))
		attrs[22]: street_address (address, string)
		attrs[23]: city (city, string)
		attrs[24]: state (state, string)
		attrs[25]: zip (zip, string or int?)
		attrs[26]: country (country, string)
		time_zone: time zone (int) (+/-, UTC)
		group_share: group sharing (bool)
		get_five: top five (bool)
		notes: "Anything else we should know?" (string)
		"""
		self.userID = attrs[1]
		if attrs[2] != '':
			self.age = int(attrs[2])
		else:
			self.age = 0
		self.gender = get_gender(attrs[3])
		self.pair_gender = get_pair_gender(attrs[4])
		self.height = get_height(attrs[5])
		self.fitness_goals = get_fitness_goals(attrs[6])
		self.competitions = get_competitions(attrs[7])
		self.workout_freq = get_workout_freq(attrs[8])
		self.training_age = get_training_age(attrs[10])
		self.athlete_type = get_athlete_type(attrs[11])
		self.diet_type = get_diet_type(attrs[12])
		self.weight_goal = get_weight_goal(attrs[14])
		self.special = get_special(attrs[15])
		self.buddy_prefs = get_buddy_prefs(attrs[16])
		self.ptype = get_personality(attrs[17])
		self.comm_freq = get_comm_freq(attrs[18])
		self.contact = get_contact(attrs[19])
		self.local = get_local_pref(attrs[20])
		self.priorities = get_priorities(attrs[21])
		# More location stuff here
		self.time_zone = get_time_zone(attrs[27])
		self.group_share = get_group_share(attrs[28])
		self.get_five = get_top_five(attrs[29])
		self.notes = attrs[31]
	
# A Survey class, which consists of a bunch of People, and the 
# rank_buddies and sort_buddies functions

# Each criteria will have its own subfunction

# Get the spreadsheet
def get_data(filename):
	# Filename needs to be TAB SEPARATED
	people = []
	with open(filename) as mylist:
		people_attrs = []
		for line in mylist:
			people.append(line.split('\t'))
	# Get rid of questions row and test answer
	people = people[2:]
	return people

def get_People(people):
	database = []
	for i in xrange(len(people)):
		database.append(Person(people[i]))
	return database