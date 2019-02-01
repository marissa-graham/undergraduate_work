# Problem 1

def summy(n):
	pad = 16 - n%15
	shape = n+pad
	to_sum = np.append(np.arange(1,n,1), np.zeros(pad)).reshape(shape/3,3)
	partial_sum = np.sum(to_sum[:,2])
	to_sum = to_sum.reshape(shape/5,5)
	return partial_sum + np.sum(to_sum[:,4]) - np.sum(to_sum[2::3,4])

# Problem 2

def fibby(n):
	"""Find the sum of all Fibonacci numbers less than n."""
	tot = 2
	f1 = 1
	f2 = 2
	nextone = 3
	while nextone <= n:
		if nextone % 2 == 0:
			tot += nextone
		f2 = nextone
		f1 = f2 - f1
	return tot

# Problem 3

def prime_factor(n, check=1000):
	"""Find the largest factor of n. Check is a crude way to
	make large numbers manageable.
	"""
	
	if n > 1000000000:
		for i in xrange(2, check):
			while n % i == 0:
				n /= i
	print n
	def find_factors(a):
		factors = np.empty(0)
		for i in xrange(int(a)-1, 1, -1):
			if a % i == 0:
				factors = np.append(factors, i)
		if len(factors) > 0:
			return factors[0]
		else:
			return None
	has_factor = find_factors(n)
	biggest = has_factor
	while has_factor:
		biggest = has_factor
		has_factor = find_factors(has_factor)
	if biggest:
		return int(biggest)
	else:
		return n

# Problem 4

def palindromes():
	"""Find the largest palindrome that's a product of three digit
	numbers.
	"""
	palindromes = np.zeros(0, dtype='int32')
	for i in xrange(1000):
		for j in xrange(1000):
			if str(i*j) == str(i*j)[::-1]:
				palindromes = np.append(palindromes, i*j)
	return np.max(palindromes)

# Problem 5
def cheater():
	"""Find the smallest positive number evenly divisible by all numbers
	from 1 to 20, but the lazy way.
	"""
	return 19*17*13*11*7*5*3*2*2*2*3*2

# Problem 6

def problem6(n):
    """Difference between sum squared and sums squared."""
    return np.sum(np.arange(1,n+1))**2 - np.sum(np.arange(1,n+1)**2)


# Problem 7

def problem7(n):
    """Return the nth prime number."""
    primes = [2]
    i = 3
    while len(primes) < n:
        is_prime = True
        for j in primes:
            if i % primes[j] == 0:
                is_prime = False
        if is_prime:
            primes.append(i)
    return primes[-1]

# Problem 8

def problem8():
    """Find the thirteen adjacent digits in a 1000-digit number with
    the greatest product.
    """
    things = []
    with open('bignumber.txt') as number:
        for word in number:
            things.append(word.strip())
    digits = np.array([], dtype='int64')
    for i in xrange(len(things)):
        for j in xrange(0,len(things[i])):
            digits = np.append(digits, int(things[i][j]))
    products = np.array([0], dtype='int64')
    for i in xrange(0,len(digits)):
        products = np.append(products, np.prod(digits[i:i+13]))
        return np.max(products)


# Problem 9

def problem9():
    foundit = False
    for i in xrange(1,998):
        a = i
        for j in xrange(1,998):
            b = j
            c = 1000 - a - b
            if a**2 + b**2 == c**2:
                foundit = True
                return a*b*c
    return foundit


# Problem 10

def problem10(n):
    """Find the sum of all the primes less than n."""
    tot = 2
    primes = np.array([2],dtype='int64')
    for i in xrange(3,n,2):
        prime = True
        sqrt = np.sqrt(i)
        j = 0
        while primes[j] <= int(sqrt) and prime:
            if i % primes[j] == 0:
                prime = False
            j += 1
        if prime:
            primes = np.append(primes,i)
            tot += i
        print '\r' + str(i)
    return tot

# Problem 11

def problem11():
    """Find the largest product of four consecutive numbers in a grid,
    diagonals, up, down, etc.
    """
    data = []
    with open('grid.csv') as grid:
        gridreader = csv.reader(grid, delimiter=' ')
        for row in gridreader:
            data.append(row)
    for i in xrange(len(data)):
        for j in xrange(len(data[i])):
            data[i][j] = int(data[i][j])
    data = np.array(data)
    prod = 0
    for i in xrange(17):
        for j in xrange(17):
            prod1 = np.prod(data[i:i+3,j])
            if prod1 >= prod:
                prod = prod1
            prod2 = np.prod(data[i,j:j+3])
            if prod2 >= prod:
                prod = prod2
            prod3 = data[i,j]*data[i+1,j+1]*data[i+2,j+2]*data[i+3,j+3]
            if prod3 >= prod:
                prod = prod3
    for i in xrange(19,3,-1):
        for j in xrange(17):
            prod1 = data[i,j]*data[i-1,j+1]*data[i-2,j+2]*data[i-3,j+3]
            if prod1 >= prod:
                prod = prod1
    return prod

# Problem 12

def problem12(n):
    """Find the first triangle number with n divisors. Naive method."""
    tri_prev = 1
    i = 1
    max_divisors = [(1,1,1)]
    biggest = 1
    while biggest <= n:
        i += 1
        tri_next = tri_prev + i
        divisors = 2
        j = 2
        to_check = n/2
        while j < to_check:
            if tri_next % j == 0:
                divisors += 2
                to_check = n/j
            j += 1
        if divisors > biggest:
            biggest = divisors
            max_divisors.append((i, tri_next, divisors))
        print '\\r' + str(tri_next) + ', the ' + str(i) +  'th triangle number, has ' + str(divisors) + ' divisors. Max divisors: ' + str(biggest),        
        tri_prev = tri_next
    return max_divisors

def alt12(n):
    max_divisors = 2
    i = 2
    while max_divisors < n:
        # Find prime factors of n
        pass

# 76576500, the 12375th triangle number, has 576 divisors. Max divisors: 576"

# Problem 13

def problem13():
    """Find the first ten digits of a long list of numbers."""
    to_sum = np.zeros(0, dtype='int64')
    with open('numbers.txt') as numbers:
        for line in numbers:
            to_sum = np.append(to_sum, int(line.strip()))
    tot = np.sum(to_sum)
    return tot


# Problem 14

def problem14(n):
    """Find the starting number <n that produces the longest
    Collatz sequence.
    """
    longest = 1
    starter = 1
    for i in xrange(1,n):
        if i % 10000 == 0:
            print i
        k = i
        chain = [k]
        
        while k != 1:
            if k % 2 == 0:
                k /= 2
            else:
                k = 3*k + 1
            chain.append(k)
        
        length = len(chain)
        if length > longest:
            starter = i
            longest = length
    return longest, starter

# Problem 15
import numpy as np
import scipy.sparse as spar

def problem15(n):
    """Find the number of paths to the bottom right corner
    in a n by n grid, moving right and down."""
    m = (n+1)**2
    
    superdiag = np.diag(np.ones(n),k=1)
    
    A = np.diag(np.ones(n*(n+1)),k=n+1)
    
    for i in xrange(n+1):
        A[i*(n+1):(i+1)*(n+1),i*(n+1):(i+1)*(n+1)] = superdiag
        
    A = spar.csr_matrix(A)
    A = A**2
    return A**n
    
# Problem 16

def problem16():
    """Find the digital root of 2^1000."""
    number = str(2**1000)
    tot = 0
    for i in xrange(len(number)):
        tot += int(number[i])
    return tot

# Problem 17

def problem17():
    """Find the number of letters it would take to write out
    all the numbers from 1 to 1000 inclusive in words."""
    
    digits = [0,3,3,5,4,4,3,5,5,4]
    teens = [3,6,6,8,8,7,7,9,8,8]
    tens = [0,3,6,6,5,5,5,7,6,6]
    
    tot = 0
    for i in xrange(1,1001):
        old = tot
        if len(str(i)) == 1:
            tot += digits[i]
        
        elif len(str(i)) == 2:
            if str(i).startswith('1'):
                tot += teens[i-10]
                
            else:
                tot += tens[i/10]
                tot += digits[i%10]
                
        elif len(str(i)) == 3:
            tot += digits[i/100]
            tot += 7
            j = i % 100
            if str(i).endswith('00'):
                tot += 0
            else:
                tot += 3
                
                if len(str(j)) == 1:
                    tot += digits[j]
                
                elif len(str(j)) == 2:
                    if str(j).startswith('1'):
                        tot += teens[j-10]

                    else:
                        tot += tens[j/10]
                        tot += digits[j%10]
        
        if len(str(i)) == 4:
            tot += 11   
            
    return tot

# Problem 18
# Code also gives problem 67

def problem18(filename):
    """Find the maximum total from top to bottom of a triangle given
    in filename.
    """
    fromtxt = np.array([])
    numrows = 0
    with open(filename) as myfile:
        for line in myfile:
            numrows += 1
            for i in line.split():
                fromtxt = np.append(fromtxt,int(i))

    A = np.zeros((numrows,numrows))
    for i in xrange(numrows):
        A[i,0:i+1] = fromtxt[i*(i+1)/2:(i+1)*(i+2)/2]

    for i in xrange(numrows):
        for j in xrange(numrows-i):
            if j == numrows-i-1:
                pass
            elif j == 0:
                A[numrows-i-2,j] += np.max(A[numrows-i-1,j:j+2])
            else:
                A[numrows-i-2,j] += np.max(A[numrows-i-1,j:j+2])

    print A[0,0]

# Problem 19

def problem19():
    """Find the number of Sundays that fell on the first of the month
    during the twentieth century.
    """
    Sundaycounter = 0
    
    daytracker = 0
    weektracker = 1
    monthcounter = 0
    monthtracker = 1
    yearcounter = 1901
    
    monthlengths = [31,28,31,30,31,30,31,31,30,31,30,31]
    
    while yearcounter <= 2000:
        
        weektracker += 1
        daytracker += 1
        monthtracker += 1
        
        if weektracker == 7:
            weektracker = 0
            
        # Leap year conditions
        if yearcounter % 4 == 0:
            if yearcounter % 100 == 0 and yearcounter % 400 != 0:
                if daytracker == 365:
                    daytracker = 0
                    yearcounter += 1

                if monthtracker == monthlengths[monthcounter]:
                    monthtracker = 1
                    monthcounter += 1
                    if monthcounter == 12:
                            monthcounter = 0

                if monthtracker == 1 and weektracker == 0:
                    Sundaycounter += 1
            else:
                if daytracker == 366:
                    daytracker = 0
                    yearcounter += 1
                
                if monthcounter == 1:
                    if monthtracker == monthlengths[monthcounter+1]:
                        monthtracker = 1
                        monthcounter += 1
                        if monthcounter == 12:
                            monthcounter = 0
                    
                if monthcounter != 1:
                    if monthtracker == monthlengths[monthcounter]:
                        monthtracker = 1
                        monthcounter += 1
                        if monthcounter == 12:
                            monthcounter = 0

                if monthtracker == 1 and weektracker == 0:
                    Sundaycounter += 1
                
        else:
            if daytracker == 365:
                daytracker = 0
                yearcounter += 1

            if monthtracker == monthlengths[monthcounter]:
                monthtracker = 1
                monthcounter += 1
                if monthcounter == 12:
                    monthcounter == 0

            if monthtracker == 1 and weektracker == 0:
                Sundaycounter += 1
                    
        return Sundaycounter


# Problem 20

def problem20():
	prod = 1
	for i in xrange(1,101):
	    prod *= i
	print prod

	tot = 0
	for i in xrange(len(str(prod))):
	    tot += int(str(prod)[i])
	print tot

# Problem 21
import numpy as np

def amicables(k):
    """Find the sum of all the amicable numbers less than k."""
    def d(n):
        tot = 1
        for i in xrange(2,n/2+1):
            if n % i == 0:
                tot += i
        return tot

    tot = 0
    d_stuff = np.array([d(i) for i in xrange(k)])
    d_stuff[0] = 0
    d_stuff[1] = 0
    is_amicable = np.array([False for i in xrange(k)])

    for i in xrange(2,k):
        print '\r' + str(i),
        if ~is_amicable[i]:
            for j in xrange(2,k):
                if i != j:
                    if d_stuff[i] == j and d_stuff[j] == i:
                        is_amicable[i] = True
                        is_amicable[j] = True

    for i in xrange(10000):
        if is_amicable[i]:
            tot += i

    print tot


# Problem 22
import csv

def namescore():
    """Get total of name scores in given file."""
    names = []
    with open('projectEuler22.txt') as data:
        namereader = csv.reader(data, delimiter=',')
        for row in namereader:
            names.append(row)

    names = names[0]
    names = sorted(names)

    def score(name):
        score = 0
        for i in xrange(len(name)):
            score += ord(name[i]) - 64
        return score

    namescore = 0
    for i in xrange(len(names)):
        print '\r' + str(i),
        namescore += (i+1)*score(names[i])

    return namescore

# Problem 23
import numpy as np

def problem23():
    def d(n):
        tot = 1
        for i in xrange(2,n/2+1):
            if n % i == 0:
                tot += i
        return tot

    tot = 0
    abundant = []
    for i in xrange(3,28124):
        print '\r' + str(i),
        if d(i) > i:
            abundant.append(i)

    print '\n'
    works = np.array([False for i in xrange(28124)])
    for i in abundant:
        print '\r' + str(i),
        for j in abundant:
            if i+j < 28124:
                works[i+j] = True

    for i in xrange(28124):
        if ~works[i]:
            tot += i

    print tot

# Problem 24
import itertools

def problem24():
    perms = []
    for i in itertools.permutations([0,1,2,3,4,5,6,7,8,9],10):
        perms.append(i)

    perms = sorted(perms)
    print perms[999999]

# Problem 25

def problem25():
    maxlength = 0
    fold = 1
    fnow = 1
    counter = 2
    while maxlength < 1000:
        counter += 1
        fnext = fold + fnow
        fold = fnow
        fnow = fnext
        maxlength = len(str(fnext))

    return counter


# Problem 26
import decimal

def problem26():
    best = 6
    best_num = 7
    for i in xrange(1,1000):
        print '\r' + str(i),
        decimal.getcontext().prec = 10000
        frac = decimal.Decimal(1)/decimal.Decimal(i)
        thing = str(frac)
        if len(thing) > 100:
            longest = thing[10]
            done = False
            j = 0
            while done == False and j < 1000:
                if longest == thing[11+j:11+j+len(longest)]:
                    done = True
                else:
                    longest += thing[11+j]
                j += 1
            if len(longest) > best:
                best_num = i
                best = len(longest)

    return best_num, best     

# Problem 27

def problem27():
    def is_prime(n):
        is_prime = True
        to_append = True
        if n < 2:
            is_prime = False 
        i = 2
        while is_prime and i <= np.sqrt(n)+1:
            if n % i == 0:
                is_prime = False
            i += 1    
        return is_prime

    maximum_primes = 0
    coeff_product = 0

    for a in xrange(-999,1000):
        print '\r', a, maximum_primes, coeff_product,
        for b in xrange(-999,1000):
            n = 0
            prime = is_prime(abs(b))
            while prime:
                n += 1
                prime = is_prime(abs(n**2+a*n+b))
            if n+1 > maximum_primes:
                maximum_primes = n+1
                coeff_product = a*b

    print '\n', maximum_primes, coeff_product
    return coeff_product

# Problem 28
import numpy as np

def problem28():
    numbers = np.arange(1001**2+1)
    tot = 1
    for i in xrange(1,501):
        newrow = numbers[4*i**2-4*i+2:4*i**2+4*i+2]
        for j in xrange(1,5):
            tot += newrow[2*i*j-1]
    print tot


# Problem 29

def problem29():
    combinations = set()
    for a in xrange(2,101):
        print '\r',a,
        for b in xrange(2,101):
            combinations = combinations.union(set([a**b]))
    return len(combinations)


# Problem 30

def problem30():
    tot = 0
    for i in xrange(2,1000000):
        stringsum = 0
        for j in xrange(len(str(i))):
            stringsum += int(str(i)[j])**5
        if i == stringsum:
            tot += i
    return tot

# Problem 31

def naive():
    numways = 1
    left = 200
    for i in xrange(3):
        for j in xrange((200-100*i)//50+2):
            for k in xrange((200-100*i-50*j)//20+2):
                for l in xrange((200-100*i-50*j-20*k)//10+2):
                    for m in xrange((200-100*i-50*j-20*k-10*l)//5+2):
                        for n in xrange((200-100*i-50*j-20*k-10*l-5*m)//2+2):
                            for p in xrange(200-100*i-50*j-20*k-10*l-5*m-2*n+2):
                                howmuch = 100*i+50*j+20*k+10*l+5*m+2*n+p
                                print '\r',i,j,k,l,
                                if howmuch == 200:
                                    numways += 1
    return numways

# Problem 32
import itertools

def problem32():
    pandigital = set()
    for i in itertools.permutations([1,2,3,4,5,6,7,8,9],9):
        prod1 = int(str(i[0])+str(i[1]))
        prod2 = int(str(i[2])+str(i[3])+str(i[4]))
        prod3 = int(str(i[5])+str(i[6])+str(i[7])+str(i[8]))

        also = i[0]
        also2 = int(str(i[1])+str(i[2])+str(i[3])+str(i[4]))
        if prod1*prod2 == prod3:
            pandigital.add(prod3)
            print '\nYay!',i,prod1,prod2,prod3,prod1*prod2,'\n'
        elif also*also2 == prod3:
            pandigital.add(prod3)
            print '\nYay!',i,also,also2,prod3,also*also2,'\n'
        print '\r',i,prod1,prod2,prod3,prod1*prod2,

    tot = 0
    for x in pandigital:
        tot += x
    return tot

# Problem 33

def problem33():
    weirdos = []
    for i in xrange(10,100):
        for j in xrange(10,100):
            i1 = int(str(i)[0])
            i2 = int(str(i)[1])
            j1 = int(str(j)[0])
            j2 = int(str(j)[1])
            if i < j and j2 != 0 and j1 != 0:
                if i1 == j1:
                    if i/float(j) == float(i2)/j2:
                        weirdos.append((i,j))
                if i1 == j2:
                    if i/float(j) == float(i2)/j1:
                        weirdos.append((i,j))
                if i2 == j1:
                    if i/float(j) == float(i1)/j2:
                        weirdos.append((i,j))
                if i2 == j2:
                    if i/float(j) == float(i1)/j1:
                        weirdos.append((i,j))
    return weirdos

# Bit of wolfram alpha to get actual answer

# Problem 34
import math

def check(n):
    factorials = [math.factorial(i) for i in xrange(10)]
    works = 0
    # Don't count 1 and 2
    tot = -3
    for i in xrange(n):
        if i % 1000 == 0:
            print '\r',n-i,
        factsum = 0
        for j in xrange(len(str(i))):
            factsum += factorials[int(str(i)[j])]
        if factsum == i:
            works += 1
            tot += i

    print '\n', works, tot
    
#check(100000)
# 40730

# Problem 35

import numpy as np

def problem35(n):
    """Find all circular primes below n."""

    def lookuptable(n):
        primebools = np.array([True for i in xrange(n)])
        primebools[0] = False
        primebools[1] = False
        done = int(np.floor(np.sqrt(n)))+1
        for i in xrange(2,done):
            primebools[2*i::i] = False
        primes = np.array(np.where(primebools == True))[0]
        primes = set(primes.flat)
        return primes
    
    primes = lookuptable(n)
    circulars = []
    for x in primes:
        circular = True
        for i in xrange(1,len(str(x))):
            z = str(x)
            rotated = int(z[i:]+z[:i])
            if rotated not in primes:
                circular = False
        if circular:
            circulars.append(x)
    return len(circulars)       
   
def lookuptable(n):
    primebools = np.array([True for i in xrange(n)])
    primebools[0] = False
    primebools[1] = False
    done = int(np.floor(np.sqrt(n)))+1
    for i in xrange(2,done):
        primebools[2*i::i] = False
    primes = np.array(np.where(primebools == True))[0]
    primes = set(primes.flat)
    return primes

# Problem 36

def problem36():
    n = 1000000
    tot = 0
    for i in xrange(n):
        if i % 1000 == 0:
            print '\r',n-i,
        str_i = str(i)
        if str_i == str_i[::-1]:
            bin_i = '{0:b}'.format(i)
            if bin_i == bin_i[::-1]:
                tot += i
    print tot

# Problem 37
import numpy as np

def problem37(n):
    
    def lookuptable(n):
        primebools = np.array([True for i in xrange(n)])
        primebools[0] = False
        primebools[1] = False
        done = int(np.floor(np.sqrt(n)))+1
        for i in xrange(2,done):
            primebools[2*i::i] = False
        primes = np.array(np.where(primebools == True))[0]
        primes = set(primes.flat)
        return primes

    primes = lookuptable(n)
    truncatables = lookuptable(n)
    
    for x in primes:
        str_x = str(x)
        truncatable = True
        for i in xrange(1,len(str_x)):
            if int(str_x[i:]) not in primes:
                truncatable = False
            if int(str_x[:-i]) not in primes:
                truncatable = False
        if truncatable == False:
            truncatables.remove(x)
            
    truncatables.remove(2)
    truncatables.remove(3)
    truncatables.remove(5)
    truncatables.remove(7)
    
    return len(truncatables),sum(truncatables),truncatables   

# Problem 38
import itertools

def problem38():
    largest = 0
    for i in itertools.permutations([1,2,3,4,5,6,7,8,9]):
        print '\r',i,
        int_i = reduce(lambda a,x:a*10+x,i)
        str_i = str(int_i)
        for j in xrange(1,5):
            k = j
            n = 2
            works = True
            while k < 9 and works:
                check1 = n*int(str_i[:j])
                length = len(str(check1))
                check2 = int(str_i[k:k+length])
                if check1 == check2:
                    n += 1
                    k += length
                else:
                    works = False
            if 9-k == 0 and works and int_i > largest:
                largest = int_i

    return largest     

# Problem 39

def problem39():
    max_sols = 0
    max_p = 0
    for p in xrange(12,1001):
        print '\r',p,
        num_sols = 0
        for a in xrange(3,p):
            for b in xrange(4,p-a):
                if a**2+b**2 == (p-a-b)**2:
                    num_sols += 1
        if num_sols > max_sols:
            max_p = p
            max_sols = num_sols

    return max_sols, max_p

# Problem 40

def problem40(n):
    def make_champer(n):
        champer = []
        for i in xrange(n):
            for j in xrange(len(str(i))):
                champer.append(int(str(i)[j]))
        return champer
    
    d = make_champer(200000)
    print d[1]*d[10]*d[100]*d[1000]*d[10000]*d[100000]*d[1000000]


# Problem 41

def problem41():
    
    def lookuptable(n):
        primebools = np.array([True for i in xrange(n)])
        primebools[0] = False
        primebools[1] = False
        done = int(np.floor(np.sqrt(n)))+1
        for i in xrange(2,done):
            print '\r',done-i,
            primebools[2*i::i] = False
        primes = np.array(np.where(primebools == True))[0]
        primes = set(primes.flat)
        return primes

    primes = lookuptable(7654322)
    print len(primes)
    largest = 0
    for i in itertools.permutations([1,2,3,4,5,6,7],7):
        int_i = reduce(lambda a,x:a*10+x,i)
        if int_i in primes:
            if int_i > largest:
                largest = int_i
    print largest

# Problem 42

def problem42():
    words = []
    with open('projectEuler42.txt') as data:
        wordreader = csv.reader(data, delimiter=',')
        for row in wordreader:
            words.append(row)
    words = words[0]

    longword = 0
    for i in xrange(len(words)):
        if len(words[i]) > longword:
            longword = len(words[i])

    triangles =set([n*(n+1)/2 for n in xrange(28)])

    def score(word):
        score = 0
        for i in xrange(len(word)):
            score += ord(word[i]) - 64
        return score

    num_triangles = 0
    for i in xrange(len(words)):
        if score(words[i]) in triangles:
            num_triangles += 1

    print num_triangles

# Problem 43

def problem43():
    primes = [1,2,3,5,7,11,13,17]
    tot = 0
    for i in itertools.permutations([j for j in xrange(10)],10):
        works = True
        print '\r',i,
        for j in xrange(1,8):
            if int(str(i[j])+str(i[j+1])+str(i[j+2])) % primes[j] != 0:
                works = False
        if works:
            tot += reduce(lambda a,x:a*10+x,i)
    print tot

# Problem 44

def problem44():
    k = 10000
    pentagonals = set([n*(3*n-1)/2 for n in xrange(1,2*k)])
    works = []
    p = lambda n:n*(3*n-1)/2
    for i in xrange(1,k):
        print '\r',i,
        pi = p(i)
        for j in xrange(1,i):
            pj = p(j)
            if pi+pj in pentagonals and pi-pj in pentagonals:
                works.append(abs(pi-pj))

    print works
    print min(works)

# Problem 45

def problem45():
    k = 56000
    works = None
    triangles = set([n*(n+1)/2 for n in xrange(286,k)])
    pentagons = set([n*(3*n-1)/2 for n in xrange(166,k)])
    hexagons = set([n*(2*n-1) for n in xrange(144,k)])
    for x in hexagons:
        print '\r',x,
        if x in triangles and x in pentagons and works == None:
            works = x
    print works


def problem46():
    def lookuptable(n):
        primebools = np.array([True for i in xrange(n)])
        primebools[0] = False
        primebools[1] = False
        done = int(np.floor(np.sqrt(n)))+1
        for i in xrange(2,done):
            primebools[2*i::i] = False
        primes = np.array(np.where(primebools == True))[0]
        primes = set(primes.flat)
        return primes

    k = 6000
    primes = lookuptable(k)
    false = []
    for i in xrange(3,k,2):
        print '\r',i,
        works = False
        if i not in primes:
            for x in primes:
                if x < i:
                    sqrt = np.sqrt((i-x)/2)
                    if sqrt - np.floor(sqrt) == 0:
                        works = True
        if works != True and i not in primes:
            false.append(i)
            print '\n',i
    print false                 

# Problem 47

def problem47():
    def lookuptable(n):
        primebools = np.array([True for i in xrange(n)])
        primebools[0] = False
        primebools[1] = False
        done = int(np.floor(np.sqrt(n)))+1
        for i in xrange(2,done):
            primebools[2*i::i] = False
        primes = np.array(np.where(primebools == True))[0]
        return primes

    def factor_table(n,primes,setprimes):
        """Return the a list with the number of factors of
        each number up to n.
        """
        table = np.array([0,0])
        for i in xrange(2,n):
            if i in setprimes:
                table = np.append(table,1)
            else:
                primedivisors = 0
                for j in xrange(len(primes[primes < n])):
                    x = primes[j]
                    divides = False
                    while i % x == 0:
                        i /= x
                        divides = True
                    if divides:
                        primedivisors += 1
                table = np.append(table,primedivisors)

        return table

    def num_divisors(n,primes,setprimes):
        """Return the number of divisors of n."""
        if n in setprimes:
            return 1
        else:
            primedivisors = 0
            i = 0
            x = primes[0]
            while x <= n and primedivisors < 6:
                divides = False
                while n % x == 0:
                    n /= x
                    divides = True
                if divides:
                    primedivisors += 1
                i += 1
                x = primes[i]
            return primedivisors


    k = 200000
    print 'Creating lookup table...'
    primes = lookuptable(k)
    setprimes = set(primes)

    works = None
    i = 0
    while i < k-5 and works == None:
        print '\r',i,
        #if num_divisors(i,primes,setprimes) == 3:
        #    if num_divisors(i+1,primes,setprimes) == 3:
        #        if num_divisors(i+2,primes,setprimes) == 3:
        #            works = i
        if num_divisors(i,primes,setprimes) == 4:
            if num_divisors(i+1,primes,setprimes) == 4:
                if num_divisors(i+2,primes,setprimes) == 4:
                    if num_divisors(i+3,primes,setprimes) == 4:
                        works = i
        i += 1
    print '\n',works

# Problem 48

def problem48():
    print sum([i**i for i in xrange(1,1001)])