
# coding: utf-8

# In[5]:

import bjCommon
import numpy as np
from matplotlib import pyplot as plt
from scipy import rand
import time

# Problem 1
def lcg(n,a=1103515245,c=12345,mod=2**31-1,seed=4329):
    """
    Return an array of random numbers using the parameter specified
    """
    X = np.empty(n)
    X[0] = float((a*seed + c) % mod)
    for i in xrange(n-1):
        X[i+1] = float((a*X[i] + c) % mod)
    return X/mod

# Problem 2
def between(n,x1,x2,a=1103515245,c=12345,mod=2**31-1,seed=4329):
    """
    Scale the output of lcg to be between x1 and x2
    """
    orig = lcg(n, a, c, mod, seed)
    scaled = abs(x2-x1)*orig + x1
    return scaled.astype(int)

def problem3():
    """
    Graph the randoms from our lcg and also from scipy.rand
    """
    a = np.random.rand(512,512)
    b = lcg(512*512).reshape((512,512))
    plt.imshow(a)
    plt.show()
    plt.imshow(b)
    plt.show()

# Problem 4
def getSweepsEasy(games,n=10000001):
    """
    Return an n x games x 52 Array that represents the shuffles for <games> for
    the first n seeds
    
    In other words, each "row" of your array will be a separate call to
    bjCommon.shuffle with the games specified and a new seed (so seeds 0 to n-1)
    for each row. Return this as a numpy array of the integer representation
    of the cards (exactly what you get from bjCommon.shuffle)
    
    Also, make the default argument for n be the number of seeds you think need
    to be checked
    """
    A = np.empty((n,games,52))
    for i in xrange(n):
        A[i] = bjCommon.shuffle(n=games, a=2521, c=13, mod=65536, seed=i)
    return A

# Problem 5
def crackBlackJack(sweeps,cardTuples):
    """
    Return the sweeps whose first three cards from the nth game match the nth
    3-tuple in cardTuples, in their string representation
    
    For example, if cardTuples is a 2 x 3 list of tuples crackBlackJack should 
    return the sweeps in "sweeps" whose first 3 cards in the first game match
    cardTuples[0] and whose first 3 cards in the second game match cardTuples[1]
    
    For each tuple in cardTuples, you will eliminate the shuffles of seeds that
    do not match
    """
    pass

# Problem 6
def getSweepsHard(games,time,approx=120):
    """
    Return an n x games x 52 Array that represents the shuffles for <games> for
    the n possible seeds from the 2*approx second interval (above and below time)
    """
    pass


# In[7]:




# In[ ]:



