
# coding: utf-8

# In[1]:

import numpy as np
import time
import gmmhmm 
import MFCC
from scipy.io.wavfile import read as wread
import glob
import pickle

def sample_gmmhmm(gmmhmm, n_sim):
    """
    Simulate sampling from a GMMHMM.
    Returns
    -------
    states : ndarray of shape (n_sim,)
        The sequence of states
    obs : ndarray of shape (n_sim, K)
        The generated observations (column vectors of length K)
    """
    A, weights, means, covars, pi = gmmhmm
    states = []
    obs = []
    for i in xrange(n_sim):
        sample_c = np.argmax(np.random.multinomial(1, weights[1,:]))
        sample = np.random.multivariate_normal(means[1,sample_c,:],covars[1,sample_c,:,:])
        states.append(sample_c)
        obs.append(sample)
    states = np.array(states)
    obs = np.array(obs)
    return states, obs

def initialize(n):    
    A = np.array([[.65, .35], [.15, .85]])
    pi = np.array([.8, .2])
    weights = np.array([[.7, .2, .1], [.1, .5, .4]])
    means1 = np.array([[0., 17., -4.], [5., -12., -8.], [-16., 22., 2.]])
    means2 = np.array([[-5., 3., 23.], [-12., -2., 14.], [15., -32., 0.]])
    means = np.array([means1, means2])
    covars1 = np.array([5*np.eye(3), 7*np.eye(3), np.eye(3)])
    covars2 = np.array([10*np.eye(3), 3*np.eye(3), 4*np.eye(3)])
    covars = np.array([covars1, covars2])
    gmmhmm = [A, weights, means, covars, pi]
    
    pi = np.random.rand(n)
    pi = pi/np.sum(pi)
    A = np.random.rand(n,n)
    for i in xrange(n):
        A[i] = A[i]/np.sum(A[i])
    return pi, A
    
# Problem 2
def get_sampledicts():
    names = ['Mathematics', 'Biology', 'Political', 'Statistics', 'Psychology']
    sampledict = {}
    for name in names:
        sampledict[name] = []
        for fname in glob.glob("Samples/"+name+" *"):
            w = wread(fname)
            sampledict[name].append(MFCC.extract(w[1]))
            #print MFCC.extract(w[1]).shape, fname
    return names, sampledict

# Problem 3
def train(names, samples):
    best_models = []
    for name in names:
        print name
        best = -np.inf
        best_model = None
        for i in xrange(1):
            print 'Test number ',i,' of 10'
            startprob, transmat = initialize(5)
            model = gmmhmm.GMMHMM(n_components=5, n_mix=3, transmat=transmat, startprob=startprob, cvtype='diag')
            # these values for covars_prior and var should work well for this problem
            model.covars_prior = 0.01
            model.fit(samples[name][:10], init_params='mc', var=0.1)
            print "Training on ", name, "result: ", model.logprob, 'at time ', time.asctime(time.localtime())
            if model.logprob > best:
                best = model.logprob
                best_model = model
                #f = open(model_name(name),'w')
                #pickle.dump(model,f)
                #print "New best prob: ", best
                #f.close()
        best_models.append(best_model)
    #print best_models
    return best_models

def test(names, samples, models):
    results = []
    for name in names:
        for i in xrange(20,30):
            best = -np.inf
            best_model = -1
            for j in xrange(len(models)):
                score = models[j].score(i)
                if score > best:
                    best = score
                    best_model = j
            results.append(best_model)
            
    results = np.array(results)
    return results


# In[2]:

print 'Getting samples'
names, samples = get_sampledicts()
print 'Samples have been got'


# In[ ]:

models = train(names, samples)
pickle.dump(models,open('best_models.p','wb'))


# In[ ]:

results = test(names, samples, models)

print results