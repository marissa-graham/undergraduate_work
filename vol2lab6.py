'''

Solution specification for Lab 6

'''
import numpy as np
import math as m
from matplotlib import pyplot as plt
from collections import Counter
from sklearn import datasets as ds
from sklearn import neighbors


def scale_data(training_set, scaling_vector):
    '''
    This function accepts a training set and a scaling vector and then returns the scaled data.
    Training set will be a tuple with two arrays:
        training_set[0] is an nxm array of features
        training_set[1] is an nx1 array of labels

    Scaling vector will scale each row of training_set[0].

    Return a tuple of scaled data and labels.
    '''
    features = np.array(training_set[0], dtype=np.float_)
    for index, x in np.ndenumerate(features):
        features[index[0],index[1]] = float(x)/float(scaling_vector[index[1]])
    return (features, labels)


def calc_distance(training_set, scaling_vector, to_classify):
    '''
    Given a training set and scaling vector, use the previous method to scale the training_set.
    Once the training set is scaled, find the distance of each row to to_classify.
    
    Once again, training_set is a tuple of features(nxm) and labels(nx1).
    to_classify is an array of size nx1.

    Return an nx2 array of labels in the first column and distances in the second column.
    '''
    def euc_metrc(v1, v2):
        '''
        Calculate the euclidean distance between two vectors.
        '''
        v1 = np.array(v1)
        v2 = np.array(v2)
        d = np.abs(v1-v2)
        return m.sqrt(np.dot(d,d))

    labels = training_set[1]
    training_set = scale_data(training_set, scaling_vector)
    to_classify = np.array(to_classify)
    n = training_set[0].shape[0]

    index = n
    closest = [index, float("inf")]
    output = []
    for i in xrange(n):
        d = euc_metrc(to_classify, training_set[i])
        if d < closest[1]:
            closest[0] = i
            closest[1] = d
        output.append([labels,d])

    return output

'''
Problem 3:  Describe the iris dataset in a few sentence

The iris dataset is stored in a 150x4 numpy array, with the rows being
the samples and the columns being various attributes of 
each sample. These are sepal length, sepal width, petal length, and 
petal width. And these are flowers, not the colored part of eyeballs.

'''

def average_integer(postal_data, integer):
    '''
    postal_data will be a tuple (points, labels) from the given postal data.
    Find every point in the postal data that is labelled with the given integer and display the average.
    '''
    labels = postal_data[1]
    given = points[1][labels[1]==integer]
    sum_given = np.sum(given, axis=0)
    avg_given = (sum_given/len(given)).reshape((28,28))
    plt.gray()
    plt.imshow(avg_given)
    plt.show()

def build_classifier(postal_data, n_neighbors, w, metric=2):
    '''
    Build and return a KNearestNeighbor classifier with sklearn
    '''
    neighbrs = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors, weights=w, p=metric)
    neighbrs.fit(postal_data[0], postal_data[1])
    return neighbrs

def classify_postal_data(classifier, test_points, test_labels):
    '''
    Classify the test_points using the classifier you built.
    test_labels contains the true labels for the test points.
    Return your error rate, i.e. num_correct/len(test_points)
    '''
    test = classifier.predict(test_points[1])
    errors = test - test_labels[1]
    wrong = np.nonzero(errors)[0].size
    total = test_labels[1].size
    return float(total-wrong)/float(total)

labels, points, testlabels, testpoints = np.load('PostalData.npz').items()
postal_data = (points[1], labels[1])

print 'Building classifiers...'
classifier1 = build_classifier(postal_data, 4, 'uniform')
classifier2 = build_classifier(postal_data, 4, 'distance')
classifier3 = build_classifier(postal_data, 10, 'uniform')
classifier4 = build_classifier(postal_data, 10, 'distance')
classifier5 = build_classifier(postal_data, 1, 'uniform')
classifier6 = build_classifier(postal_data, 1, 'distance')

print 'Calculating errors...'
error1 = classify_postal_data(classifier1, testpoints, testlabels)
print 'error1 complete'
error2 = classify_postal_data(classifier2, testpoints, testlabels)
print 'error2 complete'
error3 = classify_postal_data(classifier3, testpoints, testlabels)
print 'error3 complete'
error4 = classify_postal_data(classifier4, testpoints, testlabels)
print 'error4 complete'
error5 = classify_postal_data(classifier5, testpoints, testlabels)
print 'error5 complete'
error6 = classify_postal_data(classifier6, testpoints, testlabels)
print 'error6 complete\n'

print 'Error for 4 neighbors, uniform weights: ' + str(error1)
print 'Error for 4 neighbors, distance weights: ' + str(error2)
print 'Error for 10 neighbors, uniform weights: ' + str(error3)
print 'Error for 10 neighbors, distance weights: ' + str(error4)
print 'Error for 1 neighbor, uniform weights: ' + str(error5)
print 'Error for 1 neighbor, distance weights: ' + str(error6)  
