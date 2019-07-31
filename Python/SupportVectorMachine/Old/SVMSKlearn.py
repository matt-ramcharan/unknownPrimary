from sklearn import svm
from sklearn import metrics
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def gaussianKernelGramMatrixFull(X1, X2, sigma=1):
    """(Pre)calculates Gram Matrix K"""

    gram_matrix = np.zeros((X1.shape[0], X2.shape[0]))
    for i, x1 in enumerate(X1):
        for j, x2 in enumerate(X2):
            #x1 = x1.flatten()
            #x2 = x2.flatten()
            gram_matrix[i, j] = np.exp(- np.sum( np.power((x1 - x2),2) ) / float( 2*(sigma**2) ) )
    return gram_matrix

train_data = pd.read_csv('sonar.tr',sep='\s+', header = None, names=list(range(1,62)))
train_input = train_data.loc[:,:train_data.shape[1]-1]
train_labels = train_data.loc[:,train_data.shape[1]]

test_data = pd.read_csv('sonar.ts',sep='\s+', header = None, names=list(range(1,62)))
test_input = test_data.loc[:,:test_data.shape[1]-1]
test_labels = test_data.loc[:,test_data.shape[1]]

kernel = 'gaussian'

if kernel=='linear':
    #Create a svm Classifier

    #Linear
    clf = svm.SVC(kernel='linear')

    # Train the model using the training sets

    clf.fit(train_input, train_labels)

    # Predict the response for test dataset
    y_pred = clf.predict(test_input)

    # Model Accuracy: how often is the classifier correct?
    print("Accuracy:", metrics.accuracy_score(test_labels, y_pred))

    # Model Precision: what percentage of positive tuples are labeled as such?
    print("Precision:", metrics.precision_score(test_labels, y_pred))

    # Model Recall: what percentage of positive tuples are labelled as such?
    print("Recall:", metrics.recall_score(test_labels, y_pred))

elif kernel=='gaussian':
    # Create a svm Classifier

    #Gaussian
    C=0.1
    clf = svm.SVC(C = C, kernel="precomputed")

    #Train the model using the training sets

    model = clf.fit(gaussianKernelGramMatrixFull(train_input, train_input), train_labels)

    #Predict the response for test dataset
    y_pred = model.predict(gaussianKernelGramMatrixFull(test_input, train_input))

    # Model Accuracy: how often is the classifier correct?
    print("Accuracy:",metrics.accuracy_score(test_labels, y_pred))


    # Model Precision: what percentage of positive tuples are labeled as such?
    print("Precision:",metrics.precision_score(test_labels, y_pred))

    # Model Recall: what percentage of positive tuples are labelled as such?
    print("Recall:",metrics.recall_score(test_labels, y_pred))