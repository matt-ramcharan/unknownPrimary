from sklearn import svm, metrics
from sklearn.utils import shuffle
from sklearn.metrics import classification_report, confusion_matrix
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

train_data = np.load('/home/matt/Documents/TechnicalProject/DeepLearningGenomicMed/Python/MultitaskLearn/train_alt.npy')
test_data = np.load('/home/matt/Documents/TechnicalProject/DeepLearningGenomicMed/Python/MultitaskLearn/test_alt.npy')


train_input = train_data[1, :, :-1]
train_labels = train_data[1, :, -1]

test_input = test_data[1, :, :-1]
test_labels = test_data[1, :, -1]

test_input, test_labels = shuffle(test_input, test_labels, random_state=0)
train_input, train_labels = shuffle(train_input, train_labels, random_state=0)

# train_input, train_labels = shuffle(test_data[1, :, :-1], test_data[1, :, -1], random_state=0)
# test_input, test_labels = shuffle(train_data[1, :, :-1], train_data[1, :, -1], random_state=0)

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

# elif kernel=='gaussian':
#     # Create a svm Classifier
#
#     #Gaussian
#     C=0.1
#     clf = svm.SVC(C = C, kernel="precomputed",verbose=1)
#
#     #Train the model using the training sets
#
#     model = clf.fit(gaussianKernelGramMatrixFull(train_input, train_input), train_labels)
#
#     #Predict the response for test dataset
#     y_pred = model.predict(gaussianKernelGramMatrixFull(test_input, train_input))
#
#     # Model Accuracy: how often is the classifier correct?
#     print("Accuracy:",metrics.accuracy_score(test_labels, y_pred))
#
#
#     # Model Precision: what percentage of positive tuples are labeled as such?
#     print("Precision:",metrics.precision_score(test_labels, y_pred))
#
#     # Model Recall: what percentage of positive tuples are labelled as such?
#     print("Recall:",metrics.recall_score(test_labels, y_pred))

elif kernel == 'gaussian':
    # Create a svm Classifier
    svclassifier = svm.SVC(kernel='rbf')
    svclassifier.fit(train_input, train_labels)
    y_pred = svclassifier.predict(test_input)

    print(confusion_matrix(test_labels, y_pred))
    print(classification_report(test_labels, y_pred))
