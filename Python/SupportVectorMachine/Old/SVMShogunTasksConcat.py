# import numpy as np
# import shogun
# from shogun.Features import LibSVM
# from shogun.Kernel import *
# from shogun.Classifier import AccuracyMeasure
# from shogun.Evaluation import RealFeatures, BinaryLabels, AccuracyMeasure
import pandas as pd
from shogun.Loss import L2R_L2LOSS_SVC
import matplotlib.pyplot as plt
from shogun.Evaluation import RealFeatures, BinaryLabels, LibSVM, AccuracyMeasure, ROCEvaluation
import numpy as np
from shogun.Kernel import GaussianKernel, LibLinear, LinearKernel

#MultiTask Dataset (reduced)
# train_data = np.load('/home/matt/Documents/TechnicalProject/DeepLearningGenomicMed/Python/MultitaskLearn/train_alt.npy')
# test_data = np.load('/home/matt/Documents/TechnicalProject/DeepLearningGenomicMed/Python/MultitaskLearn/test_alt.npy')

#Ratsh Dataset
train_data = np.load('/home/matt/Documents/TechnicalProject/DeepLearningGenomicMed/Python/MultitaskLearn/train.npy')
test_data = np.load('/home/matt/Documents/TechnicalProject/DeepLearningGenomicMed/Python/MultitaskLearn/test.npy')

np.apply_along_axis(np.random.shuffle, 1, train_data)
np.apply_along_axis(np.random.shuffle, 1, test_data)

# #Tasks concatenated
test_data = np.reshape(test_data, [64000, 101])
train_data = np.reshape(train_data, [640, 101])

train_feats = train_data[:, :-1]
test_feats = test_data[:, :-1]
train_label = train_data[:, -1]
test_label = test_data[:, -1]

# #Test Train Swap
# train_feats = test_data[:, :-1]
# test_feats = train_data[:, :-1]
# train_label = test_data[:, -1]
# test_label = train_data[:, -1]

features_train = RealFeatures(train_feats.T)
features_test = RealFeatures(test_feats.T)
labels_train = BinaryLabels(train_label.T)
labels_test = BinaryLabels(test_label.T)
epsilon = 0.001
# C = 1.0
C = 100000

# Gaussian
gauss_kernel = GaussianKernel(features_train, features_train, 1)
svm = LibSVM(C, gauss_kernel, labels_train)


# #Linear
# linear_kernel = LinearKernel(features_train, features_train)
# svm = LibSVM(C, linear_kernel, labels_train)


# svm.set_epsilon(epsilon)

svm.train()

labels_predict = svm.apply(features_test)
train_pred = svm.apply(features_train)
alphas = svm.get_alphas()
b = svm.get_bias()

eval = AccuracyMeasure()
train_eval = AccuracyMeasure()
# accuracy = eval.evaluate(labels_predict.get_labels(), labels_test)
# print(accuracy)
train_eval.evaluate(train_pred, labels_train)
train_accuracy = train_eval.get_accuracy() * 100

eval.evaluate(labels_predict, labels_test)
accuracy = eval.get_accuracy() * 100

roc = ROCEvaluation()
roc_pred = roc.evaluate(labels_predict, labels_test)
test = roc.get_ROC()
plt.plot(test[0], test[1], marker='.')
plt.show()

print('Concatenated')
print('Alphas:', alphas)
print('Bias:', b)
print('AUC:', roc_pred)

print('Train Accuracy(%):', train_accuracy)
print('Accuracy(%):', accuracy)
