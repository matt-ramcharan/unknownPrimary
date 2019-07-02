#import shogun
import pandas as pd
import numpy as np
from shogun.Kernel import MultitaskKernelTreeNormalizer, GaussianKernel
from shogun.Evaluation import RealMatrixFeatures,RealFeatures, BinaryLabels, LibSVM, AccuracyMeasure, ROCEvaluation, MultitaskROCEvaluation, Task
from shogun.Classifier import MKLClassification
import matplotlib.pyplot as plt

train_data = np.load('train.npy')
test_data = np.load('test.npy')
correlations = np.load('corr.npy')


# features_train = train_data[:,:,:-1]
# features_test = train_data[:,:,:-1]
# labels_train = train_data[:,:,-1]
# labels_test = test_data[:,:,-1]

features_train = train_data[:, :, :-1]
features_test = train_data[:, :, :-1]
labels_train = train_data[:, :, -1]
labels_test = test_data[:, :, -1]

features_train = Task(RealFeatures(features_train))
features_test = Task(RealFeatures(features_test))

labels_train = BinaryLabels(labels_train.T)
labels_test = BinaryLabels(labels_test.T)

##Multitask

epsilon = 0.001
# C = 1.0
C = 100000

# Multitask Kernel

mt_kernel = MultitaskKernelTreeNormalizer()
gauss_kernel = GaussianKernel(features_train, features_train, 1.0)

libsvm  = LibSVM()
svm = MultitaskKernelTreeNormalizer(LibSVM)
svm.set_interleaved_optimization_enabled(False)
svm.set_kernel(gauss_kernel)
svm.set_labels(labels_train)

# svm = LibSVM(C, mt_kernel, labels_train)

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

print('Multitask')
print('Alphas:', alphas)
print('Bias:', b)
print('AUC:', roc_pred)

print('Train Accuracy(%):', train_accuracy)
print('Accuracy(%):', accuracy)