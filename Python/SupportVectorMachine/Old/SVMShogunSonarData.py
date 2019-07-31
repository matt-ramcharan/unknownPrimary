import pandas as pd
from shogun.Loss import L2R_L2LOSS_SVC
from shogun.Evaluation import RealFeatures, BinaryLabels, LibSVM, AccuracyMeasure, ROCEvaluation
from shogun.Kernel import GaussianKernel, LibLinear
import matplotlib.pyplot as plt

# #Sonar Dataset
train_data = pd.read_csv('sonar.tr',sep='\s+', header = None, names=list(range(1,62)))
train_feats = train_data.loc[:,:train_data.shape[1]-1].as_matrix()
train_label = train_data.loc[:,train_data.shape[1]].as_matrix()

test_data = pd.read_csv('sonar.ts',sep='\s+', header = None, names=list(range(1,62)))
test_feats = test_data.loc[:,:test_data.shape[1]-1].as_matrix()
test_label = test_data.loc[:,test_data.shape[1]].as_matrix()


features_train = RealFeatures(train_feats.T)
features_test = RealFeatures(test_feats.T)
labels_train = BinaryLabels(train_label.T)
labels_test = BinaryLabels(test_label.T)
# epsilon = 0.001
# C = 1.0
C = 100000

# Gaussian
gauss_kernel = GaussianKernel(features_train, features_train, 1)

svm = LibSVM(C, gauss_kernel, labels_train)

# #Linear
# svm = LibLinear(C, features_train, labels_train)
# svm.set_liblinear_solver_type(L2R_L2LOSS_SVC)

#svm.set_epsilon(epsilon)

svm.train()


labels_predict = svm.apply(features_test)

alphas = svm.get_alphas()
b = svm.get_bias()

eval = AccuracyMeasure()

roc = ROCEvaluation()
roc_pred = roc.evaluate(labels_predict, labels_test)
test = roc.get_ROC()
plt.plot(test[0],test[1], marker='.')
plt.show()

eval.evaluate(labels_predict,labels_test)
accuracy=eval.get_accuracy()*100



print('Alphas:', alphas)
print('Bias:', b)
print('AUROC:', roc_pred)
print('Accuracy(%):', accuracy)