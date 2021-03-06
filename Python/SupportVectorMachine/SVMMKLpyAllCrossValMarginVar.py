#load data
# print ('loading \'breast cancer\' dataset...', end='')
# from sklearn.datasets import load_breast_cancer
# ds = load_breast_cancer()
# X,Y = ds.data, ds.target
import pandas as pd
ds = pd.read_csv(r'C:\Users\matt-\Documents\Uni\TechnicalProject\unknownPrimary\Python\DataFormatting\FullDataColoRectalBreastPancreas93.csv')
Y = ds['Label']
from sklearn.preprocessing import LabelEncoder
labelencoder_y = LabelEncoder()
Y = labelencoder_y.fit_transform(Y)
# print(Y)
X = ds.drop('Label', axis=1).values


print ('done')

'''
WARNING: be sure that your matrix is not sparse! EXAMPLE:
from sklearn.datasets import load_svmlight_file
X,Y = load_svmlight_file(...)
X = X.toarray()
'''

#preprocess data
print ('preprocessing data...', end='')
from MKLpy.preprocessing import normalization, rescale_01
X = rescale_01(X)	#feature scaling in [0,1]
X = normalization(X) #||X_i||_2^2 = 1

#train/test split
from sklearn.model_selection import train_test_split
Xtr,Xte,Ytr,Yte = train_test_split(X,Y, test_size=.25, random_state=42, shuffle=True)
print ('Splitting Done')


print ('computing Linear Kernel', end='\n')
from MKLpy.metrics import pairwise
KLtr = [pairwise.homogeneous_polynomial_kernel(Xtr, degree=1)]
KLte = [pairwise.homogeneous_polynomial_kernel(Xte,Xtr, degree=1)]
# from MKLpy.metrics import pairwise
# KLtr = [pairwise.homogeneous_polynomial_kernel(Xtr, degree=d) for d in range(degrees)]
# KLte = [pairwise.homogeneous_polynomial_kernel(Xte,Xtr, degree=d) for d in range(degrees)]
# # KLtr = [pairwise.homogeneous_polynomial_kernel(Xtr, degree=1)]
# # KLte = [pairwise.homogeneous_polynomial_kernel(Xte,Xtr, degree=1)]
print ('done')

# MKL algorithms
from MKLpy.algorithms import EasyMKL, KOMD	# KOMD is not a MKL algorithm but a simple kernel machine like the SVM
# from MKLpy.model_selection import cross_val_score, cross_val_predict
from LOOCV import cross_val_predict
from sklearn.svm import SVC
import numpy as np
best_results ={}

#clf = EasyMKL(learner=SVC(C=10), lam=0, multiclass_strategy='ovo').fit(KLtr,Ytr)

for C in [1e-4,1e-3,1e-2,1e-1,1,1e2,1e3,1e4,1e5]:
    base_learner = SVC(C=C)	# simil hard-margin svm
    scores = cross_val_predict(KLtr, Ytr, EasyMKL(learner=base_learner, lam=0, multiclass_strategy='ovo'), score='accuracy')
    #print ('Validation scores are: ' + str(scores), end='\n')
    acc = np.mean(scores)
    print('Acc: %.9f with C: %i' %(acc,C))
    if not best_results or best_results['score'] < acc:
        best_results = {'C' : C, 'score' : acc}
print('Best validation accuracy: %.9f with C: %i' %(best_results['score'],best_results['C']))



# # #MKL algorithms
# from MKLpy.algorithms import AverageMKL, EasyMKL, KOMD	#KOMD is not a MKL algorithm but a simple kernel machine like the SVM
# # print ('training AverageMKL...', end='')
# # clf = AverageMKL().fit(KLtr,Ytr)	#a wrapper for averaging kernels
# # print ('done')
# # print(clf.weights)			#print the weights of the combination of base kernels
# # K_average = clf.ker_matrix	#the combined kernel matrix
#
#
# print ('training EasyMKL...', end='')
# clf = EasyMKL(lam=0.1).fit(KLtr,Ytr)		#combining kernels with the EasyMKL algorithm
# #lam is a hyper-parameter in [0,1]
# print ('done')
# print (clf.weights)
#
#evaluate the solution
from sklearn.metrics import accuracy_score, roc_auc_score, balanced_accuracy_score
clf = EasyMKL(learner=SVC(C=best_results['C']), lam=0, multiclass_strategy='ovo').fit(KLtr,Ytr)
print (clf.weights)
tr_pred = clf.predict(KLtr)
# tr_err = accuracy_score(Ytr, tr_pred)
tr_err = balanced_accuracy_score(Ytr, tr_pred)
print ('Training Error: %.3f' % tr_err)
y_pred = clf.predict(KLte)					#predictions
# y_score = clf.decision_function(KLte)		#rank

# accuracy = accuracy_score(Yte, y_pred)
accuracy = balanced_accuracy_score(Yte, y_pred)
# roc_auc = roc_auc_score(Yte, y_score)

print ('Accuracy score: %.9f' % accuracy)# , roc AUC score: % .9f') # % (accuracy))#, roc_auc))


# #select the base-learner
# #MKL algorithms use a hard-margin as base learned (or KOMD in the case of EasyMKL). It is possible to define a different base learner
# from sklearn.svm import SVC
# base_learner = SVC(C=0.1)
# clf = EasyMKL(learner=base_learner)
# clf = clf.fit(KLtr,Ytr)
# print(clf)

# import matplotlib.pyplot as plt
# from sklearn.metrics import roc_curve, auc
#
# # Compute ROC curve and ROC area for each class
# fpr = dict()
# tpr = dict()
# roc_aucc = dict()
# for i in range(2):
#     fpr[i], tpr[i], _ = roc_curve(Yte, y_score)
#     roc_aucc[i] = auc(fpr[i], tpr[i])
#
# # Compute micro-average ROC curve and ROC area
# fpr["micro"], tpr["micro"], _ = roc_curve(Yte.ravel(), y_score.ravel())
# roc_aucc["micro"] = auc(fpr["micro"], tpr["micro"])
#
# plt.figure()
# lw = 2
# plt.plot(fpr[1], tpr[1], color='darkorange',
#          lw=lw, label='ROC curve (area = %0.2f)' % roc_aucc[1
#     ])
# plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver operating characteristic: Colon, Rectum, Breast and Pancreas')
# plt.legend(loc="lower right")
# plt.savefig('Figs/ROCAll.pdf',format='pdf')
# plt.show()

f = open("All.txt","w+")
f.write('Best validation accuracy: %.3f with C: %.9f' %(best_results['score'],best_results['C']))
f.write('Training Error: %.9f' % tr_err)
f.write('Accuracy score: %.9f, roc AUC score: %.9f' % (accuracy, roc_auc))
f.close()
