#load data
# print ('loading \'breast cancer\' dataset...', end='')
# from sklearn.datasets import load_breast_cancer
# ds = load_breast_cancer()
# X,Y = ds.data, ds.target
import pandas as pd
ds = pd.read_csv(r'C:\Users\matt-\Documents\Uni\TechnicalProject\unknownPrimary\Python\DataFormatting\FullDataBreastPancreas182.csv')
Y = ds['Label']
from sklearn.preprocessing import LabelEncoder
labelencoder_y = LabelEncoder()
Y = labelencoder_y.fit_transform(Y)
# print(Y)
X = ds.drop('Label', axis=1).values

# Delete all columns with zeroes
mask = (X==0)
alt_mask = [not(any(col)) for col in mask.T]
X = X[:,alt_mask]


print ('imported data')

# from fancyimpute import KNN
# from numpy import nan
#
# X[X==0] = nan
#
#
# knnImpute = KNN(k=3)
# X_filled_knn = knnImpute.fit_transform(X)
#

# print ('Zeroes imputed')

'''
WARNING: be sure that your matrix is not sparse! EXAMPLE:
from sklearn.datasets import load_svmlight_file
X,Y = load_svmlight_file(...)
X = X.toarray()
'''

#preprocess data
print ('preprocessing data...', end='\n')
from MKLpy.preprocessing import normalization, rescale_01
X = rescale_01(X)	#feature scaling in [0,1]
X = normalization(X) #||X_i||_2^2 = 1


# #compute normalized homogeneous polynomial kernels with degrees 0,1,2,...,10.
# print ('computing monotone Conjunctive Kernels...', end='')
# from MKLpy.metrics import pairwise
# from MKLpy.preprocessing import kernel_normalization
# #WARNING: the maximum arity of the conjunctive kernel depends on the number of active variables for each example,
# # that is 4 in the case of iris dataset binarized
# KL = [kernel_normalization(pairwise.monotone_conjunctive_kernel(X, c=c)) for c in range(5)]
# print ('done')

# print ('computing Linear Kernel', end='\n')
# from MKLpy.metrics import pairwise
# KL = [pairwise.homogeneous_polynomial_kernel(X, degree=1)]

import numpy as np

print ('Computing RBF Kernel')
from sklearn.metrics.pairwise import rbf_kernel
from MKLpy.preprocessing import kernel_normalization
grange = np.linspace(10e-3,10e3,10)
KLC = [kernel_normalization(rbf_kernel(X,gamma=c)) for c in grange]

from sklearn.metrics import accuracy_score, roc_auc_score
from MKLpy.model_selection import train_test_split
from LOOCV import cross_val_predict
from sklearn.svm import SVC
import numpy as np
from MKLpy.algorithms import EasyMKL, KOMD  # KOMD is not a MKL algorithm but a simple kernel machine like the SVM
gamma_best_results = {}
#Gamma tuning
for c, KL in enumerate(KLC):
    ind = grange[c]
    KL = [KL]
    #train/test KL split (N.B. here we split a kernel list directly)

    KLtr,KLte,Ytr,Yte = train_test_split(KL, Y, test_size=.25)#, random_state=42)

    # MKL algorithms

    # from MKLpy.model_selection import cross_val_score, cross_val_predict

    print ('tuning lambda for EasyMKL...', end='\n')
    base_learner = SVC(C=10000)	# simil hard-margin svm
    best_results = {}
    for lam in [0, 0.01, 0.1, 0.2, 0.9, 1]:	# possible lambda values for the EasyMKL algorithm
        # MKLpy.model_selection.cross_val_predict performs the cross validation automatically, it optimizes the accuracy
        # the counterpart cross_val_score optimized the roc_auc_score (use score='roc_auc')
        # WARNING: these functions will change in the next version
        scores = cross_val_predict(KLtr, Ytr, EasyMKL(learner=base_learner, lam=lam), score='accuracy')
        print ('Validation scores are: ' + str(scores), end='\n')
        acc = np.mean(scores)
        if not best_results or best_results['score'] < acc:
            best_results = {'lam' : lam, 'score' : acc}
    print('Best validation accuracy: %.3f with lambda: %i' %(best_results['score'],best_results['lam']))

    #evaluation on the test set

    print ('done')
    clf = EasyMKL(learner=base_learner, lam=best_results['lam']).fit(KLtr,Ytr)
    y_pred = clf.predict(KLte)
    y_score = clf.decision_function(KLte)		#rank
    accuracy = accuracy_score(Yte, y_pred)
    roc_auc = roc_auc_score(Yte, y_score)
    tr_pred = clf.predict(KLtr)
    tr_err = accuracy_score(Ytr, tr_pred)
    print ('Training Error: %.3f' % tr_err)
    print ('accuracy on the test set: %.3f, roc AUC score: %.3f' % (accuracy, roc_auc))
    # print (clf)


    if not gamma_best_results or gamma_best_results['score'] < accuracy:
        gamma_best_results = {'KL': ind, 'score': accuracy, 'roc':roc_auc}
        Yte_best = Yte
        y_score_best = y_score

print('Best  accuracy: %.3f and ROCAUC: %.3f with gamma: %.3f' %(gamma_best_results['score'], gamma_best_results['roc'], gamma_best_results['KL']))

base_learner = SVC(C=10000)
clf = EasyMKL(learner=base_learner, lam=best_results['lam']).fit(KLtr, Ytr)
y_pred = clf.predict(KLte)
y_score = clf.decision_function(KLte)  # rank
accuracy = accuracy_score(Yte, y_pred)
roc_auc = roc_auc_score(Yte, y_score)
tr_pred = clf.predict(KLtr)
tr_err = accuracy_score(Ytr, tr_pred)
print('Training Error: %.3f' % tr_err)
print('accuracy on the test set: %.3f, roc AUC score: %.3f' % (accuracy, roc_auc))

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_aucc = dict()
for i in range(2):
    fpr[i], tpr[i], _ = roc_curve(Yte_best, y_score_best)
    roc_aucc[i] = auc(fpr[i], tpr[i])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(Yte_best.ravel(), y_score_best.ravel())
roc_aucc["micro"] = auc(fpr["micro"], tpr["micro"])

plt.figure()
lw = 2
plt.plot(fpr[1], tpr[1], color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_aucc[1
    ])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()

f= open("Figs/RBF.txt","w+")
f.write('Best validation accuracy: %.9f with C: %i' %(gamma_best_results['score'],gamma_best_results['C']))
f.write('Training Error: %.9f' % tr_err)
f.write('Accuracy score: %.9f, roc AUC score: %.9f' % (accuracy, roc_auc))
f.close()
