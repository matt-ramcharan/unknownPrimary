import pandas as pd
ds = pd.read_csv(r'C:\Users\matt-\Documents\Uni\TechnicalProject\unknownPrimary\Python\DataFormatting\FullData.csv')
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

#preprocess data
print ('preprocessing data...', end='')
from MKLpy.preprocessing import normalization, rescale_01
X = rescale_01(X)	#feature scaling in [0,1]
X = normalization(X) #||X_i||_2^2 = 1

#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics

from sklearn.model_selection import train_test_split
# Split dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.25) # 70% training and 30% test

from sklearn.svm import SVC
base_learner = SVC(C=10000, gamma='auto')

from sklearn.ensemble import AdaBoostClassifier

import numpy as np
from sklearn.model_selection import LeaveOneOut
#loo = LeaveOneOut()
#from tqdm import tqdm
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import LeaveOneOut

best_results = {}
for lr in [0.01, 0.1, 0.2, 0.9, 1]:	# possible lambda values for the EasyMKL algorithm
    # MKLpy.model_selection.cross_val_predict performs the cross validation automatically, it optimizes the accuracy
    # the counterpart cross_val_score optimized the roc_auc_score (use score='roc_auc')
    # WARNING: these functions will change in the next version
    # scores = cross_val_predict(X_train, y_train, AdaBoostClassifier(base_estimator = base_learner, learning_rate=lr, n_estimators=100, algorithm = 'SAMME', random_state=0), score='accuracy')
    #scores=[]
    clf = AdaBoostClassifier(base_estimator=base_learner, learning_rate=lr, n_estimators=100, algorithm='SAMME',
                             random_state=0)
    print('Starting Validation  ')
    scores = cross_val_score(clf, X=X, y=Y, cv=LeaveOneOut())
    print(scores)

    print ('Validation scores are: ' + str(scores), end='\n')
    acc = np.mean(scores)
    print('Validation Accuracy is: %.3f' % acc)
    if not best_results or best_results['score'] < acc:
        best_results = {'lr' : lr, 'score' : acc}
print('Best validation accuracy: %.3f with lr: %i' %(best_results['score'],best_results['lr']))



#model = clf.fit(X_train, y_train)

# AdaBoostClassifier(algorithm='SAMME.R', base_estimator=None,
#         learning_rate=1.0, n_estimators=100, random_state=0)

#
# y_pred = model.predict(X_test)
#
# # Model Accuracy, how often is the classifier correct?
# print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

#print(clf.feature_importances_)
#print(clf.predict([[0, 0, 0, 0]]))
#print(clf.score(X, Y))