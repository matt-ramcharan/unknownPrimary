import pandas as pd
ds = pd.read_csv('/home/matt/Documents/TechnicalProject/unknownPrimary/Python/DataFormatting/FullDataColoRectal93.csv')
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

from sklearn.model_selection import train_test_split

#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics

# Split dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3) # 70% training and 30% test



# from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier

clf = GradientBoostingClassifier()#n_estimators=100, random_state=0)

model = clf.fit(X_train, y_train)

# AdaBoostClassifier(algorithm='SAMME.R', base_estimator=None,
#         learning_rate=1.0, n_estimators=100, random_state=0)


y_pred = model.predict(X_test)

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import LeaveOneOut
print('Starting Validation  ')
scores = cross_val_score(clf , X = X_train , y = y_train , cv = LeaveOneOut())
import numpy as np
print(scores)
print(np.mean(scores))
# Model Accuracy, how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

print(clf.feature_importances_)
#print(clf.predict([[0, 0, 0, 0]]))
print(clf.score(X, Y))