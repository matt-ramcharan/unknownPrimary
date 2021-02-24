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

from sklearn.model_selection import train_test_split

#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics

# Split dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.25, random_state=1) # 70% training and 30% test

from sklearn.ensemble import AdaBoostClassifier

clf = AdaBoostClassifier(n_estimators=100, random_state=0)

model = clf.fit(X_train, y_train)

# AdaBoostClassifier(algorithm='SAMME.R', base_estimator=None,
#         learning_rate=1.0, n_estimators=100, random_state=0)


y_pred = model.predict(X_test)
y_pred_train = model.predict(X_train)
# Model Accuracy, how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
print("Train Accuracy:",metrics.accuracy_score(y_train, y_pred_train))

print(clf.feature_importances_)
#print(clf.predict([[0, 0, 0, 0]]))
#print(clf.score(X, Y))