#load data
# print ('loading \'breast cancer\' dataset...', end='')
# from sklearn.datasets import load_breast_cancer
# ds = load_breast_cancer()
# X,Y = ds.data, ds.target
import pandas as pd
ds = pd.read_csv('/home/matt/Documents/TechnicalProject/unknownPrimary/Python/DataFormatting/FullDataColoRectal93.csv')
Y = ds['Label']
from sklearn.preprocessing import LabelEncoder
labelencoder_y = LabelEncoder()
Y = labelencoder_y.fit_transform(Y)
# print(Y)
X = ds.drop('Label', axis=1).values
print ('imported data')

#preprocess data

print ('preprocessing data...', end='\n')

# # Delete all columns with zeroes
# mask = (X==0)
# alt_mask = [not(any(col)) for col in mask.T]
# X = X[:,alt_mask]

# from MKLpy.preprocessing import normalization, rescale_01
# X = rescale_01(X)	#feature scaling in [0,1]
# X = normalization(X) #||X_i||_2^2 = 1

from sklearn.feature_selection import SelectKBest, f_regression
#train_data = X.apply(pd.to_numeric).astype('float32')

import numpy as np
from matplotlib import pyplot as plt

kb = SelectKBest(score_func=f_regression, k=60483)
kb.fit_transform(X,Y)
indices = np.argsort(kb.scores_)[::-1]
selected_features = []
for i in range(5):
  selected_features.append(ds.columns[indices[i]])
plt.figure()
plt.bar(selected_features, kb.scores_[indices[range(5)]], color='r', align='center')
plt.xticks(rotation=45)
plt.show()
print(selected_features)