import pandas as pd
ds = pd.read_csv('/home/matt/Documents/TechnicalProject/unknownPrimary/Python/DataFormatting/FullDataColoRectal93.csv')
Y = ds['Label']
from sklearn.preprocessing import LabelEncoder
labelencoder_y = LabelEncoder()
Y = labelencoder_y.fit_transform(Y)
# print(Y)
X = pd.array(ds.drop('Label', axis=1), dtype='Int64')
from numpy import NaN
mask = X==0
# X=X.astype(pd.Int64Dtype())
X[mask] = NaN

from fancyimpute import KNN
# X is the complete data matrix
# X_incomplete has the same values as X except a subset have been replace with NaN
# Use 3 nearest rows which have a feature to fill in each row's missing features

X_filled_knn = KNN(k=3).fit_transform(X.T)
pass