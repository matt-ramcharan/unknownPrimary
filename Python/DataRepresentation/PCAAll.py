import pandas as pd
ds = pd.read_csv(r'C:\Users\matt-\Documents\Uni\TechnicalProject\unknownPrimary\Python\DataFormatting\FullDataColoRectalBreastPancreas93.csv')
Y = ds['Label']
from sklearn.preprocessing import LabelEncoder
labelencoder_y = LabelEncoder()
Y = labelencoder_y.fit_transform(Y)
# print(Y)
X = ds.drop('Label', axis=1).values

print ('done')

#preprocess data
print ('preprocessing data...', end='')
from MKLpy.preprocessing import normalization, rescale_01
X = rescale_01(X)	#feature scaling in [0,1]
X = normalization(X) #||X_i||_2^2 = 1


from sklearn.decomposition import PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(X)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

finalDf = pd.concat([principalDf, pd.Series(Y,name="target")], axis = 1)

from matplotlib import pyplot as plt

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
#targets = [3,1,2,0]#Rectum, colon, pancreas, breast, 3,1,2,0
targets = [0,1,2,3]
colors = ['r', 'b', 'g', 'm']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
#ax.legend(['Rectum', 'Colon', 'Pancreas', 'Breast'])
ax.legend(['Breast', 'Colon','Pancreas', 'Rectum' ])
ax.grid()
plt.savefig('PCAAll.pdf')
fig.show()