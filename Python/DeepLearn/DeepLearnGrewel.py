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
print ('imported data')

#preprocess data

print ('preprocessing data...', end='\n')

# Delete all columns with zeroes
mask = (X==0)
alt_mask = [not(any(col)) for col in mask.T]
X = X[:,alt_mask]

from MKLpy.preprocessing import normalization, rescale_01
X = rescale_01(X)	#feature scaling in [0,1]
X = normalization(X) #||X_i||_2^2 = 1

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=.25, shuffle=True, random_state=1)

#Model

import numpy
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout

# fix random seed for reproducibility
seed = 7
numpy.random.seed(seed)

classifier = Sequential()
#First Hidden Layer
classifier.add(Dense(2000, activation='relu', kernel_initializer='random_normal', input_dim=X.shape[1]))#Second  Hidden Layer
# classifier.add(Dense(2000, activation='relu', kernel_initializer='random_normal'))
classifier.add(Dense(66, activation='relu', kernel_initializer='random_normal'))

#classifier.add(Dense(1, activation='sigmoid', kernel_initializer='random_normal'))

# classifier.add(Dropout(0.5))

#classifier.add(Dense(5, activation='relu', kernel_initializer='random_normal'))#Output Layer
classifier.add(Dense(1, activation='sigmoid', kernel_initializer='random_normal'))

#Compiling the neural network
classifier.compile(optimizer ='adam',loss='binary_crossentropy', metrics =['accuracy'])
#Fitting the data to the training dataset
history = classifier.fit(X_train,y_train, validation_split=.2, batch_size=1, epochs=150, verbose=0)

eval_model=classifier.evaluate(X_test,y_test)

y_pred=classifier.predict(X_test)
y_pred =(y_pred>0.5)
print(eval_model)

eval_model_train=classifier.evaluate(X_train,y_train)
print(eval_model_train)

from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)
print(cm)

print('Accuracy: %.3f' % eval_model[1])

from matplotlib import pyplot as plt
# list all data in history
#print(history.history.keys())
# summarize history for accuracy
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')

plt.show()
# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.savefig('LossBreastPan.pdf', format='pdf')
plt.show()

