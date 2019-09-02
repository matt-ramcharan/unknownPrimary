#load data
# print ('loading \'breast cancer\' dataset...', end='')
# from sklearn.datasets import load_breast_cancer
# ds = load_breast_cancer()
# X,Y = ds.data, ds.target
import pandas as pd
ds = pd.read_csv(r'C:\Users\matt-\Documents\Uni\TechnicalProject\unknownPrimary\Python\DataFormatting\FullDataColoRectal93.csv')
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

#Reshape to 'images'
X_train.resize(X_train.shape[0], 120, 120,  1)
X_test.resize(X_test.shape[0], 120, 120, 1)

#Model

import numpy
from keras.models import Sequential
from keras.layers import Dense, Dropout, Conv2D, Flatten

# fix random seed for reproducibility
seed = 7
numpy.random.seed(seed)

classifier = Sequential()
#add model layers
classifier.add(Conv2D(64, kernel_size=3, activation='relu', input_shape=(120, 120, 1)))
classifier.add(Conv2D(128, kernel_size=3, activation='relu'))
classifier.add(Conv2D(256, kernel_size=3, activation='relu'))
classifier.add(Flatten())
# classifier.add(Dropout(0.25))
# classifier.add(Dense(36864, activation='relu'))
# classifier.add(Dense(1024, activation='relu'))
# classifier.add(Dense(512, activation='relu'))
classifier.add(Dense(1, activation='softmax'))

#Compiling the neural network
classifier.compile(optimizer ='adam',loss='binary_crossentropy', metrics =['accuracy'])
#Fitting the data to the training dataset
history = classifier.fit(X_train,y_train, validation_split=.2, batch_size=10, epochs=150, verbose=0)

eval_model=classifier.evaluate(X_test,y_test)

y_pred=classifier.predict(X_test)
y_pred =(y_pred>0.5)
print(eval_model)

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
plt.show()

