import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

def svm_func(x,y):
    ## Where x is training data and y is labels

    w = np.zeros(len(x[0]))
    #Learning rate
    l_rate = 1
    #Epoch
    epoch=10000
    #Output Classification
    out = []
    #Training
    for e in range(epoch):
        for i, val in enumerate(x):
            val1 = np.dot(x[i], w)
            if y[i]*val1<1:
                w = w + l_rate * ((y[i]*x[i]) - 2*(1/epoch)*w)
            else:
                w = w + l_rate * (- 2 * (1 / epoch) * w)

    for i, val in enumerate(x):
        out.append(np.dot(x[i], w))
    return w, out


train_data = pd.read_csv('sonar.tr',sep='\s+', header = None, names=list(range(1,62)))
train_input = train_data.loc[:,:train_data.shape[1]-1]
train_labels = train_data.loc[:,train_data.shape[1]]

test_data = pd.read_csv('sonar.ts',sep='\s+', header = None, names=list(range(1,62)))
test_input = test_data.loc[:,:test_data.shape[1]-1]
test_labels = test_data.loc[:,test_data.shape[1]]





#Input data
# x = np.array([
#     [-2,4,-1],
#     [4,1,-1],
#     [1, 6, -1],
#     [2, 4, -1],
#     [6, 2, -1],
#
# ])

#output label
# y = np.array([-1, -1, 1, 1, 1])

for val, inp in enumerate(x):
    if y[val] == -1:
        plt.scatter(inp[0], inp[1], s=100, marker='_', linewidths=5)
    else:
        plt.scatter(inp[0], inp[1], s=100, marker='+', linewidths=5)

plt.plot([-2,6],[6,1])
plt.show()

w, out = svm_func(x,y)

print('Calculated weights are: \n', w)

print('Output labels are: \n', out)
print('Which are \n', np.sign(out)==y)
print(str((np.sum(np.sign(out)==y)/len(out))*100) + '% training accuracy')

##Test input

u = np.array([
    [-1, 3, -1],
    [5, 5, -1],

])

plt.scatter(-1, 3, s=100, marker='_', linewidths=5)
plt.scatter(5, 5, s=100, marker='+', linewidths=5)


x1 = [w[0], w[1], -w[1], w[0]]
x2 = [w[0], w[1], w[1], -w[0]]

x1x2 = np.array([x1, x2])
X, Y, U, V = zip(*x1x2)
ax = plt.gca()
ax.quiver(X, Y, U, V, scale=1, color='blue')
plt.show()

result = []
for i, val in enumerate(u):
        result.append(np.dot(u[i], w))
print('test result')
print(result)