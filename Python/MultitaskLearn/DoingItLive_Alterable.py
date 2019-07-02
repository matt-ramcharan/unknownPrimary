#import shogun
import numpy as np
import binarytree
import scipy
from binarytree import build
import pandas as pd
from scipy.spatial.distance import euclidean, pdist, squareform
import matplotlib.pyplot as plt
import pickle


def similarity_func(u, v):
    return np.dot(u, v)


def change_repres_rand(top,node,m):
    #top[node].pprint()
    subtree = top[node].values
    repres_index = np.random.choice(range(len(top[node].represN)), m, replace=False)
    for subnode in subtree:
        cur_rep = getattr(top[subnode], 'represN')
        new_rep = cur_rep
        for j in repres_index:
            new_rep[j] = cur_rep[j]*-1

        setattr(top[subnode], 'represN', new_rep)
        #print(getattr(top[subnode],'represN'))
        #print('\n')
    return top
def rand_change_recurse(root,m):
    order = root.preorder
    for node in order:
        root = change_repres_rand(root, node.value, m)
        #print([leav.represN for leav in root.leaves])
    return root


if __name__ == "__main__":

    # Build a tree from list representation

    ##Parameters
    # tree_size = 7
    # repres_size = 5
    # m=1

    tree_size = 64-1
    m = 5
    repres_size = 100

    values = range(0,tree_size)
    repres_orig = [[1]*repres_size] * tree_size
    root = build(values,repres=repres_orig.copy())
    #print(repres_orig)



    root.pprint()
    #root = change_repres_rand(root, 1, m)
    root = rand_change_recurse(root, m)

    #root.pprint()
    represes = [rep.represN for rep in root.leaves]
    #np.array(represes).T.tolist()
    labels = root.leaves

    DF_var = pd.DataFrame(represes)
    DF_var.index = labels

    dists = pdist(DF_var, similarity_func)
    DF_euclid = pd.DataFrame(squareform(dists))
    #print(DF_euclid)
    plt.matshow(DF_euclid.corr())
    plt.xlabel('Leaves')
    plt.ylabel('Leaves')
    plt.title(r'$\mu_d$ dot products')
    plt.colorbar()
    plt.show()

    #randoms = np.random.normal(represes, 20)
    #print(randoms)
    #
    # rands_DF_var = pd.DataFrame(randoms)
    # rands_DF_var.index = labels
    #
    # rands_dists = pdist(rands_DF_var, similarity_func)
    # rands_DF_euclid = pd.DataFrame(squareform(rands_dists))
    # # print(DF_euclid)
    # plt.matshow(rands_DF_euclid.corr())
    # plt.xlabel('Leaves')
    # plt.ylabel('Leaves')
    # plt.title('Sampled from distribution dot products')
    # plt.show()

    mu_pos = DF_var.multiply(2).values
    train_points = 10
    test_points = 1000
    var = 0.5

    pos_train_dataset = np.swapaxes(np.array([np.random.normal(mu_pos, var) for td in range(train_points)]),0,1)
    pos_train_dataset = np.dstack((pos_train_dataset, np.ones((32,train_points))))

    pos_test_dataset = np.swapaxes(np.array([np.random.normal(mu_pos, var) for td in range(test_points)]),0,1)
    pos_test_dataset = np.dstack((pos_test_dataset, np.ones((32,test_points))))


    mu_neg = DF_var.multiply(-2)
    neg_train_dataset = np.swapaxes(np.array([np.random.normal(mu_neg, var) for td in range(train_points)]),0,1)
    neg_train_dataset = np.dstack((neg_train_dataset, -1*np.ones((32,train_points))))

    neg_test_dataset = np.swapaxes(np.array([np.random.normal(mu_neg, var) for td in range(test_points)]),0,1)
    neg_test_dataset = np.dstack((neg_test_dataset, -1*np.ones((32,test_points))))

    train_dataset = np.concatenate((pos_train_dataset,neg_train_dataset), axis=1)
    test_dataset = np.concatenate((pos_test_dataset, neg_test_dataset), axis=1)

    np.save('train_alt', train_dataset)
    np.save('test_alt', test_dataset)
    np.save('corr_alt', represes)
