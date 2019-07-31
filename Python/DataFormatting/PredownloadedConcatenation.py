import pandas as pd
import os
import numpy as np

base_dir = "/home/matt/Documents/TechnicalProject/Data/GDC"

##Number of samples

num = 500


## Breast Data
breast_dir = base_dir+"/Breast"

breast = []

for count,(r,_,filename) in enumerate(os.walk(breast_dir)):
    #print(filename)
    #print(r)
    if count<num:
        for file in filename:
            if file.endswith(".gz"):
                if count == 1:
                    #breast_data = pd.read_csv("/home/matt/Documents/TechnicalProject/Data/GDC/Breast/dc6315b8-b182-49c3-b372-6bd8c210358c/a18f26c8-9497-4824-a5ba-7b192f4ca987.htseq.counts.gz", compression='gzip', delimiter='\t', names=['Gene', 'Count'])
                    breast = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t', names=[str(count)], index_col=0)
                    #breast.append(breast_data)
                else:
                    breast_data = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t',
                                              names=[str(count)],index_col=0)
                    breast = pd.concat((breast, breast_data[str(count)]), axis=1)
                    continue
            else:
                continue


##COMPLETE CHECKS SAME GENES ARE USED IN FILES
#savemat('BreastData500.mat', breast)
breast.to_csv('BreastData' + str(breast.shape[1]) + '.csv')

# #Debugging
# frame = None
# breast = None

# ## Colon Data
colon_dir = base_dir+"/Colon"

colon = []

for count,(r,_,filename) in enumerate(os.walk(colon_dir)):
    if count<num:
        for file in filename:
            if file.endswith(".gz"):
                if count == 1:
                    colon = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t', names=[str(count)], index_col=0)
                else:
                    colon_data = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t',
                                              names=[str(count)], index_col=0)
                    colon = pd.concat((colon, colon_data[str(count)]), axis=1)
                    continue
            else:
                continue


##COMPLETE CHECKS SAME GENES ARE USED IN FILES
#savemat('ColonData500.mat', colon)
colon.to_csv('ColonData'+ str(colon.shape[1]) + '.csv')

## Rectum Data#
rectum_dir = base_dir+"/Rectum"

rectum = []

for count,(r,_,filename) in enumerate(os.walk(rectum_dir)):
    if count<num:
        for file in filename:
            if file.endswith(".gz"):
                if count == 1:
                    rectum = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t', names=[str(count)], index_col=0)
                else:
                    rectum_data = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t',
                                              names=[str(count)], index_col=0)
                    rectum = pd.concat((rectum, rectum_data[str(count)]), axis=1)
                    continue
            else:
                continue

rectum.to_csv('RectumData'+ str(rectum.shape[1]) + '.csv')
