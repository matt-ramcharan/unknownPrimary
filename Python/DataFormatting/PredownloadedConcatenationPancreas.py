import pandas as pd
import os
import numpy as np

base_dir = r"C:\Users\matt-\Documents\Uni\TechnicalProject\unknownPrimary\Python\DataFormatting"
num = 500
## Pancreas Data#
pancreas_dir = base_dir+"/Pancreas"

pancreas = []

for count,(r,_,filename) in enumerate(os.walk(pancreas_dir)):
    if count<num:
        for file in filename:
            if file.endswith(".gz"):
                if count == 1:
                    pancreas = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t', names=[str(count)], index_col=0)
                else:
                    pancreas_data = pd.read_csv(str(r) + '/' + str(file), compression='gzip', delimiter='\t',
                                              names=[str(count)], index_col=0)
                    pancreas = pd.concat((pancreas, pancreas_data[str(count)]), axis=1)
                    continue
            else:
                continue

pancreas.to_csv('PancreasData'+ str(pancreas.shape[1]) + '.csv')
