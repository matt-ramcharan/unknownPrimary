import pandas as pd
from pathlib import Path


breastPath = "BreastData500.csv"
pancreasPath = "PancreasData182.csv"

breast = pd.read_csv(breastPath,index_col=0).T
breast = breast.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
breast['Label'] = 'Breast'

pancreas = pd.read_csv(pancreasPath,index_col=0).T
pancreas = pancreas.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
pancreas['Label'] = 'Pancreas'

minlen = min(pancreas.shape[0], breast.shape[0])
panhead = pancreas.head(n=minlen)
brehead = breast.head(n=minlen)

dataset = pd.concat([brehead, panhead])
dataset.to_csv('FullDataBreastPancreas' + str(minlen) + '.csv', index=False)