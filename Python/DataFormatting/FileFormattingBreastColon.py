import pandas as pd
from pathlib import Path


dataDir = Path(r"C:\Users\matt-\Documents\Uni\TechnicalProject\GeneBinaryClass")
breastPath = dataDir / "BreastData500.csv"
colonPath = dataDir / "ColonData500.csv"

breast = pd.read_csv(breastPath,index_col=0).T
breast = breast.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
breast['Label'] = 'Breast'

colon = pd.read_csv(colonPath,index_col=0).T
colon = colon.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
colon['Label'] = 'Colon'


dataset = pd.concat([breast, colon])
dataset.to_csv('FullDataBreastColon.csv', index=False)
