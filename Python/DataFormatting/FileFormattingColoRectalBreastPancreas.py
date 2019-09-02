import pandas as pd
from pathlib import Path

# dataDir = Path(r"C:\Users\matt-\Documents\Uni\TechnicalProject\GeneBinaryClass")
# dataDir = Path(r"~\Documents\TechnicalProject\unknownPrimary\Python\DataFormatting")
# rectumPath = dataDir / "RectumData93.csv"
# colonPath = dataDir / "ColonData500.csv"
rectumPath = "RectumData93.csv"
colonPath = "ColonData432.csv"

breastPath = "BreastData500.csv"
pancreasPath = "PancreasData182.csv"

rectum = pd.read_csv(rectumPath,index_col=0).T
rectum = rectum.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
rectum['Label'] = 'Rectum'

colon = pd.read_csv(colonPath,index_col=0).T
colon = colon.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
colon['Label'] = 'Colon'

breast = pd.read_csv(breastPath,index_col=0).T
breast = breast.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
breast['Label'] = 'Breast'

pancreas = pd.read_csv(pancreasPath,index_col=0).T
pancreas = pancreas.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
pancreas['Label'] = 'Pancreas'


minlen = min(rectum.shape[0], colon.shape[0], breast.shape[0], pancreas.shape[0])
rechead = rectum.head(n=minlen)
colhead = colon.head(n=minlen)
panhead = pancreas.head(n=minlen)
brehead = breast.head(n=minlen)


dataset = pd.concat([rechead, colhead, panhead, brehead])
dataset.to_csv('FullDataColoRectalBreastPancreas' + str(minlen) + '.csv', index=False)
