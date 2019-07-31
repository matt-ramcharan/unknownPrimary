import pandas as pd
from pathlib import Path

# dataDir = Path(r"C:\Users\matt-\Documents\Uni\TechnicalProject\GeneBinaryClass")
# dataDir = Path(r"~\Documents\TechnicalProject\unknownPrimary\Python\DataFormatting")
# rectumPath = dataDir / "RectumData93.csv"
# colonPath = dataDir / "ColonData500.csv"
rectumPath = "RectumData93.csv"
colonPath = "ColonData432.csv"

rectum = pd.read_csv(rectumPath,index_col=0).T
rectum = rectum.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
rectum['Label'] = 'Rectum'

colon = pd.read_csv(colonPath,index_col=0).T
colon = colon.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
colon['Label'] = 'Colon'

minlen = min(rectum.shape[0], colon.shape[0])
rechead = rectum.head(n=minlen)
colhead = colon.head(n=minlen)


dataset = pd.concat([rechead, colhead])
dataset.to_csv('FullDataColoRectal' + str(minlen) + '.csv', index=False)
