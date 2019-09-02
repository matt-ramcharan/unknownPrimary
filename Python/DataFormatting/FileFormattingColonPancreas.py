import pandas as pd
from pathlib import Path

# dataDir = Path(r"C:\Users\matt-\Documents\Uni\TechnicalProject\GeneBinaryClass")
# dataDir = Path(r"~\Documents\TechnicalProject\unknownPrimary\Python\DataFormatting")
# pancreasPath = dataDir / "RectumData93.csv"
# colonPath = dataDir / "ColonData500.csv"
pancreasPath = "PancreasData182.csv"
colonPath = "ColonData432.csv"

pancreas = pd.read_csv(pancreasPath,index_col=0).T
pancreas = pancreas.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
pancreas['Label'] = 'Pancreas'

colon = pd.read_csv(colonPath,index_col=0).T
colon = colon.drop(['__no_feature', '__ambiguous', '__too_low_aQual',
       '__not_aligned', '__alignment_not_unique'],axis=1)
colon['Label'] = 'Colon'

minlen = min(pancreas.shape[0], colon.shape[0])
panhead = pancreas.head(n=minlen)
colhead = colon.head(n=minlen)


dataset = pd.concat([panhead, colhead])
dataset.to_csv('FullDataColonPancreas' + str(minlen) + '.csv', index=False)
