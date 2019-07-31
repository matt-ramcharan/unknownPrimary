import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt

dataFile1 = '~/Documents/TechnicalProject/Data/TCGA-DX-AB35-01A-21R-A41D-13/b47ba76c-1b0b-45a7-bf6f-f14865c136bd/94eaf79f-14e8-421a-ba4f-e2f1fbae2592.FPKM.txt'
dataFile2 = '~/Documents/TechnicalProject/Data/TCGA-DX-AB35-01A-21R-A41D-13/b47ba76c-1b0b-45a7-bf6f-f14865c136bd/94eaf79f-14e8-421a-ba4f-e2f1fbae2592.htseq.counts'
dataFile3 = '~/Documents/TechnicalProject/Data/TCGA-DX-AB35-01A-21R-A41D-13/b47ba76c-1b0b-45a7-bf6f-f14865c136bd/94eaf79f-14e8-421a-ba4f-e2f1fbae2592.FPKM.txt'
dataFile4 = '~/Documents/TechnicalProject/Data/TCGA-DX-AB35-01A-21R-A41D-13/b47ba76c-1b0b-45a7-bf6f-f14865c136bd/94eaf79f-14e8-421a-ba4f-e2f1fbae2592.FPKM.txt'
dataFile5 = '~/Documents/TechnicalProject/Data/TCGA-DX-AB35-01A-21R-A41D-13/b47ba76c-1b0b-45a7-bf6f-f14865c136bd/94eaf79f-14e8-421a-ba4f-e2f1fbae2592.FPKM.txt'



FPKM = pd.read_csv(dataFile,sep='\t',names=['Gene', 'Count'])
plt.plot(FPKM['Count'])
plt.show()
