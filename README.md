# unknownPrimary
Classifying Cancers of unknown primary
##Project has been finished and is available [here](https://www.overleaf.com/read/tghfbxmbbpgw)
]

# Todo
 - Collect data htseq of mRNA Breast and Colon Data.
 - KNN Impute zeroes
 - Linearise about mean of 1
 - binary classifier SVM (linear kernel, gaussian otherwise)
 - Multiclass classifier for more cancer types.
 - Colon vs Colorectal
 - CNV and MEth
 - Other data types multi kernel.

# Plan
## Data
Build up a dataset of breast and prostate expression data for training
Design an api to collect data and preprocess for future purposes

## SVM
Start with 2 cancers and do classification of each (Breast and Prostate) with expression data

# Notes

Order of number of samples GDC is
Lung (4977), Breast (3682), Colorectal (2755) (of which is colon (2305) and rectum (371)), Kidney (2201).
Try using Breast and colon sample size 200 of each.

# Resources
[https://portal.gdc.cancer.gov/]
