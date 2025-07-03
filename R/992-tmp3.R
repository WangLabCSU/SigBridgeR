library(DEGAS)
library(Rtsne)
library(ggplot2)
library(reticulate)

use_condaenv("venv_3.12", required = TRUE)
py_config()


## Warning: package 'ggplot2' was built under R version 3.5.2
scDat = read.csv('~/Data/DEGAS/scDat.csv', row.names = 1) # ?sample-gene count,
scLab = read.csv('~/Data/DEGAS/scLab.csv', row.names = 1) # ?phnotype
patDat = read.csv('~/Data/DEGAS/patDat.csv', row.names = 1) # ?TCGA-gene count,
patLab = read.csv('~/Data/DEGAS/patLab.csv', row.names = 1) # ?RFS time, phenotype?


path.data = ''
path.result = ''
initDEGAS()
set_seed_term(2)
tmpDir = 'tmp/'

ccModel1 = runCCMTLBag(
  scDat,
  scLab,
  patDat,
  patLab,
  tmpDir,
  'ClassClass',
  'DenseNet',
  3,
  5
)
