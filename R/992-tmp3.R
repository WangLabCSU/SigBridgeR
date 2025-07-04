library(DEGAS)
library(Rtsne)
library(ggplot2)
library(reticulate)

use_condaenv("venv_3.12", required = TRUE)
py_config()


## Warning: package 'ggplot2' was built under R version 3.5.2
scDat = read.csv('~/R/Project/R_code/DEGAS_data/scDat.csv', row.names = 1)
patDat = read.csv('~/Data/DEGAS/patDat.csv', row.names = 1)
patLab = read.csv('~/Data/DEGAS/patLab.csv')

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
