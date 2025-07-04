# *GBM Example

library(DEGAS)
library(Rtsne)
library(ggplot2)
library(reticulate)

use_condaenv("venv_3.12", required = TRUE)
py_config()


# *count inversed
scDat = utils::read.csv(
  '~/R/Project/R_code/DEGAS_data/GBM/scDat.csv',
  row.names = 1
)
# *bulk inversed
patDat = utils::read.csv(
  '~/R/Project/R_code/DEGAS_data/GBM/patDat.csv',
  row.names = 1
)
# ?phnotype
patLab = utils::read.csv(
  '~/R/Project/R_code/DEGAS_data/GBM/patLab.csv',
  row.names = 1
)

path.data = '~/R/Project/R_code/DEGAS_data/test_result'
path.result = '~/R/Project/R_code/DEGAS_data/test_result'
DEGAS::initDEGAS()
DEGAS::set_seed_term(2)
tmpDir = '~/R/Project/R_code/DEGAS_data/test_result/tmp/'

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

## 3-layer DenseNet BlankClass DEGAS model

## 0

## 3-layer DenseNet BlankClass DEGAS model

## 0

## 3-layer DenseNet BlankClass DEGAS model

## 0

## 3-layer DenseNet BlankClass DEGAS model

## 0

## 3-layer DenseNet BlankClass DEGAS model

## 0
