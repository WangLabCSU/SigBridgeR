library(DEGAS)
library(Rtsne)
library(ggplot2)

## Warning: package 'ggplot2' was built under R version 3.5.2
scDat = read.csv('~/Data/DEGAS/scDat.csv', row.names = 1)
patDat = read.csv('~/Data/DEGAS/patDat.csv', row.names = 1)
patLab = read.csv('~/Data/DEGAS/patLab.csv')

path.data = ''
path.result = ''
initDEGAS()
set_seed_term(2)
tmpDir = 'tmp/'

ccModel1 = runCCMTLBag(
  scDat,
  NULL,
  patDat,
  patLab,
  tmpDir,
  'BlankClass',
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

# Predicting patient outcome in cells
# ie, predicting GBM subtype association in individual cells
scpatPreds = predClassBag(ccModel1, scDat, 'pat')
colnames(scpatPreds) = colnames(patLab)


# Set seed and run tSNE
set.seed(1)
scDat_tsne = Rtsne(scDat)
colnames(scDat_tsne$Y) = c('tSNE1', 'tSNE2')
# kNN smoothing of GBM subtype association
impressions_sc_smooth = knnSmooth(scpatPreds, scDat_tsne$Y)
# Conversion of GBM subtype association to correlation
impressions_sc_smooth_cor = toCorrCoeff(impressions_sc_smooth)
tmp = data.frame(
  tSNE1 = scDat_tsne$Y[, "tSNE1"],
  tSNE2 = scDat_tsne$Y[, "tSNE2"],
  Proneural = impressions_sc_smooth_cor[, "Proneural"],
  Neural = impressions_sc_smooth_cor[, "Neural"],
  Classical = impressions_sc_smooth_cor[, "Classical"],
  Mesenchymal = impressions_sc_smooth_cor[, "Mesenchymal"],
  Patient = sub("[_].*", "", rownames(scDat))
)
# Mesenchymal
p = ggplot(
  tmp,
  aes(x = tSNE1, y = tSNE2, color = Mesenchymal, shape = Patient)
) +
  geom_point() +
  scale_color_gradient2(low = "black", mid = "lavender", high = "red")
plot(
  p +
    labs(color = 'Mesenchymal\nassociation', shape = 'Patient') +
    theme(
      legend.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1)),
      axis.text.x = element_text(size = rel(1)),
      axis.text.y = element_text(size = rel(1))
    )
)


# Proneural
p = ggplot(tmp, aes(x = tSNE1, y = tSNE2, color = Proneural, shape = Patient)) +
  geom_point() +
  scale_color_gradient2(low = "black", mid = "lavender", high = "red")
plot(
  p +
    labs(color = 'Proneural\nassociation', shape = 'Patient') +
    theme(
      legend.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1)),
      axis.text.x = element_text(size = rel(1)),
      axis.text.y = element_text(size = rel(1))
    )
)


# Neural
p = ggplot(tmp, aes(x = tSNE1, y = tSNE2, color = Neural, shape = Patient)) +
  geom_point() +
  scale_color_gradient2(low = "black", mid = "lavender", high = "red")
plot(
  p +
    labs(color = 'Neural\nassociation', shape = 'Patient') +
    theme(
      legend.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1)),
      axis.text.x = element_text(size = rel(1)),
      axis.text.y = element_text(size = rel(1))
    )
)


# Classical
p = ggplot(tmp, aes(x = tSNE1, y = tSNE2, color = Classical, shape = Patient)) +
  geom_point() +
  scale_color_gradient2(low = "black", mid = "lavender", high = "red")
plot(
  p +
    labs(color = 'Classical\nassociation', shape = 'Patient') +
    theme(
      legend.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1)),
      axis.text.x = element_text(size = rel(1)),
      axis.text.y = element_text(size = rel(1))
    )
)
