# Vasilisa A. Rudneva
# Nov 2018
# 
# This function takes as input a matrix of beta values obtained from rgSetToBetasFiltering
# > betas[1:5,1:5]
#           9297949037_R03C02 9297949037_R05C02 9305216163_R01C01 9305216153_R05C02 9297949040_R01C01
#cg00000957        0.77045408        0.84237386        0.71323971        0.78473429        0.82436785
#cg00001583        0.08450619        0.03680815        0.06126489        0.03824306        0.04667641
#cg00002028        0.08187812        0.09494766        0.14017470        0.13759215        0.10718627
#cg00002719        0.04595234        0.05079716        0.05278532        0.04935480        0.07177044
#cg00002837        0.86909894        0.80211855        0.81687188        0.86475943        0.49879025#
#
# Selects N most variable probes as SD > 0.25 then measures weighted correlations and runs t-SNE
# Returs a matrix with coordinates for t-SNE1 and t-SNE2 which can be then plotted

RuntSNE <- function(betas) {
  
  betas.sd <- apply(betas, 1, sd)
  m.sd <- names(sort(betas.sd, decreasing = TRUE))
  
  # select only most variable probes (SD > 0.25)
  probesNr <- length(which(betas.sd > 0.25))
  cor.mvprobes <- wtd.cors(betas[m.sd[1:probesNr],], y=NULL, weight=betas.sd[m.sd[1:probesNr]]-betas.sd[m.sd[probesNr]])
  
  set.seed(17112018)
  Y <- Rtsne(1-cor.mvprobes, pca = FALSE, verbose = FALSE, is_distance = TRUE, theta = 0, max_iter = 2000, dims = 2, perplexity = 10)$Y
  rownames(Y) <- rownames(cor.mvprobes)
  
  return(Y)
}
