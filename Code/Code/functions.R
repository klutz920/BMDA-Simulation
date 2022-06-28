# ========================================================
# This function is for the implementation of ANCOM-BC,
# prepared by Kevin and modified by Qiwei
# ========================================================
phylo = function(x, z){
  # matrix of counts (rows are taxa, columns are samples)
  otu_mat = as.matrix(t(x))
  
  # meta gives the group labels of the samples
  group <- factor(z, labels = c("grp1", "grp2"))
  meta = data.frame(group = group, row.names = colnames(otu_mat), stringsAsFactors = FALSE)
  
  # use phyloseq library to create the following objects needed for ANCOM
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  META = sample_data(meta)
  physeq = phyloseq(OTU, META) # this is they phyloseq object
  return(physeq)
}

# ========================================================
# This function creates the confusion table give the 
# actual binary condition (zero vs. one) and the prediction
# one.
# ========================================================
tabulate_error = function(true, pred) {
  table = matrix(0L, 2, 2);
  p <- length(true);
  for (i in 1:p) {
    table[pred[i] + 1, true[i] + 1] <- table[pred[i] + 1, true[i] + 1] + 1;
  }
  return (table);
}

# ========================================================
# This function computes the power under each pre-specified
# false discovery rate between 0 and 1 with input "step".
# "true" is the actual binary condition (zero vs. one) 
# while "pred" is the predicted continuous condition, 
# where larger values suggest condition one and 
# smaller values suggest condition zero. Thresholding always
# goes from the largest value to the smallest value.
# ========================================================
pwr_fdr <- function(true, pred, step = 0.01) {
  # Read data information
  p <- length(true)
  
  # Handle the missing data
  if (sum(is.na(pred)) > 0) {
    pred[is.na(pred)] <- min(pred, na.rm = TRUE) - 1
  }
  
  # Sort the prediction for fast computing
  index <- sort(pred, index.return = TRUE, decreasing = TRUE)$ix
  pred <- pred[index]
  true <- true[index]
  
  # Initialize the confusion table
  ct <- tabulate_error(true, pred >= pred[1] + 1) 
  fdr_temp <- 0
  pwr_temp <- ct[2, 2]/sum(ct[, 2])
  count <- 2
  
  # Compute FDR and Power under each unique cutoff
  for (j in 1:p) {
    ct[2, 1 + true[j]] <- ct[2, 1 + true[j]] + 1
    ct[1, 1 + true[j]] <- ct[1, 1 + true[j]] - 1
    if (j > 1 && pred[j] != pred[j - 1]) {
      fdr_temp[count] <- ct[2, 1]/sum(ct[2,])
      pwr_temp[count] <- ct[2, 2]/sum(ct[, 2])
      count <- count + 1
    }
  }
  fdr_temp[count] <- 1
  pwr_temp[count] <- 1
  
  # Redefine the FDR support space
  fdr <- seq(0, 1, by = step)
  pwr <- rep(NA, length(fdr))
  pwr[1] <- 0
  count <- 2
  for (ii in 2:(length(fdr_temp))) {
    while (count <= length(fdr) && fdr_temp[ii] - fdr[count] >= 0) {
      pwr[count] <- ((pwr_temp[ii] - pwr_temp[ii - 1])/(fdr_temp[ii] - fdr_temp[ii - 1]))*(fdr[count] - fdr_temp[ii - 1]) + pwr_temp[ii - 1]
      count <- count + 1
    }
  }
  pwr[length(fdr)] <- 1
  fdr[length(fdr)] <- 1
  
  return(data.frame(FDR = fdr, Power = pwr))
}



# ========================================================
# This function computes the true positive rate (power) under 
# each pre-specified false positive rate between 0 and 1 with 
# input "step". "true" is the actual binary condition 
# (zero vs. one)while "pred" is the predicted continuous 
# condition, where larger values suggest condition one and 
# smaller values suggest condition zero. Thresholding always
# goes from the largest value to the smallest value.
# ========================================================
roc <- function(true, pred, step = 0.01) {
  # Read data information
  p <- length(true)
  
  # Handle the missing data
  if (sum(is.na(pred)) > 0) {
    pred[is.na(pred)] <- min(pred, na.rm = TRUE) - 1
  }
  
  # Sort the prediction for fast computing
  index <- sort(pred, index.return = TRUE, decreasing = TRUE)$ix
  pred <- pred[index]
  true <- true[index]
  
  # Initialize the confusion table
  ct <- tabulate_error(true, pred >= pred[1] + 1) 
  fpr_temp <- 0
  pwr_temp <- ct[2, 2]/sum(ct[, 2])
  count <- 2
  
  # Compute FPR and Power under each unique cutoff
  for (j in 1:p) {
    ct[2, 1 + true[j]] <- ct[2, 1 + true[j]] + 1
    ct[1, 1 + true[j]] <- ct[1, 1 + true[j]] - 1
    if (j > 1 && pred[j] != pred[j - 1]) {
      fpr_temp[count] <- ct[2, 1]/sum(ct[, 2])
      pwr_temp[count] <- ct[2, 2]/sum(ct[, 2])
      count <- count + 1
    }
  }
  fpr_temp[count] <- 1
  pwr_temp[count] <- 1
  
  # Redefine the FPR support space
  fpr <- seq(0, 1, by = step)
  pwr <- rep(NA, length(fpr))
  pwr[1] <- 0
  count <- 2
  for (ii in 2:(length(fpr_temp))) {
    while (count <= length(fpr) && fpr_temp[ii] - fpr[count] >= 0) {
      pwr[count] <- ((pwr_temp[ii] - pwr_temp[ii - 1])/(fpr_temp[ii] - fpr_temp[ii - 1]))*(fpr[count] - fpr_temp[ii - 1]) + pwr_temp[ii - 1]
      count <- count + 1
    }
  }
  pwr[length(fpr)] <- 1
  fpr[length(fpr)] <- 1
  
  return(data.frame(FPR = fpr, TPR = pwr))
}

# ========================================================
# This function computes area under the ROC curve based on
# the direct output of the roc function. The first and second 
# columns of the data must be the false positive rates and 
# true positive rates in an non-decreasing order.
# ========================================================
auc <- function(data) {
  l <- dim(data)[1]
  fpr <- data[,1]
  pwr <- data[,2]
  temp <- 0
  for (ii in 2:l) {
    temp <- temp + (pwr[ii] + pwr[ii - 1])*(fpr[ii] - fpr[ii - 1])/2
  }
  return (temp)
}



# ========================================================
# This function calculates the posterior probability of 
# inclusion cutoff that controls the desired false discovery 
# rate.
# ========================================================
BayFDR <- function(PPI, alpha){
  PPI_sorted = sort(PPI, decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI_sorted[1:k])
    k = k + 1
  }
  return(PPI_sorted[k])
}

