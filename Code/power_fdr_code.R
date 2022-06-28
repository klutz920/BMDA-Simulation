# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_microbiome/shared/nicole")

# Load library
library(plyr)
library(ggplot2)

# Load functions
source("code/functions.R")

# Load data
input_path <- "result_nicole/"
files <- list.files(input_path)
N <- length(files) # The total number of result files

# Read results
h <- 0.005                 # The step size of false discovery rate (FDR) between 0 and 1
k <- length(seq(0, 1, h))  # The total number of FDR per result
loc <- rep(NA, k*N)        # feces or skins sample
n <- rep(NA, k*N)          # 12 or 54 samples per group
seed <- rep(NA, k*N)       # data replicate 
sigma <- rep(NA, k*N)      # effect size exp(1) or exp(2)
method <- rep(NA, k*N)     # BMDA, ANCOM-BC, LEfSe, edgeR, DESeq2, t-test, or Wilcoxon
fdr <- rep(NA, k*N)        # False discovery rate
pwr <- rep(NA, k*N)        # Power
for (i in 1:length(files)) {
  print(i)
  
  # Extract scenario information
  temp <- unlist(strsplit(files[i], split = "_"))
  loc[((i - 1)*k + 1):(i*k)] <- temp[2]
  n[((i - 1)*k + 1):(i*k)] <- as.numeric(temp[5])
  seed[((i - 1)*k + 1):(i*k)] <- as.numeric(temp[3])
  sigma[((i - 1)*k + 1):(i*k)] <- round(as.numeric(temp[4]), 2)
  method[((i - 1)*k + 1):(i*k)] <- sub(".Rdata", "", temp[6])
  
  # Load result
  load(paste0(input_path, files[i]))
  true <- grepl("TP", colnames(Y))
  if (grepl("anova|krustal|edger|deseq2", temp[6])) {
    pred <- -R$p.value 
  } 
  if (grepl("ancombc", temp[6])) {
    pred <- -R$pvalue 
  }
  if (grepl("lefse", temp[6])) {
    pred <- R$score 
  }
  if (grepl("zinb", temp[6])) {
    pred <- R$gamma_ppi
  }
  res <- pwr_fdr(true, pred, h)
  fdr[((i - 1)*k + 1):(i*k)] <- res$FDR
  pwr[((i - 1)*k + 1):(i*k)] <- res$Power
}
results <- data.frame(location = loc, sample.size.per.group = n, replicate.no = seed, effect.size = sigma, method = method, false.discovery.rate = fdr, power = pwr)
results$method <- factor(results$method, levels = c("ancombc", "anova", "deseq2", "edger", "krustal", "lefse", "zinb"), labels = c("ANCOM-BC", "t-test", "DESeq2", "edgeR", "LEfSe", "Wilcoxon", "BMDA"))
results$sample.size.per.group <- factor(results$sample.size.per.group, levels = c(12, 54), labels = c("12 samples per group", "54 samples per group"))
results$effect.size <- factor(results$effect.size, levels = c(2.72, 7.39), labels = c("Effect size = exp(1)", "Effect size = exp(2)"))

# Marginalize the results over replicates
results.s <- ddply(results, .variables = c("location", "sample.size.per.group", "effect.size", "method", "false.discovery.rate"), summarise, power = mean(power))
write.csv(results.s, file = "power_fdr_for_michael.csv")

# Generate the graphical summary
ggplot(data = subset(results.s, effect.size == "Effect size = exp(1)" & method %in% c("ANCOM-BC", "DESeq2", "edgeR", "LEfSe", "BMDA"))) +
  geom_line(mapping = aes(x = false.discovery.rate, y = power, color = method)) +
  facet_grid(sample.size.per.group ~ location) +
  labs(x = "False discovery rate", y = "Power", color = "") +
  coord_cartesian(xlim = c(0, 0.20), ylim = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))







