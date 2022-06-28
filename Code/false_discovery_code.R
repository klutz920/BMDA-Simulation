# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_microbiome/shared/nicole/")

# Load library
library(plyr)
library(ggplot2)
library(bayefdr)

# Load functions
source("code/functions.R")

# Load data
input_path <- "result_nicole_perm/"
files <- list.files(input_path)
N <- length(files) # The total number of result files

# Read results
alpha <- 0.05            # Significance level
loc <- rep(NA, N)        # Feces or skins sample
n <- rep(NA, N)          # 12 or 54 samples per group
seed <- rep(NA, N)       # Data replicate index
sigma <- rep(NA, N)      # Effect size: exp(1) or exp(2)
method <- rep(NA, N)     # Method name: BMDA, ANCOM-BC, LEfSe, edgeR, DESeq2, t-test, or Wilcoxon
fd <- rep(NA, N)         # Number of false discoveries 
for (i in 1:length(files)) {
  print(i)
  
  # Extract scenario information
  temp <- unlist(strsplit(files[i], split = "_"))
  loc[i] <- temp[2]
  n[i] <- as.numeric(temp[5])
  seed[i] <- as.numeric(temp[3])
  sigma[i] <- round(as.numeric(temp[4]), 2)
  method[i] <- sub(".Rdata", "", temp[6])
  
  # Load result
  load(paste0(input_path, files[i]))
  if (grepl("anova|krustal|edger|deseq2", temp[6])) {
    fd[i] <- sum(p.adjust(R$pvalue, method = "fdr") <= alpha, na.rm = TRUE)
  } 
  if (grepl("ancombc", temp[6])) {
    fd[i] <- sum(p.adjust(R$pvalue, method = "fdr") <= alpha, na.rm = TRUE)
  }
  if (grepl("lefse", temp[6])) {
    fd[i] <- sum(R$score >= 2, na.rm = TRUE)
  }
  if (grepl("zinb", temp[6])) {
    if (sum(R$gamma_ppi >= 0.5) == 0) {
      fd[i] <- 0
    } else {
      # fd[i] <- sum(R$gamma_ppi >= BayFDR(R$gamma_ppi, alpha))
      temp_2 <- efdr_search(probs = R$gamma_ppi, target_efdr = alpha, min_threshold = 0.5001)
      fd[i] <- sum(R$gamma_ppi >= temp_2$threshold[optimal(temp_2)])
    }
  }
}
results <- data.frame(location = loc, sample.size.per.group = n, replicate.no = seed, effect.size = sigma, method = method, false.discovery = fd)
results$method <- factor(results$method, levels = c("ancombc", "anova", "deseq2", "edger", "krustal", "lefse", "zinb"), labels = c("ANCOM-BC", "t-test", "DESeq2", "edgeR", "LEfSe", "Wilcoxon", "BMDA"))
results$sample.size.per.group <- factor(results$sample.size.per.group, levels = c(12, 54), labels = c("12 samples per group", "54 samples per group"))
results$effect.size <- factor(results$effect.size, levels = c(2.72, 7.39), labels = c("Effect size = exp(1)", "Effect size = exp(2)"))

# Marginalize the results over replicates
results.s <- ddply(results, .variables = c("location", "sample.size.per.group", "effect.size", "method"), summarise, false.discovery = mean(false.discovery))
write.csv(results, file = "false_discovery_result_for_michael.csv")

# Generate the graphical summary
ggplot(data = subset(results, effect.size == "Effect size = exp(1)" & method %in% c("ANCOM-BC", "DESeq2", "edgeR", "LEfSe", "BMDA"))) +
  geom_boxplot(mapping = aes(x = method, y = false.discovery, color = method)) +
  facet_grid(sample.size.per.group ~ location) +
  labs(x = "", y = "Number of false discoveries", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))

ggplot(data = subset(results.s, effect.size == "Effect size = exp(1)" & method %in% c("ANCOM-BC", "DESeq2", "edgeR", "LEfSe", "BMDA"))) +
  geom_bar(mapping = aes(x = method, y = false.discovery, fill = method, color = method), stat = "identity") +
  facet_grid(sample.size.per.group ~ location) +
  labs(x = "", y = "Number of false discoveries", fill = "") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_color_discrete(guide = 'none') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))







