# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_microbiome/shared/nicole/")

# Load library
library(plyr)
library(ggplot2)
library(bayefdr)

# Load functions
source("code/functions.R")

# Load data
input_path <- "result_nicole_real_history/"
files <- list.files(input_path)
N <- length(files) # The total number of result files

# Read results
alpha <- 0.05            # Significance level
seed <- rep(NA, N)       # Data replicate index
method <- rep(NA, N)     # Method name: BMDA, ANCOM-BC, LEfSe, edgeR, DESeq2, t-test, or Wilcoxon
fd <- rep(NA, N)         # Number of false discoveries 
for (i in 1:length(files)) {
  print(i)
  
  # Extract scenario information
  temp <- unlist(strsplit(files[i], split = "_"))
  seed[i] <- as.numeric(temp[2])
  method[i] <- sub(".Rdata", "", temp[3])
  
  # Load result
  load(paste0(input_path, files[i]))
  if (grepl("anova|krustal|edger|deseq2", temp[3])) {
    fd[i] <- sum(p.adjust(R$p.value, method = "fdr") <= alpha, na.rm = TRUE)
  } 
  if (grepl("ancombc", temp[3])) {
    fd[i] <- sum(p.adjust(R$pvalue, method = "fdr") <= alpha, na.rm = TRUE)
  }
  if (grepl("lefse", temp[3])) {
    fd[i] <- sum(R$score >= 2, na.rm = TRUE)
  }
  if (grepl("zinb", temp[3])) {
    if (sum(R$gamma_ppi >= 0.5) == 0) {
      fd[i] <- 0
    } else {
      # fd[i] <- sum(R$gamma_ppi >= BayFDR(R$gamma_ppi, alpha))
      temp_2 <- efdr_search(probs = R$gamma_ppi, target_efdr = alpha, min_threshold = 0.5001)
      fd[i] <- sum(R$gamma_ppi >= temp_2$threshold[optimal(temp_2)])
    }
  }
}
results <- data.frame(replicate.no = seed, method = method, false.discovery = fd)
results$method <- factor(results$method, levels = c("ancombc", "anova", "deseq2", "krustal", "lefse", "zinb"), 
                         labels = c("ANCOM-BC", "t-test", "DESeq2", "Wilcoxon", "LEfSe", "BMDA"))

# Marginalize the results over replicates
# results.s <- ddply(results, .variables = c("method"), summarise, false.discovery = mean(false.discovery))
write.csv(results, file = "false_discovery_result_for_michael_history.csv")

# Generate the graphical summary
ggplot(data = (results)) +
  geom_boxplot(mapping = aes(x = method, y = false.discovery, color = method)) +
  labs(x = "", y = "Number of false discoveries", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))

ggplot(data = subset(results, method %in% c("ANCOM-BC", "t-test", "Wilcoxon", "BMDA"))) +
  geom_boxplot(mapping = aes(x = method, y = false.discovery, color = method)) +
  labs(x = "", y = "Number of false discoveries", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))






