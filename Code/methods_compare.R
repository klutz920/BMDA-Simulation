# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_microbiome/shared/nicole/")

# Load libraries
library(SummarizedExperiment)
library(lefser)
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

# Load functions
Rcpp::sourceCpp('code/core_zinb_x5.cpp') # Core function for implementing BMDA
source('code/compare_other_method.R')
source('code/functions.R')

# Input settings
input_path <- paste0("simulate_weiss.Rdata");
output_path <- paste0("result_nicole_perm/");
methods <- c("ANOVA", "Krustal", "DESeq2", "edgeR", "ZINB-DPP", "LEfSe", "ANCOM-BC");

# Load files
load(input_path);
for (i in 1:400) {
  print(i);
  Y <- simu.data[[i]];
  z <- as.numeric(substr(rownames(Y), unlist(gregexpr("grp", rownames(Y))) + 3, nchar(rownames(Y)))) - 1;
  n <- dim(Y)[1];
  p <- dim(Y)[2];
  K <- max(z) + 1;
  
  # Run the following two lines for permutation
  # set.seed(i)
  # z <- z[sample(n, n)]
  
  if ("LEfSe" %in% methods) {
    start_time <- proc.time();
    group <- factor(z, labels = c("grp1", "grp2"))
    counts <- as.matrix(t(Y))
    colData <- DataFrame(Treatment = group, row.names = colnames(counts))
    expr <- SummarizedExperiment(assays = list(counts = counts), colData = colData)
    lefse <- lefser(
      expr = expr,
      kruskal.threshold = 1.0,
      wilcox.threshold = 1.0,
      lda.threshold = 0,  # any taxa with lda > 2 are significant
      groupCol = "Treatment",
      blockCol = NULL,
      assay = 1L,
      trim.names = FALSE
    )
    gamma <- data.frame(Names = rownames(counts))
    gamma <- left_join(gamma, lefse, by = "Names")
    score <- gamma$scores
    names(score) <- gamma$Names
    R <- list(score = score)
    end_time <- proc.time();
    time <- end_time - start_time;
    save(R, Y, time, file = paste0(output_path, names(simu.data)[[i]], "_lefse.Rdata"));
  }
  
  if ("ANCOM-BC" %in% methods) {
    start_time <- proc.time();
    dat <- Y
    physeq <- phylo(x = dat, z)
    ancom <- ancombc(physeq, formula = "group", group = "group", p_adj_method = "none")
    ancom <- data.frame(Names = rownames(ancom$res$p_val), scores = ancom$res$p_val)
    names(ancom) <- c("Names", "scores")
    gamma <- data.frame(Names = colnames(dat))
    gamma <- left_join(gamma, ancom, by = "Names")
    pvalue <- gamma$scores
    names(pvalue) <- gamma$Names
    R <- list(pvalue = pvalue)
    end_time <- proc.time();
    time <- end_time - start_time;
    save(R, Y, time, file = paste0(output_path, names(simu.data)[[i]], "_ancombc.Rdata"));
  }
  
  if ("ZINB-DPP" %in% methods) {
    iter <- 10000;
    aggregate <- FALSE;
    S <- matrix(0, nrow = 1, ncol = 1);
    G <- matrix(0, nrow = 1, ncol = 1);
    
    # s_tss <- rescale(size_factor_estimator(Y, method = "TSS"), constraint = "product");
    s_1 <- rep(1, n);
    # Y_temp <- Y/s_tss;
    # Y_temp[which(Y_temp == 0)] <- NA;
    # Y_temp <- log(Y_temp);
    # b_0 <- apply(Y_temp, 2, var, na.rm = TRUE);
    # B <- matrix(1, nrow = K, ncol = p);
    # for (k in 1:K) {
    #   B[k,] <- apply(Y_temp[which(z == k - 1),], 2, var, na.rm = TRUE);
    # }
    start_time <- proc.time();
    R <- zinb_model_estimator(Y, z, s_1, iter, TRUE, S, aggregate, FALSE, 1, 100, FALSE, G);
    end_time <- proc.time();
    time <- end_time - start_time;
    save(R, Y, time, file = paste0(output_path, names(simu.data)[[i]], "_zinb_dpp.Rdata"));
    # save(R, Y, time, file = paste0(output_path, names(simu.data)[[i]], "_a=", a, "_b=", b, "_h=", h, ".Rdata"));
  }
  
  if ("ANOVA" %in% methods) {
    if (K == 2) {
      R <- mbTest("t", t(Y), as.factor(z));
    } else {
      R <- mbTest("anova", t(Y), as.factor(z));
    }
    save(R, Y, file = paste0(output_path, names(simu.data)[[i]], "_anova.Rdata"));
  }
  
  if ("Krustal" %in% methods) {
    if (K == 2) {
      R <- mbTest("wilcox", t(Y), as.factor(z));
    } else {
      R <- mbTest("kruskal", t(Y), as.factor(z));
    }
    save(R, Y, file = paste0(output_path, names(simu.data)[[i]], "_krustal.Rdata"));
  }
  
  if ("DESeq2" %in% methods) {
    R <- mbTest("deseq2", t(Y), as.factor(z));
    save(R, Y, file = paste0(output_path, names(simu.data)[[i]], "_deseq2.Rdata"));
  }
  
  if ("edgeR" %in% methods) {
     R <- mbTest("edger", t(Y), as.factor(z));
     save(R, Y, file = paste0(output_path, names(simu.data)[[i]], "_edger.Rdata"));
  }
}