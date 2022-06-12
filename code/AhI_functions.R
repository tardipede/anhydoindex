area_UC <- function(Q, p, limits = c(0, 24), step = 0.1, prop = 0.9) {
  require(bayestestR)

  x <- seq(limits[1], limits[2], by = step)
  lambda <- -log(1 - prop) / Q
  y <- (1 - exp(-lambda * x)) * p
  AUC_abs <- area_under_curve(x, y, method = c("trapezoid"))
  AUC_max <- diff(limits)
  AUC_rel <- AUC_abs / AUC_max

  return(AUC_rel)
}

area_UC_list <- function(Q, p, prop, limits = c(0, 24)) {
  data <- data.frame(Q = Q, p = p)
  results <- as.numeric(unlist(apply(data, MARGIN = 1, FUN = function(x) {
    area_UC(Q = x[1], p = x[2], prop = prop, limits = limits)
  })))
  return(results)
}

my_summary_function <- function(x) {
  results <- c(mean(x, na.rm = T), sd(x, na.rm = T), quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  names(results) <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "9.75%")
  return(results)
}

bayes_testing <- function(x, y) {
  pval <- as.numeric(pd_to_p(p_direction(x - y)))
  mean_diff <- abs(mean(x - y))
  es <- mean_diff / sqrt(var(x) + var(y))

  results <- c(pval, es)
  names(results) <- c("pval", "es")
  return(results)
}

AhI_run_model <- function(data,
                          column_IDs = c(1, 2, 3, 4),
                          prop = 0.9,
                          n.iter = 1000000) {
  require(R2jags)

  # Rename columns to match the script below
  colnames(data)[column_IDs] <- c("population", "time", "moving", "total")

  # Set the population factors in alphabetical order
  data$population <- factor(data$population, levels = sort(unique(data$population)))
  # Save the populations order
  pops_order <- levels(data$population)

  # Prepare data for JAGS
  data.jags <- list(
    time = data$time,
    moving = data$moving,
    total = data$total,
    species_num = as.numeric(as.factor(data$population)),
    Nsp = length(unique(data$population)),
    prop = prop
  )

  f <- tempfile("model.txt")
  sink(f)
  cat("

  # Prepare JAGS modelz
  model {
    # Priors
    for (sp in 1:Nsp) {                # weak estimate Q and p separately for each species
      p[sp] ~ dunif(0,1)               # low informative prior bound to be >0
      Q[sp] ~ dnorm(0,1.0E-3)T(0,)     # low informative prior
    }
    # Likelihood
    for (i in 1:length(time)){
      lambda[i] <- -log(1-prop)/Q[species_num[i]]
      mu[i] <- (1-(exp(-1*lambda[i]*time[i]))) * p[species_num[i]]
      mu_regularized[i]  <- ifelse(mu[i] == 0, 0.00001, ifelse(mu[i] == 1, 0.99999, mu[i])) # This is to avoid the model to crash when p is exactly 0 or 1
      moving[i]~dbinom(mu_regularized[i],total[i])
    }} ")

  closeAllConnections()
  # Run JAGS model

  n.burnin=floor(n.iter/3)
  n.thin = max(1, floor((n.iter - n.burnin)/1000))

  results.jags <- jags(
    data = data.jags,
    parameters.to.save = c("p", "Q"),
    model.file = f,
    n.chains = 3,
    n.iter = n.iter,
    n.burnin=n.burnin,
    n.thin = n.thin
  )

  attr(results.jags, "data") <- data
  attr(results.jags, "params") <- list(prop = prop)
  attr(results.jags, "num_to_pop") <- data.frame(population = pops_order, number = 1:length(pops_order))

  return(results.jags)
}

AhI_analyze_model <- function(model, pairwise.test = T, limits = c(0, 24)) {
  require(R2jags)
  require(tidyverse)
  if (pairwise.test == T) {
    require(bayestestR)
  }

  # Extract chains
  chains <- data.frame(do.call(rbind, as.mcmc(model)))

  # Extract Q and p posteriors and name them with the population names
  Qs <- chains[, grep("Q", colnames(chains))]
  Qs <- Qs[, order(as.numeric(unlist(colnames(Qs) %>% str_match_all("[0-9]+"))))]
  colnames(Qs) <- attributes(model)$num_to_pop$population

  ps <- chains[, grep("p", colnames(chains))]
  ps <- ps[, order(as.numeric(unlist(colnames(ps) %>% str_match_all("[0-9]+"))))]
  colnames(ps) <- attributes(model)$num_to_pop$population

  # Calculate the AhI indexes
  prop <- attributes(model)$params$prop

  temp_AhI <- list()
  for (i in 1:length(attributes(model)$num_to_pop$population)) {
    temp_AhI[[i]] <- area_UC_list(
      Q = Qs[, attributes(model)$num_to_pop$population[i]],
      p = ps[, attributes(model)$num_to_pop$population[i]], prop = prop, limits = limits
    )
  }

  AhIs <- data.frame(do.call(cbind, temp_AhI))
  colnames(AhIs) <- attributes(model)$num_to_pop$population


  # Calculate recovery speed
  # Speed = 1-(Qs-attributes(model)$params$minQ)/(attributes(model)$params$maxQ-attributes(model)$params$minQ)


  # Calculate summary tables
  p_summary <- t(apply(ps, MARGIN = 2, FUN = my_summary_function))
  Q_summary <- t(apply(Qs, MARGIN = 2, FUN = my_summary_function))
  AhI_summary <- t(apply(AhIs, MARGIN = 2, FUN = my_summary_function))


  # Put them together in a list
  results <- list(AhI = AhI_summary, Q = Q_summary, p = p_summary,
                  p_chains = ps, Q_chains = Qs, AhI_chains = AhIs)

  if (pairwise.test == T) {

    # Create a table with the possible pairwise combinations
    populations <- attributes(model)$num_to_pop$population
    pairwise.table <- expand.grid(pop1 = populations, pop2 = populations)
    pairwise.table <- pairwise.table[pairwise.table$pop1 != pairwise.table$pop2, ]

    pairwise.table <- data.frame(t(apply(pairwise.table, MARGIN = 1, FUN = sort)))
    pairwise.table <- pairwise.table[duplicated(pairwise.table), ]

    pairwise.table <- data.frame(pairwise.table,
      p_pval = rep(NA, nrow(pairwise.table)), p_es = rep(NA, nrow(pairwise.table)),
      Q_pval = rep(NA, nrow(pairwise.table)), speed_es = rep(NA, nrow(pairwise.table)),
      AhI_pval = rep(NA, nrow(pairwise.table)), AhI_es = rep(NA, nrow(pairwise.table))
    )


    for (i in 1:nrow(pairwise.table)) {
      pairwise.table[i, 3:4] <- p_pvals <- bayes_testing(ps[, pairwise.table[i, 1]], ps[, pairwise.table[i, 2]])
      pairwise.table[i, 5:6] <- Q_pvals <- bayes_testing(Qs[, pairwise.table[i, 1]], Qs[, pairwise.table[i, 2]])
      pairwise.table[i, 7:8] <- AhI_pvals <- bayes_testing(AhIs[, pairwise.table[i, 1]], AhIs[, pairwise.table[i, 2]])
    }

    results[[7]] <- pairwise.table
    names(results)[7] <- "pairwise_p"
  }

  return(results)
}

plot_AhI <- function(analyzed_model) {
  require(ggplot2)
  require(ggrepel)
  require(ggstance)
  require(RColorBrewer)

  plotdata <- data.frame(
    treatment = rownames(analyzed_model$AhI),
    AhI = analyzed_model$AhI[, 1],
    p = analyzed_model$p[, 1],
    p_sd = analyzed_model$p[, 2],
    Q = analyzed_model$Q[, 1],
    Q_sd = analyzed_model$Q[, 2]
  )

  high_pvals <- subset(analyzed_model$pairwise_p[, c(1, 2, 7, 8)], AhI_pval > 0.05)
  high_pvals <- merge(high_pvals, plotdata[, c(1, 3, 5)], by.x = 1, by.y = 1)
  colnames(high_pvals)[5:6] <- c("p1", "Q1")
  high_pvals <- merge(high_pvals, plotdata[, c(1, 3, 5)], by.x = 2, by.y = 1)
  colnames(high_pvals)[7:8] <- c("p2", "Q2")



  ggplot(plotdata) +
    theme_bw() +
    geom_segment(
      data = high_pvals,
      aes(x = p1, y = -Q1, xend = p2, yend = -Q2, size = log(1 / AhI_es)), col = "red", alpha = 0.5
    ) +
    geom_text_repel(aes(x = p, y = -Q, label = treatment),
      min.segment.length = 0, size = 4,
      segment.colour = "black", force_pull = 1, box.padding = 1
    ) +
    geom_linerange(aes(x = p, ymin = -(Q - Q_sd), ymax = -(Q + Q_sd))) +
    geom_linerangeh(aes(y = -Q, xmin = p - p_sd, xmax = p + p_sd)) +
    geom_point(aes(x = p, y = -Q, fill = AhI), pch = 21, size = 5) +
    xlab("Max survival") +
    ylab("Recovery speed") +
    # coord_fixed(ratio=1,xlim=c(0,1),ylim=c(0,1))+
    scale_size_continuous(range = c(0.5, 2)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_fill_distiller(palette = "Spectral", direction = 1)
}
