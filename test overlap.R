AhI_analyze_model <- function(model, pairwise.test = T, limits = c(0, 24)) {
  require(R2jags)
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
  results <- list(AhI = AhI_summary, Q = Q_summary, p = p_summary)

  if (pairwise.test == T) {

    # Create a table with the possible pairwise combinations
    populations <- attributes(model)$num_to_pop$population
    pairwise.table <- expand.grid(pop1 = populations, pop2 = populations)
    pairwise.table <- pairwise.table[pairwise.table$pop1 != pairwise.table$pop2, ]

    pairwise.table <- data.frame(t(apply(pairwise.table, MARGIN = 1, FUN = sort)))
    pairwise.table <- pairwise.table[duplicated(pairwise.table), ]

    pairwise.table <- data.frame(pairwise.table,
                                overlap = rep(NA, nrow(pairwise.table)),
                                t_dist = rep(NA, nrow(pairwise.table)),
                                p_dir = rep(NA, nrow(pairwise.table)))



    for (i in 1:nrow(pairwise.table)) {
      pairwise.table[i,3] <- overlap_2d(ps[, pairwise.table[i, 1]],Qs[, pairwise.table[i, 1]],ps[, pairwise.table[i, 2]],Qs[, pairwise.table[i, 2]])
      pairwise.table[i,4] <- t_dist(list(ps[, pairwise.table[i, 1]],Qs[, pairwise.table[i, 1]],ps[, pairwise.table[i, 2]],Qs[, pairwise.table[i, 2]]))
      x_sd = sd(c(ps[, pairwise.table[i, 1]],ps[, pairwise.table[i, 2]]))
      y_sd = sd(c(Qs[, pairwise.table[i, 1]],Qs[, pairwise.table[i, 2]]))
      pairwise.table[i, 5] <- pd_2d(list((ps[, pairwise.table[i, 1]])/x_sd,(Qs[, pairwise.table[i, 1]])/y_sd,
                                         (ps[, pairwise.table[i, 2]])/x_sd,(Qs[, pairwise.table[i, 2]])/y_sd))



    }

    results[[4]] <- pairwise.table
    names(results)[4] <- "pairwise_p"
  }

  return(results)
}

analyzed_mod = AhI_analyze_model(model_Milnesium)


library(dplyr)
Qg = gather(Qs, "treatment", "Q")
pg = gather(ps, "treatment", "p")

plot_df = cbind(pg,Qg)
sum(plot_df[,1] == plot_df[,3]) == nrow(plot_df)

plot_df = plot_df[,c(1,2,4)]
library(ggplot2)
library(tidyverse)
library(ggExtra)
library(ggrepel)
require(ggstance)

high_pval = subset(analyzed_mod$pairwise_p, overlap>0.05)
high_pval$p1 = analyzed_mod$p[match(high_pval$X1,rownames(analyzed_mod$p)),1]
high_pval$p2 = analyzed_mod$p[match(high_pval$X2,rownames(analyzed_mod$p)),1]
high_pval$Q1 = analyzed_mod$Q[match(high_pval$X1,rownames(analyzed_mod$Q)),1]
high_pval$Q2 = analyzed_mod$Q[match(high_pval$X2,rownames(analyzed_mod$Q)),1]

ggplot(plot_df)+
  theme_bw()+
  stat_bag(aes(x= p, y = Q, fill = treatment),
           prop = 0.80, col="black",
           alpha=0.1, size=0.2, show.legend=F) +
  geom_point(data = aggregate(plot_df[,2:3], by=list(plot_df[,1]), FUN = mean),
             aes(x=p, y=Q))+
  scale_y_reverse()+
  geom_segment(data = high_pval, aes(x=p1, xend=p2,y=Q1,yend=Q2),
               col="red", lty="dashed", alpha=0.75)+
  geom_text_repel(data = aggregate(plot_df[,2:3], by=list(plot_df[,1]), FUN = mean),
                  aes(x = p, y = Q, label = Group.1),
                  min.segment.length = 0, size = 4,
                  segment.colour = "black", force_pull = 1, box.padding = 1
  )


