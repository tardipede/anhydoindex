
simulate_data_ahi <- function(arguments, pm = 0.9, times = c(1, 4, 24), total = 40, N.reps = 10) {
  p <- arguments[1]
  Q <- arguments[2]

  lambda <- -log(1 - pm) / Q
  moving_p <- (1 - exp(-lambda * times)) * p

  N.moving <- matrix(unlist((lapply(moving_p, FUN = function(x) {
    rbinom(N.reps, total, x)
  }))))

  data_fake <- data.frame(
    group = paste0("rep_", rep(1:N.reps, length(times))),
    time = rep(times, each = N.reps),
    moving = N.moving,
    total = rep(total, length(N.moving))
  )
  return(data_fake)
}


design_matrix <- expand.grid(p = seq(0.01, 0.99, length.out=90), Q = seq(0.5, 24, length.out=90))


run_model_function = function(x){
  mod = data.frame(AhI_run_model(simulate_data_ahi(x, total = 100, N.reps = 10, times = c(1, 4, 8, 12, 24)), n.iter = 500000)$BUGSoutput$summary)


  if("models.rds" %in% list.files()){
    temp.list = readRDS("models.rds")
    new.list = c(temp.list,list(mod))
    saveRDS(new.list, "models.rds")
  }

  if(!"models.rds" %in% list.files()){saveRDS(list(mod),"models.rds")}
  return(mod)

}


models = apply(design_matrix[1:nrow(design_matrix),], MARGIN = 1, FUN = run_model_function )



extract_values <- function(x) {
  data_temp <- data.frame(subset(x, rownames(x) != "deviance"))
  data_vector <- c(data_temp[, 1], data_temp[, 2])
  names(data_vector) <- c(paste0(rownames(data_temp), "_AVG"), paste0(rownames(data_temp), "_SD"))

  return(data_vector)
}

data_vectors <- do.call(rbind, lapply(models, extract_values))

data_all <- cbind(design_matrix, data_vectors)
data_all <- tidyr::gather(data_all, "variable", "value", 3:ncol(data_all))
data_all$variable <- stringr::str_remove(data_all$variable, "[[:punct:]][[:digit:]]+[[:punct:]]")

library(ggplot2)
library(patchwork)

ggplot(subset(data_all, variable == "Q_AVG")) +
  geom_point(aes(x = Q, y = value, col = p)) +
  geom_abline(slope = 1, intercept = 0, col = "red")


p_6 = ggplot(subset(data_all, variable == "p_AVG" & Q<6)) +
  theme_bw()+
  geom_point(aes(x = p, y = value, col = Q), alpha=0.75) +
  geom_abline(slope = 1, intercept = 0, col = "red")+
  coord_fixed(ratio = 1)

p_12 = ggplot(subset(data_all, variable == "p_AVG" & Q<12 & Q>=6)) +
  theme_bw()+
  geom_point(aes(x = p, y = value, col = Q), alpha=0.75) +
  geom_abline(slope = 1, intercept = 0, col = "red")+
  coord_fixed(ratio = 1)

p_18 = ggplot(subset(data_all, variable == "p_AVG" & Q<18 & Q>=12)) +
  theme_bw()+
  geom_point(aes(x = p, y = value, col = Q), alpha=0.75) +
  geom_abline(slope = 1, intercept = 0, col = "red")+
  coord_fixed(ratio = 1)

p_24 = ggplot(subset(data_all, variable == "p_AVG" & Q<24 & Q>18)) +
  theme_bw()+
  geom_point(aes(x = p, y = value, col = Q), alpha=0.75) +
  geom_abline(slope = 1, intercept = 0, col = "red")+
  coord_fixed(ratio = 1)



(p_6|p_12)/(p_18|p_24)


library(tidyverse)

subset(data_all, variable == "Q_SD") %>%
  group_by(p, Q) %>%
  mutate(
    p = p[1],
    Q = Q[1],
    value = mean(value)
  ) %>%
  ggplot() +
  theme_bw() +
  geom_tile(aes(x = p, y = -Q, fill = value)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 4)

subset(data_all, variable == "p_SD") %>%
  group_by(p, Q) %>%
  mutate(
    p = p[1],
    Q = Q[1],
    value = mean(value)
  ) %>%
  ggplot() +
  theme_bw() +
  geom_tile(aes(x = p, y = -Q, fill = value)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.02)



subset(data_all, variable == "Q_AVG") %>%
  group_by(p, Q) %>%
  mutate(
    p = p[1],
    Q = Q[1],
    value = mean(value)
  ) %>%
  ggplot() +
  theme_bw() +
  geom_tile(aes(x = p, y = -Q, fill = Q - value)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)


subset(data_all, variable == "p_AVG") %>%
  group_by(p, Q) %>%
  mutate(
    p = p[1],
    Q = Q[1],
    value = mean(value)
  ) %>%
  ggplot() +
  theme_bw() +
  geom_tile(aes(x = p, y = -Q, fill = p - value)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

