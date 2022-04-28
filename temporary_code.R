data = expand.grid(speed = seq(0,1,0.005), survival = seq(0,1,0.005))

data$Q = (1-data$speed + 0.25)*(24-0.25)

data$AhI = apply(data, MARGIN = 1, FUN = function(x){area_UC(x[3],x[2], limits = c(0,24))})

ggplot(data)+
  theme_bw()+
  geom_tile(aes(x=survival, y=speed, fill=AhI))+
  scale_fill_distiller(palette = "Spectral", direction = 1)

ggplot(data)+
  theme_bw()+
  geom_contour_filled(aes(x=survival, y=speed, z=AhI), binwidth = 0.1, col="black")+
  scale_fill_brewer(palette = "Spectral", direction = 1)
