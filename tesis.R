### Directorio de la carpeta tesis
setwd("~/Dropbox/tesis2")
source("utils.R")

### Carga las dependencias
load_dependencies()

datos_marg <- read.table("Datos/marg_mex2010.csv",header = TRUE, sep = ",",encoding="UTF-8",stringsAsFactors=F ) 
datos_marg$gmarg = factor(datos_marg$gmarg, levels = c("Muy bajo", "Bajo", "Medio", "Alto","Muy alto" ))
claves <- read.table("Datos/claves mun.csv",header = TRUE, sep = ",")
ent_abr = c("Ags","BC","BCS","Cam","Coah","Col","Chis","Chi","DF","Dgo","Gto","Gro","Hgo","Jal","Mex","Mich","Mor","Nay","NL","Oax","Pue","Qro","QRoo","SLP","Sin","Son","Tab","Tam","Tlax","Ver","Yuc","Zac")
datos_marg$abr_ent = as.factor(datos_marg$clave_ent)
levels(datos_marg$abr_ent) = ent_abr

### Leer polígonos
mapa_mun <-readShapePoly("Mapa/MUNICIPIOS")
mapa_mun <- mapa_mun[order(mapa_mun$CVE_ENT), ]
### Unificar id's de munciipios
mapa_mun$CVE_MUN <- (claves$id.ia)
mapa_mun <- mapa_mun[order(mapa_mun$CVE_MUN), ]

mapa_mun$imarges <- datos_marg$imarg_esc

# Permutación de valores para ejemplo
index = seq(1:dim(mapa_mun)[1])
set.seed(7)
r_index = sample(index)
mapa_mun$imarges_r <- datos_marg$imarg_esc[r_index]
#write.csv(mapa_mun[c(1:3)],"datos mapa.csv")

### Análisis exploratorio ###
nb_mun <- poly2nb(mapa_mun,mapa_mun$CVE_MUN)
listw_mun <- nb2listw(nb_mun,style="W")
lags = lag.listw(listw_mun, mapa_mun$imarges)

moran_plot_df = data.frame(marg = mapa_mun$imarges, lag = lags, name =datos_marg$nom_mun)
moran_fit = lm( lag~marg,moran_plot_df)
xtable(moran_fit)
moran_plot = ggplot(moran_plot_df, aes(x=marg, y=lag)) + 
  geom_vline(xintercept = mean(mapa_mun$imarges),color="black", size=.5, alpha=.6, linetype="longdash")+
  geom_hline(yintercept = mean(mapa_mun$imarges),color="black", size=.5, alpha=.6, linetype="longdash")+
  geom_point(color="#ff7f0e", size=1.3) +
  geom_smooth(method = "lm") + xlab("Índice de Marginación") + ylab("Índice de marginación vecinos")+
  geom_text(data = subset(moran_plot_df, abs(moran_fit$residuals)> quantile(abs(moran_fit$residuals), .995) |marg > quantile(marg, .999)  | marg < quantile(marg, .0001) |lag > quantile(lag, .999)  | lag < quantile(marg, .001)), 
            vjust=1, 
            aes(label=name), 
            size=2.7)

moran_plot
ggsave ("latex/plots/moran_plot.pdf", moran_plot, dpi=600)

datos_marg$clave_ent = as.factor(datos_marg$clave_ent) 
gmarg_plot = ggplot(datos_marg) + geom_bar(aes(x = gmarg, fill=gmarg)) +scale_fill_discrete("Grado")  +
  facet_wrap(~abr_ent, ncol=4, scales="free_y") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())

gmarg_plot

ggsave ("latex/plots/gmarg_plot.pdf", gmarg_plot,height=400,width=320, units = "mm",scale=.8,dpi=600)

### Top marginación
print(xtable(datos_marg[datos_marg$lugar %in% as.character(1:5),c(3, 15,16,17)]), include.rownames=FALSE)
print(xtable(datos_marg[datos_marg$lugar %in% c("2,452","2,453","2,454","2,455","2,456"),c(3, 15,16,17)]), include.rownames=FALSE)

moran.plot(mapa_mun$imarges, listw_mun)
           

#mapa_mun <- readOGR
gpclibPermit()
mapa_fort_ent<-fortify(mapa_mun,region="CVE_ENT")

mapa_fort_mun<-fortify(mapa_mun,region="CVE_MUN")
# plot(mapa_mun)


### Calcular centroides ###
centroids <- getSpPPolygonsLabptSlots(mapa_mun)
# plot(centroids)


#### Normalizamos datos ####
datos_marg_norm<-apply(as.matrix(datos_marg[,5:13]), 1, function(x){
  x/norm(as.matrix(x), "F")
})
datos_marg_norm <- as.data.frame(t(datos_marg_norm))

#### Determinamos el número de clusters ####
set.seed(8)
nskclust<-clusGap(datos_marg_norm,FUNcluster = skmeans,K.max=10)
print(nskclust,  method="Tibs2001SEmax")
names(nskclust)
xtable(nskclust$Tab, digits=c(1,4,4,4,5))
# k*= 5 
plot(nskclust)
gap_res = as.data.frame(nskclust$Tab)
gap_res$k = 1:dim(gap_res)[1]
gap_plot = ggplot(data = gap_res[3:10,]) + aes(x = k, y = gap) +
  geom_line(color = gg_color_n(2)) + 
  geom_point(size=3, color=gg_color_n(1)) + scale_shape(solid = FALSE) +
  geom_errorbar(aes(ymin= gap-SE.sim , ymax=gap+SE.sim), color=gg_color_n(1), width=0.1) + 
  labs(y = expression(Gap[k])) + scale_x_continuous(breaks=c(3:10))

gap_plot
ggsave("latex/plots/gap_plot.pdf", gap_plot)

##### Corremos k-medias esfércias para 5 grupos #####
set.seed(11)
grupos = skmeans(as.matrix(datos_marg_norm), 5)

cluster_count = summary(as.factor(grupos$cluster))
xtable(data.frame(count=cluster_count, percent=cluster_count/sum(cluster_count)))

clustering_plot = ggplot(data.frame(cluster = as.factor(grupos$cluster))) + geom_bar(aes(x=cluster, fill=cluster) ,alpha=.9)+
  scale_fill_manual("Grupos",values=map_colors)
ggsave("latex/plots/clustering_plot.pdf", clustering_plot)

mapa_mun$g5 <- as.factor(grupos$cluster)
datos_marg$g5 <- as.factor(grupos$cluster)
datos_marg$g5_lab <- as.factor(paste("Grupo", grupos$cluster, sep=" "))

# Permutación de etiquetas para ejemplo
mapa_mun$g5_r <- as.factor(grupos$cluster[r_index])

### Caracterización de los grupos
xtable(as.table(grupos$prototypes))

datos_marg$rank1 = apply(datos_marg_norm, 1,function(x){1-sum(grupos$prototypes[1,] * x)})
datos_marg$rank2 = apply(datos_marg_norm, 1,function(x){1-sum(grupos$prototypes[2,] * x)})
datos_marg$rank3 = apply(datos_marg_norm, 1,function(x){1-sum(grupos$prototypes[3,] * x)})
datos_marg$rank4 = apply(datos_marg_norm, 1,function(x){1-sum(grupos$prototypes[4,] * x)})
datos_marg$rank5 = apply(datos_marg_norm, 1,function(x){1-sum(grupos$prototypes[5,] * x)})

xtable(arrange(datos_marg, rank1)[1:5,c("nom_mun",  "abr_ent", "imarg_esc", "gmarg", "rank1") ], digits=c(0,0,0,2,0,5))
xtable(arrange(datos_marg, rank2)[1:5,c("nom_mun",  "abr_ent", "imarg_esc", "gmarg", "rank2") ], digits=c(0,0,0,2,0,5))
xtable(arrange(datos_marg, rank3)[1:5,c("nom_mun",  "abr_ent", "imarg_esc", "gmarg", "rank3") ], digits=c(0,0,0,2,0,5))
xtable(arrange(datos_marg, rank4)[1:5,c("nom_mun",  "abr_ent", "imarg_esc", "gmarg", "rank4") ], digits=c(0,0,0,2,0,5))
xtable(arrange(datos_marg, rank5)[1:5,c("nom_mun",  "abr_ent", "imarg_esc", "gmarg", "rank5") ], digits=c(0,0,0,2,0,5))

arrange(datos_marg, desc(imarg_esc))[1:5,c("nom_mun", "imarg_esc", "abr_ent", "gmarg", "rank5","ppob_min_5m_h", "g5")]
arrange(datos_marg, imarg_esc)[1:5,c("nom_mun", "imarg_esc", "abr_ent", "gmarg", "rank5","ppob_min_5m_h", "g5")]


gmarg_plot_grupos = ggplot(datos_marg) + geom_bar(aes(x = gmarg, fill=gmarg)) +scale_fill_discrete("Grado")  +
  facet_wrap(~g5_lab, ncol=5) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
gmarg_plot_grupos

ggsave("latex/plots/gmarg_group.pdf",gmarg_plot_grupos)


legend = g_legend(ggplot(datos_marg)+ geom_boxplot(aes(y = ppob_analf, x=g5, fill=g5), alpha=.8)+ 
  scale_fill_manual("Grupos",values=map_colors) +theme(legend.position="bottom"))

densities_by_group = arrangeGrob(
  ggplot(datos_marg)+ geom_boxplot(aes(y = ppob_analf, x=g5, fill=g5), alpha=.8)+theme(legend.position="none", plot.margin=unit(c(0,0,0,0), "cm"))+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("analf") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = ppob_sp, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("sprim") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = pocup_sd_ne, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("sdren") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = pocup_see, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("selec") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = pocup_sae, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("sagua") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = pviv_hac, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("hacina") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = pocup_cpt, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("pisot") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = ppob_min_5m_h, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("pl5khab") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")),
  ggplot(datos_marg)+ geom_boxplot(aes(y = ppob_i2sm, x=g5, fill=g5), alpha=.8)+theme(legend.position="none")+ 
    scale_fill_manual("Grupos",values=map_colors) +ylab("bingreso") + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")), 
ncol=3)

boxplot_bygroup = arrangeGrob(densities_by_group, legend, nrow=2,heights=c(10, 1))
ggsave("latex/plots/boxplot_bygroup.pdf", boxplot_bygroup, scale=1.5)


# Data frame para mapas
datos_marg_gg <- merge(mapa_fort_mun, mapa_mun, by.x = 'id',by.y="CVE_MUN")
format_coord = function(x){
  round(x/100000, 2)
}
######## Gráficas mapas #########
### Índice de marginación ###
  map_imarg <- ggplot(datos_marg_gg, aes(long, lat, group=group)) + geom_polygon(aes(fill=imarges))
  map_imarg <- map_imarg  + scale_fill_gradientn("Índice de\nmarginación\nescalado", colours=brewer.pal(9, "YlOrRd")[1:9]) +
    scale_y_continuous(labels = format_coord) + scale_x_continuous(labels = format_coord) +
    labs(x="long",y="lat") + theme(plot.margin=unit(c(0,0,0,0),"mm")) + coord_equal() #+ coord_map(project="globular") 
  map_imarg  <- map_imarg  + geom_path(data = mapa_fort_mun, colour = "white", size = .01, alpha = .1)
  map_imarg  <- map_imarg  + geom_path(data = mapa_fort_ent, colour = "white", size = .2)
  ggsave ("latex/maps/mapmarg2.png" ,map_imarg, scale=1, dpi=600)


### Mapa de grupos ###
  map_colors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  map <- ggplot(datos_marg_gg, aes(long, lat ,group=group)) + geom_polygon(aes(fill=g5))
  map <- map  + scale_fill_manual("Grupos",values=map_colors) +labs(x="long",y="lat") + 
    scale_y_continuous(labels = format_coord) + scale_x_continuous(labels = format_coord) +
    theme(plot.margin=unit(c(0,0,0,0),"mm")) + coord_equal()#+ coord_map(project="globular") 
  
  # agregar bordes por estado
  map  <- map  + geom_path(data = mapa_fort_mun, colour = "white", size = .01, alpha = .1)
  # agregar bordes por municipio
  map  <- map  + geom_path(data = mapa_fort_ent, colour = "white", size = .2)
  ggsave ("latex/maps/map5g.png", map, scale=1, dpi=600)


### Mapa de índice de marginación permutado  (ejemplo)  ###
  map_imarg_r <- ggplot(datos_marg_gg, aes(long, lat, group=group)) + geom_polygon(aes(fill= imarges_r) )+ 
    scale_fill_gradientn("Índice de\nmarginación\nescalado", colours=brewer.pal(9, "YlOrRd")[1:9]) +
    scale_y_continuous(labels = format_coord) + scale_x_continuous(labels = format_coord) +
    labs(x="long",y="lat") + theme(plot.margin=unit(c(0,0,0,0),"mm")) + coord_equal()  + 
    geom_path(data = mapa_fort_mun, colour = "white", size = .01, alpha = .1) +
    geom_path(data = mapa_fort_ent, colour = "white", size = .2)
  ggsave ("latex/maps/rmapmarg.png" ,map_imarg_r, scale=1, dpi=600)

### Mapa de grupos aleatorizado (ejemplo) ###
  map <- ggplot(datos_marg_gg, aes(long, lat ,group=group)) + geom_polygon(aes(fill=g5_r))
  map <- map  + scale_fill_manual("Grupos",values=map_colors) +labs(x="long",y="lat") + 
    scale_y_continuous(labels = format_coord) + scale_x_continuous(labels = format_coord) +
    theme(plot.margin=unit(c(0,0,0,0),"mm")) + coord_equal()#+ coord_map(project="globular") 
  
  map  <- map  + geom_path(data = mapa_fort_mun, colour = "white", size = .01, alpha = .1)
  # agregar bordes por municipio
  map  <- map  + geom_path(data = mapa_fort_ent, colour = "white", size = .2)
  ggsave ("latex/maps/rmap5g.png", map, scale=1, dpi=600)


############## Autocorrelaciones Espaciales #####################

###Lista de Vecinos###
nb_mun <- poly2nb(mapa_mun,mapa_mun$CVE_MUN)
summary(nb_mun)
listw_mun <- nb2listw(nb_mun,style="W")
names(listw_mun)
listw_mun$neighbours

##### I Moran ######
i.mc <- moran.mc(mapa_mun$imarges, listw_mun, nsim = 9999)
print(i.mc) 
i.density.plot = ggplot() + geom_density(aes(i.mc$res[1:9999]), fill="blue", alpha=".3") + 
  geom_vline(xintercept=i.mc$statistic, color="red", linetype="dashed")+ xlim(min(i.mc$res), max(i.mc$res)*1.1) +
  xlab("simulaciones de Monte Carlo")
ggsave("latex/plots/moran_density.pdf", i.density.plot)

i.r.test <- moran.test(mapa_mun$imarges, listw_mun, randomisation=T)
print(i.r.test)
i.n.test <- moran.test(mapa_mun$imarges, listw_mun, randomisation=F)
print(i.n.test)
i.density.plot = ggplot() + geom_density(aes(i.mc$res[1:9999]), fill="blue", alpha=".3") + 
  geom_vline(xintercept=i.mc$statistic, color="red")+ xlim(min(i.mc$res), max(i.mc$res)*1.1)

i.r.test.example <- moran.test(mapa_mun$imarges_r, listw_mun, randomisation=T)
i.r.test.example

x = seq(.9,1.1,.001)
i.r.test.density.plot = ggplot() + stat_bin(aes(x=i.mc$res[1:9999], y=..density..), fill="blue", colour="black", alpha=".6",binwidth=.005) + 
  geom_vline(xintercept=i.r.test.example$estimate[[1]], color="red", linetype="dashed", size=1)
i.r.test.density.plot
  geom_line(aes(x=x, y=dnorm(x, mean = c.test$estimate[[2]], sd = sqrt(c.test$estimate[[3]]))))

##### C Geary ######
c.mc <- geary.mc(mapa_mun$imarges, listw_mun, nsim = 9999)
print(c.mc)

c.mc.density.plot = ggplot() + geom_density(aes(c.mc$res[1:9999]), fill="red", alpha=".3") + 
  geom_vline(xintercept=c.mc$statistic, color="red", linetype="dashed")+ xlim(min(c.mc$res), max(c.mc$res)*1.1) +
  xlab("simulaciones de Monte Carlo")

ggsave("latex/plots/geary_density.pdf", c.mc.density.plot)
c.mc.density.plot
c.mc$statistic

c.r.test <- geary.test(mapa_mun$imarges, listw_mun, randomisation=T)
print(c.r.test)
c.n.test <- geary.test(mapa_mun$imarges, listw_mun, randomisation=F)
print(c.n.test)
##### Join Counts #####
### 5 Grupos Marginación ###
jc.multi <- joincount.multi(mapa_mun$g5,listw_mun)
print(jc.multi) 
xtable(jc.multi)

jc.test <- joincount.test(mapa_mun$g5, listw_mun)
print(jc.test)
i = 0
jc_df = ldply(jc.test, function(x){
  data.frame( "stddev" = x$statistic, "estadístico" = x$estimate[1], "esperanza" = x$estimate[2], 
              "varianza" = x$estimate[3], "valorp" = x$p.value)
#   x$estimate
#   x$p.value
  })
xtable(jc_df)


jc.mc <- joincount.mc(mapa_mun$g5,listw_mun, nsim = 9999)
dataset_jc_colors = ldply(c(1:5),function(i){
  data.frame(sim=jc.mc[[i]]$res[1:9999], res=rep(jc.mc[[i]]$statistic,9999),grupo=rep(paste("Grupo ",i, sep=""),9999))
})

jc_mc_df = ldply(jc.mc, function(x){
  data.frame( "estadístico" = x$statistic, "esperanza" = x$estimate[1], 
              "varianza" = x$estimate[2], "rank" =jc.mc[[1]]$parameter ,"valorp" = x$p.value)
  #   x$estimate
  #   x$p.value
})
xtable(jc_mc_df)
jc.mc.density.plot = ggplot(dataset_jc_colors, aes(x=sim, fill=grupo)) + geom_density(alpha=.7) +scale_fill_manual("Grupos",values=map_colors)+
  geom_vline(aes(xintercept=res, group=grupo), color="red", linetype="dashed", size=1) +#geom_vline(aes(xintercept=res*1.05, group=grupo), alpha=0) +
  facet_wrap(~grupo,ncol=2)    + xlab("muestra de Monte Carlo")+ ylab("densidad")+ theme(legend.position="none")
jc.mc.density.plot
ggsave("latex/plots/jc.mc.density.plot.pdf", jc.mc.density.plot)

