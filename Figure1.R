rm(list=ls())
library(data.table)
library(ggplot2)
library(dplyr)
source("AVGC_function.R")

colors.complete = structure(c("#3bc9db", "#f7b731","#d9d9d9"),
                            names=c("high_complete", "low_complete","others"))

dt <- fread("virus.ckv",sep="\t")
p1 <- ggplot(dt, aes(x = V2, color=type, fill = type)) +  
  geom_histogram(position = "stack", bins = 100,alpha=0.8) + 
  scale_fill_manual(values = colors.complete) +   
  scale_color_manual(values = colors.complete) +
  scale_x_log10()+
  theme_bw()
p1


pie_dt = dt %>%
  group_by(type) %>%
  count()
p2 <- pie_plot(pie_dt, value="n", fill="type", col=colors.complete)


raredt <- fread("virus.clu.count.rare")
raredt.ns <- fread("virus.clu.count.no_single.rare")

raredt1 <- raredt %>%  
  group_by(V2) %>%                       
  summarize(Average_V3 = mean(V3))
raredt1$Type <- "rare"

raredt.ns1 <- raredt.ns %>%  
  group_by(V2) %>%                       
  summarize(Average_V3 = mean(V3))
raredt.ns1$Type <- "no_single_rare"

cob <- rbind(raredt1,raredt.ns1)

p3 <- ggplot(cob, aes(x = V2, y = Average_V3, color = Type)) +  
  geom_line(size = 0.8) +            
  theme_bw()
p3

