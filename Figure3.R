rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)

load("../00.data/virus.profile.relative.RData")
map_vir = read.table("../00.data/virus.group", sep="\t", header=T, check.names=F)
load("../00.data/colors.RData")


plot_dt = data.frame(preva = rowSums(dt.virus>0), mean = rowMeans(dt.virus)) %>%
  merge(map_vir, by.x='row.names', by.y='virus')

plot_dt$preva_rate = plot_dt$preva / ncol(dt.virus)

ggpubr::show_point_shapes()




x = plot_dt %>%
  filter(subfamily.type != "") %>%
  mutate(x=mean*preva) %>%
  group_by(subfamily.type) %>%
  summarise(xx=sum(x)) %>%
  # arrange(desc(xx)) %>%
  slice_max(xx, n=20)

colors.fill = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#e31a1c", "#fdbf6f", "#ff7f00", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f")
names(colors.fill) = x$subfamily.type

map_vir %>% filter(subfamily %in% x$subfamily.type)

p1 <- ggplot(plot_dt, aes(x=preva_rate, y=log10(mean), size=mean*preva))+
  geom_point(aes(fill=subfamily.type), shape=21, color="transparent")+
  geom_vline(xintercept = c(0.1,0.5), lty="dashed")+
  theme_bw()+
  scale_fill_manual(values=colors.fill, na.value="grey90")+
  xlim(c(0,1))
p1


xp <- ggplot(plot_dt, aes(x=preva_rate))+
  geom_histogram(bins=100, color="black", fill="#a6cee3")+
  theme_bw()+
  xlim(c(0,1))

p <- ggarrange(plotlist=list(xp,p1), ncol=1, nrow=2, heights=c(1,3), common.legend = T, legend = "right")

p

# ggsave("scatter_prevalence.pdf", p, width=8, height=8)


tmp_func <- function(inf){
  # inf = "../10.virus.host/subfamily.host.phylum"
  dt1 = read.table(inf, sep="\t") %>% mutate(rate = V3 / V4)
  
  ## 删选以某一个分类为宿主的subfamily
  x1 = dt1 %>% filter(rate >= 0.5 & V2 != "multiple")
  
  ## 计算每个subfamily总共有多少vOTU有宿主
  x2 = dt1 %>%
    group_by(V1) %>%
    summarise(
      tot.host = sum(V3),
      tot.votu = mean(V4)
    )
  
  
  x3 = merge(x1, x2, by='V1', all.y=T) %>%
    mutate(
      V2 = ifelse(is.na(V2), "other", V2),
      V3 = ifelse(is.na(V3), 0, V3),
      V4 = tot.votu,
      rate = V3 / V4,
      other = tot.host - V3,
      fill = ifelse(V2 == "other","other", "current"),
      group = V2,
      no.host = tot.votu - tot.host
    )
  rownames(x3) = x3$V1
  x3
}



tmp2_func <- function(x, y){
  x = x[y$subfamily.type, ] %>% mutate(y=1, x=1:nrow(.))
  
  ggplot()+
    geom_scatterpie(data=x,
                    aes(x, y, r=0.5), cols = c("V3","other","no.host"), linewidth=0.1)+
    geom_text(data=x, mapping=aes(x,y, label=paste(V1, V2)), nudge_y =-1, angle=90, hjust=1)+
    coord_equal()+
    theme_void()+
    theme(legend.position = "none")+
    scale_fill_manual(values=c("V3"="#66ccff", "other" = "gray","no.host"="white"))+
    scale_y_continuous(limits=c(-10,5))
}

phylum = tmp_func("../10.virus.host/subfamily.host.phylum")
genus = tmp_func("../10.virus.host/subfamily.host.genus")

p1 <- tmp2_func(phylum, x)
p2 <- tmp2_func(genus, x)


