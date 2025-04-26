rm(list=ls())
load(file="adonis.Rdata")
library(dplyr)      
library(ggplot2)   
library(RColorBrewer) 

res_plot <- filter(res_new,name %in% adonis_sig)
res_p <- merge(res_plot,vartype,by.x="name",by.y="variable") %>% arrange(subset,adj2)
res_p <- res_p %>% mutate(ptext = ifelse(pvalue < 0.001, "***",
                                           ifelse(pvalue < 0.01, "**",
                                                  ifelse(pvalue < 0.05, "*", " "))))
res_p$adj2 <- ifelse(res_p$adj2<0,0,100*(res_p$adj2))
  
res_p$name <- gsub("_"," ",res_p$name)
res_p$subset <- gsub("_"," ",res_p$subset)

res_p$name <- factor(res_p$name,levels=rev(unique(res_p$name)))
name_order <- rev(unique(res_p$name))
res_p$subset <- factor(res_p$subset,levels=unique(res_p$subset))
subset_order <- unique(res_p$subset)
  
sample.color <- brewer.pal(length(unique(res_p$subset)),"Paired")
  
ado_plot <- ggplot(data=res_p, aes(x= adj2, y=name, color= subset, fill= subset))+
    geom_bar(stat='identity',width=0.8)+ 
    ylab("")+xlab(expression("Effect size (adjusted R"^2*",%)")) + theme_bw()+
    geom_text(aes(label=ptext,x=0.4),nudge_y=-0.3,size=5, color="black")+
    coord_cartesian(clip = 'off') + # 保证文本不被裁剪
    scale_color_manual(values=sample.color,limits=subset_order)+
    scale_fill_manual(values=sample.color,limits=subset_order)
  

