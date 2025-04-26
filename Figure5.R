rm(list=ls())
source("AVGC_function.R")
library(pheatmap)   
library(tidyverse)  
load("maaslin/Rdata")

op.obs <- alpha_func(dt=opvir, sample_map=op.map, group="Group", ID="Sample", 
                      index="obs",
                      sample.color=color.op, 
                      box_width=0.5, 
                      title= " ", 
                      violin = F)
op.pcoa.sd <- pcoa_se_sd(dt=opvir, sample_map=opmap, group="Group", ado_group="Visit",ID="Sample", sample.color=color.op,
                         ado_method="bray", pca_method="bray",
                         err_bar_type="se",
                         cut_rate = NA, cut_num=NA,
                         err_width=0.01,
                         title=" ", x=1, y=2)

#heatmap
res_plot<- arrange(res_plot,qval)
dt.p <- op.vir %>% .[unique(res_plot$name),]

#样本顺序
meta.p <- op.map %>% arrange(Visit,Donator)
meta.p$Visit <- factor(meta.p$Visit,levels=unique(meta.p$Visit))
dt.p <- dt.p[,meta.p$Sample]

# 显著性标记
sig <- res_plot %>%
  dplyr::select(name, qval,value) %>%
  mutate(mark = ifelse(qval < 0.001, "***", ifelse(qval < 0.01, "**", ifelse(qval < 0.05,"*",NA))))

sig.mark <- sig %>%
  dplyr::mutate(across(everything(), ~ ifelse(. < 0.001, "***", ifelse(. < 0.01,"**",ifelse(. < 0.05, "*", " ")))))

sig.mark[is.na(sig.mark)] <- " "
sig.plot <- sig.mark[rownames(dt.plot),colnames(dt.plot)]  %>% as.matrix()

#breaks,0为白色
breaks <- c(seq(-5, -0.99, length.out = 100), 0, seq(0.01, 5, length.out = 100))

# 热图颜色设置
col_fun <- colorRampPalette(c("#b5c9e2", "white", "#f4978b"))(200)

# 绘制热图
p1 <- pheatmap(
  dt.p,
  color = col_fun,
  display_numbers= sig.plot,    
  fontsize_number = 12,        
  number_color = "black",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  angle_col = "45",
  fontsize = 15,
  legend_breaks = c(-5, 0, 5), 
  breaks = breaks,
  cellwidth = 30,
  silent = T,
  heatmap_legend_param = list(
    title="log2(Fold change)"
  )
)