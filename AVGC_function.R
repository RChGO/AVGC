library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(vegan)
library(ggpubr)
library(ggrepel)

pie_plot <- function(dt, value, fill, facet_my=NULL, col=NULL){
  total_color1 = c(brewer.pal(12,"Set3"), brewer.pal(12,"Paired"))
  
  if(typeof(col) == "NULL"){
    if(length(unique(dt[,fill])) > length(total_color1)){
      message("ERROR!!!\n分类太多，请不要超过默认的24个")
      exit(1237)
    }else{
      col = total_color1[1:length(unique(dt[,fill]))]
    }
  }
  
  if (typeof(facet_my) == "NULL"){
    data = dt[, c(fill,value)]
    colnames(data) = c("fill", "value")
    data = data[order(data$value, decreasing=T), ]
    
    data$fill = factor(data$fill, levels=unique(data$fill))
    
    ss2 = sum(data$value)
    plot_dt <- data %>%
      mutate(
        Perc =  round(value/ss2, digits=4) ) %>% # 分组计算百分比,两位小数
      mutate(
        label = paste0(fill, " ,", value,", ", round((Perc)*100, digits = 4), "%")
        ,ypos = cumsum(Perc) - 0.5 * Perc
        ,wght=runif(length(fill))
        ,wght=wght/sum(wght)
        ,wght=round(wght, digits=2)
      )
  }
  ggplot(plot_dt, aes(x = "", y = Perc , fill = fill)) +
    geom_bar(stat = "identity", width = 1 , color = "white", show.legend = T) +
    geom_text_repel(aes(x = 1.5,y=1-ypos, label = label),
                    ,color="black"
                    , nudge_x = 0.3
                    ,hjust=0
                    ,size = 3
                    , segment.color = "gray50",
                    lty="dashed",force=2,
                    min.segment.length = 0,
                    box.padding = 0.5 ,
                    segment.size = .2) + 
    coord_polar("y", start = 0) +
    theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.key = element_rect(fill = 'transparent'), 
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())+
    scale_fill_manual(values=col)
  
}

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

zy_alpha = function(dt=NA, sample_map=NA, group="Group", ID="Sample", # 必须参数
                    index="shannon", # 计算参数
                    sample.color=NA, # 美化参数
                    box_width=0.5, # 箱式图宽度
                    title="alpha diversity", # 文字参数,
                    violin = F
){
  ## colors 
  if (any(is.na(sample.color))){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  message(paste(length(sample.color), "of groups to plot"))
  
  ## align dt and group
  inter <- intersect(colnames(dt),sample_map[,ID])
  if (length(inter) == 0) {  
    stop("Error: No overlapping columns found between dt and sample_map.")  
  } 
  dt <- dt[,colnames(dt) %in% inter]
  dt = dt[rowSums(dt)!=0,]
  sample_map <- filter(sample_map, sample_map[,ID] %in% inter)
  sample_map <- filter(sample_map, !is.na(!!sym(group)))
  
  #alpha
  if(tolower(index) == "obs"){
    alpha = data.frame(alpha=colSums((dt>0)+0))
  }else{
    alpha = data.frame(alpha = vegan::diversity(t(dt),index=index))
  }
  
  dm = merge(alpha,sample_map, by.x='row.names', by.y=ID)
  comp = combn(as.character(unique(dm[,group])),2,list)
  
  p = ggplot(dm, aes(x=.data[[group]], y=alpha,fill=.data[[group]]))
  if(isTRUE(violin)){
    p <- p+
      geom_violin()+
      geom_boxplot(width=box_width, fill="white",
                   position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.colour = NA)
  }else{
    p <- p+ 
      geom_boxplot(position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.color="#c1c1c1")
  }
  
  ylabs = structure(c("Number of OTUs","Shannon index", "1 - Simpson index", "Invsimpson index"),
                    names=c("obs", "shannon", "simpson","invsimpson"))
  ylab = ylabs[tolower(index)]
  
  
  p <- p+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values=sample.color)+
    #geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=sigFunc.show)+
    geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1)+
    labs(title=title, y = ylab, x=NULL)
  
  return(p)
}

pcoa_se_sd <- function(dt=NA, sample_map=NA, group=NA, ado_group=NA,ID=NA, sample.color=NULL,
                          ado_method="bray", pca_method="bray",
                          err_bar_type="se",
                          cut_rate = NA, cut_num=NA,
                          err_width=0.01,
                          title="PCoA", x=1, y=2){
  if(is.finite(cut_rate) || is.finite(cut_num)){
    sample_map = filter_group(sample_map, group=group, cut_rate = cut_rate, cut_num=cut_num)
  }
  
  ## colors 
  if ( typeof(sample.color) == "NULL" ){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  # 统计每个分组各有多少,作为新的图例
  group_summ <- sample_map %>%
    dplyr::select(all_of(group)) %>%
    dplyr::group_by(across({{group}})) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(new_label=paste(!!sym(group), " (", count, ")", sep=""))
  
  new_label <- structure(group_summ$new_label,names=as.character(unlist(group_summ[,group])))
  
  message(paste(length(unique(sample_map[,group])), "of groups to plot"))
  
  otu.dist = vegdist(t(dt), method = pca_method) # calc dist matrix
  
  if(length(unique(sample_map[,ado_group])) > 1){
    ## adonis
    ado = adonis2(otu.dist~sample_map[,ado_group], method = ado_method,na.action = na.omit)
    ado_r2 = round(ado$R2[1], digits = 4)
    ado_p = ado$`Pr(>F)`[1]
  }else{
    ado_r2 = NA
    ado_p = NA
  }
  
  ## PCoA
  pcoa = cmdscale(otu.dist, k=10, eig=T)
  eigs = signif(pcoa$eig/sum(pcoa$eig), 4)*100
  point = pcoa$points
  
  colnames(point) = paste("pcoa.", 1:ncol(point),sep="")
  
  xlab = paste("PCoA", x, " (",eigs[x],"%)", sep="")
  ylab = paste("PCoA", y, " (",eigs[y],"%)", sep="")
  
  substitle <- paste0("'R'^2~'='~'", ado_r2, "'~~italic('p')~'='~'", ado_p, "'") %>% 
    as.formula() %>% 
    eval()
  group = ifelse(length(unique(sample_map[,group])) == 1, sample.color[1], group)
  point = point[sample_map[,ID], ]
  dm = merge(point, sample_map, by.x='row.names', by.y=ID)
  tmp_mean = aggregate(point, by=list(c(sample_map[,group])), mean)
  tmp_sd  = aggregate(point, by=list(c(sample_map[,group])), sd)
  tmp_sd  = aggregate(point, by=list(c(sample_map[,group])), sd)
  
  
  plot_data = do.call(data.frame, aggregate(point, by=list(c(sample_map[,group])), 
                                            FUN = function(x){
                                              c(mean = mean(x), sd = sd(x), se = sd(x)/sqrt(length(x)))
                                            }))
  
  plot_x = paste("pcoa.", x, ".mean", sep="")
  plot_y = paste("pcoa.", y, ".mean", sep="")
  
  x_offset = paste("pcoa.", x, ".", err_bar_type, sep="")
  y_offset = paste("pcoa.", y, ".", err_bar_type, sep="")
  
  p1 <- ggplot(data=plot_data, aes(x = .data[[plot_x]], y = .data[[plot_y]], color = Group.1))+
    geom_point(size=10) +
    geom_errorbar(aes( xmin = .data[[plot_x]] - .data[[x_offset]],
                       xmax = .data[[plot_x]] + .data[[x_offset]]),
                  width=err_width
    )+
    geom_errorbar(aes( ymin = .data[[plot_y]] - .data[[y_offset]],
                       ymax = .data[[plot_y]] + .data[[y_offset]]),
                  width=err_width
    )+
    theme_bw()+
    theme(panel.grid = element_blank(),
          text = element_text(color="black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color="black", linewidth=0.25),
          panel.border = element_rect(colour="black", linewidth=0.25))+
    scale_fill_manual(values=sample.color, guide="none")+
    scale_color_manual(values=sample.color, labels=new_label)+
    labs(x=xlab, y=ylab, title=title, subtitle = substitle)
  list(plot=p1, new_label=new_label)
}

