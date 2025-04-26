rm(list=ls())
library(dplyr)
library(ComplexHeatmap)

dt = read.table("./subfamily.host.phylum-genus", sep="\t")
genecount = read.table("./subfamily.gene.count.f", sep="\t", row.names=1)
kegg = read.table("./subfamily.kegg.count.metabolism.f.f",sep="\t", quote = "", comment.char = "")
virus.map = read.table("../00.data/virus.group", sep="\t", header = T)

## lifestyle
subfamily.style = virus.map %>% filter(subfamily != "") %>% group_by(subfamily, lifestyle) %>% count() %>% dcast(formula = subfamily ~ lifestyle, value.var="n")
subfamily.style[is.na(subfamily.style)] = 0
subfamily.style <- subfamily.style %>%
  mutate(tot = rowSums(select(., Mix, Temperate, Virulent)),
         mix = Mix/tot*100,
         temp = Temperate / tot*100,
         vir = Virulent / tot * 100)
rownames(subfamily.style) = subfamily.style$subfamily


## kegg rate
genecount = genecount[kegg$V1, ]
kegg$rate = kegg$V3 / genecount

dtf = dt %>%
  mutate(rate = V3/V4) %>%
  filter(rate > 0.5, V4>80, V2 != "multiple") %>%
  extract(V2,c("phylum","genus"),"(.*)xxxxxx(.*)") %>%
  mutate(assign = V3, other = V4-V3)

dtf <- dtf %>%
  merge( dtf %>% group_by(phylum) %>% summarise(pc = sum(V3)), by='phylum') %>%
  merge( dtf %>% group_by(genus) %>% summarise(gc = sum(V3)),by='genus') %>%
  arrange(desc(pc), desc(gc), desc(V3))

top_anno = HeatmapAnnotation(
  # lifestyle = anno_customize(rep("pie",86), graphics = graphics),
  lifestyle = anno_barplot(subfamily.style[dtf$V1,c("mix","temp",'vir')], gp = gpar(fill=colors.lifestyle), which = 'row', bar_width=1),
  
  Novotu = anno_barplot(dtf[,c("V3","other")], gp = gpar(fill = colors.vir.rate ), which='row', bar_width = 1),
  
  genus = dtf$genus,
  genus_name = anno_text(dtf$genus),
  phylum = dtf$phylum,
  phylum_name = anno_text(dtf$phylum),
  
  
  annotation_height = unit(c(4, 4,0.5,4,0.5,4), "cm"),
  
  col = list(
    genus = colors.taxo,
    phylum = colors.taxo,
    lifestyle = colors.lifestyle,
    virus_rate = colors.vir.rate
  )
)

draw(top_anno)

align_heat <- function(heat_matrix, select_names){
  no = setdiff(select_names, colnames(heat_matrix))
  heat_matrix[,no] = 0
  heat_matrix = heat_matrix[,select_names]
  heat_matrix[is.na(heat_matrix)] = 0
  heat_matrix[rowSums(heat_matrix)!=0,]
}

heat_matrix = dcast(kegg, V2 ~ V1, value.var="rate")
rownames(heat_matrix) = heat_matrix$V2; heat_matrix = heat_matrix[,-1]
heat_matrix = align_heat(heat_matrix, dtf$V1)


heat_matrix = apply(heat_matrix, 2, function(x)(x/sum(x)))
heat_matrix[is.na(heat_matrix)] = 0

heat_matrix = heat_matrix[apply(heat_matrix, 1, max) > 0.5, ]

keggf = kegg %>% filter(V2 %in% rownames(heat_matrix)) %>% select(V2,V4) %>% unique()
rownames(keggf) = keggf$V2


mycol = colorRamp2(breaks =c(0, 0.2, 0.4, 0.6, 0.8, 1 ), colors=c("white", "#f0f9e8", "#a8ddb5", "#7bccc4", "#4eb3d3", "#08589e" ))

nh = nrow(heat_matrix)
nw = ncol(heat_matrix)
myunit = 5
max_nchar = max(nchar(rownames(heat_matrix)))


ha <- Heatmap(heat_matrix,
              
              row_split = as.factor(keggf[rownames(heat_matrix),"V4"]),
              row_title_rot = 0,
              show_row_dend = F,
              top_annotation = top_anno,
              row_names_gp = gpar(fontsize=8),
              
              clustering_distance_row = "binary",
              column_split = as.factor(dtf$phylum),
              
              cluster_columns = F,
              col = mycol,
              border = T,
              rect_gp = gpar(col = "grey", lwd = 0.1), # 内部线条颜色
              height = nh * unit(myunit,"mm"), width = nw * unit(myunit,"mm"), # 保持单元格是方的
              row_names_max_width = unit(max_nchar, "char"), # 有时候legend和标签文字会重叠
              
              
)


