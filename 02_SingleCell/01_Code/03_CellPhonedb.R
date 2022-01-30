
pbmc_assays_data <- as.data.frame(pbmc@assays[["RNA"]]@data)

write.table(as.matrix(pbmc@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)

All_Cell <- colnames(pbmc_assays_data)
All_Cell_Type <- result_main_hpca$HPCA_Main

meta_data <- cbind(All_Cell, All_Cell_Type)
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA

write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)


mynet <- read.delim("count_network.txt", check.names = FALSE)
table(mynet$count)
mynet %>% filter(count>0) -> mynet  # 有零会报错
head(mynet)

net<- graph_from_data_frame(mynet)


plot(net)


allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局

E(net)$width  <- E(net)$count/10  # 边点权重（粗细）


plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7)


net2 <- net  # 复制一份备用

for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
}  # 这波操作谁有更好的解决方案？ 

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 


plot(net, edge.arrow.size=.1, 
     edge.curved=0.2, # 只是调了这个参数
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7)

# dev.off()
length(unique(mynet$SOURCE)) # 查看需要绘制多少张图，以方便布局
par(mfrow=c(4,5), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1) 
  
}




length(unique(mynet$SOURCE))
# par(mfrow=c(7, 3), mar=c(.5, .5, .5, .5))
# par(mfrow=c(1, 1), mar=c(1, 1, 1, 1))
for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  # 故技重施
  
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  file <- paste0(i, "---", unique(mynet$SOURCE)[i], ".png")
  # png(file, width = 1200, height = 1200)
  # png(file)
  par(mar=c(1, 3.5, 1, 5.5))
  plot1 <-plot(net1, edge.arrow.size=1, 
               edge.curved=1,
               edge.label = E(net1)$count, # 绘制边的权重
               vertex.color=allcolour,
               vertex.frame.color="#555555",
               vertex.label.color="black",
               layout = coords,
               vertex.label.cex=2)
  # 
  print(plot1)
  # plot(net1, edge.arrow.size=1,
  #      edge.curved=1,
  #      edge.label = E(net1)$count, # 绘制边的权重
  #      vertex.color=allcolour,
  #      vertex.frame.color="#555555",
  #      vertex.label.color="black",
  #      layout = coords,
  #      vertex.label.cex=2)
  # dev.off()
}
par(mar=c(1, 3.5, 1, 5.5))
plot(net1, edge.arrow.size=1, 
     edge.curved=.5,
     edge.label = E(net1)$count, # 绘制边的权重
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=2)

