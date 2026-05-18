library(ggplot2)
library(ggrepel)  #用于添加样本标签
library(maps)  #用于载入一个世界地图模板

#读取示例的站点经纬度坐标
site <- read.delim('533sample_hqmaginfo.txt')
my_colors_ordered <- c(
  "#A7C0DE",  # 浅蓝灰
  "#729ECE",  # 蓝紫
  "#6C91C2",  # 中蓝
  "#403990",  # 深蓝紫
  "#9478B5",  # 柔和紫
  "#DBA9A9",  # 灰粉红
  "#C77D8C",  # 粉褐
  "#A4514F",  # 棕红
  "#CF3D3E",  # 鲜红
  "#F46F43",  # 橙红
  "#E56A54",  # 番茄橙
  "#FF9896",  # 粉红
  "#F9C784",  # 奶油橙
  "#D59730",  # 金黄
  "#C59D94",  # 褐粉
  "#A3C294",  # 橄榄绿
  "#3DA284",  # 深青绿
  "#1AC990",  # 青绿色
  "#88B2C4",  # 蓝灰
  "#82716D"   # 中性灰棕
)

#绘制世界地图的简图

p <- ggplot() +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +  #世界地图模板
  theme_bw() +  #以下用于修改主题、配色、坐标轴等
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(-150, -100, -50, 0, 50, 100, 150), expand = c(0, 0), 
                     labels = c('150°W', '100°W', '50°W', '0', '50°E', '100°E', '150°E')) +
  scale_y_continuous(breaks = c(-60, -30, 0, 30, 60), expand = c(0, 0), 
                     labels = c('60°S', '30°S', '0', '30°N', '60°N')) +
  labs(x = 'Longitude', y = 'Latitude', color = 'Order')

p <- ggplot() +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +  #世界地图模板
  theme_bw() +  #以下用于修改主题、配色、坐标轴等
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(-150, -100, -50, 0, 50, 100, 150), expand = c(0, 0),
                     labels = c('150°W', '100°W', '50°W', '0', '50°E', '100°E', '150°E')) +
  scale_y_continuous(breaks = c(-60, -30, 0, 30, 60), expand = c(0, 0),
                     labels = c('60°S', '30°S', '0', '30°N', '60°N')) +
  labs(x = 'Longitude', y = 'Latitude', color = 'Order') +
  theme(axis.title.x=element_text(vjust=1,size=15),axis.title.y=element_text(vjust=1,size=15), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),)

p

#绘制采样站点的地图
p1 <- p +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +  #世界地图模板
  geom_point(data = site, aes(x = Longitude, y = Latitude, color = Order), size = 4) +  #将数据中的样本点绘制在地图中
  #geom_label_repel(data = site, aes(x = Longitude, y = Latitude, label = id), size = 2, show.legend = FALSE) +  #添加样本标签
  scale_color_manual(values = c("blue",  "black", "green", "brown", "#FF1493", "pink", "red","purple",'skyblue','#8BC847'))

p1 <- p +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +  #世界地图模板
  geom_point(data = site, aes(x = Longitude, y = Latitude, color = Order), size = 4) +  #将数据中的样本点绘制在地图中
  scale_color_manual(values = my_colors_ordered) + theme(
    legend.title = element_text(size = 15, face = "plain", family = "Arial"),  # 加粗，使用衬线字体
    legend.text = element_text(size = 12, family = "Arial")                  # 使用衬线字体
  )

p1
ggsave(filename = "采样点2.pdf",plot = p1,width = 12,height = 6,)
ggsave(filename = "采样点2.png",plot = p1,width = 12,height = 6,)


p2 <- p +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +  #世界地图模板
  geom_point(data = site, aes(x = Longitude, y = Latitude), color = 'black') +  #将数据中的样本点绘制在地图中
  geom_label_repel(data = site, color = 'white', segment.colour = 'black', 
                   aes(x = Longitude, y = Latitude, label = id, fill = Order),max.overlaps = 200) +  #添加样本标签
  scale_fill_manual(values = c("blue",  "black", "green", "brown", "#FF1493", "pink", "red","purple",'skyblue','#8BC847'))+
  guides(fill = guide_legend(override.aes = list(label = "")))  # 设置图例颜色标签为空
p2 <- p +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +  #世界地图模板
  geom_point(data = site, aes(x = Longitude, y = Latitude), color = 'black') +  #将数据中的样本点绘制在地图中
  geom_label_repel(data = site, color = 'white', segment.colour = 'black', size = 5,
                   aes(x = Longitude, y = Latitude, label = id, fill = Order),max.overlaps = 200) +  #添加样本标签
  scale_fill_manual(values = my_colors_ordered)+
  guides(fill = guide_legend(title.theme = element_text(size = 20),override.aes = list(label = ""))) +
  theme(legend.key.size = unit(0.62, "cm"),
        legend.text = element_text(size = 15, face = "italic"))
p2 <- p +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +  #世界地图模板
  geom_point(data = site, aes(x = Longitude, y = Latitude), color = '#4682B4',size=4,alpha=0.4) +  #将数据中的样本点绘制在地图中
  guides(fill = guide_legend(title.theme = element_text(size = 20),override.aes = list(label = ""))) +
  theme(legend.key.size = unit(0.62, "cm"),
        legend.text = element_text(size = 15, face = "italic"))
p2 <- p +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +   # 世界地图
  geom_point(data = site, aes(x = Longitude, y = Latitude, color = factor(ifhqmag)),size = 4, alpha = 0.6) +             # 根据 ifhqmag 分组上色
  scale_color_manual(values = c("0" = "#4682B4", "1" = "#CD5C5C")) +
  guides(fill = guide_legend(title.theme = element_text(size = 20),override.aes = list(label = ""))) +
  theme(legend.key.size = unit(0.62, "cm"),
        legend.text = element_text(size = 15, face = "italic")) +
  theme(legend.position = "none") 


p2
ggsave(filename = "采样点.pdf",plot = p2,width = 12,height = 6,)
ggsave(filename = "采样点.png",plot = p2,width = 12,height = 6,)
#如要从中截取部分展示
#比如说图中没有东太平洋的采样点，因此只展示 180°W~100°E 的区域
p1 + coord_cartesian(xlim = c(-180, 100))
