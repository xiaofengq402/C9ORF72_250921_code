rm(list = ls())

rm(list = ls())
setwd("I:/Analysis/Python_code/results")

deg <- read.csv("ols_results_20250103_192328.csv",row.names="Column")
deg <- tibble::rownames_to_column(deg, "Column")

dim(deg)
deg <- deg[order(deg$p_value_age, decreasing = F),]
library(ggplot2)
library(ggrepel)
library(tidyverse)
##看一下上下调基因的数量
deg$Change = as.factor(ifelse(deg$age> 0.5 & abs(deg$p_value_age) <0.05,
                              
                              ifelse(deg$age > 0.5 ,'UP','DOWN'),'STABLE'))
## 定义校正p值<0.05和差异倍数大于动态阈值的结果
table(deg$Change) #看一下上下调基因的数量


## 定义一个标题，加上我们想显示具体的上下调的基因数。
this_tile <- paste0('Volcano plot',
                    '\nCutoff for FoldChange is ', '2',
                    '\nThe number of up genes is ',nrow(deg[deg$Change =='UP',]) ,
                    '\nThe number of down genes is ',nrow(deg[deg$Change =='DOWN',])
)
## 画一个没有基因标签的火山图
ggplot(deg, aes(x=age, y=-log10(p_value_age),color=Change)) + 
  geom_point(alpha=1, size=2) +  # 设置点的透明度和大小
  theme_bw(base_size = 12) +  #设置一个主题背景
  xlab("Log2(Fold change)") + # x轴名字
  ylab("-Log10(Pavlue)") + # y轴名字
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('steelblue','gray','brown')) + # 各自的颜色
  geom_hline(yintercept = -log10(0.05), lty = 4) + #定义p值和线形
  geom_vline(xintercept = c(-1, 1), lty = 4)+ #定义差异倍数和线形
  labs(title = this_tile) + #加上题目
  xlim(0,1)+
  ylim(0,5)

### 筛选出需要标注的基因
#rank前30
deg$rank <- rank(deg$P.value, ties.method = 'min')
demo <- subset(deg, rank <= 20 & Change != "STABLE")
dim(demo)


# ### ISG基因和趋化因子
# demo <- rbind(demo, subset(deg, SYMBOL %in% c("Gzmb", "Ifit1", "Rsad2", "Ifitm1",
#                                               "Zbtb16", "Ifi44", "Plin1", "Il11",
#                                               "CCL2", "Cxcl2", "Cxcl13"))
# )
# dim(demo)


### 选出C9ORF72
select <- c('C9orf72')
demo <- subset(deg, Column %in% select)

dim(demo)

### 去掉没有显著差异的
demo <- subset(demo, Change != "STABLE")
dim(demo)

demo_down <- subset(demo, Change == "DOWN")
demo_up <- subset(demo, Change == "UP")



library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# 创建基础图形
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# 创建基础图形
first_plot <- ggplot(deg, aes(x = age, y=-log10(p_value_age))) + 
  # 使用深灰色边框的点
  geom_point(aes(color = Change), alpha = 0.7, size = 3, shape = 21, fill = "white", stroke = 1.5, color = "darkgrey") + 
  geom_point(data = demo_down, size = 4, color = "red", alpha = 0.9) +  # 用红色标注的主要数据点
  geom_point(data = demo_up, size = 4, color = "red", alpha = 0.9) + 
  
  # 移除填充区域
  # annotate("rect", xmin = -Inf, xmax = 0.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") + 
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, alpha = 0.2, fill = "grey") + 
  # annotate("rect", xmin = 0.5, xmax = Inf, ymin = 0.5, ymax = Inf, alpha = 0.2, fill = "orange") +
  
  # 主题和标签设定
  theme_bw(base_size = 16) + 
  labs(
    x = "表达与年龄的相关性", 
    y = "表达与SASP的相关性"
  ) +
  theme(
    plot.title = element_blank(),  # 去掉标题
    legend.title = element_blank(),
    legend.position = "top",  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 1),  
    axis.ticks.length= unit(0.25, "cm")
  ) + 
  scale_color_manual(values = brewer.pal(3, "Set1")) + 
  
  coord_fixed(ratio = 0.2) + 
  geom_text_repel(data = demo, aes(label = Column),
                  size = 4,  
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = "black",
                  direction = "both",
                  color = "red",  # 标注点的文本颜色
                  show.legend = FALSE)

# 输出图形
print(first_plot)

# 第二个图的增强
library(ggplot2)

# 为数据分组
library(ggplot2)

# 为数据分组
library(ggplot2)

# 为数据分组
deg$ColorGroup <- with(deg, ifelse(age > 0.5 & sasp > 0.5, "Both High",
                                   ifelse(age > 0.5 & sasp <= 0.5, "Age High",
                                          ifelse(age <= 0.5 & sasp > 0.5, "SASP High", "Low"))))

second_plot <- ggplot(deg, aes(x = age, y = sasp, color = ColorGroup, shape = ColorGroup)) +
  geom_point(alpha = 0.5, size = 4) +
  scale_color_manual(values = c("Both High" = "red", "Age High" = "orange", "SASP High" = "blue", "Low" = "gray")) +
  scale_shape_manual(values = c(16, 16, 16, 16)) + # 为每个组使用不同的形状
  labs(title = "Figure", x = "Correlation of expression and age", y = "Correlation of expression and SASP") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "top"
  )
library(ggplot2)

# 为数据分组
deg$ColorGroup <- with(deg, ifelse(age > 0.5 & sasp > 0.5, "Both High",
                                   ifelse(age > 0.5 & sasp <= 0.5, "Age High",
                                          ifelse(age <= 0.5 & sasp > 0.5, "SASP High", "Low"))))

# 绘制散点图
second_plot <- ggplot(deg, aes(x = age, y = sasp, color = ColorGroup, shape = ColorGroup)) +
  geom_point(alpha = 0.5, size = 4) +  # 普通点
  scale_color_manual(values = c("Both High" = "red", 
                                "Age High" = "orange", 
                                "SASP High" = "blue", 
                                "Low" = "gray")) +
  scale_shape_manual(values = c(16, 16, 16, 16)) +  # 使用统一形状
  labs(title = "分组颜色的散点图", x = "Correlation of expression and age", y = "Correlation of expression and SASP_score") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

# 筛选 C9ORF72 基因数据
c9orf72_data <- deg[deg$Column == "C9orf72", ]

# 单独标注 C9ORF72 基因
second_plot <- second_plot + 
  geom_point(data = c9orf72_data, 
             aes(x = age, y = sasp), 
             color = "red",  # 更深的颜色
             size = 5,  # 更大的点
             shape = 16) +  # 自定义形状
  geom_segment(aes(x = c9orf72_data$age, y = c9orf72_data$sasp, 
                   xend = c9orf72_data$age, yend = c9orf72_data$sasp + 0.2), 
               color = "black", 
               size = 0.5) +  # 线段的宽度
  geom_text(aes(x = c9orf72_data$age, y = c9orf72_data$sasp + 0.2, 
                label = "C9ORF72"), 
            color = "red",  # 更深的颜色
            size = 5,  # 字体大小
            hjust = 0.5)  # 水平居中对齐
print(second_plot)
### 导出图片到ppt
### 如果需要矢量图，先导出为svg，再用AI做成eps
library(export)
graph2ppt(file = "output",
          append = TRUE,
          width = 5.5, height = 4)

#BiocManager::install('clusterProfiler')

#library(clusterProfiler)
#kegg_geneset2 <- clusterProfiler::read.gmt("c2.cp.kegg.v7.0.symbols.gmt")
#head(kegg_geneset2)
#dat = readLines("c2.cp.kegg.v7.0.symbols.gmt")
#kegg_geneset2

