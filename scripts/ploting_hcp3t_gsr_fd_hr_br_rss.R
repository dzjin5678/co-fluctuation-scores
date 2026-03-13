library(R.matlab)
library(ggplot2)


hcp_gar_rss <- readMat("./hcp_fd_rss_rfMRI_REST2_RL.mat");

Data1 = data.frame(Group = hcp_gar_rss$label.all, Value = hcp_gar_rss$coef.all)
Data1$Group = factor(Data1$Group,levels = c(1,2))

color_arr <- array(data = NA, dim = 7, dimnames = NULL)
color_arr[1] <- rgb(120/255, 154/255, 192/255)
color_arr[2] <- rgb(120/255, 154/255, 192/255)


#使用ggplot2包生成箱线图
P1 <- ggplot(Data1,aes(x=Group,y=Value,fill=Group))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.3,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Group), width=0.15, shape = 21, size=1.5, alpha=0.8, stroke=NA)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = color_arr)+  #设置填充的颜色
  scale_color_manual(values=c("black","black","black","black","black","black","black","black"))+ #设置散点图的圆圈的颜色为黑色
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Arial",size=11,face="plain", angle=0), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(colour="black",family="Arial",size=11,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(colour="black",family="Arial",size = 11,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(colour="black",family="Arial",size = 11,face="plain"), #设置x轴的标题的字体属性
        # plot.title = element_text(family="Arial",size=11,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour="gray", linewidth=0.2),
        axis.ticks.x = element_line(colour="gray", linewidth=0.1),
        axis.ticks.y = element_line(colour="gray", linewidth=0.1))+
  ylab("Correlation Coefficient between \n FD and Global RSS (r)") + 
  scale_x_discrete(labels = c("noGSR", "GSR")) + 
  scale_y_continuous(limits = c(-1, 1))

P1
jpeg(file = sprintf("./hcp_fd_rss_rfMRI_REST2_RL.jpg"), 
     width = 1200,height = 1200,units = "px",res = 450)
print(P1)
dev.off()