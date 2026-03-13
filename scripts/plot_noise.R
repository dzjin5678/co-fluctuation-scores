library(R.matlab)
library(ggplot2)


hcp_noise_rss <- readMat("../results_bins/noise/fd_gs_ds_avg1_global_rss_corr.mat");

noise <- "FD"
Data1 = data.frame(Value = hcp_noise_rss$fd.corr.mean.all)
noise <- "GS"
Data1 = data.frame(Value = hcp_noise_rss$gs.corr.mean.all)
noise <- "RV"
Data1 = data.frame(Value = hcp_noise_rss$ds.corr.mean.all)
noise <- "HRV"
Data1 = data.frame(Value = hcp_noise_rss$avg1.corr.mean.all)

color_arr <- array(data = NA, dim = 7, dimnames = NULL)
color_arr[1] <- rgb(120/255, 154/255, 192/255)
P1 <- ggplot(Data1,aes(x="",y=Value))+
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.3,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=""), width=0.15, shape = 21, size=1.5, alpha=0.8, stroke=NA)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
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
  ylab(sprintf("Correlation Coefficient between \n %s and Global RSS (r)", toupper(noise))) + 
  xlab("") + 
  scale_x_discrete(labels = c("Global RSS")) + 
  scale_y_continuous(limits = c(-1, 1))
P1
jpeg(file = sprintf("../results_bins/noise/%s_global_rss_corr.jpg", toupper(noise)), 
     width = 1200,height = 1200,units = "px",res = 450)
print(P1)
dev.off()