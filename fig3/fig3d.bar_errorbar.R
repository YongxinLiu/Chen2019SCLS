library(ggplot2)
library(plyr)
library(dplyr)
library("reshape2")

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))

# Figure 3d tps25/30 c3-4
average = read.table("figure3d.txt", header=T, sep="\t") 
sd = read.table("figure3d_sd.txt", header=T, sep="\t") 
melt_average = as.data.frame(melt(average, id.vars="X"))
melt_sd = as.data.frame(melt(sd, id.vars="X"))
melt_all = cbind(melt_average,melt_sd$value)
colnames(melt_all) = c("x","y","mean","sd")
melt_all$x  = factor(melt_all$x, levels=average$X)   # set group order

p=ggplot(melt_all, aes(x=x, y=mean, colour=y, group=y)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", width=.1,size=.5) +
  geom_line() + main_theme #  + geom_point(size=.5)
p

dodge <- position_dodge(width=0.9)
p=ggplot(melt_all, aes(x=x, y=mean, colour=y, group=y, fill=y)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", stat = "identity",  position = dodge, width=.5, size=.5) +
  geom_bar(stat = "identity",position=dodge, width=0.7) + main_theme +
  scale_colour_manual(values=as.character(rainbow(4)[c(3,4)])) +
  scale_fill_manual(values=as.character(rainbow(4)[c(3,4)])) +
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
p
ggsave("fig3d.bar_errorbar.pdf", p, width = 4, height = 6)

