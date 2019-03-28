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

# Figure 3C TPS tissue-specific Compound 1-4
average = read.table("fig3c.mean.txt", header=T, sep="\t") 
sd = read.table("fig3c.SD.txt", header=T, sep="\t") 
melt_average = as.data.frame(melt(average, id.vars="X"))
melt_sd = as.data.frame(melt(sd, id.vars="X"))
melt_all = cbind(melt_average,melt_sd$value)
colnames(melt_all) = c("x","y","mean","sd")
melt_all$x  = factor(melt_all$x, levels=average$X)   # set group order

p=ggplot(melt_all, aes(x=x, y=mean, colour=y, group=y)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", width=.1,size=.5) +
  geom_line() + main_theme +theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))

p
ggsave("fig3c_line_errorbar.pdf", p, width = 4, height = 3)
