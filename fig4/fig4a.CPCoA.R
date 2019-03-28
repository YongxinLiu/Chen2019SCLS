# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
library("ggplot2") # load related packages
library("grid")
library("scales")
library("vegan")

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
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

# Public file 1. "design.txt"  Design of experiment
design = read.table("design.txt", header=T, row.names= 1, sep="\t") 

# setting subset design
if (TRUE){
	sub_design = subset(design,GroupID %in% c("DM3","DM4","TPS25a","TPS25b","TPS30a","TPS30b","WT4") ) # select group1
}else{
	sub_design = design
}
if (TRUE){
	sub_design = subset(sub_design,batch %in% c("3") ) # select group2
}

# Set group style, single or combine
if (FALSE){
	sub_design$group=paste(sub_design$GroupID,sub_design$batch,sep = "")
}else{
	sub_design$group=sub_design$GroupID
}

# Set group order
if ("TRUE" == "TRUE") {
    sub_design$group  = factor(sub_design$group, levels=c("DM3","DM4","TPS25a","TPS25b","TPS30a","TPS30b","WT4"))   # set group order
}

print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers

# 按组上色
colors = data.frame(group=unique(sub_design$group), 
                    color=rainbow(length(unique(sub_design$group)))) 
#shapes = data.frame(group=unique(sub_design$group),shape=c(1:length(unique(sub_design$group))))
# 按批次设置形状
shapes = data.frame(group=unique(sub_design$time),shape=c(15:(length(unique(sub_design$time))+14)))


# Function for analysis CPCoA/CCA result
variability_table = function(cca){
  chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table = cbind(chi, chi/chi[1])
  colnames(variability_table) = c("inertia", "proportion")
  rownames(variability_table) = c("total", "constrained", "unconstrained")
  return(variability_table)
}

# Load css OTU table, shape color same with beta diversity
otu_table = read.table("otu_table_css.txt", sep="\t", header=T, row.names= 1) # CSS norm otu table
idx = rownames(sub_design) %in% colnames(otu_table) 
sub_design = sub_design[idx,]
sub_otu_table = otu_table[, rownames(sub_design)] 

# Constrained analysis OTU table by genotype
capscale.gen = capscale(t(sub_otu_table) ~ group, data=sub_design, add=F, sqrt.dist=T, distance="bray") 

# ANOVA-like permutation analysis
perm_anova.gen = anova.cca(capscale.gen)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen = variability_table(capscale.gen)
eig = capscale.gen$CCA$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

# extract the weighted average (sample) scores
points = capscale.gen$CCA$wa[, 1:3]
points = as.data.frame(points)
colnames(points) = c("x", "y","z")
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])

# plot CPCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=time)) +
  geom_point(alpha=.7, size=1.5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape)+
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) + main_theme
p
ggsave(paste( "fig4a.CPCoA.pdf", sep=""), p, width = 5, height = 3)


