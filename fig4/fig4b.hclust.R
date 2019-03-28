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

# Hclust
x <- cbind(as.matrix(points$x),as.matrix(points$y))
colnames(x) <- c("x", "y")

rownames(x) = rownames(points)
dat_dist <- vegdist(x, method="euclidean")

dat_dist_clu <- hclust(dat_dist, "average")

plot(dat_dist_clu, main = "Hcluste with euclidean", lwd=1.5, cex=.5)

# We can now create a basic dendrogram:
hc=hclust(dat_dist, "average")
dhc=as.dendrogram(hc)

# Actually, each leaf of the tree has several attributes, like the color, the shape … You can take back these attributes and look at them like that
# 提取叶的属性，观察树的基本信息
specific_leaf=dhc[[1]][[1]][[1]]
specific_leaf
attributes(specific_leaf)

# 准备注释文件
anno = cbind(rownames(points), points[,c("time","group")])
colnames(anno)[1]="sample"
head(anno)

# So if I Want to color each leaf of the Tree, I have to change the attributes of each leaf. This can be done using the dendrapply function. So I create a function that add 3 attributes to the leaf : one for the color (“lab.col”) ,one for the font “lab.font” and one for the size (“lab.cex”).
# 函数，设置叶形状为time 1，2，3分别为方、圆和三角，对应15/16/17；颜色按组彩虹色
i=0
colLab<<-function(n){
  if(is.leaf(n)){
    #I take the current attributes
    a=attributes(n)
    # 匹配样本名位置
    ligne=match(attributes(n)$label,anno[,1])
    # 按位置设置颜色和形状
    # 二列time设置形状
    time=anno[ligne,2];
    if(time=="T1"){shape_time=15};if(time=="T2"){shape_time=16};if(time=="T3"){shape_time=17}
    group=anno[ligne,3];
    if(group=="DM3"){col_group=rainbow(7)[1]};if(group=="DM4"){col_group=rainbow(7)[2]};if(group=="TPS25a"){col_group=rainbow(7)[3]}
    if(group=="TPS25b"){col_group=rainbow(7)[4]};if(group=="TPS30a"){col_group=rainbow(7)[5]};if(group=="TPS30b"){col_group=rainbow(7)[6]}
    if(group=="WT4"){col_group=rainbow(7)[7]}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=shape_time,col=col_group,lab.col=col_group,lab.font=1,lab.cex=1))
  }
  return(n)
}

# Finally I just have to apply this to my dendrogram
dL <- dendrapply(dhc, colLab)
plot(dL , main="structure of the population")


pdf("fig4b.hclust.pdf", width=16, height=10)
plot(dL, main = "Hclust with euclidean")
dev.off()





# Reference http://www.r-graph-gallery.com/31-custom-colors-in-dendrogram/

# # Dataset 
# sample=paste(rep("sample_",24) , seq(1,24) , sep="")
# specie=c(rep("dicoccoides" , 8) , rep("dicoccum" , 8) , rep("durum" , 8))
# treatment=rep(c(rep("High",4 ) , rep("Low",4)),3)
# data=data.frame(sample,specie,treatment)
# for (i in seq(1:5)){
#   gene=sample(c(1:40) , 24 )
#   data=cbind(data , gene)
#   colnames(data)[ncol(data)]=paste("gene_",i,sep="")
# }
# data[data$treatment=="High" , c(4:8)]=data[data$treatment=="High" , c(4:8)]+100
# data[data$specie=="durum" , c(4:8)]=data[data$specie=="durum" , c(4:8)]-30
# data
# 
# # Euclidean distance
# rownames(data)=data[,1]
# dist=dist(data[ , c(4:8)] , diag=TRUE)
# 
# # We can now create a basic dendrogram:
# hc=hclust(dist)
# dhc=as.dendrogram(hc)
# 
# # Actually, each leaf of the tree has several attributes, like the color, the shape … You can take back these attributes and look at them like that
# # 提取叶的属性，观察树的基本信息
# specific_leaf=dhc[[1]][[1]][[1]]
# specific_leaf
# attributes(specific_leaf)
# 
# # So if I Want to color each leaf of the Tree, I have to change the attributes of each leaf. This can be done using the dendrapply function. So I create a function that add 3 attributes to the leaf : one for the color (“lab.col”) ,one for the font “lab.font” and one for the size (“lab.cex”).
# # 函数，设置叶标签颜色、字体和大小
# i=0
# colLab<<-function(n){
#   if(is.leaf(n)){
#     
#     #I take the current attributes
#     a=attributes(n)
#     
#     #I deduce the line in the original data, and so the treatment and the specie.
#     ligne=match(attributes(n)$label,data[,1])
#     treatment=data[ligne,3];
#     if(treatment=="Low"){col_treatment="blue"};if(treatment=="High"){col_treatment="red"}
#     specie=data[ligne,2];
#     if(specie=="dicoccoides"){col_specie="red"};if(specie=="dicoccum"){col_specie="Darkgreen"};if(specie=="durum"){col_specie="blue"}
#     
#     #Modification of leaf attribute
#     attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,col=col_treatment,lab.col=col_specie,lab.font=1,lab.cex=1))
#   }
#   return(n)
# }
# 
# # Finally I just have to apply this to my dendrogram
# dL <- dendrapply(dhc, colLab)
# 
# plot(dL , main="structure of the population")
# legend("topright", 
#        legend = c("High Nitrogen" , "Low Nitrogen" , "Durum" , "Dicoccoides" , "Dicoccum"), 
#        col = c("red", "blue" , "blue" , "red" , "Darkgreen"), 
#        pch = c(20,20,4,4,4), bty = "n", pt.cex = 1.5, cex = 0.8 , 
#        text.col = "black", horiz = FALSE, inset = c(0.1, 0.1))
