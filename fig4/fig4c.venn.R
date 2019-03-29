library(VennDiagram)
pdf(file="tps25a_tps25b_enriched.pdf", onefile=FALSE,paper="special", width=4, height=4, pointsize=8)

if (! FALSE) {
	data <- read.table(file="otu.txt", sep="\t", quote="")
	num <- 0
	if("TPS25avsWT4_enriched" != "A"){
		TPS25avsWT4_enriched <- data[grepl("\\<TPS25avsWT4_enriched\\>",data[,2]),1]
		num <- num + 1
	}
	if("TPS25bvsWT4_enriched" != "B"){
		TPS25bvsWT4_enriched <- data[grepl("\\<TPS25bvsWT4_enriched\\>",data[,2]),1]
		num <- num + 1
	}
	if("C" != "C"){
		C <- data[grepl("\\<C\\>",data[,2]),1]
		num <- num + 1
	}
	if("D" != "D"){
		D <- data[grepl("\\<D\\>",data[,2]),1]
		num <- num + 1
	}
	if("E" != "E"){
		E <- data[grepl("\\<E\\>",data[,2]),1]
		num <- num + 1
	}
	color_v <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:num]
	#label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
	if(num == 5){
		p <- venn.diagram( 
			x = list(TPS25avsWT4_enriched=TPS25avsWT4_enriched, D=D,
			E=E, TPS25bvsWT4_enriched=TPS25bvsWT4_enriched,
			C=C),
			filename = NULL, col = "black", lwd = 1, 
			fill = color_v,
			alpha = 0.50,
			label.col = c("black"),
			cex = 1, fontfamily = "Helvetica",
#			cat.col = c("black"),cat.cex = 1.1, margin=0.1, 
			cat.col = color_v, cat.cex = 1, margin=0.1, 
			cat.fontfamily = "Helvetica"
		)
	}else if(num == 4){
		p <- venn.diagram( 
			x = list(TPS25avsWT4_enriched=TPS25avsWT4_enriched, D=D, TPS25bvsWT4_enriched=TPS25bvsWT4_enriched,
			C=C),
			filename = NULL, col = "black", lwd = 1, 
			fill = color_v,
			alpha = 0.50,
			label.col = c("black"),
			cex = 1, fontfamily = "Helvetica",
#			cat.col = c("black"), cat.cex = 1.1, margin=0.05, 
			cat.col = color_v, cat.cex = 1, margin=0.05, 
			cat.fontfamily = "Helvetica"
		)
	} else if (num==3) {
		p <- venn.diagram( 
			x = list(TPS25avsWT4_enriched=TPS25avsWT4_enriched, TPS25bvsWT4_enriched=TPS25bvsWT4_enriched, C=C),
			filename = NULL, col = "transparent", 
			fill = color_v,
			alpha = 0.50,
			label.col = c("black", "black", "black", "black", "black", "black", "black"),
			cex = 1, fontfamily = "Helvetica", cat.default.pos="text",
			cat.pos=0,  magrin=0.1, 
#			cat.col = c("black", "black", "black"),cat.cex = 1,cat.fontfamily = "Helvetica"
			cat.col = color_v, cat.cex = 1, cat.fontfamily = "Helvetica"
		)
	} else if (num==2) {
		p <- venn.diagram( 
			x = list(TPS25avsWT4_enriched=TPS25avsWT4_enriched, TPS25bvsWT4_enriched=TPS25bvsWT4_enriched),
			filename = NULL, col = "transparent", 
			fill = color_v,
			alpha = 0.50,
			label.col = c("black"),
			cex = 1, fontfamily = "Helvetica",
			cat.default.pos="outer",
			cat.pos=0, margin=0.1,  
			cat.col = color_v,cat.cex = 1, cat.fontfamily = "Helvetica"
		)
	}
	grid.draw(p)
} else {
#---venn plot for given numbers---------
	numList <- c()
	labelList <- c()
	num <- length(labelList)
	color_v <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:num]
	
	if (num==2) {
		draw.pairwise.venn(area1=numList[1], area2=numList[2],
		cross.area=numList[3], category=labelList, lwd=rep(1,1),
		lty=rep(2,2), col="transparent", fill=color_v,
		cat.col=color_v)
	} else if (num==3) {
		draw.triple.venn(area1=numList[1], area2=numList[2],
		area3=numList[3], n12=numList[4], n23=numList[5],
		n13=numList[6], n123=numList[7], 
		category=labelList, col="transparent", fill=color_v,
		cat.col=color_v, reverse=FALSE)
	}else if (num==4) {
		draw.quad.venn(area1=numList[1], area2=numList[2],
		area3=numList[3], area4=numList[4], n12=numList[5], 
		n13=numList[6], n14=numList[7], n23=numList[8],
		n24=numList[9], n34=numList[10], n123=numList[11], 
	    n124=numList[12], n134=numList[13], n234=numList[14], 
   		n1234=numList[15], 	   
		category=labelList, col="transparent", fill=color_v,
		cat.col=color_v, reverse=FALSE)
	}else if (num==5){
		#draw.quintuple.venn()
	}
}
dev.off()
system("rm *.log")
