#マイクロアレイのスタンダードな解析
#正規化、アノテーション、X倍発現、cmap、gsea、ヒートマップ、倍率検討　まで

#--------------------------------設定が必要--------------------------------------------------------

#基本設定
#作業ディレクトリの設定
#使用するCELファイル、このプログラム、gsym-U133-na33.txt、U133A-psid.txtがそのディレクトリに入っているようにする
setwd("/Users/haruka/Desktop/WORK/2014/141001_POC解析まとめ")
title <- "POC"

#正規化
#mas5で正規化した後のファイル名を決める
mas5_filnm <- "POC20-24_mas5"

#X倍発現
#50未満→50にしている,logはとっていない
cont <- 6　#コントロールが含まれている列番号
row_samp <- 1:5
th.rat <- 2　#何倍以上発現したプローブを出力したいのか
th.val <- 100　#比を取ったときに分子の発現量が何以上であればよいか
samp_name <- c(20:24)

#ヒートマップ
#ソートする時の基準となるサンプルの列番号
sortnum <- 1

#--------------------------------設定が必要--------------------------------------------------------


#結果フォルダの作成
res_dirnm <- "result/"
dir.create(res_dirnm)
rat_dirnm <- "result/rat"
dir.create(rat_dirnm)
cmap_dirnm <- "result/cmap"
dir.create(cmap_dirnm)
gsea_dirnm <- "result/gsea"
dir.create(gsea_dirnm)
heatmap_dirnm <- "result/heatmap"
dir.create(heatmap_dirnm)
plot_dirnm <- "result/plot"
dir.create(plot_dirnm)


#mas5で正規化
library(affy)
af <- ReadAffy()
rset <- mas5(af, sc=100)

write.exprs(rset, file=paste(mas5_filnm,".txt",sep=""), sep="\t")

#遺伝子名付与

data <- read.table(paste(mas5_filnm,".txt",sep=""), header=T, sep="\t", row.names=1)
gid <- rownames(data)
sname <- substr(colnames(data), 13, nchar(colnames(data))-21)

gsym <- read.table("gsym-U133-na33.txt",sep="\t",header=T,quote="\"",row.names=1)
gid.gsym <- as.vector(gsym[gid,3:4])

out <- cbind(gid.gsym,data)

mas5_gene_filnm <- paste(mas5_filnm,"_gene.txt",sep="")

write.table(t(c("Probeset_ID","Symbol","Title",sname)),mas5_gene_filnm,append=T,col.names=F,row.names=F,sep="\t",quote=F)
write.table(out,mas5_gene_filnm,append=T,col.names=F,row.names=T,sep="\t",quote=F)


#140930_X倍以上発現or1/X以下発現したプローブと遺伝子名とそのレシオを出力　ngeneでまわさない方法

data_anno <- read.table(paste(mas5_filnm,"_gene.txt",sep=""),row.names=1,header=T,sep="\t",quote="\"")
nsamp <- length(row_samp)
colname <- c("Probeset_ID","Symbol","Title","Ratio")

for(i in 1:nsamp){
	
	outputname_up <- paste(title,samp_name[i],"_rat_up_",th.rat,".txt",sep="")
	outputname_down <- paste(title,samp_name[i],"_rat_down_",th.rat,".txt",sep="")
	
	rat <- pmax(data[,row_samp[i]]/data[,cont])	
	up_sign <- rat >=th.rat & data[,row_samp[i]] >= th.val
	down_sign <- rat <= 1/th.rat & data[,cont] >= th.val
	gid_list_up <- c()
	gid_list_down <- c()
	gid_list_up <- c(gid_list_up,gid[up_sign])
	gid_list_down <- c(gid_list_down,gid[down_sign])
	
	ct_up <- pmax(data[gid_list_up,cont],50)
	x_up <- pmax(data[gid_list_up,row_samp[i]],50)
	rat_up <- x_up/ct_up
	ct_down <- pmax(data[gid_list_down,cont],50)
	x_down <- pmax(data[gid_list_down,row_samp[i]],50)
	rat_down <- x_down/ct_down
	
	up_matrix <- matrix(0,length(gid_list_up),length(colname))
	down_matrix <- matrix(0,length(gid_list_down),length(colname))
	
	up_matrix[,1] <- gid_list_up
	down_matrix[,1] <- gid_list_down
	up_matrix[,2] <- as.character(data_anno[gid_list_up,1])
	down_matrix[,2] <- as.character(data_anno[gid_list_down,1])
	up_matrix[,3] <- as.character(data_anno[gid_list_up,2])
	down_matrix[,3] <- as.character(data_anno[gid_list_down,2])
	up_matrix[,4] <- rat_up
	down_matrix[,4] <- rat_down
	
	colnames(up_matrix) <- colname
	colnames(down_matrix) <- colname	
	
	write.table(up_matrix,paste("result/rat/",outputname_up,sep=""),append=T,quote=F,sep="\t",row.names=F,col.names=T)
	write.table(down_matrix,paste("result/rat/",outputname_down,sep=""),append=T,quote=F,sep="\t",row.names=F,col.names=T)		    
			
}

#Cmap
#Cmapに含まれるプローブidのみを抜き出したい

filename_up <- list(rep(NA,nsamp))
filename_down <- list(rep(NA,nsamp))
outputname_up <- list(rep(NA,nsamp))
outputname_down <- list(rep(NA,nsamp))


for(i in 1:nsamp){
	
	filename_up[i] <- paste(title,samp_name[i],"_rat_up_2.txt",sep="")
	filename_down[i] <- paste(title,samp_name[i],"_rat_down_2.txt",sep="")
	outputname_up[i] <- paste(title,samp_name[i],"_rat_up_2_for_cmap.txt",sep="")
	outputname_down[i] <- paste(title,samp_name[i],"_rat_down_2_for_cmap.txt",sep="")
	
}


filename <- list(filename_up,filename_down)
outputname <- list(outputname_up,outputname_down)

for(i in 1:2){		
for(j in 1:nsamp){
		
data_rat <- read.table(paste("result/rat/",filename[[i]][[j]],sep=""),header=T,sep="\t",row.names=1,quote="\"")	
cmap_gene_data <- read.table("U133A-psid.txt")
cmap_gene_list <- c()
cmap_gene_list <- as.character(cmap_gene_data[,1])
pick_up_gene <- rownames(data_rat[cmap_gene_list,])
NA_id <- grep("NA",pick_up_gene)
pick_up_gene_2 <- pick_up_gene[-NA_id]
write.table(t(pick_up_gene_2),paste("result/cmap/",outputname[[i]][[j]],sep=""),append=T,sep="\n",row.names=F,col.names=F,quote=F)

}
}

#GSEA
description <- as.character(data_anno[,1])
output <- c()

for(i in 1:nsamp){
	output <- c(output,paste(title,samp_name[i],"_for_gsea.txt",sep=""))
	write.table(t(c("NAME","DESCRIPTION",paste(title,samp_name[i],sep=""),"Control")),paste("result/gsea/",output[i],sep=""),append=T,quote=F,sep="\t",row.names=F,col.names=F)
}

for(i in 1:nsamp){

	Chem_sign <- data[,i] >= th.val & data[,cont] >= th.val
	#&&にすると、TRUE/FALSEが一個しかでてこなくなる

	NAME <- rownames(data[Chem_sign,])
	DESCRIPTION <- description[Chem_sign]
	CHEM <- data[Chem_sign,i]
	CONTROL <- data[Chem_sign,cont]

	ngene <- length(NAME)
	gsea_matrix <- matrix(0,ngene,4)	
	gsea_matrix[,1] <- NAME
	gsea_matrix[,2] <- DESCRIPTION
	gsea_matrix[,3] <- CHEM
	gsea_matrix[,4] <- CONTROL
	
	write.table(gsea_matrix,paste("result/gsea/",output[i],sep=""),append=T,quote=F,sep="\t",row.names=F,col.names=F)
		
}


#ヒートマップ作成

gid_list <- c()

for(i in 1:nsamp){
	
	rat <- pmax(data[,row_samp[i]]/data[,cont],50)	
	up_sign <- rat >=th.rat & data[,row_samp[i]] >= th.val
	down_sign <- rat <= 1/th.rat & data[,cont] >= th.val
	gid_list <- c(gid_list,gid[up_sign],gid[down_sign])
		
}

gid_list_uniq <- unique(gid_list)
data_2 <- data[gid_list_uniq,]

ngene <- nrow(data_2)
ratio_mat <- matrix(0,ngene,nsamp)

for(i in 1:nsamp){
	
	ct <- pmax(data_2[,cont],50)
	x <- pmax(data_2[,row_samp[i]],50)
	ratio_mat[,i] <- log(x/ct)
	
}

colnames(ratio_mat) <- sname[row_samp]
rownames(ratio_mat) <- gid_list_uniq

#sort 
od <- order(ratio_mat[,sortnum],decreasing=T)
ratio_mat_sort <- ratio_mat[od,]

# ヒートマップ表示
#install.packages("gplots")
library(gplots)
dist.p <- function(x) as.dist((1-cor(t(x)))/2)
dist.s <- function(x) as.dist((1-cor(t(x), method="spearman"))/2)
hclust.ward <- function(x) hclust(x, "ward.D")
#ウォード法という距離計算の仕方らしい。ward→ward.Dになぜか変更された

col1 <- greenred(60)

png(file=paste("result/heatmap/",title,"_heatmap.png",sep=""), width=600, height=600)
hv <- heatmap.2(ratio_mat_sort, Rowv=F, distfun=dist.p, hclustfun=hclust.ward, dendrogram="col", margin=c(5,2),
           breaks=seq(-3,3,0.1), col=col1, scale="none", #ColSideColors=patientcolors,
           key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none", labRow=F, cexCol=1)

#保存する        
dev.off()     


#倍率の検討のためのプロット

ncol <- ncol(data)

poc <- as.list(rep(NA,ncol))

for(i in 1:ncol){
	
	poc[[i]] <- pmax(50,data[,i])
		
}

for(i in 1:nsamp){
	
	outputname <- paste("result/plot/",title,samp_name[i],"_plot.pdf",sep="")
	pdf(outputname)
	x <- poc[[cont]]
	y <- poc[[i]]
	plot(x,y)
	abline(0,1,col=1)
	abline(0,2,col=2)
	abline(0,3,col=3)
	dev.off()

}
