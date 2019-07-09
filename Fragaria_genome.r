rm(list = ls())
pdf('/path/to/file/PhylogenyWindowsFragaria.pdf', width = 12, height = 6)

plot(-1000,-1000,xlim=c(0,310000000),ylim=c(-0.002,0.032),ylab="Subgenome",xlab="Haploid Genome Position",main="",yaxt="n",xaxt="n",bty = "l")
axis(side=1,at=c(12126511.5,43928638,82766404,123882480.5,160551478.5,200164166,237176575.5),lab=c(1,2,3,4,5,6,7))
axis(side=2,at=0,lab="Camarosa\nvesca",col.axis = "#e41a1c")
axis(side=2,at=0.01,lab="Camarosa\nviridis",col.axis = "#4daf4a")
axis(side=2,at=0.02,lab="Camarosa\nnipponica",col.axis = "#984ea3")
axis(side=2,at=0.03,lab="Camarosa\niinumae",col.axis = "#377eb8")

camvestrees<-read.table('/path/to/file/Camarosa_vesca_tree_windows_closest_diploid.txt', header=T)
attach(camvestrees)
lines(camvestrees$Position,camvestrees$Position*0,type="p",pch="|",col=as.character(camvestrees$Color),cex=3)

camvirtrees<-read.table('/path/to/file/Camarosa_viridis_tree_windows_closest_diploid.txt', header=T)
attach(camvirtrees)
lines(camvirtrees$Position,camvirtrees$Position*0+0.01,type="p",pch="|",col=as.character(camvirtrees$Color),cex=3)

camniptrees<-read.table('/path/to/file/Camarosa_nipponica_tree_windows_closest_diploid.txt', header=T)
attach(camniptrees)
lines(camniptrees$Position,camniptrees$Position*0+0.02,type="p",pch="|",col=as.character(camniptrees$Color),cex=3)

camiintrees<-read.table('/path/to/file/Camarosa_iinumae_tree_windows_closest_diploid.txt', header=T)
attach(camiintrees)
lines(camiintrees$Position,camiintrees$Position*0+0.03,type="p",pch="|",col=as.character(camiintrees$Color),cex=3)

legend(255000000,0.030,c("F vesca sister","F vesca clade","F viridis sister","F viridis clade","F nipponica sister","F nipponica clade","F iinumae sister","F iinumae clade","F nilgerrensis sister","F nilgerrensis clade","Ambiguous","Absent"),col=c("#e41a1c","#f8c1c2","#4daf4a","#cce9cb","#984ea3","#e0c4e4","#377eb8","#bfd6ec","#ff7f00","#ffdab4","black","grey70"),lty=-1,pch=15,box.lty = 0,cex=1.2)

dev.off()

pdf('/path/to/file/PhylogenyPropotionsFragaria.pdf', width = 7, height = 6)

plot(-1000,-1000,xlim=c(0.8,5.72),ylim=c(0,100),xaxt="n",xlab="",ylab="",bty = "l")
axis(side=1,at=1,lab="Camarosa\nvesca",col.axis = "#e41a1c",cex.lab=0.5, tck = 0)
axis(side=1,at=2,lab="Camarosa\nviridis",col.axis = "#4daf4a",cex.lab=0.5, tck = 0)
axis(side=1,at=3,lab="Camarosa\nnipponica",col.axis = "#984ea3",cex.lab=0.5, tck = 0)
axis(side=1,at=4,lab="Camarosa\niinumae",col.axis = "#377eb8",cex.lab=0.5, tck = 0)

#Examining Camarosa_vesca
rect(0.8,0,1.2,64.6514,col="#e41a1c",border=NA)
rect(0.8,64.6514,1.2,98.0497,col="#f8c1c2",border=NA)
rect(0.8,98.0497,1.2,98.2935,col="#ffdab4",border=NA)
rect(0.8,98.2935,1.2,98.586,col="#4daf4a",border=NA)
rect(0.8,98.586,1.2,98.781,col="#cce9cb",border=NA)
rect(0.8,98.781,1.2,99.1223,col="#984ea3",border=NA)
rect(0.8,99.1223,1.2,99.3173,col="#e0c4e4",border=NA)
rect(0.8,99.3173,1.2,99.9511,col="#bfd6ec",border=NA)
rect(0.8,99.9511,1.2,99.9999,col="black",border=NA)

#Examining Camarosa_viridis
rect(1.8,0,2.2,0.7583,col="#e41a1c",border=NA)
rect(1.8,0.7583,2.2,3.9337,col="#f8c1c2",border=NA)
rect(1.8,3.9337,2.2,4.5498,col="#ff7f00",border=NA)
rect(1.8,4.5498,2.2,8.8626,col="#ffdab4",border=NA)
rect(1.8,8.8626,2.2,9.0048,col="#4daf4a",border=NA)
rect(1.8,9.0048,2.2,10.1422,col="#cce9cb",border=NA)
rect(1.8,10.1422,2.2,10.2844,col="#984ea3",border=NA)
rect(1.8,10.2844,2.2,10.9479,col="#e0c4e4",border=NA)
rect(1.8,10.9479,2.2,12.5593,col="#377eb8",border=NA)
rect(1.8,12.5593,2.2,99.7631,col="#bfd6ec",border=NA)
rect(1.8,99.7631,2.2,100.0001,col="black",border=NA)

#Examining Camarosa_nipponica
rect(2.8,0,3.2,0.5155,col="#e41a1c",border=NA)
rect(2.8,0.5155,3.2,5.1078,col="#f8c1c2",border=NA)
rect(2.8,5.1078,3.2,5.4827,col="#ff7f00",border=NA)
rect(2.8,5.4827,3.2,9.1378,col="#ffdab4",border=NA)
rect(2.8,9.1378,3.2,9.3252,col="#4daf4a",border=NA)
rect(2.8,9.3252,3.2,11.0122,col="#cce9cb",border=NA)
rect(2.8,11.0122,3.2,11.1996,col="#984ea3",border=NA)
rect(2.8,11.1996,3.2,12.1837,col="#e0c4e4",border=NA)
rect(2.8,12.1837,3.2,13.9175,col="#377eb8",border=NA)
rect(2.8,13.9175,3.2,99.7657,col="#bfd6ec",border=NA)
rect(2.8,99.7657,3.2,100,col="black",border=NA)

#Examining Camarosa_iinumae
rect(3.8,0,4.2,0.2322,col="#e41a1c",border=NA)
rect(3.8,0.2322,4.2,2.7868,col="#f8c1c2",border=NA)
rect(3.8,2.7868,4.2,2.8797,col="#ff7f00",border=NA)
rect(3.8,2.8797,4.2,3.2977,col="#ffdab4",border=NA)
rect(3.8,3.2977,4.2,3.6228,col="#cce9cb",border=NA)
rect(3.8,3.6228,4.2,3.7621,col="#e0c4e4",border=NA)
rect(3.8,3.7621,4.2,85.3692,col="#377eb8",border=NA)
rect(3.8,85.3692,4.2,99.8142,col="#bfd6ec",border=NA)
rect(3.8,99.8142,4.2,100,col="black",border=NA)

legend(4.3,100,c("F vesca sister","F vesca clade","F viridis sister","F viridis clade","F nipponica sister","F nipponica clade","F iinumae sister","F iinumae clade","F nilgerrensis sister","F nilgerrensis clade","Ambiguous"),col=c("#e41a1c","#f8c1c2","#4daf4a","#cce9cb","#984ea3","#e0c4e4","#377eb8","#bfd6ec","#ff7f00","#ffdab4","black"),lty=-1,pch=15,box.lty = 0)

dev.off()

