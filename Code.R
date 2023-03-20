#https://github.com/YinLiLin/CMplot

SNP <- read.csv(file = 'pig.csv')
head(SNP)

#run each section separately and save as each figure (figures overwrite)

CMplot(SNP,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("blue", "yellow", "red"),file="jpg",memo="",dpi=300,
    main="SNP density",file.output=TRUE,verbose=TRUE,width=9,height=6)


CMplot(SNP,type="p",plot.type="c",chr.labels=paste("Chr",c(1:18),sep=""),r=0.4,cir.legend=TRUE,
        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
        memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)


CMplot(SNP,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18),sep=""),
      threshold=c(1e-6,1e-4),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red",
      "blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
      bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)


CMplot(SNP,type="p",plot.type="m",LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,
    file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)


CMplot(SNP, plot.type="m", col=c("grey30","grey60"), LOG10=TRUE, ylim=c(2,12), threshold=c(1e-6,1e-4),
        threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
        chr.den.col=NULL, signal.col=c("red","green"), signal.cex=c(1.5,1.5),signal.pch=c(19,19),
        file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)


CMplot(SNP, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),
        threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1.5,1.5),
        signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
        width=14,height=6)



signal <- SNP$Position[which.min(SNP$trait2)]
SNPs <- SNP$SNP[SNP$Chromosome==13 & 
        SNP$Position<(signal+1000000)&SNP$Position>(signal-1000000)]
CMplot(SNP, plot.type="m",LOG10=TRUE,col=c("grey30","grey60"),highlight=SNPs,
        highlight.col="green",highlight.cex=1,highlight.pch=19,file="jpg",memo="",
        chr.border=TRUE,dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)



SNPs <-  SNP[SNP$trait2 < 1e-4, 1]
CMplot(SNP,type="h",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="p",
        highlight.col=NULL,highlight.cex=1.2,highlight.pch=19,file="jpg",memo="",
        dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6,band=0.6)


SNPs <-  SNP[SNP$trait2 < 1e-4, 1]
CMplot(SNP,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="h",
        col=c("grey30","grey60"),highlight.col="darkgreen",highlight.cex=1.2,highlight.pch=19,
        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)




SNPs <-  SNP[
	SNP$trait1 < 1e-4 |
	SNP$trait2 < 1e-4 |
	SNP$trait3 < 1e-4, 1]
CMplot(SNP,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
        threshold=1e-4,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
        signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)


CMplot(SNP[SNP$Chromosome==13, ], plot.type="m",LOG10=TRUE,col=c("grey60"),highlight=SNPs,
        highlight.col="green",highlight.cex=1,highlight.pch=19,file="jpg",memo="", 
        threshold=c(1e-6,1e-4),threshold.lty=c(1,2),threshold.lwd=c(1,2), width=9,height=6,
        threshold.col=c("red","blue"),amplify=FALSE,dpi=300,file.output=TRUE,verbose=TRUE)


SNPs <- SNP[SNP[,5] < (0.05 / nrow(SNP)), 1]
genes <- paste("GENE", 1:length(SNPs), sep="_")
set.seed(666666)
CMplot(SNP[,c(1:3,5)], plot.type="m",LOG10=TRUE,col=c("grey30","grey60"),highlight=SNPs,
        highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=c(15:17), highlight.text=genes,      
        highlight.text.col=c("red","blue","green"),threshold=0.05/nrow(pig60K),threshold.lty=2,   
        amplify=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)


SNPs <- list(
	SNP$SNP[SNP$trait1<1e-6],
	SNP$SNP[SNP$trait2<1e-6],
	SNP$SNP[SNP$trait3<1e-6]
)
CMplot(SNP, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
        threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
        signal.cex=1, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
        highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4)



CMplot(SNP,plot.type="q",box=FALSE,file="jpg",memo="",dpi=300,
    conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
    file.output=TRUE,verbose=TRUE,width=5,height=5)





SNP$trait1[sample(1:nrow(SNP), round(nrow(SNP)*0.80))] <- NA
SNP$trait2[sample(1:nrow(SNP), round(nrow(SNP)*0.25))] <- NA
CMplot(SNP,plot.type="q",col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),threshold=1e-6,
        ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",conf.int=TRUE,box=FALSE,multracks=
        TRUE,cex.axis=2,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,ylim=c(0,8),width=5,height=5)




