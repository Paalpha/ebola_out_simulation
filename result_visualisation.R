# visualization of simulation experiments
tiff(file = "SampleVisualisation.tiff", width = (2300*3.7), height = (2800*3), units = "px", res = 600)
par(mfcol=c(4,3))
par(mar = c(0.5,2.5,1,0), oma = c(3,0.5,1,1))
par(mgp = c(1, 0.0, 0))
par(las=1)
#10% Missingness
plot(Cnconfirmed_01MARfinal$p,Cnconfirmed_01MARfinal$sens.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)],
     xlab="",ylab="",title(main = "MAR missingness=0·1"),xaxt="n", cex=1, yaxt='n', ylim=c(0.52,0.72),
     pch = c(16, 17, 4)[factor(Cnconfirmed_01MARfinal$pw)])
title(ylab="Sensitivity",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.52,0.72,by=0.04)
axis(2,at=ylabels,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$sens.Valid-(Cnconfirmed_01MARfinal$sens.Valid.1*1.96),
        Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$sens.Valid+(Cnconfirmed_01MARfinal$sens.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
legend("topright", title = "", bty = "n",
       legend = c("p=0·5","p=0·75","p=1"),cex=1,
       pch = c(16, 17, 4)) 


plot(Cnconfirmed_01MARfinal$p,Cnconfirmed_01MARfinal$spec.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.42,0.74),
     pch = c(16, 17, 4)[factor(Cnconfirmed_01MARfinal$pw)])
title(ylab="specificity",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.42,0.74,by=0.04)
axis(2,at=ylabels,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$spec.Valid-(Cnconfirmed_01MARfinal$spec.Valid.1*1.96),
        Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$spec.Valid+(Cnconfirmed_01MARfinal$spec.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 


plot(Cnconfirmed_01MARfinal$p,Cnconfirmed_01MARfinal$PCC.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.44,0.72),
     pch = c(16, 17, 4)[factor(Cnconfirmed_01MARfinal$pw)])
title(ylab="PCC",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.44,0.72,by=0.04)
axis(2,at=ylabels,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$PCC.Valid-(Cnconfirmed_01MARfinal$PCC.Valid.1*1.96),
        Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$PCC.Valid+(Cnconfirmed_01MARfinal$PCC.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 


plot(Cnconfirmed_01MARfinal$p,Cnconfirmed_01MARfinal$ROC.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.52,0.76),
     pch = c(16, 17, 4)[factor(Cnconfirmed_01MARfinal$pw)])
title(ylab="AUC",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=TRUE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.52,0.76,by=0.04)
axis(2,at=ylabels,tck=-0.01,cex.axis=1)
abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$ROC.Valid-(Cnconfirmed_01MARfinal$ROC.Valid.1*1.96),
        Cnconfirmed_01MARfinal$p, Cnconfirmed_01MARfinal$ROC.Valid+(Cnconfirmed_01MARfinal$ROC.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_01MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 

#20% Missingness
plot(Cnconfirmed_02MARfinal$p,Cnconfirmed_02MARfinal$sens.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)],
     xlab="",ylab="",title(main = "MAR missingness=0·2"),xaxt="n", cex=1, yaxt='n', ylim=c(0.52,0.72),
     pch = c(16, 17, 4)[factor(Cnconfirmed_02MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.52,0.72,by=0.04)
axis(2,at=ylabels, labels=FALSE ,tck=-0.01,cex.axis=1)
abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$sens.Valid-(Cnconfirmed_02MARfinal$sens.Valid.1*1.96),
        Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$sens.Valid+(Cnconfirmed_02MARfinal$sens.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 
legend("topright", title = "", bty="n",
       legend = c("LR","RF","BRT","BART","ANN"),
       col =c("cyan4","chartreuse4","darkorchid","black","brown"),lwd =1.5,cex = 0.9,
       lty = "solid") 

plot(Cnconfirmed_02MARfinal$p,Cnconfirmed_02MARfinal$spec.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.42,0.74),
     pch = c(16, 17, 4)[factor(Cnconfirmed_02MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.42,0.74,by=0.04)
axis(2,at=ylabels,labels=FALSE,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$spec.Valid-(Cnconfirmed_02MARfinal$spec.Valid.1*1.96),
        Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$spec.Valid+(Cnconfirmed_02MARfinal$spec.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 


plot(Cnconfirmed_02MARfinal$p,Cnconfirmed_02MARfinal$PCC.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.44,0.72),
     pch = c(16, 17, 4)[factor(Cnconfirmed_02MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.44,0.72,by=0.04)
axis(2,at=ylabels,labels=FALSE, tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$PCC.Valid-(Cnconfirmed_02MARfinal$PCC.Valid.1*1.96),
        Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$PCC.Valid+(Cnconfirmed_02MARfinal$PCC.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 


plot(Cnconfirmed_02MARfinal$p,Cnconfirmed_02MARfinal$ROC.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.52,0.76),
     pch = c(16, 17, 4)[factor(Cnconfirmed_02MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=TRUE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.52,0.76,by=0.04)
axis(2,at=ylabels,labels=FALSE,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$ROC.Valid-(Cnconfirmed_02MARfinal$ROC.Valid.1*1.96),
        Cnconfirmed_02MARfinal$p, Cnconfirmed_02MARfinal$ROC.Valid+(Cnconfirmed_02MARfinal$ROC.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_02MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 
#40% Missingness
plot(Cnconfirmed_04MARfinal$p,Cnconfirmed_04MARfinal$sens.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)],
     xlab="",ylab="",title(main = "MAR missingness=0·4"),xaxt="n", cex=1, yaxt='n', ylim=c(0.52,0.72),
     pch = c(16, 17, 4)[factor(Cnconfirmed_04MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.52,0.72,by=0.04)
axis(2,at=ylabels, labels=FALSE ,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$sens.Valid-(Cnconfirmed_04MARfinal$sens.Valid.1*1.96),
        Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$sens.Valid+(Cnconfirmed_04MARfinal$sens.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 

#legend("topright", title = "", bty="n",
#       legend = c("BART","ANN"),
#      col =c("black","brown"),lwd =1.5,cex = 1,
#       lty = "solid") 
#legend("topleft", title = "", bty="n",
#legend = c("Glm","RF","BRT","BART","ANN"),
#col =c("cyan4","chartreuse4","darkorchid","black","brown"),lwd =1.5,cex = 0.9,
#lty = "solid") 
plot(Cnconfirmed_04MARfinal$p,Cnconfirmed_04MARfinal$spec.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.42,0.74),
     pch = c(16, 17, 4)[factor(Cnconfirmed_04MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.42,0.74,by=0.04)
axis(2,at=ylabels,labels=FALSE,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$spec.Valid-(Cnconfirmed_04MARfinal$spec.Valid.1*1.96),
        Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$spec.Valid+(Cnconfirmed_04MARfinal$spec.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 


plot(Cnconfirmed_04MARfinal$p,Cnconfirmed_04MARfinal$PCC.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.44,0.72),
     pch = c(16, 17, 4)[factor(Cnconfirmed_04MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=FALSE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.44,0.72,by=0.04)
axis(2,at=ylabels,labels=FALSE, tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$PCC.Valid-(Cnconfirmed_04MARfinal$PCC.Valid.1*1.96),
        Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$PCC.Valid+(Cnconfirmed_04MARfinal$PCC.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 
plot(Cnconfirmed_04MARfinal$p,Cnconfirmed_04MARfinal$ROC.Valid,col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)],
     xlab="",ylab="",title(main = ""),xaxt="n", cex=1, yaxt='n', ylim=c(0.52,0.76),
     pch = c(16, 17, 4)[factor(Cnconfirmed_04MARfinal$pw)])
title(ylab="",line=1.7)
xlabels<-c(50,65,80)
axis(1,at=xlabels,labels=TRUE,tck=-0.01,cex.axis=1)
ylabels<-seq(0.52,0.76,by=0.04)
axis(2,at=ylabels,labels=FALSE,tck=-0.01,cex.axis=1)
#abline(h = 0.5,lty=2,col="black")
arrows( Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$ROC.Valid-(Cnconfirmed_04MARfinal$ROC.Valid.1*1.96),
        Cnconfirmed_04MARfinal$p, Cnconfirmed_04MARfinal$ROC.Valid+(Cnconfirmed_04MARfinal$ROC.Valid.1*1.96), 
        length=0, angle=90, code=3,lwd=1,
        col=c("cyan4","chartreuse4","darkorchid","black","brown")[factor(Cnconfirmed_04MARfinal$Motype)])
#legend(35,0.759, title = "bf",
# legend = c("0.5","0.75","1"),
#  col = c(1:3),lty = "solid")
#legend("topleft", title = "", bty = "n",
#legend = c("p=0·5","p=0·75","p=1"),cex=0.9,
#pch = c(16, 17, 4)) 
mtext("Proportion of training data (%)",side=1,line=1,outer=TRUE,cex=1)
dev.off()
##################################################################################################################