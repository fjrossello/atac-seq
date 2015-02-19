## Load CENTIPEDE library
library(CENTIPEDE)

## Load arguments
argv <- commandArgs(TRUE)
runmode<-argv[1]

## Load file locations, depending on "nocons" or "withcons" mode
if (grepl("^nocons",runmode)) {
   cutmatFn<-argv[2]
   annoFn<-argv[3]
   outBedFn<-argv[4]
   outPdfFn<-argv[5]

} else if (grepl("^withcons",runmode)) {

   cutmatFn<-argv[2]
   consFn<-argv[3]
   annoFn<-argv[4]
   outBedFn<-argv[5]
   outPdfFn<-argv[6]

} else {
   print("ERROR: run mode neither matches nocons nor withcons")
   stop()
}


## Read in files common to both run modes
cutmat<-read.table(cutmatFn)

annomat<-read.table(annoFn)


## setup scoring matrix
if (grepl("^nocons",runmode)) {
   scoremat<-cbind( rep(1,dim(annomat)[1]), annomat[,5])
} else if (grepl("^withcons",runmode)) {
   consScores<-read.table(consFn)
   scoremat<-cbind( rep(1,dim(annomat)[1]), annomat[,5], consScores$V1)
}


motifLen<-(annomat[1,"V3"] - annomat[1,"V2"])



##fpd140725_1210  If motif is on (-) strand have to switch first half of row with second half
matWidth<-ncol(cutmat)
cutmat[annomat$V6=="-",]<-cutmat[annomat$V6 == "-" ,c((matWidth/2 + 1):matWidth, 1:(matWidth/2))]
print(paste("Flipped fw and rv halves of the matrix for ",nrow(cutmat[annomat$V6 == "-",])," hits on negative strand"))


## Fit CENTIPEDE model
centFit<-fitCentipede(Xlist = list(DNase=as.matrix(cutmat)), Y=scoremat)

## add PWM, phylop (if applicable), and CENTIPEDE score to output file.
if (grepl("^withcons",runmode)) {
   annomat$V5<-paste("pwm=",annomat$V5,";phylop=",consScores$V1,";centipede=",centFit$PostPr,sep="")
} else {
   annomat$V5<-paste("pwm=",annomat$V5,";centipede=",centFit$PostPr,sep="")
}


## Print out sites with > 0.95 posterior probability are predicted sites.
predSites<-annomat[centFit$PostPr > 0.95,]
write.table(predSites,file=outBedFn,quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)


## Create images of cutsites and footprint profile.
pdf(outPdfFn)
imageCutSites(cutmat[order(centFit$PostPr),][c(1:100, (dim(cutmat)[1]-100):(dim(cutmat)[1])),])
plotProfile(centFit$LambdaParList[[1]],Mlen=motifLen)
dev.off()

print(paste("Finished properly: ", outBedFn))
