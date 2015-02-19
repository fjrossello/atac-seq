


loadFragLens<-function(fn, nLinesSubsamp, maxLen, subSamplerFn){
## head results/ATACseq/viz/131217_SN651_0223_BH0Y9TADXX/pvalb_ATAC_rep1.bed
##   chr1    3000147 3000215 .       0       +

# Figure out how many lines in the whole BED file
   nLines<-system2("wc",paste(" -l ",fn,sep=""),stdout=TRUE)
   nLines<-as.numeric(strsplit(nLines, " ")[[1]][1])
   print(paste("nLines is ",nLines))

# Figure out what fraction of lines to keep
   fracKeep<-c(nLinesSubsamp/nLines)
   print(paste("frac keep: ",fracKeep))

# Subsample BED file and compute fragment lengths
   readCom<-paste("perl ",subSamplerFn," ",fracKeep," < ",fn)
   print(paste("Reading subsample of fragment lengths from ",fn))
   a<-scan(pipe(readCom))

   a<-a[a<=maxLen]

   return(a)
}


fitLengthDistr<-function(a, outFnBase, maxLen){

   nComp<-6
   maxIter<-5000
   initialGuessMeans<-c(50,200,380,560,740,920)
   # 1. sub-nucleosomal
   # 2. mono-
   # 3. di-
   # 4. tri-
   # 5. tetra-
   # 6. penta-
   initialGuessSD<-rep(30,nComp)

# Fit distribution with gaussian mixture model
   print("Fitting gaussian mixture components")
   library(mixtools)
   mixmdl<-normalmixEM(a, k=nComp, mu=initialGuessMeans,
                       sigma=initialGuessSD,maxit=maxIter)


   if (0) {
   i<-1
   while (i <= nComp) {
      print(paste(i,". lambda=",mixmdl$lambda[i],
                  " mu=",mixmdl$mu[i],
                  " sigma=",mixmdl$sigma[i]))
      i<-i+1
   }
   }


   print("Determine P(label|length) >= 90% cutoffs")
   posteriors<-matrix(nrow=maxLen,ncol=nComp,data=NA)

   len<-1
   while (len <= maxLen) {
      posteriors[len,]<-c(mixmdl$lambda * dnorm(len,mean=mixmdl$mu,sd=mixmdl$sigma))
      posteriors[len,]<-posteriors[len,]/sum(posteriors[len,])
      len<-len+1
   }
   posteriors<-cbind(posteriors, c(1:maxLen))


   outFn<-paste(outFnBase,".components.txt",sep="")
   sink(outFn)
   cat("n", "lambda", "mu", "sigma", "start_pr90", "end_pr90",sep="\t")
   cat("\n")

   i<-1
   starts<-vector(length=nComp,mode="integer")
   ends<-vector(length=nComp,mode="integer")
   while (i <= nComp) {
      starts[i]<-min(posteriors[posteriors[,i] > 0.9, (nComp + 1)])
      ends[i]<-max(posteriors[posteriors[,i] > 0.9, (nComp + 1)])
      cat(i, mixmdl$lambda[i], mixmdl$mu[i], mixmdl$sigma[i],
                  starts[i], ends[i],sep="\t")
      cat("\n")
      i<-i+1
   }
   sink()

   outFn<-paste(outFnBase,".length_distr_fit.pdf",sep="")
   pdf(outFn)
   plot(mixmdl, which=2)
   lines(density(a,bandwith=1), lwd=3)
   abline(v=starts,lwd=2,lty=1)
   abline(v=ends,lwd=2,lty=2,col="darkgrey")
   dev.off()
}


args <- commandArgs(trailingOnly = TRUE)
subSamplerFn<-args[3]

a<-loadFragLens(fn=args[1], nLinesSubsamp=100000, maxLen=1000, subSamplerFn=subSamplerFn)
fitLengthDistr(a,outFnBase=args[2], maxLen=1000)
