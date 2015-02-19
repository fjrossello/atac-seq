argv <- commandArgs(TRUE)
a<-read.table(argv[1])

pdf(argv[2])
hist(a$V1, breaks=2000, lty=1, main="BOWTIE2, transposon-trimmed ATAC", xlab="fragment length (nt)", lwd=2, xaxs="i", xlim=c(0,1000))
axis(side=1,at=seq(0,1000,by=100),tcl=-0.2,label=FALSE)
dev.off()

pdf(argv[3])
plot(density(a$V1, bw=1), lty=1, main="BOWTIE2, transposon-trimmed ATAC", xlab="fragment length (nt)", lwd=2, xaxs="i", log="y", xlim=c(0,1400), ylim=c(1E-5, 1E-2))
axis(side=1,at=seq(0,1400,by=100),tcl=-0.2,label=FALSE)
dev.off()
