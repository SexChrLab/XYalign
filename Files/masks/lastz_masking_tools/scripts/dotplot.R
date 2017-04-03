args = commandArgs(trailingOnly=TRUE)

dots = read.table(args[1], header=T)
jpeg(args[2])
plot(dots, type="l", xlab="UCSC chrY", ylab="UCSC chrY")
dev.off()
