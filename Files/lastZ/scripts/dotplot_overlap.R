args = commandArgs(trailingOnly=TRUE)

dots = read.table(args[1], header=T)
masked_dots = read.table(args[2], header=T)
jpeg(args[3])
plot(dots, type="l", xlab="UCSC chrY", ylab="UCSC chrY")
lines(masked_dots, col="red")
dev.off()
