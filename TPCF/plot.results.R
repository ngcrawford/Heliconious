results$b_no_log <- 10^results$Bin
results <- subset(results, Bin != 0)
results <- results[complete.cases(results), ]
results$cat <- do.call(paste, c(results[c("Comparison", "Chrm")], sep = ""))
p <- ggplot(results, aes(x=b_no_log, y=Mean, color=cat))
p <- p + geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.1)
p <- p + geom_line() + geom_point()
p <- p + scale_x_log10( name='Bin Size (base-pairs)', breaks=c(1e+04, 1e+05, 1e+06, 1e+07))
p <- p + theme_bw()
p <- p + theme(legend.position="none") # no legend
p <- p + scale_y_continuous(name= expression(1 + bar(xi) %+-% SEM)) 
p <- p + ggtitle("Cydno, Melpo, Pachi") + theme(plot.title = element_text(lineheight=.8, face="bold"))
show(p)
