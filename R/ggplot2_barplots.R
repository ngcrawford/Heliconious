plot_bar <- function(gatkdata) { 
  gatkdata.sub <- gatkdata$TiTvVariantEvaluator[c(4, 13, 16),]
  print(gatkdata.sub)
  p <- ggplot(gatkdata.sub, aes(x=JexlExpression, y=tiTvRatio)) 
  p <- p + geom_bar(aes(fill=factor(JexlExpression))) 
  p <- p + opts(legend.position = "none", axis.text.x = theme_blank())
  p
  }