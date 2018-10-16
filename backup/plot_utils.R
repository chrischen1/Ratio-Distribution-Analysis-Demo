multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
    }
  }
}

plot_distribution = function(ratio , title_name ,left_line,right_line,max_ratio = 6,text_size=12){
  ratio[ratio>max_ratio]=max_ratio
  library(ggplot2)
  df = data.frame(ratio)
  g1<-ggplot(df, aes(x=ratio)) + geom_histogram(binwidth=0.05, color="black", fill="cornflowerblue") +
    geom_freqpoly(binwidth=0.05, size=1,col="gray24")+ theme_bw()  + xlim(0, max_ratio)+
    # scale_fill_gradient("Frequency",low = "green", high = "red") +
    labs(title=title_name, x="Ratio", y ="Frequency",family="sans")+
    theme(plot.title = element_text(color="black", size=text_size,face="bold")) +
    theme(legend.text = element_text( size=text_size),legend.title = element_text( size=text_size), legend.key.size = unit(1,"cm")) +
    theme(axis.title = element_text(color="black", size=text_size),axis.text=element_text(size=text_size))+
    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted"),panel.grid.minor.y = element_blank())+
    geom_vline(xintercept = c(left_line,1,right_line),color = "black", size=0.8)+theme(plot.title = element_text(hjust = 0.5))
  g1 = g1 + annotate("text", x = c(0.4,1.2,2.4), 
                     y = ggplot_build(g1)$layout$panel_params[[1]]$y.range[2]*0.98, 
                     label = sprintf("%.2f",round(c(left_line,1,right_line),2)),size=text_size*0.3,family="sans")
  return(g1)
}

plot_distribution_pairs <- function(r,cis_genes,title_name = '',left_line = 0.67,right_line = 1.5,max_ratio = 6){
  g <- multiplot(plotlist = list(plot_distribution(r[intersect(cis_genes,names(r))],title_name = 'Cis genes',left_line = left_line,right_line = right_line,max_ratio = max_ratio),
                                 plot_distribution(r[!names(r)%in%cis_genes],title_name = 'Trans genes',left_line = left_line,right_line = right_line,max_ratio = max_ratio)))
  return(g)
}

plot_scatter <- function(Fold_Change,Expre,P_Value,plot_title='',left_line=0.67,right_line=1.5,hide_black_dots =F,show_lines=T,cex=0.72,xlim=6,ylim=20){
  # Red indicates P_Value<0.05 and log2Fold_Change<-1, green is P_Value<0.05 and log2Fold_Change>1)
  # red indicates up regulated, green is downregulated

  if(hide_black_dots){
    plot(log2(Fold_Change), log2(Expre),col='white',
         cex.main=2, main=plot_title,cex.main=1.35*cex,family="sans", pch=20, xlab= '', ylab="",cex.lab=1.35*cex, xlim=c(-xlim,xlim), ylim=c(0,xlim),cex.axis=1.35*cex)
    
  }else{
    plot(log2(Fold_Change), log2(Expre),cex=cex,main=plot_title,cex.main=1.35*cex,family="sans", pch=20, xlab= '', ylab="",cex.lab=1.35*cex, xlim=c(-xlim,xlim), ylim=c(0,20),cex.axis=1.35*cex)
  }
  title(xlab="log"["2"]~"(Fold Change)", cex.lab=1.35*cex,line=2.5,family="sans")
  title(ylab=("log"["2"]~"(Average Expression)"), cex.lab=1.35*cex,line=2.5,family="sans")
  #lines
  if(show_lines){
    abline(v=log2(right_line),lty=3,lwd=2,col='black')
    text(log2(2)+0.5,0, paste0("",right_line), col = "black", adj = c(0, -.1),cex=1.5*cex)
    abline(v=log2(left_line),lty=3,lwd=2,col='black')
    text(log2(0.5)-0.9,0, paste0("",left_line), col = "black", adj = c(0, -.1),cex=1.5*cex)
  }
  abline(v=log2(1),lty=3,lwd=2,col='black')
  # text(0,0, paste("",0), col = "black", adj = c(0, -.1),cex=1.5*cex)
  
  
  up_ind <- P_Value<.05 & log2(Fold_Change) > 0
  down_ind <- P_Value<.05 & log2(Fold_Change) < 0
  
  points(log2(Fold_Change)[up_ind], log2(Expre)[up_ind], pch=20, col="red",cex=cex)
  points(log2(Fold_Change)[down_ind], log2(Expre)[down_ind], pch=20, col="green",cex=cex)
}

plot_scatter_pairs <- function(Fold_Change,Expre,P_Value,cis_genes,plot_title='',left_line=0.67,right_line=1.5,hide_black_dots =F,show_lines=T,cex=1.2,xlim=6,ylim=20){
  par(mfrow=c(2,1))
  trans_genes <- names(Fold_Change)[!names(Fold_Change)%in%cis_genes]
  plot_scatter(Fold_Change = Fold_Change[cis_genes],plot_title='Cis genes',Expre = Expre[cis_genes],P_Value = P_Value[cis_genes],left_line = left_line,right_line = right_line)
  plot_scatter(Fold_Change = Fold_Change[trans_genes],plot_title='Trans genes',Expre = Expre[trans_genes],P_Value = P_Value[trans_genes],left_line = left_line,right_line = right_line)
  par(mfrow=c(1,1))
  
}

