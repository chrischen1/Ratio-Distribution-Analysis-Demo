#' wrapper for getting fold change, pvalue and FDR, by per cell line per time point
#' 
#' @param cnt p by n matrix for p genes across n samples
#' @param grp_table dataframe with 3 columns: group, condition and control
#'  group: contains information which treatment samples will be compared against control cases in each group
#'  condition: indicates type of treatment, replicates have same condition, do NOT use numbers only.
#'  control: TRUE for controls and FALSE for treatments
#'  order of well in samples annotation must be the same as the columns in count table
#' @param norm_method: Method of normalization, one if 'TMM','UQ','RLE','none'
#' @param combine_fdr T for combine FDR and p-values with group and F for compute pairwisely
#' @param w n by p matrix for n samples and p factors for batch effect correction from RUVSeq
#' @return a list, each element is a matrix containing p-value, LR,logCPM, LR and FDR for each group.treatment
edgeR_wrapper <- function(cnt,grp_table,norm_method = 'none',w = NULL,combine_fdr = F){
  library(edgeR)
  return_list <- list()
  library(edgeR)
  if(sum(rownames(grp_table) %in% colnames(cnt)) < nrow(grp_table)){
    warning('rownames in group table not compatible with colnames in count table')
    return(0)
  }
  grp_table$condition <- as.character(grp_table$condition)
  grp_table$group <- as.character(grp_table$group)
  
  cnt <- cnt[,rownames(grp_table)]
  design <- model.matrix(~condition,data = grp_table)
  # add RUV batch effect correction when w exists
  if(!is.null(w))  design <- cbind(design,w)
  y <- DGEList(counts=cnt, group=grp_table$condition)
  y <- calcNormFactors(y,method=norm_method)
  # Calculate overall dispersions when called first time
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  if(combine_fdr){
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2:(ncol(design)))
    lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt),]
    colnames(lrt_tab) <- gsub('logFC.condition','',colnames(lrt_tab))
    return(lrt_tab)
  }
  CommonDisp <- y$common.dispersion
  TagwiseDisp <- y$tagwise.dispersion
  for(group_slt in unique(grp_table$group)){
    grp_table_slt <- grp_table[grp_table$group==group_slt,]
    colname_control <- rownames(grp_table_slt)[grp_table_slt$control]
    for(condition_slt in unique(grp_table_slt$condition[!grp_table_slt$control])){
      colname_treatment <- rownames(grp_table_slt)[grp_table_slt$condition==condition_slt]
      grp_slt <- c(rep('control',length(colname_control)),rep('treatment',length(colname_treatment)))
      y_slt <- DGEList(counts = cnt[,c(colname_control,colname_treatment)],group  = grp_slt)
      design_slt <- model.matrix(~grp_slt,data=data.frame(grp_slt,stringsAsFactors = F))
      y_slt <- calcNormFactors(y_slt,method=norm_method)
      if(length(colname_control)==1 & length(colname_treatment)==1){
        # When both control and treatment lacking replicates, use overall dispersion instead
        y_slt$common.dispersion <- CommonDisp
        y_slt$tagwise.dispersion <- TagwiseDisp
      }else{
        y_slt <- estimateGLMCommonDisp(y_slt, design_slt)
        y_slt <- estimateGLMTagwiseDisp(y_slt, design_slt)
      }
      fit <- glmFit(y_slt, design_slt)
      lrt <- glmLRT(fit, coef=2:(ncol(design_slt)))
      lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt),]
      return_list[[paste(group_slt,condition_slt,sep = '.')]] <- lrt_tab
    }
  }
  return(return_list)
}

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
  library(ggplot2)
  ratio[ratio>max_ratio]=max_ratio
  xlabel = c(sprintf("%.1f", round(seq(0,5.85,0.5),2)),paste(">",max_ratio,sep=""))
  df = data.frame(ratio)
  g1<-ggplot(df, aes(x=ratio)) + geom_histogram(binwidth=0.05, color="black", fill="cornflowerblue") +
    geom_freqpoly(binwidth=0.05, size=0.5,col="gray24")+ theme_bw()  + 
    labs(title=title_name, x="Ratio", y ="Frequency",family="sans")+
    theme(plot.title = element_text(color="black", size=text_size,face="bold")) +
    theme(legend.text = element_text( size=text_size),legend.title = element_text( size=text_size), legend.key.size = unit(1,"cm")) +
    theme(axis.title = element_text(color="black", size=text_size),axis.text=element_text(size=text_size))+
    theme(panel.grid.major = element_line(colour = "grey"),panel.grid.minor.y = element_blank())+
    geom_vline(xintercept = left_line,color = "forestgreen", size=0.8)+
    geom_vline(xintercept = 1,color = "black", size=0.8)+
    geom_vline(xintercept = right_line,color = "deeppink2", size=0.8)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous("Ratio",limits = c(0,max_ratio),breaks = seq(0,max_ratio,0.5),  labels=xlabel,minor_breaks = seq(-0.3,max_ratio,0.1)) 
  return(g1)
}

plot_distribution_pairs <- function(r,cis_genes,title_name = '',left_line = 0.67,right_line = 1.5,max_ratio = 6){
  g <- multiplot(plotlist = list(plot_distribution(r[intersect(cis_genes,names(r))],title_name = 'Cis genes',left_line = left_line,right_line = right_line,max_ratio = max_ratio),
                                 plot_distribution(r[!names(r)%in%cis_genes],title_name = 'Trans genes',left_line = left_line,right_line = right_line,max_ratio = max_ratio)))
  return(g)
}

plot_scatter <- function(Fold_Change,Expre,P_Value,plot_title='',left_line=0.67,right_line=1.5,hide_black_dots =F,show_lines=T,cex=0.72,xlim=6,ylim=20){

  if(hide_black_dots){
    plot(log2(Fold_Change), log2(Expre),col='white',
         cex.main=2, main=plot_title,cex.main=1.35*cex,family="sans", pch=20, xlab= '', ylab="",cex.lab=1.35*cex, xlim=c(-xlim,xlim), ylim=c(0,xlim),cex.axis=1.35*cex)
    
  }else{
    plot(log2(Fold_Change), log2(Expre),cex=cex,main=plot_title,cex.main=1.35*cex,family="sans", pch=20, xlab= '', ylab="",cex.lab=1.35*cex, xlim=c(-xlim,xlim), ylim=c(0,20),cex.axis=1.35*cex)
  }
  title(xlab="log"["2"]~"(Fold Change)", cex.lab=1.35*cex,line=2.5,family="sans")
  title(ylab=("log"["2"]~"(Avg. EXP)"), cex.lab=1.35*cex,line=2.5,family="sans")

    if(show_lines){
    abline(v=log2(right_line),lwd=1.5,col='deeppink2')
    abline(v=log2(left_line),lwd=1.5,col='forestgreen')
  }
  abline(v=log2(1),lty=3,lwd=2,col='black')

  
  up_ind <- P_Value<.05 & log2(Fold_Change) > 0
  down_ind <- P_Value<.05 & log2(Fold_Change) < 0
  
  points(log2(Fold_Change)[up_ind], log2(Expre)[up_ind], pch=20, col="red",cex=cex)
  points(log2(Fold_Change)[down_ind], log2(Expre)[down_ind], pch=20, col="green",cex=cex)
}

plot_scatter_pairs <- function(Fold_Change,Expre,P_Value,cis_genes,plot_title='',left_line=0.67,right_line=1.5,hide_black_dots =F,show_lines=T,cex=1.2,xlim=6,ylim=20){
  par(mfrow=c(2,1),mar=c(4, 4, 2, 2), mgp=c(0, 1, 0))
  trans_genes <- names(Fold_Change)[!names(Fold_Change)%in%cis_genes]
  plot_scatter(Fold_Change = Fold_Change[cis_genes],plot_title='Cis genes',Expre = Expre[cis_genes],P_Value = P_Value[cis_genes],left_line = left_line,right_line = right_line)
  plot_scatter(Fold_Change = Fold_Change[trans_genes],plot_title='Trans genes',Expre = Expre[trans_genes],P_Value = P_Value[trans_genes],left_line = left_line,right_line = right_line)
 
  
}

