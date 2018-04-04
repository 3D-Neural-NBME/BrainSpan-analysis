library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

# This script should recreate plots in the manuscript.

final_annotated = read.delim("~/Dropbox (Personal)/3D_Neural_NBME:BrainSpan_analysis/final.merged.brainspan.organoid.correlation.matrix.txt")
orders_timepoint = read.delim("~/Dropbox (Personal)/3D_Neural_NBME:BrainSpan_analysis/orders_timepoint.txt")

# Select tissues and conditions to plot here
tissues = c("V1C","DFC","A1C","M1C")
conditions =  c("ConJ","Con14X","Con4X","Con1X","ConT")
#conditions =  c("Con7X","Con3X","Con6X","Con10X","Con16X","Con11X")
#conditions =  c("ConC","ConD","ConI","ConJ")
#conditions =  c("ConC","ConD","ConI","ConJ")
#conditions =  c("ConO","Con18X")
#conditions =  c("Con17X","Con3X","Con2X","ConP","Con12X","Con10X","Con5X","Con8X","Con13X","Con11X","Con15X","Con9X")
#conditions =  c("Con2X","Con13W")
#conditions =  c("ConI","ConG","ConJ","ConH","ConK","ConO","Con18X","Con17X","Con3X","Con7X","Con2X","ConP",'Con12X',"Con10X","Con6X","Con5X","Con8X")

final_annotated_0105 = subset(final_annotated, condition %in% conditions)

condition_orders <-data.frame(condition=conditions,condition_order = letters[1:length(conditions)],stringsAsFactors = F)

final_annotated_0105_conditions = merge(final_annotated_0105,condition_orders,by="condition",all.x=T)

final_annotated_0105_conditions$timepoint<-gsub("\\."," ",final_annotated_0105_conditions$timepoint)
orders_timepoint$timepoint<- gsub("\\."," ",orders_timepoint$timepoint)

# Point + line plots, e.g. Fig2D, Fig4B, and several supplemental figures 
for(k in tissues){
  tmp = subset(final_annotated_0105_conditions, tissue ==k)
  tmp = subset(tmp, timepoint %in% c("12 pcw", "16 pcw", "19 pcw", "37 pcw" ))
  tmp_orders_timepoint <- subset(orders_timepoint, timepoint %in% tmp$timepoint )
  limits<-aes(x=condition_order,ymax=Mean+SEM,ymin=Mean-SEM)
  
  p = ggplot(tmp,aes(y = Mean,x = condition_order,colour = timepoint, group = timepoint)) + geom_point(size = 4)+geom_line()+
    scale_x_discrete(labels = condition_orders$condition) +
    xlab("") + ylab("Mean correlation") + theme_bw()+
    ggtitle(k) + 
    theme(axis.text = element_text(size=13),legend.text=element_text(size=13),axis.title.y = element_text(size = 15), plot.title=element_text(hjust =0.5))+
    geom_errorbar(limits, width = 0.1)
  print(p)
  }


# Bar plots, e.g. Fig5C. SuppFig15A

summary_file = subset(final_annotated, final_annotated$timepoint %in% c( "12.pcw", "19.pcw", "37.pcw") & final_annotated$tissue %in%  c("V1C","DFC","A1C","M1C")
)

summary_file<-subset(summary_file, condition %in% c("ConI","ConG","ConJ","ConH","ConK","ConO","Con18X","Con17X","Con3X","Con7X","Con2X","ConP","Con12X","Con10X","Con6X",'Con5X',"Con8X") )
summary_file$condition_ordered = factor(summary_file$condition, levels =c("ConI","ConG","ConJ","ConH","ConK","ConO","Con18X","Con17X","Con3X","Con7X","Con2X","ConP","Con12X","Con10X","Con6X",'Con5X',"Con8X") )

limits<-aes(x=condition_ordered,ymax=Mean+SEM,ymin=Mean-SEM)

for(tiss in tissues){
  tmp = subset(summary_file, summary_file$tissue == tiss)
  p <- ggplot(tmp,aes(x=condition_ordered,y=Mean,fill=timepoint))+
    geom_bar(stat = 'identity',position = 'dodge', width = 0.7) +
    xlab("")+ylab("Mean correlation") +
    ggtitle(tiss) +
    theme_bw()+
    theme(axis.text = element_text(size=16),legend.text=element_text(size=13),panel.grid.major=element_blank(),axis.text.x = element_text(size=13 ),axis.title.y=element_text(size=18), plot.title=element_text(hjust = 0.5))+
    geom_errorbar(limits,width=0.3,position=position_dodge(width=0.7))+
    coord_cartesian(ylim=c(0.4, 0.65)) +
    scale_fill_manual(values = c( "#1B9E77","#7570B3" ,"#E7298A" )) 
  
  print(p)
}

# Heatmaps, not featured in manuscript
for(k in tissues){
  tmp = subset(final_annotated_0105_conditions, tissue ==k)
  tmp_orders_timepoint <- subset(orders_timepoint, timepoint %in% tmp$timepoint )
  
  p = ggplot(tmp,aes(y=timepoint_order,x=condition_order,fill=Mean))+ geom_tile() +
    scale_y_discrete(labels=tmp_orders_timepoint$timepoint) + 
    scale_x_discrete(labels=condition_orders$condition)+xlab("") +
    ylab("")+theme_classic() +
    ggtitle(k) +
    scale_fill_gradientn(colours=c("blue","dodgerblue","lightgreen","yellow","orange","red"))
  
  print(p)
}