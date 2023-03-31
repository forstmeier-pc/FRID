setwd('~/R')


library(R2easyR)
library(ggplot2)
#library(tidyverse)
#library(extrafront)

data = read.csv('NC_0455122Severeacut.FINAL.csv')

p <- ggplot(data, aes(x=data$X5..Boundary, color = as.character(as.integer(as.logical(Number.of.PKs))))) + 
  geom_segment(x=0,
               xend=29903,
               y=0,
               yend=0, 
               color='black',
               size=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color='black',size=20),
        #axis.line.y = element_line(),
        #axis.ticks.y = element_blank(),
        axis.title.y = element_text(color='black',size=20),
        axis.text.x = element_text(vjust = 17, color='black',size=20),
        #axis.title.x = element_text(size = 30, vjust=35),
        axis.ticks.x = element_line(),
        axis.title.x=element_blank(),
        legend.text = element_text(size=30),
        legend.position = "none")+ 
        #legend.position = c(0.5,0.15) )+
        #axis.ticks.y = element_blank().
        #axis.text.y=element_blank(),
        #axis.text.x = element_text(size = 10))+
  ylim(-1,2.26)+
  xlim(0,30000)+
  scale_color_manual(name = '',
                     values=c("0"="red", "1"="blue","2"='green', '3'='green'),
                     labels = c('ScanFold Structures without Predicted Pseudoknots',"
                                              ScanFold Structures with Predicted Pseudoknots (non-conserved)",
                                              'ScanFold Structures with Predicted Pseudoknots (conserved)', 
                                              "sgRNA 3' Ends")) +
  scale_y_continuous(name=paste('\U0394','G (kcal/mol)         ',sep=''), limits = c(-.3,1.26), breaks=c(0,.125,.25), labels = c('0','-25','-50'))

p = p +
  geom_segment(x=13496, xend=13476,
               y=.48,yend=0.005, linetype='dashed',
               arrow=arrow(length = unit(2, "mm"), type='closed'))

mfes = read.delim('middle_pos.txt',sep='\t')
for (i in 1:length(mfes$dg)){
  p = p+
    geom_segment(x=mfes$middle[i], xend=mfes$middle[i],
                 y=0,
                 yend=mfes$dg[i]/(-200),
                 size=.8,
                 color='black')
}

pfps = c(922,2818,9331,14774, 16043, 17829,19164,20814,21706,22541, 23027,23980,24518,26478)
for (i in 1:length(data$X5..Boundary)){

  p = p + geom_segment(x=data$X5..Boundary[i],
                       xend=data$X3..Boundary[i],
                       y=0,yend=0, size = 3)
  if (as.character(as.integer(as.logical(data$Number.of.PKs[i]))) != "0"){
    p = p + geom_segment(x=data$X5..Boundary[i],
                         xend=data$X3..Boundary[i],
                         y=0,yend=0, size = 3, color='blue')
  }
  if (data$X5..Boundary[i] %in% pfps){
    p = p + geom_segment(x=data$X5..Boundary[i],
                         xend=data$X3..Boundary[i],
                         y=0,yend=0, size = 3, color='green')
  }
}


gff = read.delim('R_script_gff_SARS_CoV_2.txt', sep='\t')
for (i in 1:length(gff$name)){
  print(gff$start[i])
  print(gff$stop[i])
  print(gff$name[i])
  p = p +
    geom_segment(x=gff$start[i],
                 xend=gff$stop[i],
                 y=(.1+(gff$level[i]/10))+.08,
                 yend=(.1+(gff$level[i]/10))+.08,
                 size = 1.25, color='black') +
    geom_text(mapping = aes(), 
              x=((as.integer(gff$start[i])+as.integer(gff$stop[i])))/2, 
              y=(.1+(gff$level[i]/10))+.08 +.025,
              color='black', 
              label=gff$name[i],
              size=5)
}

p = p +
  geom_text(mapping = aes(),
            x=133, 
            y=(.2+.08 +.025),
            color='black', 
            label='  ORF1ab',
            size=5) +
  geom_text(mapping = aes(), 
            x=266/2, 
            y=(.1+(2/10))+.08 +.025,
            color='black', 
            label='ORF1a',
            size=5)

sgRNA = read.delim("sgRNAs_SARS_CoV_2_kim.txt",sep='\t')
for (i in 1:length(sgRNA$start)){
  p = p + geom_segment(x=sgRNA$start[i],
               xend=sgRNA$start[i],
               y=-.0425,
               yend=-0.01,
               color='black',
               size=1,
               arrow=arrow(length = unit(2, "mm"), type='closed'))
}

p = p + 
  geom_segment(x=500,
                     xend=500,
                     y=-.321,
                     yend=-0.285,
                     color='black',
                     size=1,
                     arrow=arrow(length = unit(2, "mm"), type='closed')) +
  geom_segment(x=300,
               xend=700,
               y=-.15,
               yend=-0.15,
               color='red',
               size=5)+
  geom_segment(x=300,
               xend=700,
               y=-.2,
               yend=-.2,
               color='blue',
               size=5)+
  geom_segment(x=300,
               xend=700,
               y=-.25,
               yend=-.25,
               color='green',
               size=5) +
  geom_text(mapping=aes(x=1000, y=-0.145, hjust=0), label='ScanFold Structures without Predicted Pseudoknot', color = "black",size=5) +
  geom_text(mapping=aes(x=1000, y=-0.195, hjust=0), label='ScanFold Structures with Predicted Pseudoknot (non-conserved)', color = "black",size=5) +
  geom_text(mapping=aes(x=1000, y=-0.245, hjust=0), label='ScanFold Structures with Predicted Pseudoknot (conserved)', color = "black",size=5) +
  geom_text(mapping=aes(x=1000, y=-0.299, hjust=0), label="sgRNA TRS-L Attachment Site", color = "black",size=5)


prot = read.delim("protein_labels_SARS_CoV_2.txt", sep='\t')
for (i in 1:length(prot$start)){
  p = p +
    geom_segment(x=prot$start[i],
                 xend =prot$stop[i],
                 y=(.1+(as.integer(prot$level[i])/10))+.08,
                 yend=(.1+(as.integer(prot$level[i])/10))+.08,
                 size=3, color = "black")
  p = p +
    geom_text(mapping = aes(), 
              x=((as.integer(prot$start[i])+as.integer(prot$stop[i])))/2, 
              y=(.1+(as.integer(prot$level[i])/10))+.08 +0.025,
              color='black', 
              label=prot$name[i],
              size=5)
}


p=p+geom_segment(x=-1500,xend=-1500,y=-.001,yend=.25,color='black',size=1)


wavy_data = read.csv("SARS_genome_wavy_graph_background_data.csv")
Y_adj = 0.8
p=p+
  geom_segment(x=0,xend=1000,
               y=Y_adj,
               yend=(Y_adj+((.25/8)*3)),
               color = 'grey',
               size = 2) +
  geom_segment(x=1000,xend=3000,
               y=(Y_adj+((.25/8)*3)),
               yend=(Y_adj+((.25/8)*8)),
               color = 'grey',
               size = 2) +
  geom_segment(x=3000,xend=5000,
               y=(Y_adj+((.25/8)*8)),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=5000,xend=7000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=7000,xend=9000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=9000,xend=11000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=11000,xend=13000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=13000,xend=15000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*2)),
               color = 'grey',
               size = 2) +
  geom_segment(x=15000,xend=17000,
               y=Y_adj+((.25/8)*2),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=17000,xend=19000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*2)),
               color = 'grey',
               size = 2) +
  geom_segment(x=19000,xend=21000,
               y=(Y_adj+((.25/8)*2)),
               yend=(Y_adj+((.25/8)*5)),
               color = 'grey',
               size = 2) +
  geom_segment(x=21000,xend=23000,
               y=(Y_adj+((.25/8)*5)),
               yend=(Y_adj+((.25/8)*5)),
               color = 'grey',
               size = 2) +
  geom_segment(x=23000,xend=25000,
               y=(Y_adj+((.25/8)*5)),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=25000,xend=27000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*4)),
               color = 'grey',
               size = 2) +
  geom_segment(x=27000,xend=29000,
               y=(Y_adj+((.25/8)*4)),
               yend=(Y_adj+((.25/8)*3)),
               color = 'grey',
               size = 2) +
  geom_segment(x=29000,xend=30000,
               y=(Y_adj+((.25/8)*3)),
               yend=(Y_adj+((.25/8)*0)),
               color = 'grey',
               size = 2) 
  
  
p
