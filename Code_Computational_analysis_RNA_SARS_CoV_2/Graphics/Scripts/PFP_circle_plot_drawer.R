setwd('~/R')


#peter.ct.reformat(paste(13476, '.ct',sep=''))


library(R2easyR)
library(ggplot2)
library(ggrepel)
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + (r * sin(tt)/4)
  return(data.frame(x = xx, y = yy))
}
peter.ct.reformat = function(file.ct){
  
  con = file(file.ct)
  
  Lines = readLines(con)
  print(Lines)
  Lines = Lines[-which(Lines == "")]
  print(Lines)
  for (i in 1:length(Lines)){
    Lines[i] <- gsub("\t\t", "   ", Lines[i])
    #print(i)
    #print(gsub("\t\t", "   ", Lines[i]))
  }
  Lines[1] <- "N   Nucleotide  N-1   N+1   BP   N"
  #Lines
  writeLines(Lines, file(paste(file.ct,'_edited.ct',sep='')))
  
  close(con)
  close(file(paste(file.ct,'_edited.ct',sep='')))
  return (paste(file.ct,'_edited.ct',sep=''))
}
drawLoops <- function(begin, index, pairs, p, dl, pk=FALSE){
  i = pairs[index]
  j = pairs[index + 1]
  diameter = abs(j-i)
  center = (diameter/2) + begin - 1 + as.numeric(i)
  dat <- circleFun(c(center,.5),diameter,npoints = dl)
  shade = 'black'
  line = 'solid'
  if (pk == TRUE){
    #shade = 'blue4'
    shade = 'red'
  }
  if (pk == 'maybe'){
    shade = 'red'
    line = 'dashed'
  }
  p <- p + 
    geom_line(mapping = aes(x=dat$x,y=dat$y),alpha = 0.8, size = 2, color = shade, linetype = line ) 
  #    scale_color_discrete('Plot')
  return(p)
}

drawLoopsAny <- function(begin,index, pairs, p, dl, pk=FALSE){
  i = pairs[index]
  j = pairs[index + 1]
  diameter = abs(j-i)
  center = (diameter/2) + begin - 1 + as.numeric(i)
  dat <- circleFun(c(center,.5),diameter,npoints = dl)
  shade = 'black'
  line = 'solid'
  if (pk == TRUE){
    #shade = 'blue4'
    shade = 'red'
  }
  if (pk == 'maybe'){
    shade = 'red'
    line = 'dashed'
  }
  p <- p + 
    geom_line(mapping = aes(x=dat$x,y=dat$y),alpha = 0.8, size = 2, color = shade, linetype = line ) 
  #    scale_color_discrete('Plot')
  return(p)
}

addWTLabel <- function(i,cols,p){
  xs = cols[[i]][[1]]
  labels = cols[[i]][[2]]
  indecies = cols[[i]][[3]]
  #print(length(xs))
  #print(length(labels))
  p <- p +
    geom_text(mapping=aes(x=xs),y=0,label=labels, fontface='bold',size=15)#+geom_text(mapping=aes(x=xs),y=-1, label=indecies)
  
  return (p)
}

make_figures <- function(begin, endd){
  #begin = 13476
  #endd = 13544
  #begin = 9331
  #endd = 9415
  #begin = 1
  #endd = 147
  #begin = 14774
  #endd = 14889
  #begin=20814
  
  #begin =922
  #endd = 1048
  #being = 19164
  #endd = 19256
  
  #begin = 22541
  #endd = 22843
  begin = 2818
  endd = 2914
  
  #data <- read.delim('ROI_'+as.character(begin)+'_info_w_strains_CERTAIN_AAs.txt')
  holy_path = paste('new_ROI_directory/ROI_',as.character((begin)),'_Pass.drct/',sep='')
  hp2 = paste(holy_path,'ROI_',as.character(begin),'_info_CERTAIN_AAs_stats_wobbles_excluded_NON_SYNOM_INCLUDED.txt',sep='')
  data <- read.delim(hp2)
  
  #Need to get y-axis level formatted and added to the dataframe
  
  #data=subset(data, type!='')
  #data=subset(data,type!='ORF1ab polyprotein')
  if (begin == 13476){
    ct_file_name = peter.ct.reformat(paste('PFPs_test/',paste(as.character(begin),'_SS_WT_Laederach.ct',sep=''),sep=''))
    #ct_file_name = peter.ct.reformat(paste('PFPs_test/',paste(as.character(begin),'_SS_ALT_Laederach.ct',sep=''),sep=''))
    begin=13462
    #ct_file_name = peter.ct.reformat(paste('PFPs_test/',paste(as.character(begin),'_SS.ct',sep=''),sep=''))
    
    
  } else {
    ct_file_name = peter.ct.reformat(paste('PFPs_test/',paste(as.character(begin),'.ct',sep=''),sep=''))
  }
  ct_data = read.ct(ct_file_name)
  #ct_data = read.ct('PFPs/1.ct_edited.ct')
  ct_data <- add.dot.bracket(ct_data)
  
  pk_pairs2 = r2easyR.pk_finder(ct_data)
  sequence = paste(ct_data$Nucleotide,collapse='')
  dot_struc = paste(pk_pairs2$r2easyR.dataframe$Dotbracket,collapse='')
  dot_struc = gsub('<','(',dot_struc)
  dot_struc = gsub('>',')',dot_struc)
  #sequence =  "TGCTGTAAATTTACTTACTAATATGTTTACACCACTAATTCAACCTATTGGTGCTTTGGACATATCAGCATCTATAGTAGCTGG"
  #dot_struc = '.<.<<<<<............((((((>>>>>.(((((((.......)))))))......))))))(((>...........))).'
  brackets = strsplit(dot_struc, "")
  
  pairs = c()
  
  
  for (ii in 1:length(brackets[[1]])){
    par_start = 0
    par_end = 0
    for (i in 1:length(brackets[[1]])){
      j = brackets[[1]][i]
      if (j == '('){
        par_start = i
      }
      if (j == ')'){
        par_end = i
      }
      if (par_start != 0 && par_end != 0){
        pairs <- c(pairs, c(par_start, par_end))
        brackets[[1]][par_start] <- '.'
        brackets[[1]][par_end] <- '.'
        break
      }
    }
  }
  
  pk_pairs = c()
  
  for (ii in 1:length(brackets[[1]])){
    par_start = 0
    par_end = 0
    for (i in 1:length(brackets[[1]])){
      j = brackets[[1]][i]
      if (j == '<'){
        par_start = i
      }
      if (j == '>'){
        par_end = i
      }
      if (par_start != 0 && par_end != 0){
        pk_pairs <- c(pk_pairs, c(par_start, par_end))
        brackets[[1]][par_start] <- '.'
        brackets[[1]][par_end] <- '.'
        break
      }
    }
  }
  
  sequence = strsplit(sequence, "")
  
  ddGs = c()
  for (i in data$ddG){
    j = strsplit(i," ")
    ddGs = c(ddGs, as.numeric(j[[1]][1]))
  }
  
  Loci = c()
  for (i in data$loci){
    Loci = c(Loci, as.numeric(i))
  }
  
  colors_chosen = c()
  
  for (i in 1:length(data$loci)){
    if (data$ALT[i] == 'T'){
      data$ALT[i] = 'U'
    }
    
    
    if (data$type[i] == 'Silent'){
      colors_chosen = c(colors_chosen, 'blue')
    } else if (data$type[i] == 'Missense Conservative'){
      colors_chosen = c(colors_chosen, 'green4')
    } else if (data$type[i] == 'Missense Non-conservative'){
      colors_chosen = c(colors_chosen, 'purple3')
    } else if (data$type[i] == 'Nonsense'){
      colors_chosen = c(colors_chosen, 'red')
    }
  }
  
  data$colors = colors_chosen
    
  p <- ggplot(data, aes(x=data$loci)) 
  #pk_pairs2 = unlist(pk_pairs2$pknot.list)
  
  pk_pairs2$pknot.list
  
  pk_pairs3 = c()
  
  for (i in pk_pairs2$pknot.list){
    print(i)
    for (ii in 1:(length(i)/2)){
      pk_pairs3 = c(pk_pairs3, i[ii])
      pk_pairs3 = c(pk_pairs3, i[length(i)-ii+1])
    }
  }
  
  pk_pairs3
  dl = length(data$WT)
  dlm = dl - 1
  
  for (i in pairs){
    index = which(pairs==i)
    if (index%%2 == 1){
      p = drawLoops(begin,index, pairs, p,dl)
      
      p = drawLoops(begin,index, pk_pairs3, p,dl, pk=TRUE)
      if (length(alt_pk)!= 0){
        p = drawLoopsAny(begin, index, alt_pk,p,dl, 'maybe')
      }
    }
  }
  
  lst = c(2,3)
  colors = rainbow(length(lst))
  #p=p+
  #  geom_text(mapping = aes(x=data$loci, y=-1), label = data$ALT, colour=data$type)+
  #  scale_color_manual(name = 'amino acids',values = c('Missense Conservative'='blue','Missense Non-conservative'='green4','Silent'='purple3', 'Nonsense'='red'))
  p=p+
    geom_text(mapping = aes(x=(as.integer(data$loci)), y=-1),
              position=position_stack(vjust=-1.95),
              size=15, 
              label=data$ALT, 
              color=data$colors)+
    geom_segment(x=as.integer(data$loci), 
                 xend=as.integer(data$loci), 
                 y=-2,yend=-.8, color='black',size=1,arrow=arrow(length = unit(2, "mm"), type='closed'))+ 
    geom_segment(x=5+begin,
                 xend=7+begin,
                 y=-5,
                 yend=-5,
                 color='blue',
                 size=5)+
    geom_segment(x=5+begin,
                 xend=7+begin,
                 y=-6.1,
                 yend=-6.1,
                 color='red',
                 size=5)+
    geom_segment(x=5+begin,
                 xend=7+begin,
                 y=-7.2,
                 yend=-7.2,
                 color='green4',
                 size=5)+
    geom_segment(x=5+begin,
                 xend=7+begin,
                 y=-8.3,
                 yend=-8.3,
                 color='purple3',
                 size=5)+
    expand_limits(y=c(NA,-8))+
    
    geom_text(mapping=aes(x=8+begin, y=-4.9, hjust=0), label='Silent',size=15)+
    geom_text(mapping=aes(x=8+begin, y=-6.0, hjust=0), label='Nonsense',size=15)+
    geom_text(mapping=aes(x=8+begin, y=-7.1, hjust=0), label='Missense Conservative',size=15)+
    geom_text(mapping=aes(x=8+begin, y=-8.2, hjust=0), label='Missense Non-conservative',size=15)

  
  if (begin==13462){
    p = p +
      geom_segment(x=begin-.4,
                   xend=begin+.4,
                   y=-1.3,
                   yend=-1.3,
                   color='black',
                   size=5, alpha = .03)+
      geom_segment(x=21.4+begin,
                   xend=18.2+begin+.4,
                   y=-1.3,
                   yend=-1.3,
                   color='red',
                   size=5, alpha = .03) +
      geom_segment(x=begin-.4,
                   xend=6+begin+.4,
                   y=0,
                   yend=0,
                   color='yellow',
                   size=15, alpha = .04)+
      geom_text(mapping = aes(x=begin,y=.7),label = 'Slippery Sequence', size = 15)
    idkidc<-0
    while (idkidc < 6){
      p = p+
        geom_segment(x=1+begin+(3*idkidc)-.4,
                     xend=3+begin+(3*idkidc)+.4,
                     y=-1.3,
                     yend=-1.3,
                     color='black',
                     size=5, alpha = .03)
      
      idkidc = idkidc+1
      
    }
    idkidc<-0
    while (idkidc < 27){
      p = p+
        geom_segment(x=begin+(3*idkidc)-.4,
                   xend=2+begin+(3*idkidc)+.4,
                   y=-1.75,
                   yend=-1.75,
                   color='black',
                   size=5, alpha = .03)
        
      idkidc = idkidc+1
      
    }

    p = p +geom_text(mapping = aes(x=begin,y=-1.25),label = 'ORF1a', size = 15, hjust=1.35)
    p = p +geom_text(mapping = aes(x=begin,y=-1.7),label = 'ORF1ab', size = 15, hjust=1.32)
    
  }
  
  p
  

  if (begin == 13462){
    num_blocks = length(sequence[[1]])%/%length(data$WT)+1
  } else {
    num_blocks = length(sequence[[1]])%/%length(data$WT)+2
  }
  
  
  
  cols = list()
  for (i in 1:num_blocks){
    cols <- c(cols, list(i))
  }
  
  cols
  
  i=1
  
  while (i <= (num_blocks)){
    if (i == num_blocks){
      x = (endd-dlm-1):(endd-1)
      label=sequence[[1]][(length(sequence[[1]])-dlm):length(sequence[[1]])]
    }
    if (i != num_blocks){
      if (i==1){
        x = (begin+(dlm*(i-1))):(begin+(dlm*(i)))
        label=sequence[[1]][(dlm*(i-1)+1):((dlm*(i))+1)]
        indecies=seq((dlm*(i-1)+1),((dlm*(i))+1))
      } else {
        x = (begin+(dlm*(i-1))):(begin+(dlm*(i)))
        label=sequence[[1]][(dlm*(i-1)+1):((dlm*(i))+1)]
        indecies=seq((dlm*(i-1)+1),((dlm*(i))+1))
      }
      
    }
    
    cols[[i]] = list(x,label, indecies)
    
    i=i+1
  }
  cols
  
  i = 0
  
  while(i<(length(cols))){
    i = i+1
    #p<-p+geom_text(mapping=aes(x=cols[[i]][[1]]),y=0,label=cols[[i]][[2]])
    p = addWTLabel(i,cols,p)
  }
  p
  
  
  loci_labels = seq(10,endd-begin+1, by =10)
  loci_labels_pos = c()
  for (l in loci_labels){
    loci_labels_pos <- c(loci_labels_pos, l+begin-1)
  }
  loci_labels_pos
  master_loci_labels = data.frame(loci_labels, loci_labels_pos)
  
  
  while (length(master_loci_labels$loci_labels)<length(Loci)){
    master_loci_labels[nrow(master_loci_labels)+1,] = c(master_loci_labels[1],master_loci_labels$loci_labels_pos[1])
  }
  master_loci_labels
  p=p+geom_text(mapping = aes(master_loci_labels$loci_labels_pos, y=-1.15,label=master_loci_labels$loci_labels),color='black', size = 15)
  #grey55
  
  #p=p+
  #  geom_rect(xmin=begin+2,xmax=begin+3,ymin=2,ymax=3,fill='red')
  
  #p
  #p=p+ geom_text(mapping = aes(x=seq(9331,9362,by=1), y=0), label=sequence[[1]][0:32])+geom_text(mapping = aes(x=seq(9362,9393,by=1), y=0), label=sequence[[1]][32:63])+geom_text(mapping = aes(x=seq(9384,9415,by=1), y=0), label=sequence[[1]][54:85])
  #p=p+scale_y_continuous(breaks = seq(-30,30, by=1)) +ylim(-10,10)\
  
  #for FSE make limits -8 to 8
  p=p+scale_y_continuous(limits = c(-8.5,16), breaks = c(-1.5,-2.5,-3.3,-4.5), labels=c("SNP", "PseudoSNitches", 'RiboSNitches', 'PseudoSNitches & \nRiboSNitches'))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none',
          #axis.title.y = element_text(size = 55),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  #axis.line = element_line())
  #p + scale_color_discrete("Amino Acid Mutation Types")
  
  p
  
  #out_filename = paste(as.character(begin),'_PFP_SPOT-RNA.png',sep='')
  #out_filename = paste(as.character(begin),'_PFP_WT_Laederach.png',sep='')
  out_filename = paste(as.character(begin),'_PFP_4.svg',sep='')
  
  ggsave(out_filename,width=(endd-begin)/.9,height=30, units = 'cm', limitsize = FALSE)
  return (p)
}  


#make_figures(13476, 13544)


make_figures(13476, 13544)
alt_pk = c()


alt_pk = c(6,33,7,32,8,31,9,30,10,29)

start_stops = list(c(13476,13544), 
                   c(9331, 9415), 
                   c(14774,14889),
                   c(16043,16188),
                   c(17829,17914), 
                   c(20814,20928),
                   c(23027,23176),
                   c(922, 1048),
                   c(19164,19256),
                   c(21706,21821),
                   c(22541, 22843),
                   c(24518,24621),
                   c(26478, 26682))

for (i in start_stops){
  print(i)
  print(i[1])
  print(i[2])
  p = make_figures(i[1],i[2])
  p
}

#26478
start_stops = list(c(26478, 26682))
alt_pk = c(124,176,125,175,126,174,127,173,129,171,130,170)

#24518
start_stops = list(c(24518,24621))
alt_pk = c(12,93,21,84,24,81 )

#22541
start_stops = list(c(22541, 22843))
alt_pk = c(162,257,163,256,164,255,165,254,168,251)

#21706
start_stops = list(c(21706, 21821))
alt_pk = c(83,109, 77,113,78,114)

#922
start_stops = list(c(922, 1048))
alt_pk = c(14,25,15,24)

#19164
start_stops = list(c(19164,19256))
alt_pk = c(13,80, 15,84,16,83)

#9331
start_stops = list(c(9331, 9415))
alt_pk = c(3,68, 9, 26,4,67,5,66,1,70)

#13476
alt_pk = c()
start_stops = list(c(13476,13544))

#14774
start_stops = list(c(14774,14889))
alt_pk = c(11,18, 63,98, 64, 97, 65,96, 66,95)

#16043
start_stops = list(c(16043,16188))
alt_pk = c(19,35,22,32,23,31)

#17829
start_stops = list(c(17829,17914))
alt_pk = c(5,63,4,64,3,65,7,61)

#20814
start_stops = list(c(20814,20928))
alt_pk = c(5,35)

#23027
start_stops = list(c(23027,23176))
alt_pk = c(47,97,46,98,45,99)

#1
start_stops = list(c(1,147))
alt_pk = c(6,17)

#2818
start_stops = list(c(2818, 2914))
alt_pk = c(1,20,2,19,3,18,5,16,6,15,7,14)

#23980
start_stops = list(c(23980,24068))
alt_pk = c(7,42,8,41,5,44,4,45)

#start = 13476
#stopp = 13544


#make_figures(9331, 9415)
