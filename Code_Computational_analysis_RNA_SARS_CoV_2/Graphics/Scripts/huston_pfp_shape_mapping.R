library(ggplot2)

setwd('~/R/PFPs_test')
library(R2easyR)


peter.ct.reformat = function(file.ct){
  
  con = file(file.ct)
  
  Lines = readLines(con)
  print(Lines)
  if (Lines[2]!=""){
    return(file.ct)
  }
  Lines = Lines[-which(Lines == "")]
  #print(Lines)
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



draw_Huston <- function(start, len){
  #start = 9331
  #len = 9415-9331-1
  #start = 13476
  #len = 13544-13476-1c(21706,21821)
  #start = 21706
  #len = 21821-21706-1
  filename = paste('SHAPE_Huston', start, sep='')
  #new_ct = peter.ct.reformat(paste(start, '_rstruc_Fold.ct',sep=''))
  #filename = paste('SHAPE_Huston_rstruc', start, sep='')

  if (start == 13476){
    new_ct = peter.ct.reformat(paste(start, '_SS_WT_Laederach.ct',sep=''))
    start = 13462
    len = 13544-13462-2
  } else {
    new_ct = peter.ct.reformat(paste(start, '.ct',sep=''))
    
  }
  Reactivity <- read.shape('SARS-CoV-2_SHAPE_Reactivity.txt')
  
  p9331_df<- read.ct(new_ct)
  
  p9331_df <-add.dot.bracket(p9331_df)
  p9331_df
  pk = r2easyR.pk_finder(p9331_df)
  pk
  p9331_df =  pk$r2easyR.dataframe
  
  p9331_df$Reactivity <- Reactivity[(start):(start+len)]
  p9331_df
  for (base in 1:length(p9331_df$Reactivity)){
    if (is.na(p9331_df$Reactivity[base]) == TRUE){
      p9331_df$Reactivity[base] <- 0
    }
    else if (p9331_df$Reactivity[base] < 0){
      p9331_df$Reactivity[base] <- 0
    }
  }
  p9331_df
  palettes <- r2easyR.palettes()
  
  p9331_df <- r2easyR.color(p9331_df, palettes$Reds.c, abs_reactivity_threshold = 0.2)#, manual.scale = c(0.2,2))
  
  r2easyR.write(filename,p9331_df, colors='circles')
  
  r2easyR.stem_editor(paste(filename,'.sto',sep=''))
  
  if (is.null(pk$pknot.list) == FALSE){
    r2easyR.pknot_drawer(pknot = pk$pknot.list, R2R.sto = paste(filename,'.sto',sep=''))
  }
  (p9331_df)
}

start_stops = list(c(9331, 9415), c(14774,14889),c(16043,16188),c(17829,17914), c(20814,20928),c(23027,23176))
start_stops = list(c(13476,13544))
start_stops = list(c(922, 1048),
                   c(19164,19256),
                   c(22541, 22843),
                   c(24518,24621),
                   c(26478, 26682))

start_stops = list(c(21706,21821))

start_stops = list(c(2818, 2914))

start_stops = list(c(23980,24068))

for (i in start_stops){
  print(i[2]-i[1]-1)
  draw_Huston(i[1],i[2]-i[1]-1)
}


