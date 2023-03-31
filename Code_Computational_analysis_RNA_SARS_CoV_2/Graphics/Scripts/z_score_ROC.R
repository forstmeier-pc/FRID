setwd('~/R')
library('ggplot2')

#data = read.delim('../../Desktop/Research Lab/Bevilacqua Research Lab/work/z_score_output_ribozymes.txt', sep ='\t')

data = read.delim('../../Desktop/Research Lab/Bevilacqua Research Lab/work/FINAL_youden_scale_all_RNA_structures.txt', sep ='\t')

sens = c()
spec = c()
spec_one_minus = c()
youden = c()

#Y = sens+spec-1

for (i in 1:length(data$z.score)){
  sens <- c(sens, data$TP[i]/(data$FN[i] + data$TP[i]))
  spec <- c(spec, data$TN[i]/(data$TN[i]+data$FP[i]))
  spec_one_minus <- c(spec_one_minus, data$FP[i]/(data$TN[i]+data$FP[i]))
  youden <- c(youden, sens[i]+spec[i]-1)
}
print(sens)
print(spec)
print(spec_one_minus)

data['sens'] = sens
data['spec'] = spec
data['one_minus_spec'] = spec_one_minus
data['youden'] = youden

subset(subset(data, sens > 0), one_minus_spec<1)

good_data = (subset(subset(data, sens > 0), one_minus_spec<1))
good_data$z.score[1] <- '-2.0'

for (i in 1:length(good_data$z.score)){
  if (((-100*as.numeric(good_data$z.score[i])))%%10 != 0 && (good_data$z.score != -2.3)){
    print(good_data$z.score[i])
    print((-100*as.numeric(good_data$z.score[i])))
    print(((-100*as.numeric(good_data$z.score[i])))%%10)
    if (good_data$z.score[i] == '-2.3'){
      good_data$z.score[i]= '-2.3'
      
    } else if (good_data$z.score[i] == '-2.2'){
      good_data$z.score[i]= '-2.2'
    } else {
      good_data$z.score[i]= ''
    }
  }
}

p = ggplot(good_data) + 
  geom_point(mapping=aes(x=one_minus_spec, y=sens)) +
  geom_abline(size=1)+
  geom_path(mapping=aes(x=one_minus_spec, y=sens), size = 1) +
  geom_segment(x=data$one_minus_spec[which.max(data$youden)], 
               xend =data$one_minus_spec[which.max(data$youden)],
               y = data$one_minus_spec[which.max(data$youden)],
               yend = data$sens[which.max(data$youden)],
               size=1, linetype=2)+
  #geom_vline(xintercept = data$one_minus_spec[which.max(data$youden)], size =1) + 
  geom_text(mapping=aes(x=one_minus_spec, y=sens+.03, label=z.score, angle=90), size=6) +
  geom_text(mapping=aes(x=one_minus_spec, y=sens+.03, label=z.score, angle=90), size=6) +
  
  xlab("1-Specificity        FP/(FP+TN)")+
  ylab('Sensitvity        TP/(TP+FN)') +
  xlim(0.1,1.1)+
  ylim(0.1,1.1)

p=p+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text=element_text(size=18, color = 'black'),
        axis.title=element_text(size=24))+
  geom_label(mapping=aes(x=.8,y=.8, label = 'y=x'), size=8)+
  geom_label(mapping=aes(x = data$one_minus_spec[which.max(data$youden)]), 
             y=1.15*data$one_minus_spec[which.max(data$youden)], 
             label = 'J',
             size=8)+
  scale_y_continuous(breaks=seq(.2,1.1,.1))+
  scale_x_continuous(breaks=seq(.2,1.1,.1))
  
p
max(data$youden)
