library(ggplot2)
library(dplyr)
library(rootSolve)


to_tenth_power_labels <- function(x_labels){
  return(parse(text=paste0("10^",x_labels)))
}

expected_number_of_barcodes <- function(n_reads,n_barcodes){
  return((1-((n_barcodes-1)/n_barcodes)^n_reads)*n_barcodes)
}

text_size=10
custom_theme <- theme_bw()+#starts black white theme as starting point
  theme(axis.title.y=element_text(hjust=0.5,vjust=0.5,size=text_size),#Sets y axis title to be horizontal and appropriate size
        axis.title.x=element_text(size=text_size,angle=0),#sets x axis title
        axis.text.x = element_text(size=text_size),axis.text.y = element_text(size=text_size),#sets x axis text size
        legend.text = element_text(size=text_size),#sets legend text size
        panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),#removes the minor panel grid for a cleaner look
        legend.title = element_text(size=text_size),#sets legend title
        legend.key.size = unit(0.5,'cm'),#sets the sie of the legend keys to be a half centimeter
        strip.background.x=element_rect(fill='#FFFFFF'),strip.background.y = element_rect(fill="#FFFFFF"),#sets facet backgrouns to white
        strip.text.x=element_text(size=text_size),strip.text.y = element_text(size=text_size))#sets facet title square backgrounds to be white


#we set the directory for the data
#due to the size limitations of github we cannot include the sequencing data with this repo, so we would advise you to 
#create your won data folder and place the seuqncing data in that folder. 
#instructions on how to access sequencing data can be found within the article itself.
data_directory <- "your/data/folder/here/"

#and the directory for the outputs
figure_output_directory <- "your/output/folder/here/"

#we load the table of barcode frequencies, which was output by the barcode_extractor script  
barcodes_frequencies <- read.table(paste0(data_directory,"barcode_table.tsv"),header = T)
barcodes_frequencies[["frequency"]] <- as.numeric(barcodes_frequencies[["frequency"]])

#we set the number of phages in as detected by CFU
n_phages <- 1.8*10^7

#we check how many of the barcodes meet the correct format
#the format of the barcode is A/T/G/C for 5, then TTT, then A/T/C/G for five, then TTT, then A/T/C/G for five. So
ATCG <- '[A/T/C/G/]'
spacer <- 'TTT'
pattern_to_grep <- paste0(c(rep(ATCG,5),spacer,rep(ATCG,5),spacer,rep(ATCG,5)),collapse = '')


#for all the sequences where we got multiple probably barcodes we check if they have one that meets our criteria, then include that one if possible
indeces_for_multi_barcodes <- grep(barcodes_frequencies[["barcode"]],pattern = '/')
for(i in indeces_for_multi_barcodes){
  potential_barcodes <- strsplit(barcodes_frequencies[i,"barcode"],split = '/')[[1]]
  if(sum(nchar(potential_barcodes)==21)>0){
    if(sum(grepl(pattern_to_grep,potential_barcodes))>0){
      barcodes_frequencies[i,"barcode"] <- potential_barcodes[grepl(pattern_to_grep,potential_barcodes)][1]
    }else{
      barcodes_frequencies[i,"barcode"] <- potential_barcodes[1]
    }
  }
}

#first we check how many are length 21 and then we remove the ones that aren't
barcode_length_distribution <- table(nchar(barcodes_frequencies[,"barcode"]))
print(paste0(barcode_length_distribution['21']/sum(barcode_length_distribution)*100,' % of the barcodes extracted are of length 21'))
barcodes_frequencies <- barcodes_frequencies[nchar(barcodes_frequencies[["barcode"]])==21,]

#we grep for the pattern to ensure the correct format is found for each barcode
pattern_match <- grepl(pattern = pattern_to_grep,barcodes_frequencies[["barcode"]])
print(paste0(mean(pattern_match)*100,' % of the length 21 barcodes meet the structure we want'))

#we subset to only the ones that meet our pattern
barcodes_frequencies <- barcodes_frequencies[pattern_match,]

#now let's make a plot of often we see things repeated
plotdat <- as.data.frame(table(barcodes_frequencies[,"frequency"]))
plotdat[["Var1"]] <- as.numeric(as.character(plotdat[["Var1"]]))


hidden_vars <- data.frame()
#we make a waterfall plot
for(row in seq(1,nrow(plotdat))){
  if(row>1){
    hidden_freq <- sum(plotdat[seq(1,row-1),c("Freq")])
    hidden_vars[row,c('Var1','Freq')] <- c('Var1'=plotdat[row,"Var1"],'Freq'=hidden_freq)
  }else{
    hidden_vars[row,c('Var1','Freq')] <- c('Var1'=plotdat[row,"Var1"],'Freq'=0)
  }
}

hidden_vars <- cbind(hidden_vars,'color'=rep('Hidden',nrow(hidden_vars)))

plotdat <- cbind(plotdat,'color'=rep('Shown',nrow(plotdat)))

plotdat <- rbind(plotdat,hidden_vars)

plotdat[["color"]] <- factor(plotdat[["color"]],levels=c('Shown','Hidden'))

#we set the break to the frequency where 99% of barcodes are included
cumsum_table <- table(barcodes_frequencies[,"frequency"])
#we calculate the ways reads are distributed - just for fun
read_distribution <- as.numeric(names(cumsum_table))*cumsum_table
read_distribution_percent <- read_distribution/sum(barcodes_frequencies[["frequency"]])*100

cumsum_table <- cumsum_table[as.character(sort(as.numeric(names(cumsum_table))))]
cum_freq <- cumsum(cumsum_table)/sum(cumsum_table)
x_upper_lim <- max(as.numeric(names(cum_freq)[cum_freq<0.99]))

y_breaks <- seq(0,100,10)

#we convert to percentages
plotdat <- plotdat %>% group_by(Var1,color) %>% summarise('Percentage'=Freq/nrow(barcodes_frequencies)*100)
plotdat <- as.data.frame(plotdat)

p <- ggplot()+geom_bar(plotdat,mapping = aes(x=Var1,y=Percentage,fill=color,alpha=color),stat='identity',width=1,position='stack')
p <- p + scale_x_continuous(breaks = seq(1,x_upper_lim),limits = c(0.3,x_upper_lim+0.5),expand = c(0, 0))+scale_y_continuous(breaks=y_breaks,labels = paste0(y_breaks,' %'),limits=c(-0.3,102),expand = c(0, 0))
p <- p + custom_theme + theme(legend.position = 'none',panel.grid.major.x = element_blank())
p <- p + scale_fill_manual(values = c('Hidden'=NA,'Shown'="#8c8c8c")) + scale_alpha_manual(values = c('Hidden'=0,'Shown'=1))
p <- p + xlab('Barcode occurrence') + ylab('Cumulative frequency')
ggsave(p,filename =  paste0(figure_output_directory,'waterfall_plot.pdf'),width=unit(4,units = 'inches'),height=unit(4,units = 'inches'))


#we would like to show a plot of the the barcodes that are found most frequently
#first we sort barcodes by frequency
barcodes_frequencies <- barcodes_frequencies[order(barcodes_frequencies[,"frequency"],decreasing = T),]

#we plot the ones that are above the limiting frequency
plotdat <- barcodes_frequencies[1:20,]
plotdat <- as.data.frame(plotdat)
#we order the barcodes according to frequency
plotdat[["barcode"]] <- factor(as.character(plotdat[["barcode"]]),levels = plotdat[["barcode"]])

p <- ggplot(plotdat,aes(x=barcode,y=frequency))+geom_bar(stat='identity')+custom_theme
p <- p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1))+ylab('Frequency')+theme(axis.title.x=element_blank(),axis.text.x=element_text(hjust=0.5))
p <- p + ggtitle(paste0(paste0(nrow(plotdat),' most commonly sequenced barcodes')))+theme(plot.title = element_text(hjust=0.5)) + scale_y_continuous(trans='log10',breaks=10^seq(-10,10,1),labels = to_tenth_power_labels(seq(-10,10,1)))

ggsave(p,filename = paste0(figure_output_directory,'most_common_barcodes.pdf'),height=4.5,width=6)


#finally we try to project the real number of total barcodes to generate some stats for the article
observed_different=nrow(barcodes_frequencies)
total_reads=sum(barcodes_frequencies[,"frequency"])

#we set the function we wsh to solve
func <- function(B0) {B0*(1-((B0-1)/B0)^total_reads) - observed_different} 
#and then we find solutions between the two possible bounds
B0_solution <- uniroot(func, c(observed_different, n_phages))
B0 <- B0_solution[["root"]]

print(paste0('We estimate the total number of different barcodes to be: ',round(B0)))

#this would be equivalent to an "efficiency" of 
print(paste0('This means that there are ',round(B0/n_phages,2),' different barcodes per phage'))

print(paste0('Given the difference between the number of different barcodes we observe and the projected number, we have found ',round(observed_different/B0*100,2),' % of the total barcodes'))

