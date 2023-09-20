library("stringdist")
library('R.utils')
library('DescTools')
library('R.utils')

#We set the directory for the
data_dir <- "/Users/aoesniprbiome.com/Library/CloudStorage/OneDrive-SniprBiome/projects/barcode_prediction/submission_safe/data/"
barcode_fasta_dir <- paste0(data_dir,'barcodes.fasta')
#we evaluate all the combinations of acceptable matches and return the best one
#first if there is a perfect match of a start or stop sequence, then the surrounding nts will have a match distance of 2 - also acceptable
#so we remove any 2's that come in pairs sorrounding a zero
remove_sorrounding_perfect <- function(start_seq_matches){
  if(0 %in% start_seq_matches){
    remove_index <- c()
    for(perfect in seq(1,length(start_seq_matches))[start_seq_matches==0]){
      remove_index <- c(remove_index,perfect+1,perfect-1)
    }
    start_seq_matches <- start_seq_matches[-remove_index]
    return(start_seq_matches)
    
  }else{
    return(start_seq_matches)
  }
}

#a function that generates all k-mers of a given size in a sequence
generate_kmers <- function(sequence,size){
  #sets the size to be 1 lower for eaase of later use
  size=size-1
  #we get the start index for every kmer
  start_indeces <- seq(1,nchar(sequence)-size)
  #then we generate all the k-mers of the start indeces we found above and return those k-mers
  kmers <- lapply(start_indeces,function(x,sequence,size){
    substr(x = sequence,start = x,stop = x+size)
    },'sequence'=sequence,'size'=size)
  return(unlist(kmers))
}

#a function to extract the barcode that are flanked by a start and ending sequence, within certain thresholds of error of those start and end sequences
#takes argumnts
#x:                the sequence that we wish to extract data for
#start_seq:        the sequence that appear immediately before the barcode
#end_seq:          the sequence that appear immediately after the barcode
#intended_length:  the length of the barcodes if everything went well

#it returns the extracted barcode or NA if no barcode could be extracted
extract_barcode <- function(x,start_seq,end_seq,intended_length){
  #We check for exact unambiguous start and end sequences
  start_index <- gregexpr(start_seq, x)
  end_index <- gregexpr(end_seq, x)
  
  #if both are found only once as exact matches that are only found once and the start sequence is found before the ending sequence
  #(95% of cases) we'll just return the portion between the two.
  #first checks if only one match was found for each
  if(length(start_index[[1]])==1 & start_index[[1]]!=-1 &
     length(end_index[[1]])==1 & end_index[[1]]!=-1){
    #also checks if the end index is after the start index
    if(start_index[[1]][1]<end_index[[1]][1]){
      #returns the sequence between the two
      return(substr(x,start_index[[1]]+nchar(start_seq),end_index[[1]]-1))
    }
  }
  
  
  #if we did not find an exact match maybe the match is a bit off
  #in order to check the sequence for flanking regions that are just a bit off we would like to find close matches and extract between those
  #this is easier to do if we generate k-mers first
  #We generate k-mers of the same size as the start and end sequence 
  #if the sequence is longer than the length of the end and start sequences combined, then there might be a barcode in there
  if((nchar(start_seq)+nchar(end_seq))<nchar(x)){
    #we generate all the barcodes
    kmers <- generate_kmers(sequence = x,size = max(c(nchar(start_seq),nchar(end_seq))))
    
    #if we don't have a sequence that is long enough to contain both a start and end sequence, we return an NA
  }else{
    return(NA)
  }
  
  #if there is a close match to the flanking regions then 
  #we generate the number of mismatches to the end and start sequences
  start_distances <- lapply(kmers,FUN = DescTools::StrDist,'method'='hamming','y'=start_seq)
  end_distances <- lapply(kmers,FUN = DescTools::StrDist,'method'='hamming','y'=end_seq)
  
  #we combine the results to make matching a bit easier
  edit_distances <- rbind('start'=unlist(start_distances),'end'=unlist(end_distances))

  #if either start or end doesn't have an acceptable match we return an NA
  if(max(matrixStats::rowMins(edit_distances))>3){
    return(NA)
  }

  #For this extraction we prioritize extracting a barcode of the appropriate length
  #between two sets of plausible flanking regions. We define this as an edit distance of below 3.
  #Of course we would prefer to use slightly mismatching flanking regions so later we will prioritize close matches to distant matches.
  #However if we do have a perfect match for either the end or the start sequence we 
  #use that sequence even if slightly off matches exist
  matching_indeces <- list()
  
  #we go through the starts and ends separately
  for(row in rownames(edit_distances)){
    #we get the relevant data
    distances_row <- unlist(edit_distances[row,])
    #if there is a perfect match use the perfect match
    if(0%in%distances_row){
      matching_indeces[[row]][["0"]] <- which(edit_distances[row,]==0)
      
    #we get the distances that we consider to be good enough if there are no perfect matches
    }else if(min(distances_row)<=3){
      for(i in unique(edit_distances[row,edit_distances[row,]<=3])){
        matching_indeces[[row]][[as.character(i)]] <- which(edit_distances[row,]==i)
      }
      #if there is no match that satisfies our requirements, we return an NA since no barcode could be ectracted
    }else{
      return(NA)
    }
  }
  
  #we incrementally allow more errors:
  #if we can find a set of barcodes that are between a perfect match will force the use of that flanking region
  #but we will incrementally increase the error threshold until we reach an acceptable barcode and then return that
  #IF multiple barcodes are available, well damn it all we'll just have to return both
  #we incrementally increase the the threshold to whatever is the largest distance that is found, but never above 3
  
  #there is no need to search zero, 1 and 2 if the one either the start or the end index only has a match to 3 for example, so we get where we star the search from
  index_start <- max(unlist(lapply(matching_indeces,FUN = function(x){min(as.numeric(names(x)))})))
  for(max_distance in seq(index_start,3)){
    #we get the relevant start and end indeces that lie within our error threshold
    start_indeces <- unlist(matching_indeces[["start"]][as.character(seq(0,max_distance))])
    end_indeces <- unlist(matching_indeces[["end"]][as.character(seq(0,max_distance))])
    
    #then we generate every combination of start and end index
    start_end_combinations <- expand.grid(list(start_indeces,end_indeces))
    colnames(start_end_combinations) <- c('start','end')
    
    #if there are any combinations to go through
    if(nrow(start_end_combinations)>0){
      #we filter out the ones where the start is too before ahead of the beginning
      start_end_combinations <- start_end_combinations[start_end_combinations[["end"]]>start_end_combinations[["start"]]+nchar(start_seq),]
      
      #and then we extract all potential barcodes
      potential_barcodes <- apply(start_end_combinations,1,FUN = function(coords,sequence,len_start){
        #we return the potential sequences
        return(substr(sequence,coords['start']+len_start,coords['end']-1))
      },
      'sequence'=x,
      'len_start'=nchar(start_seq))
      
      #if we get a barcode of the exact length we need we want to stop now 
      #use that as a barcode, since there is no need to use another barcode of 
      #the correct length but found using starting and ending sequences that have more mismatches
      #than nessecary
      if(sum(nchar(potential_barcodes)==intended_length)>0){
        length_match_barcodes <- potential_barcodes[nchar(potential_barcodes)==intended_length]
        #if we get multiple barcodes, we collapse them around a dash
        return(paste0(length_match_barcodes,collapse = '/'))
        
      }
    }
  }

  
  #if we've made it to there then the last set of potential barcodes was generated
  #at max edit distance but none were the length we needed them to be so
  #we'll just return the ones that are best
  
  #we generate every combination of likely barcodes
  length_diffs <- abs(nchar(potential_barcodes)-intended_length)
  
  #we get the ones that are the closest match to the intended length
  best_barcodes <- potential_barcodes[length_diffs==min(length_diffs)]
  
  if(length(best_barcodes)>0){
    #if we identify multiple barcodes of the optimal length then we just return both I guess
    return(paste0(best_barcodes,collapse = '/'))
    #if there still isn't nothing return an NA
  }else{
    return(NA)
  }
}

add_two_tables <- function(a,b){
  if(length(a)==0){
    return(b)
  }else if(length(b)==0){
    return(a)
  }
  
  temp <- wrMisc::mergeVectors(a,b)
  return(colSums(temp,na.rm = T))
}

#we set the sequence that comes before the barcode
start_seq <- "ATTATCACTT"
#and the sequence that comes after
end_seq <- "CTTATGAGGG"
#we define the intended length of the barcode. In case of multiple potential barcodes within a sequence, we will extract the one whose length is the closest to what it is supposed to
intended_length <- 21

#we start a vector for all the sequences that we could not extract a barcode for
failed_sequneces <- c()
#then we start an empty named vector for 
barcodes_table <- c()

#we get the total number of lines for a progress print
total_lines <- R.utils::countLines(barcode_fasta_dir)

#we set the umber of lines we read at a time
n_lines_per_run <- 100000

#we open the connection to the fasta file
con = file(barcode_fasta_dir, "r")

#we set a line progress counter
i=0


#we read the line and continue to load lines until we do not get any more lines
while(TRUE){
  #we read n lines of the file until there are no more lines left, then we break
  line = readLines(con, n = n_lines_per_run)
  if ( length(line) == 0 ) {
    break
  }
  
  #we get the lines that are actually sequences
  sequences <- line[!grepl('>',line)]
  
  #we run the function ectract barcode on each sequence, setting the start and end sequence, as well as the intended length
  barcodes_temp <- unlist(lapply(sequences,FUN = extract_barcode,"start_seq"=start_seq,"end_seq"=end_seq,'intended_length'=21))
  
  #we get the full sequences in which no barcode could be identified
  failed_sequneces <- c(failed_sequneces,sequences[is.na(barcodes_temp)])

  #and we get the barcodes that were obtained successfully
  barcodes_temp <- barcodes_temp[!is.na(barcodes_temp)]
  
  #since certain barcodes occour multiple times and we are interested in saving space as well as formatting things in a more reasonable format, we
  #transform obtaine barcodes into a table of times a given barcode has occoured
  barcodes_tmp_table <- table(barcodes_temp)
  
  #we add the barcodes we just obtained to the barcodes we obtained previously
  barcodes_table <- add_two_tables(barcodes_tmp_table,barcodes_table)
  
  #we print a percentage of the total file processed
  i=i+n_lines_per_run
  print(paste0(round(i/total_lines*100,2),' %'))
}
#now that the file has been read we close the connection
close(con)

#for clarity we set the times we obtained no barcode - so a start sequence immediately followed by an end sequence, as [EMPTY]
names(barcodes_table)[names(barcodes_table)==''] <- '[EMPTY]'

#we write the table as a regular .tsv file
write.table(cbind('barcode'=names(barcodes_table),'frequency'=unname(barcodes_table)),file = paste0(data_dir,'barcode_table.tsv'),sep='\t',quote = F,col.names = T,row.names=F)
#and also return the sequences we could not find a barcode in
write.table(failed_sequneces,file = paste0(data_dir,'failed_fetch_sequences.tsv'),sep = '\t',quote=F)



