

# change to "paper"
#setwd("C:/Users/boelle/Desktop/Nanopore_ONT/paper")

setwd("C:/Users/boelle/nextcloud_SU/paper")

#liste_DIR = dir()
liste_DIR=c("OPDM_pilote03")
#
root_dir = "C:/Users/boelle/nextcloud_SU/paper"

# this function adds column 'col' to data.frame 'df'
# it pads (fills) with NA if the column is longer than pre-existing df
# example : 
# df = 1 2      col = 5
#      3 4            6
#                     7
# then result is 
# df = 1  2  5
#      3  4  6
#      NA NA 7

pad_df <- function(df, col) {
  require(dplyr)
  if(is.null(df)) {
    df=data.frame(col)
  } else if (nrow(df) < length(col)) { # pad df
    extra_rows <- length(col) - nrow(df)
    df1 <- bind_rows(df, data.frame(matrix(NA, ncol = ncol(df), nrow = extra_rows)))[1:ncol(df)]
    df <- cbind(df1,col)
  } else { # pad col
    extra_rows = nrow(df) -length(col)
    df <- cbind(df,c(col,rep(NA,extra_rows)))
  }
  return(df)
}

# written by ChatGPT 
# we used reverseComplement from Biostrings instead
reverse_complement <- function(sequence) {
  # Définir la table de complément
  complement_table <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  # Convertir la séquence en vecteur de caractères, inverser et trouver le complément pour chaque base
  complement_sequence <- sapply(strsplit(sequence, NULL)[[1]], function(base) complement_table[base])
  # Rejoindre les bases pour former la séquence complémentaire
  reverse_complement_sequence <- paste(rev(complement_sequence), collapse = "")
  return(reverse_complement_sequence)
}

# summarise all files of STR analyses in an excel file
# loops over length of flanking sequences
# extract sequences having the largest coverage
# filter on average quality and identity of flanking sequences
# align sequences starting on 'first.codon' found the closest to the beginning flanking sequence
# stops sequence on the last 'repeat.codon' found the closest to the closing flanking sequence
# writes an excel file with summary information and individual sheets 
# then calls waterfall 
process.one.gene<- function(DIR,root_dir, gene, query=NULL, ncommon=NULL, caller="guppy",
                            QUALITY=20,IDENTITY=0.9, rev.seq=FALSE, first.codon=NULL,repeat.codon=NULL) {
  #DIR : repertoire de base contenant les callers
  #CALLER : le caller (bonito, dorado..)
  #query: liste des codons a retrouver, dans l'ordre de recherche.query[1] is the repeated sequence - should be in the desired sense
  #ncommon : nombre de codons communs a montrer (apres la query), 0 si aucun 
  #QUALITY : avg quality du read
  #IDENTITY : % of identity in the left+right sequences
  #rev.seq : should sequences be reversed
  #first.codon : first codon before repeats that will be searched to cut inital sequences 
  #repeat.codon : repeat codon that will be looked for to find stop in sequences 
  
  base_dir=paste0(root_dir,"/",DIR,"/",caller)
  setwd(paste0(base_dir,"/",gene))
  
  
    liste_df_sum = list()
    liste_df_sum[[1]]=data.frame()
    names(liste_df_sum)[1]="summary"
    liste_df_sum[[2]]=data.frame()
    names(liste_df_sum)[2]="distribution"
    liste_df_sum[[3]]=data.frame()
    names(liste_df_sum)[3]="longueur"
    liste_df_sum[[4]]=data.frame()
    names(liste_df_sum)[4]="sequence"
    idx.liste=5
    longueur=NULL
    title.longueur=NULL
    read.names=NULL
    read.quality=NULL
    
    print(paste("working on",gene))
    # get _all.summary with difgferent length
    liste_summary_files = dir(pattern = paste0("overlap_",gene,"_all_[0-9]+.summary"))
    if (length(liste_summary_files)==0) {
      print(paste("no files for ",gene))
    } else {
      for (sum_file in liste_summary_files) {
      # get read length
      lr = as.numeric(gsub(".summary","",gsub(paste0("overlap_",gene,"_all_"),"",sum_file)))
      tab_lr = read.table(sum_file,sep=" ",dec=".", header=F)
      nb.read.tot=dim(tab_lr)[1]
      idx.avgQ = grep("avgQ:",tab_lr[1,])
      tab_lr = tab_lr[tab_lr[,idx.avgQ+1]>QUALITY & !is.infinite(tab_lr[,idx.avgQ+1]),]
      if(dim(tab_lr)[1]==0) {next}
      idx.colX = grep("Id:",tab_lr[1,])
      if(is.na(idx.colX)) stop("Id: column not found")
      tab_lr = tab_lr[tab_lr[idx.colX+1]>lr * 0.9,]
      tab_lr = tab_lr[tab_lr[idx.colX+2]>lr * 0.9,]
      liste_df_sum[[idx.liste]]=tab_lr
      names(liste_df_sum)[idx.liste] = paste0("overlap_",gene,"_all_",lr)
      idx.liste=idx.liste+1
      # extract la 
      if(dim(tab_lr)[1]==0) {next}
      idx.length=grep("N:",tab_lr[1,])
      # mettre a niveau le nombre de lignes
      # read lengths and sequences
      tab_lr = tab_lr[tab_lr[,idx.length]>0,]
      #ajoute longueur
      longueur=pad_df(longueur,c(tab_lr[,idx.length+1],NA))
      #ajoute read;names
      read.names=pad_df(read.names,c(tab_lr[,1],NA))
      read.quality=pad_df(read.quality,c(tab_lr[,idx.avgQ+1],NA))
      # ajoute longueur du read
      title.longueur=c(title.longueur,lr) 
    }
      
    if (!is.null(longueur)) {
        
      names(longueur)=NULL
      longueur=longueur[,order(title.longueur)]
      read.names=read.names[,order(title.longueur)]
      read.quality=read.quality[,order(title.longueur)]
    
      title.longueur=title.longueur[order(title.longueur)]
    
      longueur = setNames(longueur,paste0("L_",title.longueur))
      liste_df_sum[[3]] = longueur
  
      # we select the number with the most reads
      nb.read.align = apply(!is.na(longueur[-1,]),2,sum)
      idx.max.reads = which.max(nb.read.align)
  
  # we compute distribution of lengths
      max.read.length = which.max(apply(longueur[-1,],2,max,na.rm=T))
      grid.repeat = c(seq(0,20,2),seq(30,200,10),seq(250,500,50),seq(600,1500,100))
      repeat.length.grid = data.frame(len=grid.repeat[-1])
      for (i in 1:length(title.longueur)) {
        repeat.length.grid = cbind( repeat.length.grid, data.frame(table(cut(longueur[-1,i,drop=T], breaks=grid.repeat)))[,2])
      }
      names(repeat.length.grid)=NULL
      #summary 
      #rajouter nb reads 
      ligne.max = c("nb_read_tot",rep(nb.read.tot,length(title.longueur)))
      names(ligne.max)=NULL
      ligne.max = rbind(ligne.max,c("nb_read_rep",nb.read.align))
      names(ligne.max)=NULL
      # compute median lengths
      med.lo=c("read.short")
      med.hi=c("read.long")
      # we first get the values where it is 0
      idx.to.split.min = max(1,which(cumsum(repeat.length.grid[,2,drop=T])==0))
      idx.to.split.max = length(repeat.length.grid[,2,drop=T]) - max(which(cumsum(rev(repeat.length.grid[,2,drop=T]))==0))
    
      split.read.length = median(repeat.length.grid[idx.to.split.min:idx.to.split.max,1,drop=T])
    
      for (i in 1:length(title.longueur)) {
        read.lengths = longueur[-1,i,drop=T]
        # compute median of short reads
        med.lo=c(med.lo,median(read.lengths[read.lengths<split.read.length],na.rm=T))
        # compute the median of long reads
        med.hi=c(med.hi,median(read.lengths[read.lengths>split.read.length],na.rm=T))
      }
      ligne.max = rbind(ligne.max,med.lo)
      ligne.max = rbind(ligne.max,med.hi)
      ligne.max = data.frame(ligne.max)
      ligne.max = setNames(ligne.max, c("variable",paste0("L_",title.longueur)))
      liste_df_sum[[1]] = ligne.max
      
      ####################
      ###READ SEQUENCES
      ####################
      names.max.read = read.names[,idx.max.reads]
      flank.seq.length = title.longueur[idx.max.reads]
      file.max.read=paste0("overlap_",gene,"_all_",flank.seq.length,".fastq")
      require(ShortRead)
      fastq_file = readFastq(file.max.read)
      # get read names only
      all.read.names=unlist(lapply(as.character(ShortRead::id(fastq_file)),function(x) {y=strsplit(x," "); paste0("@",y[[1]][1])}))
      # match names in list
      positions=which(all.read.names %in% as.character(read.names[,idx.max.reads]))
      filtered_read=fastq_file[positions]
      # get sequences as character
      seq = sread(filtered_read)
      
      if (rev.seq) seq =  reverseComplement(seq)
      #transform to a list of character strings
      seq = as.data.frame(seq)[,1]
      require(tidyr)
      # we assume that if first.codon is provided, we know what should be done
      df_seq = data.frame(read=as.character(ShortRead::id(filtered_read)))
      df_seq = separate(df_seq, "read",sep=" ", into=c("name",NA,NA,NA,NA,NA,"avgQ",NA,"minQ",NA,"maxQ",NA,NA,NA,NA,NA,NA,NA,"R",NA,NA,NA,NA,NA,NA,NA,"IdL","IdR"))
      
      if (!is.null(first.codon)) {
        # find first.codon closest to length of primer
        df_seq$start=unlist(lapply(seq, function(x,len,first.codon){
                                                pos = gregexpr(first.codon,x)[[1]]; 
                                                if (pos[1]== -1) {return(NA)}
                                                else {pos[which.min(abs(pos-len))]}},
                                      len=flank.seq.length,first.codon=first.codon))
      # find last repeat codon close to length of ending primer
        df_seq$stop=unlist(lapply(seq, function(x,len,repeat.codon) {
                                                 pos=gregexpr(repeat.codon,x)[[1]];
                                                 if (pos[1]==-1) {return(NA)}
                                                 else {pos[which.min(abs(pos - (nchar(x)-len)))]+nchar(repeat.codon)-1}},
                                   len=flank.seq.length, repeat.codon=repeat.codon))

        df_seq$all.seq=seq
        # now extract from start to stop
        df_seq$seq = unlist(apply(df_seq[,c("all.seq","start","stop")],1,function(x) {if (any(is.na(c(x[2],x[3])))) return(NA) else {substr(x[1],start=x[2],stop=x[3])}}))
        df_seq$R = nchar(df_seq$seq)/nchar(repeat.codon)
        df_seq=df_seq[!is.na(df_seq$seq),]
        if (dim(df_seq)[1]==0) {
          warning("all sequences have been removed. check first.codon, rev.seq, repeat.codon")
        }
        
        names(df_seq)=c("name","avgQ","minQ","maxQ","R","IdL","IdR","start","stop","all.seq","seq")
      } else {
      #
        sub.seq = unlist(lapply(seq,function(x){as.character(subseq(x,start=as.numeric(title.longueur[idx.max.reads])+1, end=length(x)-title.longueur[idx.max.reads]))}))
        df_seq$seq=sub.seq
        names(df_seq)=c("name","avgQ","minQ","maxQ","R","IdL","IdR","seq")
      }
      #order from longest to shortest for waterfall
      df_seq=df_seq[rev(order(as.numeric(df_seq$R))),]
    
      #add title
      repeat.length.grid = setNames(repeat.length.grid, c("len.STR",paste0("L_",title.longueur)))
      liste_df_sum[[2]] = repeat.length.grid

      liste_df_sum[[4]]=df_seq

      # then write
      require(writexl)
      write_xlsx(liste_df_sum, path=paste0(base_dir,"/",gene,"_summary.xlsx"),col_names=TRUE)
      print(paste("wrote ",paste0(base_dir,"/",gene,"_summary.xlsx")))
      plot.seq.waterfall(df_seq,file_name=paste0(base_dir,"/",gene,"_plot.pdf"),query = query, ncommon = ncommon)
    }
    }
    setwd(paste0(root_dir))
  }


# this loops over a directory/caller
# finds all genes present
# then call process.one.gene over each

process.gene.files <- function(DIR,root_dir, query=NULL, ncommon=NULL, caller="guppy",QUALITY=20,IDENTITY=0.9) {
  #DIR : repertoire de base contenant les callers
  #CALLER : le caller (bonito, dorado..)
  #query: liste des codons a retrouver, dans l'ordre de recherche.query[1] is the repeated sequence - should be in the desired sense
  #ncommon : nombre de codons communs a montrer (apres la query), 0 si aucun 
  #QUALITY : avg quality du read
  #IDENTITY : % of identity in the left+right sequences
  #rev.seq : should sequences be reversed
  #first.codon : first codon before repeats that will be searched to cut inital sequences 
  
  setwd(root_dir)
  setwd(paste0(root_dir,"/",DIR,"/",caller))
# get gene dirs
  genes=  list.dirs(recursive = F,full.names = F)
# all but overlap BAMS
  genes = genes[genes!="OVERLAPS"]
# do genes
  print("found")
  print(genes)

  for (gene in genes) {
    process.one.gene(DIR,root_dir, gene, query=NULL, ncommon=NULL, caller=caller,QUALITY=QUALITY,IDENTITY=IDENTITY, rev.seq=FALSE, first.codon=NULL)
  }
  print("done.")
  setwd(root_dir)
  #  return(df_seq)
}


# plot a waterfall 
# query codons are matched in order. maybe best to first match "rare" alleles then common, otherwise
# rare may be masked if overlaping with common

plot.seq.waterfall <- function(seq, query=NULL, file_name, ncommon=NULL) {
  require(ggplot2)
  require(ggpubr)
  require(ggsci)
  
#this builds a table with column x1/x2 signalling begin/end of sequence 'seq'
# sequences in query are matched in order  
  split.sequence.in.repeats <- function(x, seq, query) {
    # take the sequence
    seq_single <- seq[x]
    # match plotted sequences in order, then replace by NNN those matched
    res <- NULL
    for (q in query) {
      coords <- gregexpr(q, seq_single) # get coords of q in seq_single
      # if no match return -1
      x1 <- coords[[1]]
      if (x1[1]!= -1) { #there are matching sequences
        x2 <- x1 + attributes(coords[[1]])$match.length - 1
        res <- rbind(res, data.frame(idx = x, x1 = x1, x2 = x2, seq = q))
        # replace with NNN
        seq_single <- gsub(q,paste0(rep("N",nchar(q)),collapse=""),seq_single)
      }
    }
    res = res[order(res$x1),]
    
    # now we look for what remains 
    if (dim(res)[1]>1) {
      res.other = data.frame(idx=x, x1=(res$x2+1)[-length(res$x2)],x2=res$x1[-1], seq="OTHER")
    } else {
      if (dim(res)[1]==0) {
        res.other=data.frame(idx=x,x1=1,x2=nchar(seq_single),seq="OTHER")
      } else {
        res.other=data.frame(idx=x,x1=c(1,min(res$x1+1,nchar(seq_single))),x2=c(res$x1-1,nchar(seq_single)),seq="OTHER")
      } 
    }
    res.other = res.other[res.other$x1 < res.other$x2,]
    res <- rbind(res,res.other)
    res <- res[order(res$x1),]
    if (res[1,"x1"]!=1) {
      res=rbind(data.frame(idx=x,x1=1,x2=res[1,"x1"]-1,seq="OTHER"),res)
    }    
    res
  }
  
  if (seq$seq[1] == "seq") seq <- seq[-1, ]
  # get rid of sequences of length less than 3 (less than a codon)
  seq = seq[nchar(seq$seq)>3,]
  
  # find common patterns of lenght 3 not listed in q to add.
  # if q is empty we discover
  seq.tomatch = seq$seq
  for (q in query) {
    seq.tomatch = lapply(seq.tomatch,gsub,pattern=q,replacement=paste0(rep("N",nchar(q)),collapse=""))
  }
  # now seq.tomatch is NNN at q, we make #3 chuncks and add them to q
  common.codons <- c()
  count.codons = c()
  while(length(common.codons) < 10) {
    #get most commmon codon
    codons = unlist(lapply(seq.tomatch,function(x) {unlist(strsplit(x, "(?<=.{3})", perl = TRUE))}))
    if (length(grep("N",codons))>0) {
      codons = codons[-grep("N",codons)]
    }
    codons = codons[nchar(codons)>2]
    if (length(codons)==0) break
    # determine most common codon
    tmp = rev(sort(table(codons)))
    most.common.codon = names(tmp)[1]
    count.codons = c(count.codons,tmp[1])
    common.codons = c(common.codons,most.common.codon)
    #substitute and loop
    seq.tomatch = lapply(seq.tomatch,gsub,pattern=most.common.codon,replacement=paste0(rep("N",nchar(most.common.codon)),collapse=""))
  }
  
  if (is.null(ncommon) & length(common.codons)>0) {
    # select codons making up 80%
    ncommon=sum(cumsum(count.codons)<0.99 * sum(count.codons))
    ncommonsup=sum(cumsum(count.codons)<0.99 * sum(count.codons))
  }
  
  query <- c(query,common.codons[0:min(ncommon,length(common.codons))])
  if (is.null(query)) {query = common.codons[0:min(ncommonsup,length(common.codons))]}
  print(query)
  
  res <- dplyr::bind_rows(lapply(seq_along(seq$seq), split.sequence.in.repeats, seq$seq, query=query))
  res$seq <- factor(res$seq,levels=c(rev(query),"OTHER"))
  
  seq.plot = ggplot(res) + geom_rect(aes(xmin=x1,xmax=x2+1,ymin=idx,ymax=idx+1,fill=seq,color=seq)) + 
    theme_bw() + ylab("Reads") + xlab("Bases") + scale_color_d3() + scale_fill_d3() 
  ggsave(filename = file_name, plot = seq.plot)
  res
}


