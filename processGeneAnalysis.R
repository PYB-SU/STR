

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


liste_DIR = dir()


liste_DIR.DM1 <-
  c("2022_GEN_DM1_pilote05_02-15235", 
    "2022_GEN_DM1_pilote05_03", "2022_GEN_DM1_pilote05_04",
    "2022_GEN_DM1_pilote05_06", "2022_GEN_DM1_pilote05_07",  
    "2022_GEN_DM1_pilote05_09", "2022_GEn_DM1_pilote06_01", "2022_GEN_DM1_pilote06_02", 
    "2022_GEN_DM1_pilote07_01-15235-homemade", "2022_GEN_DM1_pilote07_02-15235-ONT", 
    "2022_GEN_DM1_WG_pilote01",  
    "DM1_pilote03_BC1", "DM1_pilote03_BC2", "DM1_pilote04")


for (d in liste_DIR.DM1) {
  process.gene.files(d,root_dir,"guppy",QUALITY = 20, IDENTITY = 0.9, query=c("CTG", "CAG", "CCG", "CGG", "CTC" ),ncommon=1)
  process.gene.files(d,root_dir,"dorado_hac",QUALITY = 20, IDENTITY = 0.9, query=c("CTG", "CAG", "CCG", "CGG", "CTC" ),ncommon=1)
  process.gene.files(d,root_dir,"bonito_hac",QUALITY = 20, IDENTITY = 0.9, query=c("CTG", "CAG", "CCG", "CGG", "CTC" ),ncommon=1)
}

liste_DIR.OPDM <-
  c("2022_GEN_OPDM_pilote03-PAM21431", 
    "2022_GEN_OPDM_pilote06_Lib1_BC3", "2022_GEN_OPDM_pilote06_Lib1_BC4", 
    "2022_GEN_OPDM_pilote06_Lib1_BC5", "2022_GEN_OPDM_pilote06_Lib2_BC6", 
    "2022_GEN_OPDM_pilote06_Lib2_BC7", "2022_GEN_OPDM_pilote06_Lib2_BC8", 
    "2022_GEN_OPDM_pilote07_lib1", "2022_GEN_OPDM_pilote07_lib2", 
    "2022_GEN_OPDM_pilote08_Lib1_BC5", "2022_GEN_OPDM_pilote08_Lib1_BC6", 
    "2022_GEN_OPDM_pilote08_Lib1_BC7", "2022_GEN_OPDM_pilote08_Lib2_BC8", 
    "2022_GEN_OPDM_pilote08_Lib2_BC9", "2022_GEN_OPDM_pilote09_BC1", 
    "2022_GEN_OPDM_pilote09_BC2", "2022_GEN_OPDM_pilote10", "2023_GEN_OPDM_pilote05_BC1", 
    "2023_GEN_OPDM_pilote05_BC2", "2023_GEN_OPDM_pilote05_BC3", "2023_GEN_OPDM_pilote05_BC4", 
    "OPDM_pilote03", "OPDM_pilote04_BC1", "OPDM_pilote04_BC2")


for (d in liste_DIR.OPDM) {
  process.gene.files(d,root_dir,"guppy",QUALITY = 20, IDENTITY = 0.9, query=c("CGG", "CCG", "CCT"),ncommon=1)
  process.gene.files(d,root_dir,"dorado_hac",QUALITY = 20, IDENTITY = 0.9 , query=c("CGG", "CCG", "CCT"),ncommon=1)
  process.gene.files(d,root_dir,"bonito_hac",QUALITY = 20, IDENTITY = 0.9 , query=c("CGG", "CCG", "CCT"),ncommon=1)
}


# certain I know : 
# OPDM9BC2 gène Notch2NLC: Triplet mettre en 5'-3' codon GGC bleu, AGG en orange et other. Faire l'alignement à partir du 1 codon qui juxtapose la repeat. (modifié) 
process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC2",root_dir = root_dir,gene = "NOTCH2NLC", caller = "guppy",QUALITY = 20, IDENTITY = 0.9, 
                   query=c("AGG", "GGC"),ncommon=0,first.codon="CCC",rev.seq = TRUE, repeat.codon="GGC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC2",root_dir = root_dir,gene = "NOTCH2NLC", caller = "dorado_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("AGG", "GGC"),ncommon=0,first.codon="CCC",rev.seq = TRUE, repeat.codon="GGC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC2",root_dir = root_dir,gene = "NOTCH2NLC", caller = "bonito_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("AGG", "GGC"),ncommon=0,first.codon="CCC",rev.seq = TRUE, repeat.codon="GGC")

# OPDM9BC1 gène GIPC1: Triplet mettre en 3'-5' codon CCG bleu, GGC en orange et other. Faire l'alignement à partir du 1 codon qui juxtapose la repeat.
process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC1",root_dir = root_dir,gene = "GIPC1", caller = "guppy",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GGC", "CCG"),ncommon=0,first.codon="CTCC",rev.seq = TRUE, repeat.codon="CCG")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC1",root_dir = root_dir,gene = "GIPC1", caller = "dorado_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GGC", "CCG"),ncommon=0,first.codon="CTCC",rev.seq = TRUE, repeat.codon="CCG")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC1",root_dir = root_dir,gene = "GIPC1", caller = "bonito_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GGC", "CCG"),ncommon=0,first.codon="CTCC",rev.seq = TRUE, repeat.codon="CCG")

# OPDM10 gène LRP12: Triplet mettre en 3'-5' codon CCG bleu, CAG en orange et other. Faire l'alignement à partir du 1 codon qui juxtapose la repeat.
process.one.gene(DIR = "2022_GEN_OPDM_pilote10",root_dir = root_dir,gene = "LRP12", caller = "guppy",
                 QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GAC", "GCC"),ncommon=0,first.codon="GAC",rev.seq = TRUE, repeat.codon="GCC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote10",root_dir = root_dir,gene = "LRP12", caller = "dorado_hac",
                 QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GAC", "GCC"),ncommon=0,first.codon="GAC",rev.seq = TRUE, repeat.codon="GCC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote10",root_dir = root_dir,gene = "LRP12", caller = "bonito_hac",
                 QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GAC", "GCC"),ncommon=0,first.codon="GAC",rev.seq = TRUE, repeat.codon="GCC")


######################

fastq_file = read.table(paste0(root_dir,"/","2022_GEN_OPDM_pilote09_BC2","/","guppy","/","NOTCH2NLC","/","overlap_NOTCH2NLC_all_100.fastq"),sep=" ",quote="",fill = T)
desc_data=fastq_file[seq(1,dim(fastq_file)[1],4),]
seq_data=fastq_file[seq(2,dim(fastq_file)[1],4),1,drop=F]

keep_rows = which((desc_data$V7>20) & (desc_data$V27>90) & (desc_data$V28>90))
desc_data = cbind(setNames(desc_data[keep_rows,c("V1","V7")],c("read.name","qual")),setNames(seq_data[keep_rows,, drop=F],"seq"))

desc_data$rev_seq = lapply(desc_data[,"seq"], reverse_complement)
desc_data$rev_seq = unlist(lapply(desc_data[,"seq"], reverse_complement))

# find first.codon closest to length of primer
desc_data$start=unlist(lapply(desc_data$rev_seq, function(x,len){pos = gregexpr("CCC",x)[[1]]; pos[which.min(abs(pos-len))]},len=100))
# find las repeat codon close to length of ending primer
desc_data$stop=unlist(gregexpr("GGC",desc_data$rev_seq[1])[[1]][which.min(abs(gregexpr("GGC",desc_data$rev_seq[1])[[1]] - (nchar(desc_data$rev_seq[1])-100)))]+2)
# now extract from start to stop
desc_data$excerpt = unlist(apply(desc_data[,c("rev_seq","start","stop")],1,function(x) {substr(x[1],start=x[2],stop=x[3])}))

# find the CCC that is 

###################
library(BSgenome.Hsapiens.UCSC.hs1)
#repeat.coords





library(Gviz)
# make plot of methylation files.
# annotation in C:/Users/boelle/nextcloud_SU\pape/TxDb.T2T-CHM13v2.0
library(GenomicFeatures)
txdb <- loadDb("C:/Users/boelle/nextcloud_SU/paper/TxDb.T2T-CHM13v2.0")
GeneRegionTrack(txdb)


# construction of the TxDb file :
library(rtracklayer)

## Takes < 1 min, consumes about 7Gb of RAM
gff <- import("C:/Users/boelle/nextcloud_SU/paper/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz")


library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("T2T-CHM13v2.0")
seqlevels(gff) <- setNames(chrominfo$SequenceName, chrominfo$RefSeqAccn)

seqinfo(gff) <- Seqinfo(genome="T2T-CHM13v2.0")


gff = gff[gff$type!= "pseudogene",]
gff = gff[gff$type!= "miRNA",]
gff = gff[gff$type!= "snoRNA",]
gff = gff[gff$type!= "snRNA",]
gff = gff[gff$type!= "tRNA",]
gff = gff[gff$type!= "V_gene_segment",]
gff = gff[gff$type!= "C_gene_segment",]
gff = gff[gff$type!= "J_gene_segment",]
gff = gff[gff$type!= "D_gene_segment",]


#library(txdbmaker)
txdb <- makeTxDbFromGRanges(gff, taxonomyId=9606)

seqlevelsStyle(txdb) <- "NCBI"

saveDb(txdb, "C:/Users/boelle/nextcloud_SU/paper/TxDb.T2T-CHM13v2.0")

save(repeat.coords, file = "C:/Users/boelle/nextcloud_SU/paper/repeats.rda")

plot.methyl.fancy <- function(DIR,CALLER,gene,methyl, root="C:/Users/boelle/nextcloud_SU/paper/",
                              from=NULL, to=NULL, type="s", reads=c("short","long")) {
  #methyl is type of methylation, found in files
  # type.methyl 1 or 2 depends on what will be read in modkit
  if (type==1) {type.methyl=1}
  else if (type==2) {type.methyl=2}
  else if (type=="s") {type.methyl=3}
  else if (type=="b") {type.methyl=4}
  else {stop("type must be 1,2, (s)um or (b)oth")}
  
  options(ucscChromosomeNames=FALSE) # necessary when working with txdb as chromosome are not chr natively
  library(Gviz)
  library(rtracklayer)
  
  gtrack <- GenomeAxisTrack()
  # get coordinates of gene
  library(GenomicFeatures)
  
  GENE = genes(txdb, filter=list(gene_id=gene), columns=c("gene_id","tx_chrom"))
  if (length(GENE) != 1) stop("check GENE step")
  start.gene = start(GENE)
  end.gene = end(GENE)
  
  if (is.null(from)) from=start.gene
  if (is.null(to)) to=end.gene
  
  grtrack <- GeneRegionTrack(txdb,
      chromosome = as.character(chrom(GENE)), 
      start = from, end = to, 
      stacking="dense",name="gene")
  
  setwd(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/"))
  library(bsseq)
  # verifier que les fichiers existent
  for (xx in reads) {
    if (!file.exists(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/",xx,"_",gene,"_",methyl,".bedmethyl")))
      stop(paste0("file ",xx,"_",gene,"_",methyl,".bedmethyl not found"))
  }
  # quoi faire sur les methylation
  bs.methyl <- read.modkit(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/",reads,"_",gene,"_",methyl,".bedmethyl"), rmZeroCov = TRUE)

  if (type.methyl==1) {
    bs.methyl=list(bs.methyl[[1]])
    message("using 1st methylation")
  }else if (type.methyl==2) {
      if (length(bs.methyl)>1) {
        bs.methyl=list(bs.methyl[[2]])
        message("using 2nd methylation")
      } else {
        warning("requested 2nd methylation but only 1 is available - defaulting to 1st")
        bs.methyl=list(bs.methyl[[1]])
    }
  } else if (type.methyl==3) {
    if (length(bs.methyl)>1) {
      #add M of 2 to 1
      bs.methyl[[1]]@assays@data@listData$M = bs.methyl[[1]]@assays@data@listData$M +bs.methyl[[2]]@assays@data@listData$M 
      message("added 1st and 2nd methylation together")
      bs.methyl=list(bs.methyl[[1]])
    } else {
      warning("requested add methylation but only 1 methylation available - defaulting to 1st")
    }
  } else if(type.methyl=="b") {
    if (length(bs.methyl)>1) {
      #keep the 2 methylations
      bs.methyl=list(bs.methyl[[1]], bs.methyl[[2]])
    } else {
      warning("requested two methylations but only 1 methylation available - defaulting to 1st")
    }
  }
  
  mat.methyl <- NULL
  mat.cov <- NULL
  names.methyl=NULL
  methyl.sm<-list()
  for (i in 1:length(bs.methyl)) {
    # filter based on coverage - use methyl1 so that get the same ranges
    idx.cov = which(apply(getCoverage(bs.methyl[[1]]),1,sum)>ceiling((max(apply(getCoverage(bs.methyl[[1]]),1,sum))/10)))
    methyl.sm[[i]] <- BSmooth(bs.methyl[[i]][idx.cov,])
    seqlevels(methyl.sm[[i]]@rowRanges) <- as.character(chrom(GENE))
    mat.methyl = rbind( mat.methyl,
                        t(getMeth(methyl.sm[[i]], type="smooth"))  )
    names.methyl = c(names.methyl, paste0(reads,".",i))
    mat.cov = rbind(mat.cov, t(getCoverage(methyl.sm[[i]], type="Cov")))
  }
  # depending on 
  rownames(mat.methyl)=names.methyl
  rownames(mat.cov)=names.methyl

  # restrict plot to coverage of reads if set to gene length
  if ((from==start.gene) & (from < min(start(methyl.sm[[1]]@rowRanges)))) { from=min(start(methyl.sm[[1]]@rowRanges)) }
  if ((to==end.gene) & (max(end(methyl.sm[[1]]@rowRanges))<to)) {to = max(end(methyl.sm[[1]]@rowRanges)) }
  
  ctrack <- DataTrack(
    range = methyl.sm[[1]]@rowRanges,
    data = mat.cov ,
    type = "l",
    name = "coverage",
    chromosome = as.character(chrom(GENE)), 
    genome="T2T-CHM13v2.0"
  )
  
  mtrack <- DataTrack(
    range = methyl.sm[[1]]@rowRanges,
    data =  mat.methyl,
    type = "l",
    name = "methylation",
    chromosome = as.character(chrom(GENE)), 
    genome="T2T-CHM13v2.0"
  )
  
  repeat.coords = repeat.coords[repeat.coords$chrom == paste0("chr",as.character(chrom(GENE))) & repeat.coords$start > start.gene & repeat.coords$end < end.gene, ]
  
  rtrack = AnnotationTrack(
    start=repeat.coords$start, width=repeat.coords$end-repeat.coords$start,
    chromosome= as.character(chrom(GENE)),id=repeat.coords$rep,
    groupAnnotation="id", 
    name="repeats" )
  
  pdf(paste0(root,"/",DIR,"/",CALLER,"/methylation_",gene,"_",methyl,"_",type,".pdf"))
  plotTracks(list(gtrack,ctrack,mtrack,grtrack,rtrack), 
             groups=names.methyl,  
             from = from, to = to,
             fontcolor.feature=1, cex.feature=0.5)
  dev.off()  
  
}



txdb = loadDb("C:/Users/boelle/nextcloud_SU/paper/TxDb.T2T-CHM13v2.0") # -> txdb
load("C:/Users/boelle/nextcloud_SU/paper/repeats.rda") # repeat.coords


fix.collapse <- function (DIR,CALLER, gene,methyl,methyl_mod,root_dir) {
  # make sure that bedmethyl files are first combined in strands
  # this is needed when using unfiltered as these are present in +/-
  for (x in c("short","long")) {

    data=read.table(file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                x,"_",gene,"_",methyl,".bedmethyl"))
    data2 = aggregate(data[,c("V5",paste0("V",10:18))], 
                      by=list(V1=data$V1,V2=data$V2,V3=data$V3,V4=data$V4,V7=data$V7,V8=data$V8,V9=data$V9), sum)
    data2$V6="-"
    data2=data2[,paste0("V",1:18)]
    write.table(data2,file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                  x,"_",gene,"_",methyl_mod,".bedmethyl"), row.names=F,col.names=F, sep="\t",quote = FALSE)
  }
  
}

# in bedmethyl
# V10 : coverage
# V12 : mod (-> m)
# V13 : canonical
# V13 : other_mod (-> h)


split.and.collapse <- function (DIR,CALLER, gene,methyl,context,root_dir) {
  # split the bedmethyl files in different files according to h/m/a
  # collapse strands at same positions +/-
  # METHYL IS THE TYPE OF METHYL CALLED found in the name (5mC, 5mCG, 5mC_5hmC,5mCG_5hmCG,)
  # context is unfiltered / CPGfiltered
  # files are changed to m/h/a + context
  for (x in c("short","long","overlap")) {
    data=read.table(file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                x,"_",gene,"_",methyl,"_",context,".bedmethyl"))
    for (mod in unique(data$V4)) {
      # this should be h/m/a
      data2 = data[data$V4 == mod,]
      # aggregate over strands 
      data2 = aggregate(data2[,c("V5",paste0("V",10:18))], 
                        by=list(V1=data2$V1,V2=data2$V2,V3=data2$V3,V4=data2$V4,V7=data2$V7,V8=data2$V8,V9=data2$V9), sum)
      
      data2$V6="-"
      # put canonical to canonical + mod_other
      data2$V13=data2$V13 + data2$V14
      # put mod_other to 0
      data2$V14=0
      # valid_coverage
      data$V10 = data$V12 + data$V13
      # remove coverage=0
      data2=data2[data2$V10>0,]
      #reorder
      data2=data2[,paste0("V",1:18)]
      write.table(data2,file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                  x,"_",gene,"_",mod,"_",context,".bedmethyl"), row.names=F,col.names=F, sep="\t",quote = FALSE)
    }
  }
}


check.short.long <- function(DIR,CALLER, root_dir) {
  # find the "ids" files in directories containing read-id and repeat-length
  # make a table showing how many reads in each and median length
  setwd(paste0(root_dir,"/",DIR,"/",CALLER,"/"))
#  get gene dirs
  genes = list.dirs(recursive = F)
  genes = genes [genes != "./OVERLAPS"]
  res=c()
  for (gene in genes) {
    setwd(paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/"))
    liste.files = dir(pattern="(long|short|overlap)[a-zA-Z0-9_]*.ids")
    for (file in liste.files) {
      if (file.size(file)>0) {
        tab <- read.table(file)
        nb.reads = dim(tab)[1]
        len.repeat=median(tab[,2])
      } else {
        nb.reads=0
        len.repeat=NA
      }
      
      res <- rbind(res,data.frame(gene=gene, nb.reads=nb.reads, len.repeat=len.repeat,file=file,dir=DIR))
    }
    setwd(paste0(root_dir,"/",DIR,"/",CALLER,"/"))
  }
  setwd(root_dir)
  res
}


liste_DIR.DM1 <-
  c("2022_GEN_DM1_pilote05_02-15235", 
    "2022_GEN_DM1_pilote05_03", "2022_GEN_DM1_pilote05_04",
    "2022_GEN_DM1_pilote05_06", "2022_GEN_DM1_pilote05_07",  
    "2022_GEN_DM1_pilote05_09", "2022_GEN_DM1_pilote06_01", "2022_GEN_DM1_pilote06_02", 
    "2022_GEN_DM1_pilote07_01-15235-homemade", "2022_GEN_DM1_pilote07_02-15235-ONT", 
    "2022_GEN_DM1_WG_pilote01",  
    "DM1_pilote03_BC1", "DM1_pilote03_BC2", "DM1_pilote04")


#DM1
res.short.long.DM1 = c()
for (dir in liste_DIR.DM1) {
  res = check.short.long(dir,"dorado_hac", root_dir)
  res.short.long.DM1 = rbind(res.short.long.DM1,res)
}

res.short.long.OPDM =c()
for (dir in liste_DIR.OPDM) {
  res = check.short.long(dir,"dorado_hac", root_dir)
  res.short.long.OPDM = rbind(res.short.long.OPDM,res)
}

  
split.and.collapse("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "m_unfiltered",
                    root_dir, from=48594000, to= 48600000, type="1")
plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000, type="1")

split.and.collapse("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                    root_dir, from=48594000, to=48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to=48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "a_unfiltered",
                  root_dir, from=48594000, to= 48600000)

split.and.collapse("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "a_unfiltered",
                  root_dir, from=48594000, to= 48600000)

split.and.collapse("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                    root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "a_unfiltered",
                  root_dir, from=48594000, to= 48600000)

split.and.collapse("DM1_pilote03_BC1","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("DM1_pilote03_BC1","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000)

#05-06
split.and.collapse("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)


plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "m_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "h_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)


#05-07
split.and.collapse("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "m_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "h_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)






# OPDM

split.and.collapse("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                    root_dir,from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "h_CPGfiltered",
                  root_dir,from=105714000, to = 105718000)

plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir,from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "h_unfiltered",
                  root_dir,from=105714000, to = 105718000)

split.and.collapse("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12","5mCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12","5mCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                  root_dir,from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir,from=105714000, to = 105718000)

split.and.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "m_CPGfiltered",
                  root_dir, from=14620000, to=14625000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "h_CPGfiltered",
                  root_dir,from=14620000, to = 14625000)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "m_unfiltered",
                  root_dir, from=14620000, to=14625000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "h_unfiltered",
                  root_dir,from=14620000, to = 14625000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "a_unfiltered",
                  root_dir,from=14620000, to = 14625000)

split.and.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "m_CPGfiltered",
                  root_dir, from= 148518000, to=148522000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "h_CPGfiltered",
                  root_dir,from= 148518000, to=148522000)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "m_unfiltered",
                  root_dir, from= 148518000, to=148522000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "h_unfiltered",
                  root_dir,from= 148518000, to=148522000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "a_unfiltered",
                  root_dir,from= 148518000, to=148522000)


split.and.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "h_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)

plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "h_unfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "a_unfiltered",
                  root_dir, from=105714000, to = 105718000)



split.and.collapse("OPDM_pilote03","dorado_hac","LRP12","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("OPDM_pilote03","dorado_hac","LRP12","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "h_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)

plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "h_unfiltered",
                  root_dir, from=105714000, to = 105718000)





split.and.collapse("OPDM_pilote04_BC1","dorado_hac","GIPC1","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("OPDM_pilote04_BC1","dorado_hac","GIPC1","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "m_CPGfiltered",
                  root_dir, from=14620000, to = 14625000)
plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "h_CPGfiltered",
                  root_dir, from=14620000, to = 14625000)

plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "m_unfiltered",
                  root_dir, from=14620000, to = 14625000)
plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "h_unfiltered",
                  root_dir, from=14620000, to = 14625000)


split.and.collapse("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "m_CPGfiltered",
                  root_dir,from= 148518000, to=148522000)
plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "h_CPGfiltered",
                  root_dir, from= 148518000, to=148522000)

plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "m_unfiltered",
                  root_dir, from= 148518000, to=148522000)
plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "h_unfiltered",
                  root_dir, from= 148518000, to=148522000)



fix.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "6mA_collapsed",
                    root_dir)
  
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","LRP12",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","LRP12",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","LRP12",methyl = "6mA_collapsed",
                    root_dir)
  

  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","NOTCH2NLC",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","NOTCH2NLC",methyl = "6mA_collapsed",
                    root_dir)
  
  ## 9/BC2
  
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","GIPC1",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","GIPC1",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","GIPC1",methyl = "6mA_collapsed",
                    root_dir)
  
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","LRP12",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","LRP12",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","LRP12",methyl = "6mA_collapsed",
                    root_dir)
  
  
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "6mA_collapsed",
                    root_dir)
  
  # OPDM10
  plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","GIPC1",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote10","dorado_hac","GIPC1",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","GIPC1",methyl = "6mA_collapsed",
                    root_dir)
  
  plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "6mA_collapsed",
                    root_dir)
  
  
  plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  
  fix.collapse("2022_GEN_OPDM_pilote10","dorado_hac","NOTCH2NLC",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
  plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","NOTCH2NLC",methyl = "6mA_collapsed",
                    root_dir)
  
    
   #OPDM5
  
  plot.methyl.fancy("2023_GEN_OPDM_pilote05_BC4","dorado_hac","GIPC1",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)

  
  plot.methyl.fancy("2023_GEN_OPDM_pilote05_BC4","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                    root_dir)
  