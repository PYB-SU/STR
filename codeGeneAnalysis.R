

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

#' pad_df
#'
#' add a column to a data frame and pad columns if required
#'
#' @param df Existing 'data.frame'
#' @param col A column to add to 'df'
#'
#' @returns a data.frame with one column more

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


#' step_geno
#' list ref and alt alleles in the 0|1 order
#' @param i
#' @param geno indicator of ref/alt alleles position
#' @param ref reference allele
#' @param alt other allele
#' @returns a list with alleles in order ref|alt
#'
step_geno <-function(i,geno,ref,alt) {
  haps = c(as.character(ref[i]),as.character(alt[[i]]))
  if (geno[i] == "0|1") {
    return(haps)
  } else if (geno[i] == "1|0") {
    return(rev(haps))
  }
}


#' process files from analyses
#' summarise all files of STR analyses in an excel file
#' loops over length of flanking sequences
#' extract sequences having the largest coverage
#' filter on average quality and identity of flanking sequences
#' align sequences starting on 'first.codon' found the closest to the beginning flanking sequence
#' stops sequence on the last 'repeat.codon' found the closest to the closing flanking sequence
#' writes an excel file with summary information and individual sheets
#' returns a data frame for plot.seq.waterfall
#' @param DIR directory where files are located
#' @param root_dir where DIR is located
#' @param gene name of the gene directory within DIR
#' @param query list of motifs to be looked for in the repeated region. Processed in order provided.
#' @param ncommon number of motifs to be added from list of commonly found entities in repeated region
#' @param caller subdirectory
#' @param average QUALITY for read to be processed
#' @param minimum IDENTITY in flanking sequences to keep reads
#' @param reverse do the reads need to be reversed - this is for setting direction
#' @param complement do the reeds need to be complemented - and strand
#' @param first.pattern pattern to be matched corresponding to just before the repeat
#' @param last.pattern pattern to be matched corresponding to just after the repeat
#' @param by.haplotype split reads by haplotype
#'
#' @returns a list with data for plot.seq.waterfall
#'

process.one.gene<- function(DIR,root_dir, gene, caller="guppy",
                            QUALITY=20,IDENTITY=0.9, reverse=FALSE, complement=FALSE,
                            first.pattern=NULL, last.pattern=NULL,
                            by.haplotype=FALSE, method="encode",max_length=7,motifs=NULL,
                            save="yes", file_name=NULL, tol.match=10) {

  if (method=="encode") {
    if (is.null(motifs)) {
      stop("please provide motifs if using method encode for CIGARs")
    } else {
      message("motifs are looked for in order provided")
    }
  } else if (method == "tree") {
    if (max_length>10) {
      max_length=10
      message("max search tree fix at 10^10")
    }
  } else if (method=="dyn") {
    if (max_length>20) {
      max_length=20
      message("max search time fixed at 20s")
    }
  }

  base_dir=paste0(root_dir,"/",DIR,"/",caller)
  setwd(paste0(base_dir,"/",gene))

  # each df in the list will be written as a separate sheet in the XLS
    liste_df_sum = list()
    liste_df_sum[[1]]=data.frame()
    names(liste_df_sum)[1]="summary"
    liste_df_sum[[2]]=data.frame()
    names(liste_df_sum)[2]="haplotypes"
    liste_df_sum[[3]]=data.frame()
    names(liste_df_sum)[3]="distribution"
    liste_df_sum[[4]]=data.frame()
    names(liste_df_sum)[4]="longueur"
    liste_df_sum[[5]]=data.frame()
    names(liste_df_sum)[5]="sequence"
    # where we are in the list to add extra sheets
    idx.liste=6
    longueur=NULL
    title.longueur=NULL
    read.names=NULL
    read.quality=NULL

    nb.read.tot=NULL
    nb.read.match=NULL
    nb.read.qual=NULL

    message(paste("working on",gene))
    message("going through files")
    # get _all.summary with different flanking length
    liste_summary_files = dir(pattern = paste0("overlap_",gene,"_all_[0-9]+.summary"))
    if (length(liste_summary_files)==0) {
      message(paste("no files for ",gene))
    } else {
      message(paste("found ",length(liste_summary_files), "files for ",gene))
      for (sum_file in liste_summary_files) {
      # get read length
      lr = as.numeric(gsub(".summary","",gsub(paste0("overlap_",gene,"_all_"),"",sum_file)))
      tab_lr = read.table(sum_file,sep=" ",dec=".", header=F)
      nb.read.tot=c(nb.read.tot, dim(tab_lr)[1])

      #
      idx.to.keep <- c(1) # read.name
      names.to.keep <- c("name")
      #faire une table avec ce qui est utile
      idx.allele = grep("allele:",tab_lr[1,])
      if(length(idx.allele)==0) {message("allele: column not found")}
      else {
        idx.to.keep <- c(idx.to.keep, idx.allele+1)
        names.to.keep <- c(names.to.keep, "allele")
      }
      idx.avgQ = grep("avgQ:",tab_lr[1,]) # avg quality repeat
      if(length(idx.avgQ)==0) {
        message("avgQ: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.avgQ+1)
        names.to.keep <- c(names.to.keep,"avgQ")
      }
      idx.minQ = grep("minQ:",tab_lr[1,]) # min quality repeat
      if(length(idx.minQ)==0) {
        message("minQ: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.minQ+1)
        names.to.keep <- c(names.to.keep,"minQ")
      }
      idx.maxQ = grep("maxQ:",tab_lr[1,]) # max quality repeat
      if(length(idx.maxQ)==0) {
        message("maxQ: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.maxQ+1)
        names.to.keep <- c(names.to.keep,"maxQ")
      }
      idx.leftQ = grep("leftQ:",tab_lr[1,]) # avg quality left
      if(length(idx.leftQ)==0) {
        message("leftQ: column not found")
      }else {
        idx.to.keep <- c(idx.to.keep, idx.leftQ+1)
        names.to.keep <- c(names.to.keep,"leftQ")
      }
      idx.rightQ = grep("rightQ:",tab_lr[1,]) # avg quality right
      if(length(idx.rightQ)==0) {
        message("rightQ: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep,  idx.rightQ+1)
        names.to.keep <- c(names.to.keep,"rightQ")
      }
      idx.pos=grep("P:",tab_lr[1,]) # position begin/end of repeat
      if(length(idx.pos)==0) {
        message("P: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep,  idx.pos+1, idx.pos+2)
        names.to.keep <- c(names.to.keep,c("leftP","rightP"))
      }
      idx.length=grep("L:",tab_lr[1,]) # length of repeat in BP
      if(length(idx.length)==0) {
        message("L: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.length+1)
        names.to.keep <- c(names.to.keep,"length")
      }
      idx.score=grep("S:",tab_lr[1,]) # score matching alignment left/right
      if(length(idx.score)==0) {
        message("S: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep,  idx.score+1, idx.score+2)
        names.to.keep <- c(names.to.keep,c("leftS","rightS"))
      }
      idx.ALN = grep("A:",tab_lr[1,]) # score alignment left/right
      if(length(idx.ALN)==0) {
        message("Id: column not found")
      }else {
        idx.to.keep <- c(idx.to.keep,  idx.ALN+1, idx.ALN+2)
        names.to.keep <- c(names.to.keep,c("leftA","rightA"))
      }
      idx.Id = grep("Id:",tab_lr[1,]) # nombre de bases identiques left/right
      if(length(idx.Id)==0) {
        message("Id: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.Id+1,idx.Id+2)
        names.to.keep <- c(names.to.keep,c("leftI","rightI"))
      }
      idx.G = grep("G:",tab_lr[1,]) # nombre de gaps left/right
      if(length(idx.G)==0) {
        message("G: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.G+1,idx.G+2)
        names.to.keep <- c(names.to.keep,c("leftG","rightG"))
      }
      idx.T = grep("T:",tab_lr[1,]) # nombre de gaps left/right
      if(length(idx.T)==0) {
        message("T: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.T+1)
        names.to.keep <- c(names.to.keep,c("read.length"))
      }
      idx.F = grep("F:",tab_lr[1,]) # position des flanks
      if(length(idx.F)==0) {
        message("F: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.F+1,idx.F+2)
        names.to.keep <- c(names.to.keep,c("leftF","rightF"))
      }
      idx.D = grep("D:",tab_lr[1,]) # direction du read
      if(length(idx.D)==0) {
        message("D: column not found")
      } else {
        idx.to.keep <- c(idx.to.keep, idx.D+1)
        names.to.keep <- c(names.to.keep,c("direction"))
      }

      tab_lr_orig = data.frame(tab_lr[-1,idx.to.keep])
      names(tab_lr_orig)=names.to.keep
      tab_lr_orig$name=substring(tab_lr_orig$name,2)
      # filter
      tab_lr = tab_lr_orig[tab_lr_orig$leftI > lr * 0.9  & tab_lr_orig$rightI>lr * 0.9,]
      nb.read.match=c(nb.read.match,dim(tab_lr)[1])

      tab_lr = tab_lr[tab_lr$avgQ >QUALITY & !is.infinite(tab_lr$avgQ),]
      nb.read.qual=c(nb.read.qual,dim(tab_lr)[1])

      if(dim(tab_lr)[1]==0) {next}
      liste_df_sum[[idx.liste]]=tab_lr_orig
      names(liste_df_sum)[idx.liste] = paste0("overlap_",gene,"_all_",lr)
      idx.liste=idx.liste+1
      # mettre a niveau le nombre de lignes
      # read lengths and sequences
      tab_lr = tab_lr[tab_lr$length>0, ]
      #ajoute longueur
      tmp.longueur=tab_lr[,"length",drop=F]
      rownames(tmp.longueur)=tab_lr$name
      names(tmp.longueur)=paste0("length.",lr)
      tmp.quality=tab_lr[,"avgQ",drop=F]
      rownames(tmp.quality)=tab_lr$name
      names(tmp.quality)=paste0("avgQ.",lr)

      if (is.null(longueur)) {longueur=tmp.longueur}
      else {
        longueur=merge(longueur,tmp.longueur,by=0,all=T)
        rownames(longueur)=longueur$Row.names
        longueur$Row.names=NULL
      }
      if (is.null(read.quality)) {read.quality=tmp.quality}
      else {
        read.quality=merge(read.quality,tmp.quality,by=0,all=T)
        rownames(read.quality)=read.quality$Row.names
        read.quality$Row.names=NULL
      }
      #ajoute read;names
      # ajoute longueur du read
      title.longueur=c(title.longueur, lr)
    }

    message("reading haplotypes")
    #read haplotype information - if necessary
    haplotypes=NULL
    snps=data.frame("no_haplotype defined")
    if (by.haplotype) {
      if (file.exists(paste0("phased_SNPs_",gene,".vcf.gz"))) {
        require(VariantAnnotation)
        snps = readVcf(paste0("phased_SNPs_",gene,".vcf.gz"),genome = "hs1")
        GT = geno(snps)$GT
        REF = ref(snps)
        ALT=alt(snps)
        alleles_per_sample <- sapply(1:length(GT), step_geno,geno=GT,ref=REF,alt=ALT)
        snps = data.frame(ranges(snps))
        snps$H1=alleles_per_sample[1,]
        snps$H2=alleles_per_sample[2,]
        snps$nH1=NA
        snps$nH2=NA

        reads.hap = NULL
        reads.hap = read.table(paste0("haplo_",gene,".tsv"))
        names(reads.hap)=c("name","haplotype","phaseset","chr")
        haplotypes=reads.hap
        snps$nH1 = sum(reads.hap$haplotype=="H1")
        snps$nH2 = sum(reads.hap$haplotype=="H2")
      }
    }

    liste_df_sum[[2]] = snps

    message("choosing best flanking length")
    if (is.null(longueur)) {
      if (save=="yes") {
        write.csv("",file=paste0(base_dir,"/no_reads_with_quality_above_",QUALITY,"_for",gene,".csv"))
        warning(paste("no_reads_with_quality_above_",QUALITY,"_for",gene))
      }
    } else {
      longueur=longueur[,order(title.longueur)]
      read.quality=read.quality[,order(title.longueur)]

      liste_df_sum[[4]] = longueur
      # we select the number with the most reads
      nb.read.align = apply(!is.na(longueur),2,sum)
      idx.max.reads = which.max(nb.read.align)

  # we compute distribution of lengths
      max.read.length = which.max(apply(longueur,2,max,na.rm=T))
      grid.repeat = c(seq(0,20,2),seq(30,200,10),seq(250,500,50),seq(600,1500,100))
      repeat.length.grid = data.frame(len=grid.repeat[-1])
      for (i in 1:length(title.longueur)) {
        repeat.length.grid = cbind( repeat.length.grid, data.frame(table(cut(longueur[,i,drop=T], breaks=grid.repeat)))[,2])
      }
      names(repeat.length.grid)=NULL
      #summary
      library(mclust)
      model <- Mclust(unlist(as.vector(longueur))[!is.na(unlist(as.vector(longueur)))], G = 2)

      split.read.length <- exp(mean(log(model$parameters$mean)))
      #rajouter nb reads
      summary.reads=data.frame(lr=title.longueur,nb.total = nb.read.tot,nb.align = nb.read.match,nb.Q20 = nb.read.qual)
      summary.reads = summary.reads[order(title.longueur),]

      #annot_hap = read.table(file = paste0("annotation_haplo_",gene,".txt"))
      #names(annot_hap) = c("name","haplo.SNP","haplo.size")
      # compute median lengths
      med.lo=NULL
      med.hi=NULL

      rownames(haplotypes)=haplotypes$name

      med.H1=NULL
      med.H2=NULL

      for (i in 1:length(title.longueur)) {
        read.lengths = longueur[,i,drop=T]
        # compute median of short reads
        med.lo=c(med.lo,median(read.lengths[read.lengths<split.read.length],na.rm=T))
        # compute the median of long reads
        med.hi=c(med.hi,median(read.lengths[read.lengths>split.read.length],na.rm=T))

        read.lengths = longueur[,i,drop=F]
        len.hap = merge(read.lengths, haplotypes[,c("haplotype"),drop=F],by=0)
        med.H1=c(med.H1,median(len.hap[len.hap$haplotype=="H1",2],na.rm=T))
        med.H2=c(med.H2,median(len.hap[len.hap$haplotype=="H2",2],na.rm=T))
      }
      summary.reads$median.short=med.lo
      summary.reads$median.long=med.hi
      summary.reads$median.H1=med.H1
      summary.reads$median.H2=med.H2
      liste_df_sum[[1]] = summary.reads

      title.longueur=title.longueur[order(title.longueur)]

      message("computing sequences")
      ####################
      ###READ SEQUENCES
      ####################
      names.max.read = rownames(read.names[!is.na(read.names[,idx.max.reads]),])
      flank.seq.length = title.longueur[idx.max.reads]
      file.max.read=paste0("overlap_",gene,"_all_",flank.seq.length,".fastq")
      require(ShortRead)
      fastq_file = readFastq(file.max.read)
      # get read names only
      all.read.names=unlist(lapply(as.character(ShortRead::id(fastq_file)),function(x) {y=strsplit(x," "); paste0("@",y[[1]][1])}))
      all.read.names=substring(all.read.names,2)
      # match names in list
      read.names.max = rownames(longueur[!is.na(longueur[,idx.max.reads,drop=F]),])
      positions=which(all.read.names %in% read.names.max)
      filtered_read=fastq_file[positions]
      # get sequences as character
      seq = sread(filtered_read)


      #transform to a list of character strings
      seq = as.data.frame(seq)[,1]
      require(tidyr)
      # we assume that if first.codon is provided, we know what should be done
      df_seq = data.frame(read=as.character(ShortRead::id(filtered_read)))
      df_seq = separate(df_seq, "read",sep=" ", into=c("name",NA,"found",NA,"allele.charONT",NA,"qual.avg",NA,"qual.min",NA,"qual.max",NA,"qual.left",NA,"qual.right",NA,"pos.left","pos.right",NA, "repeat.length",NA,"score.left","score.right",NA,"align.left","align.right",NA,"id.left","id.right",NA,"gap.left","gap.right", NA, "read.length",NA,"flank.left","flank.right",NA,"direction"))
      df_seq$qual.avg=as.numeric(df_seq$qual.avg)
      df_seq$qual.left=as.numeric(df_seq$qual.left)
      df_seq$qual.right=as.numeric(df_seq$qual.right)
      df_seq$pos.left=as.numeric(df_seq$pos.left)
      df_seq$pos.right=as.numeric(df_seq$pos.right)
      df_seq$id.left=as.numeric(df_seq$id.left)
      df_seq$id.right=as.numeric(df_seq$id.right)
      df_seq$read.length=as.numeric(df_seq$read.length)
      df_seq$repeat.length=as.numeric(df_seq$repeat.length)
      # find

      # if we reverse, positions should be changed as well
      if (reverse) {
        seq =  as.character(reverse(DNAStringSet(seq)))
        # permute pos.
        tmp=df_seq$pos.left
        df_seq$pos.left = df_seq$read.length +1 - df_seq$pos.right
        df_seq$pos.right = df_seq$read.length +1 - df_seq$pos.left
      }
      if (complement) {
        seq =  as.character(complement(DNAStringSet(seq)))
      }

      if (!is.null(first.pattern)) {
        require(pwalign)
        # we match repeatedly for different length of primer (4,8,12)
        # then choose the start the closest to  flank.seq.length
        coords.repeat.start = NULL
        for (i in 4:nchar(first.pattern)) {
          pattern=substr(first.pattern,nchar(first.pattern)-i+1,nchar(first.pattern))
          # pairwisealignment  to seq
          pair=pairwiseAlignment(rep(pattern,length(seq)),seq, type="local")@subject@range
          coords.repeat.start=cbind(coords.repeat.start,pair@start+pair@width-1)
        }
        df_seq$start=apply(coords.repeat.start,1,
                           function(x,len){return(x[which.min(abs(x-len))])},
                           len=flank.seq.length)+1-4
        # find first.pattern closest to length of primer
        #old function
        #df_seq$start=unlist(lapply(seq, function(x,len,first.pattern){
        #                              pos = gregexpr(first.pattern,x)[[1]];
        #                                        if (pos[1]== -1) {return(NA)}
        #                                        else {pos[which.min(abs(pos-len))]}},
        #                              len=flank.seq.length,first.pattern=first.pattern))
      # find last pattern close to length of ending primer

        coords.repeat.stop = NULL
        for (i in 0:(nchar(first.pattern)-4)) {
          pattern=substr(last.pattern,1,nchar(first.pattern)-i)
          # pairwisealignment  to seq
          pair=pairwiseAlignment(rep(pattern,length(seq)),seq, type="local")@subject@range
          coords.repeat.stop=cbind(coords.repeat.stop,nchar(seq) - pair@start)
        }
        #df_seq$stop=unlist(lapply(seq, function(x,len,last.pattern) {
        #                                         pos=gregexpr(last.pattern,x)[[1]]+nchar(last.pattern)-1;
        #                                         if (pos[1]==-1) {return(NA)}
        #                                         else {pos[which.min(abs(pos - (nchar(x)-len)))]}},
        #                           len=flank.seq.length, last.pattern=last.pattern))
        df_seq$stop=nchar(seq)-apply(coords.repeat.stop,1,function(x,len){return(x[which.min(abs(x-len))])},len=flank.seq.length)+4

        df_seq$all.seq=seq
        # we get rid of sequences for which the match is not too good
        nbefore=dim(df_seq)[1]
        df_seq=df_seq[(abs(df_seq$start-flank.seq.length)<tol.match) &
                        (abs(df_seq$stop-nchar(seq)+flank.seq.length)<tol.match) &
                        (df_seq$start < df_seq$stop),]

        message(paste("removed ",nbefore-dim(df_seq)[1]," sequences with poorly localized subsequence"))
        # now extract from start to stop
        df_seq$seq = unlist(apply(df_seq[,c("all.seq","start","stop")],1,function(x) {if (any(is.na(c(x[2],x[3])))) return(NA) else {substr(x[1],start=as.numeric(x[2]),stop=as.numeric(x[3]))}}))
        # remove 4 bases on each side
        df_seq$R = df_seq$stop-df_seq$start-8+1
        df_seq=df_seq[!is.na(df_seq$seq),]
        if (dim(df_seq)[1]==0) {
          warning("all sequences have been removed. check first.codon, rev.seq, repeat.codon")
        }
      } else {
      #
        sub.seq = unlist(lapply(seq,function(x){as.character(subseq(x,start=as.numeric(title.longueur[idx.max.reads])+1, end=length(x)-title.longueur[idx.max.reads]))}))
        df_seq$seq=sub.seq
      }
      if (!is.null(haplotypes)) {
        df_seq=merge(df_seq, haplotypes ,by="name")
      } else {
        df_seq$haplotype="H0"
      }
      # remove haplotype="none"
      df_seq=df_seq[df_seq$haplotype!="none",]
      #order from longest to shortest for waterfall
      df_seq=df_seq[order(df_seq$haplotype,as.numeric(df_seq$repeat.length),decreasing=c(FALSE,TRUE),method="radix"),]
      # add haplotype longueur
      df_seq$haplotype.l = factor(df_seq$repeat.length < split.read.length,
                                  levels=c(TRUE,FALSE),
                                  labels=c("short","long"))
      #add title
      repeat.length.grid = setNames(repeat.length.grid, c("len.STR",paste0("L_",title.longueur)))
      liste_df_sum[[3]] = repeat.length.grid



#      res.seq = plot.seq.waterfall(df_seq,first.pattern=first.pattern,last.pattern=last.pattern,
      # file_name=paste0(base_dir,"/",gene,"_plot.pdf"),query = query, ncommon = ncommon, plot.order=plot.order,
#                                   first.pattern.label=first.pattern.label, last.pattern.label=last.pattern.label,
      #length.common.codons=length.common.codons)

      # merge df_seq with res.seq$qual
      message("computing CIGAR")
      df_seq$CIGAR = unlist(lapply(df_seq$seq, shortest_cigar, first.pattern=substr(first.pattern,nchar(first.pattern)-3,nchar(first.pattern)),
                                   last.pattern=substr(last.pattern,1,4),max_length=max_length, method=method, motifs=motifs))
#      df_seq$CIGAR = unlist(lapply(df_seq$seq, rle_tree_search, first.pattern=first.pattern, last.pattern=last.pattern))
      liste_df_sum[[5]]=df_seq

      # then write
      if (save=="yes") {
        require(writexl)
        if (is.null(file_name)) {
          file_name=paste0(base_dir,"/",gene,"_summary.xlsx")
        } else {
          file_name=paste0(base_dir,"/",file_name)
        }
        write_xlsx(liste_df_sum, path=file_name,col_names=TRUE)
        message(paste("wrote ",file_name))
      }
    }
    }
    setwd(paste0(root_dir))
    list(df_seq=df_seq, first.pattern=first.pattern, last.pattern=last.pattern)
  }


#' this loops over a directory/caller
#' finds all genes present
#' then call process.one.gene over each

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
    process.one.gene(DIR,root_dir, gene, caller=caller,QUALITY=QUALITY,IDENTITY=IDENTITY, reverse=FALSE, complement=FALSE)
  }
  print("done.")
  setwd(root_dir)
  #  return(df_seq)
}

#' a custom color scale for plotting STR
scale_color_custom <- function(pal, palette, n, order, alpha = 1) {
  pal <- getFromNamespace(paste0("pal_", pal), "ggsci")

  colors <- if (missing(palette)) {
    pal(alpha = alpha)(n)
  } else {
    pal(palette = palette, alpha = alpha)(n)
  }

  if (length(order) > length(colors)) {
    stop("The length of order exceeds the number of colors.", call. = FALSE)
  }
  colors <- if (!missing(order)) colors[order]

  ggplot2::scale_color_manual(values = colors)
}

#' function to split a sequence in chunks
#' @param x numero de la sequence
#' @param seq table of sequences
#' @param query the list of
#' @param first.pattern to match
#' @param last.pattern to match
#' @param first.pattern.label what to use instead of first.pattern
#' @param last.pattern.label what to use instead of last.pattern
#' @returns a data.frame with sequence chopped in small bits with coordinates
#'
#'
split.sequence.in.repeats <- function(x, seq, query,first.pattern,last.pattern, first.pattern.label, last.pattern.label) {
  # take the sequence
  seq_single <- seq[x,"seq"]
  # match first.pattern andlast.pattern
  res.first.last=data.frame(idx = c(x,x),
                            x1 = c(1,nchar(seq_single)-nchar(last.pattern)+1),
                            x2 = c(nchar(first.pattern),nchar(seq_single)),
                            seq = c(first.pattern.label,last.pattern.label))
  seq_single=paste0(paste0(rep("N",nchar(first.pattern)),collapse=""),
                    substr(seq_single,nchar(first.pattern)+1,nchar(seq_single)-nchar(last.pattern)),
                    paste0(rep("N",nchar(last.pattern)),collapse="")
  )
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
  res=rbind(res.first.last,res)
  res = res[order(res$x1),]
  res$seq_plot=res$seq

  # now we look for what remains
  if (dim(res)[1]>1) {
    res.other = data.frame(idx=x, x1=(res$x2+1)[-length(res$x2)],x2=res$x1[-1], seq=NA,seq_plot="OTHER")
    res.other$seq=unlist(apply(res.other,1,function(ro,seq){substr(seq,as.numeric(ro[2]),as.numeric(ro[3])-1)}, seq=seq_single))
  } else {
    if (dim(res)[1]==0) {
      res.other=data.frame(idx=x,x1=1,x2=nchar(seq_single),seq=seq_single,seq_plot="OTHER")
    } else {
      res.other=data.frame(idx=x,x1=c(1,min(res$x1+1,nchar(seq_single))),x2=c(res$x1-1,nchar(seq_single)),
                           seq=c(substr(seq_single,1,res$x1-1),
                                 substr(seq_single,min(res$x1+1,nchar(seq_single)),nchar(seq_single))),
                           seq_plot="OTHER")
    }
  }
  res.other = res.other[res.other$x1 < res.other$x2,]
  res <- rbind(res,res.other)
  res <- res[order(res$x1),]
  if (res[1,"x1"]!=1) {
    res=rbind(data.frame(idx=x,x1=1,x2=res[1,"x1"]-1,seq=substr(seq_single,1,res[1,"x1"]-1),seq_plot="OTHER"),res)
  }
  res$haplotype=seq[x,"haplotype"]
  res$avgQ=seq[x,"avgQ"]
  res$name=seq[x,"name"]

  res$seq_plot[is.na(res$seq)]="OTHER"
  res$seq[is.na(res$seq)]="X"

  return(res)
}

#' plot a waterfall
#' query codons are matched in order. maybe best to first match "rare" alleles then common, otherwise
#' rare may be masked if overlaping with common
#' @param seq
#' @param first.pattern
#' @param last.pattern
#' @param DIR name to
#' @param gene
#' @param root_dir
#' @param query
#' @param ncommon
#' @param plot.order
#' @param first.pattern.label
#' @param last.pattern.label
#' @param length.common.codons
#' @param save
#' @filename
#'

waterfall <- function(obj, DIR, root_dir, gene, caller, QUAL=0, query=NULL, ncommon=NULL,plot.order=NULL,
                               first.pattern.label=NULL,last.pattern.label=NULL,length.common.codons=NULL,
                      save="yes", file_name=NULL) {
  require(ggplot2)
  require(ggpubr)
  require(ggsci)

  base_dir=paste0(root_dir,"/",DIR,"/",caller)
  setwd(paste0(base_dir,"/",gene))


  seq = obj$df_seq[obj$df_seq$qual.avg > QUAL,]
  if (dim(seq)[1]==0) {
    message(paste("no reads with quality above",QUAL))
    return(NULL)
  }
  first.pattern=substr(obj$first.pattern,nchar(obj$first.pattern)-3,nchar(obj$first.pattern))
  last.pattern=substr(obj$last.pattern,1,4)

#this builds a table with column x1/x2 signalling begin/end of sequence 'seq'
# each sequence begins with a "first.pattern" and ends with a "last.pattern"
# seq is the original seq table with (name / avgQ / minQ / maxQ / R / IdL / IdR / start / stop / all.seq / seq / haplotype / phaseset/ chr)
# sequences in query are matched in order

 # get rid of sequences of length less than 3 (less than a codon)
  seq = seq[nchar(seq$seq)>nchar(first.pattern)+nchar(last.pattern),]

  nb.seq.H1 = sum(seq$haplotype=="H1")
  nb.seq.H2 = sum(seq$haplotype=="H2")

  # find common patterns of length lenght.common.codons not listed in q to add.
  # if q is empty we discover
  seq.tomatch = seq$seq
  for (q in query) {
    seq.tomatch = lapply(seq.tomatch,gsub,pattern=q,replacement=paste0(rep("N",nchar(q)),collapse=""))
  }
  if (is.null(query)) {
    length.common.codons=3
  } else {
    length.common.codons=max(nchar(query))
  }
  # now seq.tomatch is NNN at q, we make #3 chuncks and add them to q
  common.codons <- c()
  count.codons = c()
  while(length(common.codons) < 10) {
    #get most commmon codon
    codons = unlist(lapply(seq.tomatch,function(x) {unlist(strsplit(x, paste0("(?<=.{",length.common.codons,"})"), perl = TRUE))}))
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
    # select codons making up 99%
    ncommon=sum(cumsum(count.codons)<0.99 * sum(count.codons))
    ncommonsup=sum(cumsum(count.codons)<0.99 * sum(count.codons))
  }

  query <- c(query,common.codons[0:min(ncommon,length(common.codons))])
  if (is.null(query)) {query = common.codons[0:min(ncommonsup,length(common.codons))]}
  print(query)

  res <- dplyr::bind_rows(lapply(seq_along(seq$seq), split.sequence.in.repeats, seq, query=query,first.pattern=first.pattern,last.pattern=last.pattern,first.pattern.label=first.pattern.label, last.pattern.label=last.pattern.label))
  if (is.null(plot.order)) plot.order=query
  if (any(!(query %in% plot.order))) plot.order=c(plot.order, query[!(query %in% plot.order)])
  res$seq_plot <- factor(res$seq_plot,levels = c(plot.order,"START","END","OTHER"))
  # add an extra level "OTHER" to plot START/END and other

  # real plot.order
  # remove from plot.order the sequences that are not in res$seq_plot
  plot.order = plot.order[plot.order %in% res$seq_plot]

  avgQ = data.frame(unique(res[,c("name","idx","haplotype")]))
  avgQ = merge(avgQ, seq[,c("name","qual.avg")], by="name")

  # renumerote par haplotype
  min.idx.H1 = min(res$idx[res$haplotype=="H1"])
  min.idx.H2 = min(res$idx[res$haplotype=="H2"])

  res$idx2=res$idx
  res$idx2[res$haplotype=="H1"] =  res$idx2[res$haplotype=="H1"] - min.idx.H1 +1
  res$idx2[res$haplotype=="H2"] =  res$idx2[res$haplotype=="H2"] - min.idx.H2 +1

  avgQ$idx2=avgQ$idx
  avgQ$idx2[avgQ$haplotype=="H1"] =  avgQ$idx2[avgQ$haplotype=="H1"] - min.idx.H1 +1
  avgQ$idx2[avgQ$haplotype=="H2"] =  avgQ$idx2[avgQ$haplotype=="H2"] - min.idx.H2 +1


  #make subplot by haplotype
  library(ggside)

  scale_fill_custom <- function(len) {
    # add black/black/gray for start/end/other
    val=c(pal_d3(palette = c("category10", "category20", "category20b", "category20c"), alpha = 1)(len),"black","black","gray80")
    scale_fill_manual('fill', values=val )
  }

  seq.plot = ggplot(res) +geom_tile(aes(x=(x1+x2+1)/2,y=idx2+0.5,height=1,width=(x2+1-x1),fill=seq_plot)) +
    theme_bw(base_size=20) + ylab("") + xlab("Read length (bp) (+avg quality)") + scale_fill_custom(length(plot.order)) +
    facet_grid(rows=~haplotype) +
    geom_ysideline(aes(x=qual.avg,y=idx2),data=avgQ,show.legend=FALSE) +
    scale_ysidex_continuous(limits=c(0,60),breaks=seq(0,60,20)) +
    theme_ggside_bw(base_size = 8) +
    theme(ggside.axis.text.x=element_text(angle=90))
  if (is.null(file_name)) {
    file_name=paste0(base_dir,"/",gene,"Q_",QUAL,"_plot.pdf")
  }
  if (save=="yes") {
    ggsave(filename = file_name, plot = seq.plot)
  }
  return(list(seq=res,qual=avgQ))
}


run_length_encoding <- function(sequence, max_motif_length = 4) {
  encoded_sequence <- ""  # Chaîne pour stocker la séquence encodée
  i <- 1  # Indice de parcours de la séquence

  while (i <= nchar(sequence)) {
    # Déterminer la longueur maximale du motif à tester, jusqu'à 4 ou moins si en fin de séquence
    max_length <- min(max_motif_length, nchar(sequence) - i + 1)

    # Essayer les motifs de longueur 4 à 1 et vérifier les répétitions
    motif_found <- FALSE
    for (motif_length in max_length:1) {
      motif <- substr(sequence, i, i + motif_length - 1)
      j <- i + motif_length  # Début de la vérification des répétitions

      # Compter le nombre de répétitions du motif
      repeat_count <- 1
      while (j + motif_length - 1 <= nchar(sequence) && substr(sequence, j, j + motif_length - 1) == motif) {
        repeat_count <- repeat_count + 1
        j <- j + motif_length
      }

      # Si on a trouvé au moins une répétition, on encode le motif et son nombre de répétitions
      if (repeat_count > 1 || motif_length == 1) {
        encoded_sequence <- paste0(encoded_sequence, motif, repeat_count)
        i <- j  # Avancer l'indice
        motif_found <- TRUE
        break
      }
    }

    # Si aucun motif n'a été trouvé, avancer d'un caractère
    if (!motif_found) {
      encoded_sequence <- paste0(encoded_sequence, substr(sequence, i, i), "1")
      i <- i + 1
    }
  }
  if (nchar(encoded_sequence)<nchar(sequence)) {
    return(encoded_sequence)
  }  else {
    return(sequence)
  }
}


split_string <- function(string, length = 2) {
  # Generate starting positions for each substring of the given length
  starts <- seq(1, nchar(string), by = length)

  # Extract substrings using `substring()`
  substrings <- sapply(starts, function(i) substring(string, i, i + length - 1))

  return(substrings)
}


shortest_rle <- function(sequence) {
  # Function to perform RLE for a given pattern length
  rle_with_pattern_length <- function(seq, pattern_length) {
    result <- ""
    i <- 1
    while (i <= nchar(seq)) {
      # Extract the pattern of the given length
      pattern <- substring(seq, i, i + pattern_length - 1)

      # Count how many times the pattern repeats consecutively
      count <- 1
      while (i + pattern_length <= nchar(seq) && substring(seq, i + pattern_length, i + 2*pattern_length - 1) == pattern) {
        count <- count + 1
        i <- i + pattern_length
      }

      # Add the pattern and its count (if more than 1) to the result
      result <- paste0(result, pattern, ifelse(count > 1, count, "1"))
      i <- i + pattern_length
    }
    return(result)
  }

  # Try pattern lengths from 4 to 1 and choose the shortest result
#  encoded_strings <- sapply(4:1, function(len) rle_with_pattern_length(sequence, len))
  encoded_strings <- sapply(4, function(len) rle_with_pattern_length(sequence, len))

  # Return the shortest encoded string
  return(encoded_strings[which.min(nchar(encoded_strings))])
}


make.shortest.CIGAR <- function(seq,first.pattern,last.pattern) {
  n <- nchar(seq)-nchar(last.pattern)
  i <- nchar(first.pattern)+1
  result <- first.pattern

  while (i <= n) {
    motif=rep(NA,6) # quel pattern
    count=rep(NA,6) # combien de repetitions
    span=rep(NA,6)
    for (length in 1:6) {
      if (i + length - 1 <= n) {
        motif[length] <- substr(seq, i, i + length - 1)
        count[length] <- 1
        span[length] <- length

        while (i + length * count[length] + length - 1 <= n &&
               substr(seq, i + length * count[length], i + length * count[length] + length - 1) == motif[length]) {
          count[length] <- count[length] + 1
          span[length] = count[length] * length
        }
      }
    }
    #now what does the best compression ?
    od = order(span/(1:6),decreasing = T,na.last = T)
    #take the best and iterate
    result = paste0(result,"+",motif[od[1]],count[od[1]])
    #set i to after
    i <- i + span[od[1]]
  }
  result <- paste0(result,"+",last.pattern)

  return(result)
}


# Fonction pour créer un arbre de recherche pour RLE
rle_tree_search <- function  (seq,first.pattern,last.pattern) {
  require(data.tree)
  n <- nchar(seq)-nchar(last.pattern) # maximum length of seq

  i <- nchar(first.pattern)+1
  result <- first.pattern

  # Fonction pour trouver la meilleure représentation avec un motif donné
  encode_motif <- function(motif, seq, start, motif_len, motif_before) {
    count <- 1
    while (start + motif_len * count + motif_len - 1 <= n &&
           substr(seq, start + motif_len * count, start + motif_len * count + motif_len - 1) == motif) {
      count <- count + 1
    }
    encoding=paste0(motif,ifelse(count>1,count,""))
    motif = paste0(motif_before,"+",encoding)
    return(list(rep.count = count,
                encoding=encoding,
                motif=motif,
                runlength=nchar(motif),
                len=motif_len*count))
  }

  # Créer la racine de l'arbre
  rootNode <- Node$new("rootNode")
  rootNode$motif=first.pattern
  # Fonction récursive pour construire l'arbre à partir de la position i
  build_tree <- function(node, seq, i) {
    if (i > n) return()  # Si on dépasse la séquence, on arrête
    rep.count <- rep(NA,6)
    span <- rep(NA,6)
    name <- rep(NA,6)
    for (length in 1:6) {
      if (i + length - 1 <= n) {
        motif <- substr(seq, i, i + length - 1)
        enc <- encode_motif(motif, seq, i, length,node$motif)
        # Créer un nouveau nœud pour cet encodage
        child <- node$AddChild(enc$encoding)
        child$span <- enc$len
        child$rep.count <- enc$rep.count
        child$i <- i + enc$len
        child$motif=enc$motif
        child$runlength=enc$runlength
        rep.count[length] <-  enc$rep.count
        span[length] <-  enc$len+i
        name[length] <-  enc$encoding
      }
    }
    od <- order(rep.count,decreasing = T,na.last = T)
    if (min(rep.count,na.rm = T)!=max(rep.count,na.rm = T)) {
      # one sequence is repeated more than others
      node.name.to.remove=name[-od[1]]
      for (nn in node.name.to.remove) {
        node$RemoveChild(nn)
      }
    }
    #now recurse
    for (nod in node$children) {
      build_tree(nod, seq, nod$i)
    }

  }

  # Construire l'arbre à partir de la racine
  build_tree(rootNode, seq, nchar(first.pattern)+1)

  # take the motif with the minimal length
  df.rl = ToDataFrameTable(rootNode, "runlength", "motif")
  df.rl = df.rl[order(df.rl$runlength),]

  final.motif = paste0(df.rl$motif[1],"+",last.pattern)

  return(final.motif)
}


rle_tree_search2 <- function(seq, first.pattern, last.pattern) {
  require(data.tree)

  n <- nchar(seq) - nchar(last.pattern)  # Maximum length of seq

  # Function to encode a motif
  encode_motif <- function(motif, seq, start, motif_len, motif_before) {
    count <- 0
    while (start + motif_len * count <= n &&
           substr(seq, start + motif_len * count, start + motif_len * count + motif_len - 1) == motif) {
      count <- count + 1
    }
    encoding <- paste0(motif, ifelse(count > 1, count, ""))
    motif_combined <- paste0(motif_before, "+", encoding)
    list(rep.count = count, encoding = encoding, motif = motif_combined, runlength = nchar(motif_combined), len = motif_len * count)
  }

  # Create the root of the tree
  rootNode <- Node$new("rootNode", motif = first.pattern)
  # Recursive function to build the tree
  build_tree <- function(node, seq, i) {
    if (i > n) return()  # Stop if we exceed the sequence length
    rep.counts <- integer(6)
    names <- character(6)

    for (length in 1:6) {
      if (i + length - 1 <= n) {
        motif <- substr(seq, i, i + length - 1)
        enc <- encode_motif(motif, seq, i, length, node$motif)

        child <- node$AddChild(enc$encoding)
        child$span <- enc$len
        child$rep.count <- enc$rep.count
        child$i <- i + enc$len
        child$motif <- enc$motif
        child$runlength <- enc$runlength

        rep.counts[length] <- enc$rep.count
        names[length] <- enc$encoding
      }
    }
    # Remove nodes with lower repeat counts
    if (any(rep.counts < max(rep.counts, na.rm = TRUE))) {
      lapply(names[rep.counts < max(rep.counts, na.rm = TRUE)], node$RemoveChild)
    }
    # Recurse on child nodes
    lapply(node$children, build_tree, seq = seq, i = node$i)
  }

  # Build the tree and find the minimal runlength
  build_tree(rootNode, seq, nchar(first.pattern)+1)
  final.motif <- paste0(ToDataFrameTable(rootNode, "runlength", "motif")$motif[which.min(df.rl$runlength)], "+", last.pattern)

  return(final.motif)
}


shortest_cigar <- function(seq,first.pattern,last.pattern,max_length=20,method="encode",motifs) {
  seq=trimws(seq)
  if (nchar(seq)>nchar(first.pattern)+nchar(last.pattern)) {
    if (method=="encode") {
      return(encode_repeat_region(seq,first.pattern,last.pattern,motifs))
    } else if (method=="tree") {
        return(encode_rle_tree(seq,first.pattern,last.pattern,max_length))
    } else if (method=="dyn") {
      return(encode_rle_dyn(seq,first.pattern,last.pattern,max_length))
    }
  } else {
    return("")
  }
}
